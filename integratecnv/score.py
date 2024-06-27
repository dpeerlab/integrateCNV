import numpy as np

def filter_low_genes(dense_tab, exp_cells=10):
    """
    Filters out genes (columns) from the dense_tab dataframe where the number of
    cells with expression values greater than 0 is less than or equal to exp_cells.

    Parameters:
    dense_tab (pd.DataFrame): A dataframe where rows are cells and columns are genes.
    exp_cells (int): The minimum number of cells that must express a gene (expression > 0)
                     for the gene to be retained in the dataframe.

    Returns:
    pd.DataFrame: A dataframe with columns filtered based on the exp_cells threshold.
    """
    # Apply a function along the columns (axis=0) to count the number of non-zero entries
    # and filter the columns based on the condition.
    filtered_columns = dense_tab.loc[:, (dense_tab > 0).sum(axis=0) > exp_cells]

    return filtered_columns

def format_chr(chr_str):
    """
    Formats chromosome names to have a consistent format (e.g., 'chr1' -> 'chr01').
    Parameters:
        :param chr_str (str): Chromosome name.
    Returns:
        str: Formatted chromosome name.
    """
    prefix = 'chr'
    if isinstance(chr_str, int):
        number = chr_str
    elif not chr_str.startswith(prefix):
        try:
            number = int(chr_str)
        except ValueError:
            number = str(chr_str)
    else:
        number = int(chr_str[len(prefix):])
    return f'{prefix}{number:02}'

def aggregate_genes_over_regions(cna_gene_map, counts, pseudocount=0.1, verbose=False):
    """
    Aggregate gene counts over regions defined by the CNA gene map.
    Parameters:
        :param cna_gene_map: DataFrame containing regions of interest and genes present therein. Must contain columns (chr, start, end, gene_name)
        :param counts: DataFrame containing gene counts. Rows are cells and columns are genes.
    Returns:
        :return: DataFrame containing counts aggregated over regions. Rows are cells and columns are regions.
    """
    counts = counts.copy()
    counts += pseudocount
    assert 'gene_name' in cna_gene_map.columns, 'Gene name column not found in cna_gene_map. Please ensure your gene names are stored in the gene_name column.'
    assert 'chr' in cna_gene_map.columns, 'Chromosome column not found in cna_gene_map. Please ensure your chromosome names are stored in the chr column.'
    assert 'start' in cna_gene_map.columns, 'Start column not found in cna_gene_map. Please ensure your start positions are stored in the start column.'
    assert 'end' in cna_gene_map.columns, 'End column not found in cna_gene_map. Please ensure your end positions are stored in the end column.'

    cna_gene_map = cna_gene_map[['chr', 'start', 'end', 'gene_name']].set_index('gene_name')
    genes = np.intersect1d(counts.columns, cna_gene_map.index)
    if verbose:
        print(
            f'Found {len(genes)} genes in both the count matrix ({counts.shape[1]}) and the CNA gene map ({cna_gene_map.shape[0]}).')
    counts = counts[genes].T
    cna_gene_map['region'] = cna_gene_map['chr'].apply(format_chr) + ':' + cna_gene_map['start'].astype(str) + '-' + \
                             cna_gene_map['end'].astype(str)
    counts = counts.join(cna_gene_map['region'], how='inner')
    if verbose:
        print(f'Final count matrix shape contains {counts.shape[1]-1} cells and {counts.region.nunique()} regions.')
    return counts.groupby('region').sum().T


def lognorm_counts(counts_in_altered_regions, counts_in_neutral_regions):
    """
    Compute the log2 ratio of counts in altered regions to counts in neutral regions.
    :param counts_in_altered_regions: (pd.DataFrame) Counts of genes in altered regions. (cells x regions)
    :param counts_in_neutral_regions: (pd.DataFrame) Counts of genes in neutral regions. (cells x regions)
    :return: (pd.DataFrame) Log2 ratio of counts in altered regions to counts in neutral regions.
    """
    # Check that counts_in_altered_regions and counts_in_neutral_regions have the same cells (index)
    assert counts_in_altered_regions.index.equals(
        counts_in_neutral_regions.index), 'Indices of counts_in_altered_regions and counts_in_neutral_regions do not match.'

    # Compute the log2 ratio of counts in altered regions to counts in neutral regions
    ratio_log_normed = np.log2(counts_in_altered_regions.div(counts_in_neutral_regions.values, axis=0))

    return ratio_log_normed


def get_parameters(ratio_log_normed, normal_cells, verbose=True):
    """
    Compute the null mean and standard deviation for each region.
    :param ratio_log_normed: (pd.DataFrame) Log2 ratio of counts in altered regions to counts in neutral regions.(cells x regions)
    :param verbose: (bool) Whether to print progress messages.
    :return: (pd.Series, pd.Series) Null mean and standard deviation for each region.
    """
    if verbose:
        print(f'Using {len(normal_cells)} normal cells to compute the null mean and standard deviation for each region.')
    # Subset to normal cells to determine the reference mean and standard deviation for each region
    try:
        null_mean = ratio_log_normed.loc[normal_cells].mean()
        null_std = ratio_log_normed.loc[normal_cells].std()
    except KeyError:
        raise KeyError('Normal cells not found in the ratio_log_normed DataFrame index.')

    # Compute the loglikelihood for each cell and region under the null reference/normal distributions
    if verbose:
        print('Computing loglikelihood for each cell and region under the null reference/normal distributions.')

    return null_mean, null_std

def score_cnas(cell_gene_matrix, gene_means, gene_stds):
    """
    Compute the log likelihood and z-scores of each cell's gene expression profile under the assumption
    that each gene's expression follows a normal distribution.

    Parameters:
        cell_gene_matrix (pd.DataFrame): DataFrame where rows are cells and columns are genes.
        gene_means (pd.Series): Series containing the mean expression for each gene.
        gene_stds (pd.Series): Series containing the standard deviation of expression for each gene.
    Returns:
        pd.Series: A Series containing the log likelihood for each cell.
    """
    # Ensure that the gene_means and gene_stds are aligned with the cell_gene_matrix columns
    gene_means = gene_means.reindex(cell_gene_matrix.columns)
    gene_stds = gene_stds.reindex(cell_gene_matrix.columns)

    # Calculate the z-scores for each cell's gene expression
    z_scores = (cell_gene_matrix - gene_means) / gene_stds

    # Compute the log of the probability density function of the normal distribution
    log_probs = -0.5 * np.log(2 * np.pi * (gene_stds ** 2)) - (z_scores ** 2) / 2

    return log_probs, z_scores


def run_integrateCNV(ad, normal_celltypes, cna_gene_map_alt, cna_gene_map_nm, verbose=True):
    """
    Pipeline to run integrateCNV on an AnnData object.
    :param ad: (AnnData) AnnData object containing the scRNA-seq data.
    :param normal_celltypes: (list) List of celltypes to use as normal references.
    :param cna_gene_map_alt: (pd.DataFrame) DataFrame containing regions of interest and genes present therein for altered regions.
    :param cna_gene_map_nm: (pd.DataFrame) DataFrame containing regions of interest and genes present therein for neutral regions.
    :return: (pd.DataFrame, pd.DataFrame) Log likelihood scores and z-scores for each cell and region.
    """

    # We want to make sure we are using raw counts; we want to normalize only against neutral regions
    assert 'counts' in ad.layers, 'Counts layer not found in AnnData object. Please ensure your counts are stored in the .layers attribute.'
    assert 'celltype' in ad.obs, 'Celltype annotation not found in AnnData object. Please ensure your celltype annotations are stored in the .obs attribute.'

    normal_cells = ad.obs_names[ad.obs['celltype'].isin(normal_celltypes)]
    if verbose:
        print(f'Loaded RNA-seq data with {ad.n_vars} genes and {ad.n_obs} cells ({len(normal_cells)} normal cells).')

    counts = ad.to_df(layer='counts')
    # Filter lowly expressed genes
    counts = filter_low_genes(counts, exp_cells=10)

    # This step may be time-consuming; be sure the save the filtered counts to a file for future use
    counts_in_altered_regions = aggregate_genes_over_regions(cna_gene_map_alt,
                                                             counts,
                                                             pseudocount=0.01,
                                                             verbose=verbose)
    # We ignore the differences across neutral regions for now and aggregate across all neutral regions
    cna_gene_map_all_neutral = cna_gene_map_nm.copy()
    cna_gene_map_all_neutral['chr'] = 'chr00'
    cna_gene_map_all_neutral['start'] = 'neutral'
    cna_gene_map_all_neutral['end'] = 'neutral'

    counts_in_neutral_regions = aggregate_genes_over_regions(cna_gene_map_all_neutral,
                                                                       counts,
                                                                       pseudocount=0.01,
                                                                       verbose=verbose)

    ratio_log_normed = lognorm_counts(counts_in_altered_regions, counts_in_neutral_regions)
    null_mean, null_std = get_parameters(ratio_log_normed, normal_cells, verbose=verbose)

    # Compute the loglikelihood for each cell and region under the null reference/normal distributions
    lls, z_scores = score_cnas(ratio_log_normed, null_mean, null_std)
    return lls, z_scores