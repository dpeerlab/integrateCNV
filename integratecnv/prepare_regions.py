from pybedtools import BedTool
import tempfile
import os
import re
import pandas as pd

from matplotlib import pyplot as plt
import matplotlib.patches as patches

def find_files(directory, extension):
    """
    Find all files with a specific extension in a directory and its subdirectories.
    Parameters:
    :param directory: (str) The directory to search in.
    :param extension: (str) The extension to match.
    """
    matched_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(extension):
                matched_files.append(os.path.join(root, file))
    return matched_files

def multiintersect_cna_altered(cna_paths, filter_inconsistent=True, col_feature="tcn.em", verbose=True):
    """
    Processes multiple CNA files to identify overlapping regions that reliably contain copy number alterations.

    Parameters:
    :param cna_paths (list of str): List of file paths to CNA data.
        Each file should be a tab-separated file with columns ['chrom', 'loc.start', 'loc.end', 'tcn.em', 'lcn.em'].
        (As described by the FACETS documentation)
    :param filter_inconsistent (bool): Whether to filter out regions with inconsistent copy number alterations across
        samples.
    :param col_feature (str): Specifies the column to use for determining alterations ('tcn.em' is default). Options:
                            - tcn: This stands for "total copy number". It represents the estimated total copy number of
                                    both alleles at a given segment of the genome. This estimation integrates the tumor
                                    purity and ploidy assessed during the analysis. Essentially, it's a measure of how
                                    many copies of a particular region of DNA are present, accounting for both normal
                                    and tumor cells.
                            - tcn.em: This column, typically called "total copy number estimated by EM algorithm",
                                    represents the total copy  number estimated specifically by the Expectation-
                                    Maximization (EM) algorithm used within FACETS. The EM algorithm is used to fine-tune
                                    the estimates of tumor purity and ploidy and to assess the copy number states more
                                    accurately. The values in this column are the result of this algorithmic assessment
                                    and should be very similar or identical to those in the tcn column.
                            - cnlr: The column, "copy number log ration" is a log-ratio of observed to expected read
                                    depths at a particular genomic location. Observed read depth is the actual count of
                                    reads mapped to a genomic region in the tumor sample. Expected read depth is
                                    calculated based on a neutral copy number state, normalized for factors like
                                    GC-content and other systematic biases.

    Returns:
        DataFrame: Processed data of intersected regions with additional CNA information from each file.
    """

    # Create temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp()
    all_filtered_cna_files = []

    for idx, cna_path in enumerate(cna_paths):
        cna_data = pd.read_csv(cna_path, sep="\t")
        assert col_feature in cna_data.columns, "Unknown column feature specified"

        # Select regions with copy number alterations (excluding neutral regions)
        cna_data = cna_data[cna_data['tcn.em'] != 2]
        # Center copy number alterations around 0; cnlr is already centered
        if col_feature != 'cnlr':
            cna_data['alteration'] = cna_data[col_feature] - 2

        cna_data['chr'] = 'chr' + cna_data['chrom'].astype(str)
        cna_data.rename(columns={'loc.start': 'start', 'loc.end': 'end'}, inplace=True)
        cna_data = cna_data[['chr', 'start', 'end', 'alteration']]
        cna_data.sort_values(by=['chr', 'start', 'end'], inplace=True)
        # Remove duplicates
        cna_data.drop_duplicates(subset=['chr', 'start', 'end'], inplace=True)

        # Create BedTool object from filtered data
        bed_file = os.path.join(temp_dir, f"cna_{idx}.bed")
        cna_data.to_csv(bed_file, sep="\t", header=False, index=False,
                        columns=['chr', 'start', 'end', 'alteration'])
        all_filtered_cna_files.append(bed_file)

    # Use BedTool to perform multi-intersection across regions with CNAs
    bedtool_objects = [BedTool(f) for f in all_filtered_cna_files]
    intersection = BedTool().multi_intersect(i=[b.fn for b in bedtool_objects])
    intersected_df = intersection.to_dataframe(
        names=['chr', 'start', 'end', 'num_contributing_bed_files', 'overlapping_file_indices'] +
              [f'sample_{i}' for i in range(len(cna_paths))])
    intersected_df = intersected_df[['chr', 'start', 'end']]

    # Save intersected regions to a file
    intersected_df.to_csv(os.path.join(temp_dir, 'cna_intersections.bed'), sep="\t", index=False, header=False)
    intersection = BedTool(os.path.join(temp_dir, 'cna_intersections.bed'))

    # Get the actual alteration in each potential region from each FACETS file
    alterations_list = []
    for idx, cna_path in enumerate(cna_paths):
        # Perform intersection between the global regions and each sample FACETS file
        intersect_bed = intersection.intersect(b=f'{temp_dir}/cna_{idx}.bed', loj=True, wa=True, wb=True)
        intersect_df = intersect_bed.to_dataframe(
            names=['chrom_a', 'start_a', 'end_a', 'chrom_b', 'start_b', 'end_b', 'alteration'])
        # Remove duplicates
        intersect_df = intersect_df.drop_duplicates(subset=['chrom_a', 'start_a', 'end_a'])
        # Replace missing modifications with 0
        intersect_df['alteration'] = intersect_df['alteration'].replace({'.': 0}).astype(int)
        intersect_df.set_index(['chrom_a', 'start_a', 'end_a'], inplace=True)
        # Rename the modification column to the sample name
        intersect_df.rename(columns={'alteration': cna_path.split('/')[-1].rstrip('.cncf.txt')}, inplace=True)
        alterations_list.append(intersect_df[cna_path.split('/')[-1].rstrip('.cncf.txt')])

    alterations = pd.concat(alterations_list, axis=1)

    if verbose:
        print(f"Found {alterations.shape[0]} regions with copy number alterations across samples.")

    # Filter inconsistent data if specified; we would like all samples to agree on the direction of the alteration
    if filter_inconsistent:
        contains_amplifications = (alterations > 0).sum(axis=1) > 0  # Check if any sample has an amplification
        contains_deletions = (alterations < 0).sum(axis=1) > 0  # Check if any sample has a deletion
        inconsistent_regions = contains_amplifications & contains_deletions
        alterations = alterations[~inconsistent_regions]
        if verbose:
            print(f"Filtered out {inconsistent_regions.sum()} inconsistent regions.")
            print(f"Found {alterations.shape[0]} regions after filtering for consistent copy number alterations.")
    # Clean up temporary files and directory
    for f in all_filtered_cna_files:
        os.remove(f)
    os.remove(os.path.join(temp_dir, 'cna_intersections.bed'))
    os.listdir(temp_dir)
    os.rmdir(temp_dir)

    return alterations

def count_genes_per_region(cna_tab, gene_annot_tab):
    """
    Counts the number of genes in each region specified in the CNA table by intersecting with a gene annotation table.

    Parameters:
    :param cna_tab (pd.DataFrame): DataFrame containing the CNA data - columns corresponding to 'chr', 'start', 'end', followed by sample columns.
                                This should be the output of `multiintersect_cna_alt` - see documentation for details on expected content.
    :param gene_annot_tab (str): Path to the gene annotation table file.

    Returns:
        pd.DataFrame: DataFrame with the original CNA regions and an additional column for gene counts.
    """

    assert isinstance(cna_tab, pd.DataFrame), f"cna_tab must be a DataFrame, not of type {type(cna_tab)}"

    # Move chrom, start, end to index if they are columns
    if 'chr' in cna_tab.columns:
        cna_tab = cna_tab.set_index(['chr', 'start', 'end'])

    # Write the CNA table data to a temporary file if it's a DataFrame
    temp_cna_file = tempfile.mktemp(suffix=".bed")
    cna_tab.to_csv(temp_cna_file, sep="\t", index=True, header=False)
    cna_bed = BedTool(temp_cna_file)

    gene_annot_bed = BedTool(gene_annot_tab)
    result = cna_bed.intersect(gene_annot_bed, c=True, wa=True)
    count_tab = result.to_dataframe(names=list(cna_tab.columns) + ['count'])

    # Clean up the temporary file
    os.remove(temp_cna_file)

    return count_tab

def get_altered_regions(cna_files,
                        gene_annotations,
                        occurrence_cutoff=1,
                        gene_count_cutoff=20,
                        filter_inconsistent=True,
                        verbose=True):
    """
    Identifies regions showing CNAs in at least one sample by default, filtering based on gene count and occurrence across samples.

    Parameters:
        :param cna_files (list of str): List of paths to CNA files.
        :param gene_annotations (str): Path to the gene annotation table file.
        :param occurrence_cutoff (int, optional): Minimum number of samples in which a CNA must occur to be considered.
        :param gene_count_cutoff (int): Minimum number of genes that must be present in a region.
        :param filter_inconsistent (bool): Whether to filter out regions with inconsistent copy number alterations across
            samples.

    Returns:
        :return: pd.DataFrame: Filtered DataFrame containing CNA regions.
    """
    assert occurrence_cutoff >= 1, "Occurrence cutoff must be at least 1."

    cna_alterations_bed = multiintersect_cna_altered(cna_files, filter_inconsistent=filter_inconsistent, verbose=verbose)
    if verbose:
        print(f"Found {cna_alterations_bed.shape[0]} regions with copy number alterations across samples.")
    gene_count = count_genes_per_region(cna_alterations_bed, gene_annotations)
    contained_by_sufficient_samples = (cna_alterations_bed != 0).sum(axis=1) >= occurrence_cutoff
    contains_sufficient_genes = gene_count['count'] >= gene_count_cutoff

    # Select only regions that meet both criteria
    cna_alterations_bed_filtered = cna_alterations_bed[contained_by_sufficient_samples & contains_sufficient_genes]
    cna_alterations_bed_filtered =  cna_alterations_bed_filtered.reset_index()

    if verbose:
        print(f"Found {cna_alterations_bed_filtered.shape[0]} regions after filtering for sufficient occurrence number "
              f"and sufficient gene count.")

    # Convert the columns to the correct names
    rename_dict = dict(zip(cna_alterations_bed_filtered.columns[:3], ['chr', 'start', 'end']))
    cna_alterations_bed_filtered.rename(columns = rename_dict, inplace = True)

    return cna_alterations_bed_filtered

def get_neutral_regions(cna_paths):
    """
    Intersects multiple CNA files to identify neutral mutation regions and their overlap with CNA regions.

    Parameters:
        :param cna_paths (list of str): List of paths to CNA data files.
            Each file should be a tab-separated file with columns ['chrom', 'loc.start', 'loc.end', 'tcn.em', 'lcn.em'].
            (As described by the FACETS documentation)

    Returns:
        DataFrame: DataFrame containing the intersections with counts of overlaps for regions without CNAs in any sample
    """

    # Temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp()
    all_filtered_nm_files = []
    all_filtered_cna_files = []

    # Processing CNA files to save neutral mutation regions
    for idx, cna_path in enumerate(cna_paths):
        cna_data = pd.read_csv(cna_path, sep="\t")
        # Select only neutral mutations
        cna_data = cna_data[(cna_data['tcn.em'] == 2) & (cna_data['lcn.em'] == 1)]
        cna_data['alteration'] = cna_data['tcn.em'] - 2
        cna_data['chr'] = 'chr' + cna_data['chrom'].astype(str)
        cna_data = cna_data[['chr', 'loc.start', 'loc.end', 'alteration']]
        cna_data.rename(columns={'loc.start': 'start', 'loc.end': 'end'}, inplace=True)
        cna_data.sort_values(by=['chr', 'start', 'end'], inplace=True)
        # Remove duplicates
        cna_data.drop_duplicates(subset=['chr', 'start', 'end'], inplace=True)

        nm_bed_file = os.path.join(temp_dir, f"nm_{idx}.bed")
        cna_data.to_csv(nm_bed_file, sep="\t", index=False, header=False, columns=['chr', 'start', 'end', 'alteration'])
        all_filtered_nm_files.append(nm_bed_file)

    # Use BedTool to perform multi-intersection across regions without CNAs
    bedtool_objects = [BedTool(f) for f in all_filtered_nm_files]
    intersection = BedTool().multi_intersect(i=[b.fn for b in bedtool_objects], wa=True, wb=True)
    intersected_df = intersection.to_dataframe(
        names=['chrom', 'start', 'end', 'num_contributing_bed_files', 'overlapping_file_indices'] +
              [f'sample_{i}' for i in range(len(cna_paths))])

    # Save intersected regions to a file
    intersected_df.to_csv(os.path.join(temp_dir, 'cna_intersections.bed'), sep="\t", index=False, header=False)
    intersection = BedTool(os.path.join(temp_dir, 'cna_intersections.bed'))  # its.bed

    # Processing CNA files to save altered regions
    cna_frames = []
    for idx, cna_path in enumerate(cna_paths):
        cna_data = pd.read_csv(cna_path, sep="\t")
        # Filter for CNAs that are not neutral mutations
        cna_data = cna_data[cna_data['tcn.em'] != 2]
        cna_data['alteration'] = cna_data['tcn.em'] - 2
        cna_data['chr'] = 'chr' + cna_data['chrom'].astype(str)
        cna_data = cna_data[['chr', 'loc.start', 'loc.end', 'alteration']]
        cna_data.rename(columns={'loc.start': 'start', 'loc.end': 'end'}, inplace=True)
        cna_data.sort_values(by=['chr', 'start', 'end'], inplace=True)
        # Remove duplicates
        cna_data.drop_duplicates(subset=['chr', 'start', 'end'], inplace=True)
        cna_frames.append(cna_data)

    # Concatenate all CNA dataframes, sort and remove duplicates
    concatenated_cna_df = pd.concat(cna_frames)

    concatenated_cna_df.sort_values(by=['chr', 'start', 'end'], inplace=True)
    concatenated_cna_df.drop_duplicates(subset=['chr', 'start', 'end'], inplace=True)

    # Save concatenated and sorted CNA data to a file
    cna_bed_all = os.path.join(temp_dir, 'cna_all.bed')
    concatenated_cna_df.to_csv(cna_bed_all, sep="\t", index=False, header=False)

    # Perform intersection of the intersection results with the concatenated CNA file
    bed_with_count = BedTool(intersection).intersect(b=cna_bed_all, c=True)

    # Clean up temporary files and directory
    for f in all_filtered_cna_files:
        os.remove(f)
    for f in all_filtered_nm_files:
        os.remove(f)
    os.remove(os.path.join(temp_dir, 'cna_intersections.bed'))
    os.remove(cna_bed_all)
    os.rmdir(temp_dir)

    bed_with_count = bed_with_count.to_dataframe(
        names=['chr', 'start', 'end', 'num_contributing_bed_files', 'overlapping_file_indices'] +
              [f'sample_{i}' for i in range(len(cna_paths))] + ['count']).sort_values(by=['chr', 'start', 'end'])

    return bed_with_count[bed_with_count['num_contributing_bed_files'] == len(cna_paths)]

def map_gene_to_cna(cna_tab, gene_annot_tab):
    """
    Maps genes to CNA regions by intersecting CNA data with gene annotations,
    potentially outputting both intersecting entries from CNA and gene annotation tables.

    Parameters:
    :param cna_tab (pd.DataFrame): DataFrame containing the CNA data - columns corresponding to 'chrom', 'start', 'end', followed by sample columns.
                                This should be the output of `multiintersect_cna_alt` - see documentation for details on expected content.
    :param gene_annot_tab (str): Path to the gene annotation table file.

    Returns:
        pd.DataFrame: DataFrame containing the intersections with both CNA and gene annotation data.
    """
    assert isinstance(cna_tab, pd.DataFrame), f"cna_tab must be a DataFrame, not of type {type(cna_tab)}"
    assert 'chr' in cna_tab.columns, '`chr` column not found in cna_tab'
    assert 'start' in cna_tab.columns, '`start` column not found in cna_tab'
    assert 'end' in cna_tab.columns, '`end` column not found in cna_tab'

    # Write the CNA table data to a temporary file if it's a DataFrame
    temp_cna_file = tempfile.mktemp(suffix=".bed")
    cna_tab.to_csv(temp_cna_file, sep="\t", index=False, header=False)
    cna_bed = BedTool(temp_cna_file)

    gene_annot_bed = BedTool(gene_annot_tab)
    result = cna_bed.intersect(gene_annot_bed, wa=True, wb=True)
    gene_cna_tab = result.to_dataframe(
        names=list(cna_tab.columns) + ['chr_gene', 'start_gene', 'end_gene', 'gene_id', 'gene_name', 'strand'])

    # Clean up the temporary file
    os.remove(temp_cna_file)

    return gene_cna_tab

def plot_selected_regions(df, color='mediumseagreen'):
    """
    Plots the selected regions from the DataFrame as rectangles on a chromosome ideogram.
    Parameters:
        :param df (pd.DataFrame): DataFrame containing the selected regions with columns 'chr', 'start', 'end'.
        :param color (str): Color of the rectangles representing the selected regions.
    Returns:
        None
    """
    df = df.copy()

    # Ensure the DataFrame contains the required columns
    required_columns = {'chr', 'start', 'end'}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_columns}")

    # Remove 'chr' prefix if present and sort chromosomes
    df['chr'] = df['chr'].apply(lambda x: re.sub(r'^chr', '', str(x)))
    df['chr'] = pd.Categorical(df['chr'],
                               categories=[str(i) for i in range(1, 23)] + ['X', 'Y'],
                               ordered=True)
    df = df.sort_values('chr')

    # Get unique chromosomes
    chromosomes = df['chr'].unique()

    # Set up a figure for the plots
    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(10, len(chromosomes)))
    if len(chromosomes) == 1:
        axes = [axes]

    # Plot each chromosome
    for i, chrom in enumerate(chromosomes):
        chrom_data = df[df['chr'] == chrom]

        # Get the length of the chromosome
        chrom_length = max(chrom_data['end'])

        # Create a subplot for the chromosome
        ax = axes[i]
        ax.set_xlim(0, chrom_length)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_title(f'Chromosome {chrom}')

        # Plot the positions
        for _, row in chrom_data.iterrows():
            start = row['start']
            end = row['end']
            rect = patches.Rectangle((start, 0), end - start, 1,
                                     linewidth=1,
                                     edgecolor='none',
                                     facecolor=color)
            ax.add_patch(rect)

    plt.xlabel('Position')
    plt.tight_layout()
    plt.show()

def plot_paired_cnv_regions(file_list, highlighted_regions, color='mediumseagreen'):
    """
    Generates paired plots for specified chromosomes, showing copy number profile on top
    and highlighted regions on the bottom. These may be used to visualize altered regions
    or neutral regions inferred from copy number data.

    Only chromosomes present in the highlighted regions DataFrame will be plotted.

    Parameters:
    :param file_list (list of str): List of file paths to read the data from.
    :param df (pd.DataFrame): DataFrame containing the selected regions with columns 'chr', 'start', 'end'.
    :param chromosomes_to_plot (list of str): List of chromosomes to plot (e.g., ['1', '2', ..., '22', 'X', 'Y']).
    :param color (str): Color of the rectangles representing the selected regions.
    """

    # Ensure the DataFrame contains the required columns
    required_columns = {'chr', 'start', 'end'}
    if not required_columns.issubset(highlighted_regions.columns):
        raise ValueError(f"DataFrame must contain columns: {required_columns}")
    highlighted_regions = highlighted_regions.copy()

    # Remove 'chr' prefix if present and sort chromosomes
    highlighted_regions['chr'] = highlighted_regions['chr'].apply(lambda x: re.sub(r'^chr', '', str(x)))
    highlighted_regions['chr'] = pd.Categorical(highlighted_regions['chr'],
                               categories=[str(i) for i in range(1, 23)] + ['X', 'Y'],
                               ordered=True)
    highlighted_regions = highlighted_regions.sort_values('chr')

    # Function to clean chromosome labels
    def clean_chromosome_label(chrom):
        return re.sub(r'^chr', '', str(chrom))

    # Read data from each file and clean chromosome labels
    data_frames = []
    max_cnv = 0
    for file_path in file_list:
        df = pd.read_csv(file_path, sep="\t")
        df['chrom'] = df['chrom'].apply(clean_chromosome_label)
        data_frames.append(df)
        max_cnv = max(max_cnv, df['tcn.em'].max())

    # Plotting
    for chrom in highlighted_regions['chr'].cat.categories:
        chrom_data = highlighted_regions[highlighted_regions['chr'] == chrom]
        if chrom_data.empty:
            continue

        fig, (ax1, ax2) = plt.subplots(2, 1,
                                       figsize=(12, 4),
                                       sharex=True,
                                       gridspec_kw={'height_ratios': [6, 1]})

        # Plot chromosome positions

        # Plot the positions
        for _, row in chrom_data.iterrows():
            start = row['start']
            end = row['end']
            rect = patches.Rectangle((start, 0), end - start, 1,
                                     linewidth=1,
                                     edgecolor='none',
                                     facecolor=color)
            ax2.add_patch(rect)
        ax2.set_title(f'Chromosome {chrom} - Positions')
        ax2.set_xlabel('Position')
        ax2.grid(True)

        for i, df in enumerate(data_frames):
            # Filter data for the current chromosome
            chrom_data = df[df['chrom'] == chrom]
            if chrom_data.empty:
                continue

            # Prepare data for step plot
            positions = chrom_data[['loc.start', 'loc.end']].values.flatten()
            copy_numbers = chrom_data['tcn.em'].repeat(2)

            # Create step plot
            ax1.step(positions, copy_numbers, where='post', label=f'{file_list[i].split("/")[-1]}')

        ax1.set_title(f'Chromosome {chrom} - Copy Number')
        ax1.set_ylabel('Copy Number')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax1.grid(True)

        plt.tight_layout()
        plt.show()


