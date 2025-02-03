# integrateCNV

integrateCNV is a Python package designed to integrate copy number variation (CNV) data from bulk DNA sequencing with single-cell RNA sequencing data to infer single-cell copy number profiles in targeted regions of the genome likely to harbor alterations.
## Key Features

- Integration of bulk DNA-seq CNV calls with scRNA-seq data
- Support for FACETS CNV output format
- Cell-level copy number inference
- Statistical analysis of CNV-expression relationships
- Visualization tools for integrated analysis

## Installation

You can install integrateCNV using pip:

```bash
pip install git+https://github.com/dpeerlab/integrateCNV.git
```

or poetry:

```bash
poetry add git+https://github.com/dpeerlab/integrateCNV.git
```

## Usage

integrateCNV requires three main types of input data:

1. **Bulk DNA Sequencing CNV Calls**
   - Format: FACETS output (_hisens.cncf.txt files)
   - Required columns:
     ```
     ID          - Sample identifier (e.g., s_RA19_10_3_s_RA19_10_11_1_hisens)
     chrom       - Chromosome number (e.g., 1)
     loc.start   - Start position (e.g., 13118)
     loc.end     - End position (e.g., 16817418)
     tcn         - Total copy number
     lcn         - Lesser copy number
     cf          - Cellular fraction
     ```
   - Additional FACETS columns (optional):
     ```
     seg         - Segment identifier
     num.mark    - Number of markers
     nhet        - Number of heterozygous positions
     cnlr.median - Copy number log-ratio median
     mafR        - Minor allele frequency ratio
     segclust    - Segment cluster
     ```

2. **Gene Annotations**
   - Format: Tab-separated BED file
   - Required columns:
     ```
     chromosome  - Chromosome identifier (e.g., chr1)
     start       - Gene start position (e.g., 29554)
     end         - Gene end position (e.g., 31109)
     gene_id     - Ensembl gene ID (e.g., ENSG00000243485)
     gene_name   - Gene symbol (e.g., MIR1302-2)
     strand      - Strand direction (+ or -)
     ```

3. **Single-cell RNA Sequencing Data**
   - Format: AnnData object (.h5ad)
   - Requirements:
     - Gene expression matrix
     - Gene annotations matching the BED file
     - Cell metadata (optional)



## Analysis Workflow

1. **Initialize and Load Data**
   ```python
   import integratecnv as cnv
   import pandas as pd
   
   # Set paths
   gene_annot_tab = "path/to/annotations.gtf.bed"
   cna_dir = "path/to/facets/output/"
   ```

2. **Process CNV Files**
   ```python
   # Find FACETS output files
   cna_paths = cnv.prepare_regions.find_files(cna_dir, "_hisens.cncf.txt")
   
   # Determing regions that are neutral in all samples from WGS data
    cna_neutral_bed = cnv.prepare_regions.get_neutral_regions(cna_paths)

    # Determine regions that contain alterations from WGS data 
    cna_alterations_bed_filtered = cnv.prepare_regions.get_altered_regions(cna_paths, gene_annot_tab, filter_inconsistent=True, gene_count_cutoff=20)
    n_altered = cna_alterations_bed_filtered.shape[0]
    print(f'Found {n_altered} altered regions.')
   ```

3. **Map CNVs to Genes**
   ```python
    # Get genes that fall in altered and neutral regions  
    cna_gene_map_alt = cnv.prepare_regions.map_gene_to_cna(cna_alterations_bed_filtered, gene_annot_tab)
    cna_gene_map_nm = cnv.prepare_regions.map_gene_to_cna(cna_neutral_bed, gene_annot_tab)
   ```

You may then run the entire pipeline by running the function `cnv.score.run_integrateCNV(ad, normal_celltypes, cna_gene_map_alt, cna_gene_map_nm)`.


The function takes the following arguments:

- `ad`: AnnData object
- `normal_celltypes`: List of cell types to use as normal cells
- `cna_gene_map_alt`: Gene-to-CNV mapping for altered regions
- `cna_gene_map_nm`: Gene-to-CNV mapping for non-altered regions

Alternatively, you can run the pipeline step by step by running the functions described in our [example notebook](https://github.com/dpeerlab/integrateCNV/blob/master/example_notebook.ipynb).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use integrateCNV in your research, please cite our paper.


