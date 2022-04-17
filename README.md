[![Build Status](https://app.travis-ci.com/Barski-lab/scRNA-Seq-Analysis.svg?branch=main)](https://app.travis-ci.com/Barski-lab/scRNA-Seq-Analysis)
[![Python 3.8](https://img.shields.io/badge/python-3.8-green.svg)](https://www.python.org/downloads/release/python-38/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5315021.svg)](https://doi.org/10.5281/zenodo.5315021)

# CWL toolkit for single-cell data analysis

**Notes:**
- For details on how to use the published workflows version ***[v1.0.1](https://github.com/Barski-lab/scRNA-Seq-Analysis/tree/v1.0.1)*** in ***[SciDAP](https://scidap.com/)*** refer to the ***[Tutorials](https://barski-lab.github.io/sc-seq-analysis/)*** page.
- For up to date workflow description see Wiki page.

**Publications:**

- *Aizhan Surumbayeva, Michael Kotliar, Linara Gabitova-Cornell, Andrey Kartashov, Suraj Peri, Nathan Salomonis, Artem Barski, Igor Astsaturov, Preparation of mouse pancreatic tumor for single-cell RNA sequencing and analysis of the data, STAR Protocols, Volume 2, Issue 4, 2021, 100989, ISSN 2666-1667,
https://doi.org/10.1016/j.xpro.2021.100989*

--------

This repository contains R scripts and CWL tools for single-cell RNA-Seq and Multiome data analyses. The CWL tools can be chained together into the workflows as shown on the figure below.

![](./docs/images/readme/figure_1.png)

**Complete list of the available CWL tools**
| Name                                 | Description  |
|:-------------------------------------|:-------------|
| cellranger-mkref.cwl                 | Builds Cell Ranger compatible reference folder from the custom genome FASTA and gene GTF annotation files |
| cellranger-count.cwl                 | Quantifies gene expression from a single-cell RNA-Seq library |
| cellranger-aggr.cwl                  | Aggregates outputs from multiple runs of Cell Ranger Count Gene Expression |
| cellbrowser-build-cellranger.cwl     | Exports clustering results from Cell Ranger Count Gene Expression and Cell Ranger Aggregate experiments into compatible with UCSC Cell Browser format |
| cellranger-arc-mkref.cwl             | Builds Cell Ranger ARC compatible reference folder from the custom genome FASTA and gene GTF annotation files |
| cellranger-arc-count.cwl             | Quantifies chromatin accessibility and gene expression from a single-cell Multiome ATAC/RNA-Seq library |
| cellranger-arc-aggr.cwl              | Aggregates outputs from multiple runs of Cell Ranger ARC Count Chromatin Accessibility and Gene Expression |
| cellbrowser-build-cellranger-arc.cwl | Exports clustering results from Cell Ranger ARC Count Chromatin Accessibility and Gene Expression or Cell Ranger ARC Aggregate experiments into compatible with UCSC Cell Browser format |
| sc-rna-filter.cwl                    | Filters single-cell RNA-Seq datasets based on the common QC metrics |
| sc-rna-reduce.cwl                    | Integrates multiple single-cell RNA-Seq datasets, reduces dimensionality using PCA |
| sc-rna-cluster.cwl                   | Clusters single-cell RNA-Seq datasets, identifies gene markers |
| sc-multiome-filter.cwl               | Filters single-cell multiome ATAC and RNA-Seq datasets based on the common QC metrics |
| sc-atac-reduce.cwl                   | Integrates multiple single-cell ATAC-Seq datasets, reduces dimensionality using LSI |
| sc-atac-cluster.cwl                  | Clusters single-cell ATAC-Seq datasets, identifies differentially accessible peaks |
| sc-wnn-cluster.cwl                   | Clusters multiome ATAC and RNA-Seq datasets, identifies gene markers and differentially accessible peaks |
| tar-extract.cwl                      | Extracts the content of TAR file into a folder |
| tar-compress.cwl                     | Creates compressed TAR file from a folder |

**Running from command line**

All CWL files from this repository are compatible with any workflow management system or runner that implements CWL v1.0 standard (see the list [here](https://www.commonwl.org/#Implementations)). As an example, we will use [cwltool](https://github.com/common-workflow-language/cwltool) (the reference implementation) to show how to get the list of input parameters for any CWL file in order to run it from the command line. Additionally, we will show how to generate a template job definition file to be used as an alternative way of setting workflow input parameters. Please note, for a better portability and reprocibility all the tools used in our workflows are wrapped into Docker containers, thus a properly configured [Docker](https://www.docker.com/) installation is recommended.

1. To get a list of all input parameters for a CWL workflow see example below.

    ```
    cwltool sc-rna-analyze-wf.cwl --help

    usage: sc-rna-analyze-wf.cwl [-h] [--aggregation_metadata AGGREGATION_METADATA] [--barcodes_data BARCODES_DATA] [--cell_cycle_data CELL_CYCLE_DATA]
                                            --feature_bc_matrices_folder FEATURE_BC_MATRICES_FOLDER [--grouping_data GROUPING_DATA] [--highly_var_genes_count HIGHLY_VAR_GENES_COUNT]
                                            [--maximum_mito_perc MAXIMUM_MITO_PERC] [--minimum_logfc MINIMUM_LOGFC] [--minimum_pct MINIMUM_PCT] [--mito_pattern MITO_PATTERN]
                                            [--parallel_memory_limit PARALLEL_MEMORY_LIMIT] [--regress_cellcycle] [--regress_genes] [--regress_mito_perc] [--regress_rna_umi]
                                            [--rna_minimum_cells RNA_MINIMUM_CELLS] [--threads THREADS] [--vector_memory_limit VECTOR_MEMORY_LIMIT]
                                            [job_order]

    Single-cell RNA-Seq Analyze Runs filtering, normalization, scaling, integration (optionally) and clustering for a single or aggregated single-cell RNA-Seq datasets.

    positional arguments:
    job_order             Job input json file

    optional arguments:
    -h, --help            show this help message and exit
    --aggregation_metadata AGGREGATION_METADATA
                            Path to the metadata TSV/CSV file to set the datasets identities. If '--mex' points to the Cell Ranger Aggregate outputs, the aggregation.csv file can be
                            used. If input is not provided, the default dummy_metadata.csv will be used instead.
    --barcodes_data BARCODES_DATA
                            Path to the headerless TSV/CSV file with the list of barcodes to select cells of interest (one barcode per line). Prefilters input feature-barcode matrix
                            to include only selected cells. Default: use all cells.
    --cell_cycle_data CELL_CYCLE_DATA
                            Path to the TSV/CSV file with the information for cell cycle score assignment. First column - 'phase', second column 'gene_id'. If loaded Seurat object
                            already includes cell cycle scores in 'S.Score' and 'G2M.Score' metatada columns they will be removed. Default: skip cell cycle score assignment.
    --feature_bc_matrices_folder FEATURE_BC_MATRICES_FOLDER
                            Path to the compressed folder with feature-barcode matrix from Cell Ranger Count/Aggregate experiment in MEX format.
    --grouping_data GROUPING_DATA
                            Path to the TSV/CSV file to define datasets grouping. First column - 'library_id' with the values and order that correspond to the 'library_id' column
                            from the '--identity' file, second column 'condition'. Default: each dataset is assigned to its own group.
    --highly_var_genes_count HIGHLY_VAR_GENES_COUNT
                            Number of highly variable genes used in datasets integration, scaling and dimensionality reduction. Default: 3000
    --maximum_mito_perc MAXIMUM_MITO_PERC
                            Include cells with the percentage of transcripts mapped to mitochondrial genes not bigger than this value. Default: 5 (applied to all datasets)
    --minimum_logfc MINIMUM_LOGFC
                            For putative gene markers identification include only those genes that on average have log fold change difference in expression between every tested pair
                            of clusters not lower than this value. Ignored if '--diffgenes' is not set. Default: 0.25
    --minimum_pct MINIMUM_PCT
                            For putative gene markers identification include only those genes that are detected in not lower than this fraction of cells in either of the two tested
                            clusters. Ignored if '--diffgenes' is not set. Default: 0.1
    --mito_pattern MITO_PATTERN
                            Regex pattern to identify mitochondrial genes. Default: '^Mt-'
    --parallel_memory_limit PARALLEL_MEMORY_LIMIT
                            Maximum memory in GB allowed to be shared between the workers when using multiple --cpus. Default: 32
    --regress_cellcycle   Regress cell cycle scores as a confounding source of variation. Ignored if --cellcycle is not provided. Default: false
    --regress_genes       Regress genes per cell counts as a confounding source of variation. Default: false
    --regress_mito_perc   Regress the percentage of transcripts mapped to mitochondrial genes as a confounding source of variation. Default: false
    --regress_rna_umi     Regress UMI per cell counts as a confounding source of variation. Default: false
    --rna_minimum_cells RNA_MINIMUM_CELLS
                            Include only genes detected in at least this many cells. Default: 5 (applied to all datasets)
    --threads THREADS     Number of cores/cpus to use. Default: 1
    --vector_memory_limit VECTOR_MEMORY_LIMIT
                            Maximum vector memory in GB allowed to be used by R. Default: 128
    ```
2. To create a template for the job definition file see example below.

    ```
    cwltool --make-template sc-rna-analyze-wf.cwl > sc-rna-analyze-wf.yaml
    ```
    Open `sc-rna-analyze-wf.yaml` in a text editor, update input parameters, and run as follows.
    ```
    cwltool sc-rna-analyze-wf.cwl sc-rna-analyze-wf.yaml
    ```