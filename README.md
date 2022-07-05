[![Build Status](https://app.travis-ci.com/Barski-lab/sc-seq-analysis.svg?branch=main)](https://app.travis-ci.com/Barski-lab/sc-seq-analysis)
[![Python 3.8](https://img.shields.io/badge/python-3.8-green.svg)](https://www.python.org/downloads/release/python-38/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5315021.svg)](https://doi.org/10.5281/zenodo.5315021)

# CWL toolkit for single-cell sequencing data analysis

**Notes:**
- For details on how to use the published version ***[v1.0.1](https://github.com/Barski-lab/scRNA-Seq-Analysis/tree/v1.0.1)*** of workflows for scRNA-Seq data analysis in ***[SciDAP](https://scidap.com/)*** refer to the ***[Tutorials](https://barski-lab.github.io/sc-seq-analysis/)*** page.
- For up to date workflow description see ***[Wiki](https://github.com/Barski-lab/sc-seq-analysis/wiki)*** page.

**Publications:**

- *Aizhan Surumbayeva, Michael Kotliar, Linara Gabitova-Cornell, Andrey Kartashov, Suraj Peri, Nathan Salomonis, Artem Barski, Igor Astsaturov, Preparation of mouse pancreatic tumor for single-cell RNA sequencing and analysis of the data, STAR Protocols, Volume 2, Issue 4, 2021, 100989, ISSN 2666-1667,
https://doi.org/10.1016/j.xpro.2021.100989*

**Minimum software requirements:**
- [cwltool](https://github.com/common-workflow-language/cwltool) or alternative CWL runner supporting v1.0
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/) container runtime environment

--------
**How to use it**

This repository contains R [scripts](./scripts/sc_tools), CWL [tools](./tools) and examples of CWL [workflows](./workflows) for single-cell RNA-Seq and Multiome data analyses.

Each R script can be run directly from the command line following the `--help` message instructions. However, to guarantee results reproducibility we containerized them and wrapped in CWL format.

CWL tools can be combined into the workflows depending on the type of input datasets and required complexity of the analysis. For example, for single-cell RNA-Seq use **1(a) – 2(a) – 3(a)** and optionally **4(a) – 5(a,b)**; for Multiome ATAC-Seq and RNA-Seq use **1(b) – 2(b) – 2(a) – 3(a) - 3(b) - 3(c)** and optionally **4(a) – 5(a,b)**.

![](https://raw.githubusercontent.com/michael-kotliar/sc-seq-analysis-wiki-data/main/readme/scheme.gif)

All CWL tools are divided into groups to cover the major steps of data analysis. For integrity reasons we recommend starting from the raw FASTQ files and use one of the Cell Ranger based pipelines from the *Data preprocessing* group. The results of these pipelines can be optionally exported into [UCSC Cell Browser](https://cellbrowser.readthedocs.io/en/master/) (see *Visualization* group).

Both *sc-rna-filter.cwl* and *sc-multiome-filter.cwl* tools use feature-barcode matrices as the main inputs. All other tools from the *scRNA-Seq*, *scATAC-Seq and Multiome*, and *Secondary analyses* groups exchange data through [RDS](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/readRDS) files.

**Data preprocessing**
| Name                                 | Description  |
|:-------------------------------------|:-------------|
| [cellranger-mkref.cwl](./tools/cellranger-mkref.cwl)                                | Builds Cell Ranger compatible reference folder from the custom genome FASTA and gene GTF annotation files |
| [cellranger-count.cwl](./tools/cellranger-count.cwl)                                | Quantifies gene expression from a single-cell RNA-Seq library |
| [cellranger-aggr.cwl](./tools/cellranger-aggr.cwl)                                  | Aggregates outputs from multiple runs of Cell Ranger Count Gene Expression |
| [cellranger-arc-mkref.cwl](./tools/cellranger-arc-mkref.cwl)                        | Builds Cell Ranger ARC compatible reference folder from the custom genome FASTA and gene GTF annotation files |
| [cellranger-arc-count.cwl](./tools/cellranger-arc-count.cwl)                        | Quantifies chromatin accessibility and gene expression from a single-cell Multiome ATAC/RNA-Seq library |
| [cellranger-arc-aggr.cwl](./tools/cellranger-arc-aggr.cwl)                          | Aggregates outputs from multiple runs of Cell Ranger ARC Count Chromatin Accessibility and Gene Expression |

**Visualization**
| Name                                 | Description  |
|:-------------------------------------|:-------------|
| [cellbrowser-build-cellranger.cwl](./tools/cellbrowser-build-cellranger.cwl)        | Exports clustering results from Cell Ranger Count Gene Expression and Cell Ranger Aggregate experiments into compatible with UCSC Cell Browser format |
| [cellbrowser-build-cellranger-arc.cwl](./tools/cellbrowser-build-cellranger-arc.cwl) | Exports clustering results from Cell Ranger ARC Count Chromatin Accessibility and Gene Expression or Cell Ranger ARC Aggregate experiments into compatible with UCSC Cell Browser format |

**scRNA-Seq**
| Name                                 | Description  |
|:-------------------------------------|:-------------|
| [sc-rna-filter.cwl](./tools/sc-rna-filter.cwl)                                      | Filters single-cell RNA-Seq datasets based on the common QC metrics |
| [sc-rna-reduce.cwl](./tools/sc-rna-reduce.cwl)                                      | Integrates multiple single-cell RNA-Seq datasets, reduces dimensionality using PCA |
| [sc-rna-cluster.cwl](./tools/sc-rna-cluster.cwl)                                    | Clusters single-cell RNA-Seq datasets, identifies gene markers |

**scATAC-Seq and Multiome**
| Name                                 | Description  |
|:-------------------------------------|:-------------|
| [sc-multiome-filter.cwl](./tools/sc-multiome-filter.cwl)                            | Filters single-cell multiome ATAC-Seq and RNA-Seq datasets based on the common QC metrics |
| [sc-atac-reduce.cwl](./tools/sc-atac-reduce.cwl)                                    | Integrates multiple single-cell ATAC-Seq datasets, reduces dimensionality using LSI |
| [sc-atac-cluster.cwl](./tools/sc-atac-cluster.cwl)                                  | Clusters single-cell ATAC-Seq datasets, identifies differentially accessible peaks |
| [sc-wnn-cluster.cwl](./tools/sc-wnn-cluster.cwl)                                    | Clusters multiome ATAC-Seq and RNA-Seq datasets, identifies gene markers and differentially accessible peaks |

**Secondary analyses**
| Name                                 | Description  |
|:-------------------------------------|:-------------|
| [sc-ctype-assign.cwl](./tools/sc-ctype-assign.cwl)                                  | Assigns cell types for clusters based on the provided metadata file |
| [sc-rna-de-pseudobulk.cwl](./tools/sc-rna-de-pseudobulk.cwl)                        | Identifies differentially expressed genes between groups of cells coerced to pseudobulk datasets |
| [sc-rna-da-cells.cwl](./tools/sc-rna-da-cells.cwl)                                  | Detects cell subpopulations with differential abundance between datasets split by biological condition |

**Utilities**
| Name                                 | Description  |
|:-------------------------------------|:-------------|
| [tar-extract.cwl](./tools/tar-extract.cwl)                                          | Extracts the content of TAR file into a folder |
| [tar-compress.cwl](./tools/tar-compress.cwl)                                        | Creates compressed TAR file from a folder |

**Workflow examples for scRNA-Seq analysis**
| Name                                 | Description  |
|:-------------------------------------|:-------------|
| [sc-ref-indices-wf.cwl](./workflows/sc-ref-indices-wf.cwl) | Builds a Cell Ranger and Cell Ranger ARC compatible reference folders from the custom genome FASTA and gene GTF annotation files |
| [sc-rna-align-wf.cwl](./workflows/sc-rna-align-wf.cwl) | Runs Cell Ranger Count to quantify gene expression from a single-cell RNA-Seq library |
| [sc-rna-aggregate-wf.cwl](./workflows/sc-rna-aggregate-wf.cwl) | Aggregates gene expression data from multiple Single-cell RNA-Seq Alignment experiments |
| [sc-rna-analyze-wf.cwl](./workflows/sc-rna-analyze-wf.cwl) | Runs filtering, normalization, scaling, integration (optionally) and clustering for a single or aggregated single-cell RNA-Seq datasets |

**Workflows examples for Multiome analysis**
| Name                                 | Description  |
|:-------------------------------------|:-------------|
| [sc-multiome-align-wf.cwl](./workflows/sc-multiome-align-wf.cwl) | Runs Cell Ranger ARC Count to quantifies chromatin accessibility and gene expression from a single-cell Multiome ATAC and RNA-Seq library |
| [sc-multiome-aggregate-wf.cwl](./workflows/sc-multiome-aggregate-wf.cwl) | Aggregates data from multiple Single-cell Multiome ATAC and RNA-Seq Alignment experiments |
| [sc-multiome-analyze-wf.cwl](./workflows/sc-multiome-analyze-wf.cwl) | Runs filtering, normalization, scaling, integration (optionally) and clustering for a single or aggregated single-cell Multiome ATAC-Seq and RNA-Seq datasets |