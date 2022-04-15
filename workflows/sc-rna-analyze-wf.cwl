cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement


inputs:

  feature_bc_matrices_folder:
    type: File
    doc: |
      Path to the compressed folder with feature-barcode matrix from Cell Ranger Count/Aggregate
      experiment in MEX format.

  aggregation_metadata:
    type: File?
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities. If '--mex' points to
      the Cell Ranger Aggregate outputs, the aggregation.csv file can be used. If input is not
      provided, the default dummy_metadata.csv will be used instead.

  grouping_data:
    type: File?
    doc: |
      Path to the TSV/CSV file to define datasets grouping. First column - 'library_id'
      with the values and order that correspond to the 'library_id' column from the
      '--identity' file, second column 'condition'.
      Default: each dataset is assigned to its own group.

  barcodes_data:
    type: File?
    doc: |
      Path to the headerless TSV/CSV file with the list of barcodes to select
      cells of interest (one barcode per line). Prefilters input feature-barcode
      matrix to include only selected cells.
      Default: use all cells.

  rna_minimum_cells:
    type: int?
    doc: |
      Include only genes detected in at least this many cells.
      Default: 5 (applied to all datasets)

  minimum_genes:
    type:
    - "null"
    - int
    - int[]
    doc: |
      Include cells where at least this many genes are detected. If multiple values
      provided, each of them will be applied to the correspondent dataset from the
      '--mex' input based on the '--identity' file.
      Default: 250 (applied to all datasets)

  maximum_genes:
    type:
    - "null"
    - int
    - int[]
    doc: |
      Include cells with the number of genes not bigger than this value. If multiple
      values provided, each of them will be applied to the correspondent dataset from
      the '--mex' input based on the '--identity' file.
      Default: 5000 (applied to all datasets)

  rna_minimum_umi:
    type:
    - "null"
    - int
    - int[]
    doc: |
      Include cells where at least this many UMI (transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file.
      Default: 500 (applied to all datasets)

  minimum_novelty_score:
    type:
    - "null"
    - float
    - float[]
    doc: |
      Include cells with the novelty score not lower than this value, calculated for
      as log10(genes)/log10(UMI). If multiple values provided, each of them will
      be applied to the correspondent dataset from the '--mex' input based on the
      '--identity' file.
      Default: 0.8 (applied to all datasets)

  mito_pattern:
    type: string?
    doc: |
      Regex pattern to identify mitochondrial genes.
      Default: '^Mt-'

  maximum_mito_perc:
    type: float?
    doc: |
      Include cells with the percentage of transcripts mapped to mitochondrial
      genes not bigger than this value.
      Default: 5 (applied to all datasets)

  cell_cycle_data:
    type: File?
    doc: |
      Path to the TSV/CSV file with the information for cell cycle score assignment.
      First column - 'phase', second column 'gene_id'. If loaded Seurat object already
      includes cell cycle scores in 'S.Score' and 'G2M.Score' metatada columns they will
      be removed.
      Default: skip cell cycle score assignment.

  highly_var_genes_count:
    type: int?
    doc: |
      Number of highly variable genes used in datasets integration, scaling and
      dimensionality reduction.
      Default: 3000

  regress_mito_perc:
    type: boolean?
    doc: |
      Regress the percentage of transcripts mapped to mitochondrial genes as a
      confounding source of variation.
      Default: false

  regress_rna_umi:
    type: boolean?
    doc: |
      Regress UMI per cell counts as a confounding source of variation.
      Default: false

  regress_genes:
    type: boolean?
    doc: |
      Regress genes per cell counts as a confounding source of variation.
      Default: false

  regress_cellcycle:
    type: boolean?
    doc: |
      Regress cell cycle scores as a confounding source of variation.
      Ignored if --cellcycle is not provided.
      Default: false

  dimensions:
    type:
    - "null"
    - int
    - int[]
    doc: |
      Dimensionality to use in UMAP projection and when constructing nearest-neighbor
      graph before clustering (from 1 to 50). If single value N is provided, use from
      1 to N dimensions. If multiple values are provided, subset to only selected
      dimensions.
      Default: from 1 to 10

  resolution:
    type:
    - "null"
    - float
    - float[]
    doc: |
      Clustering resolution applied to the constructed nearest-neighbor graph.
      Can be set as an array.
      Default: 0.3, 0.5, 1.0

  genes_of_interest:
    type:
    - "null"
    - string
    - string[]
    doc: |
      Genes of interest to build genes expression plots.
      Default: None

  minimum_logfc:
    type: float?
    doc: |
      For putative gene markers identification include only those genes that
      on average have log fold change difference in expression between every
      tested pair of clusters not lower than this value. Ignored if '--diffgenes'
      is not set.
      Default: 0.25

  minimum_pct:
    type: float?
    doc: |
      For putative gene markers identification include only those genes that
      are detected in not lower than this fraction of cells in either of the
      two tested clusters. Ignored if '--diffgenes' is not set.
      Default: 0.1

  parallel_memory_limit:
    type: int?
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32

  vector_memory_limit:
    type: int?
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Default: 128

  threads:
    type: int?
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  raw_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_1_2_qc_mtrcs_pca_plot_png
    doc: |
      PC1 and PC2 from the QC metrics PCA (not filtered).
      PNG format

  raw_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_2_3_qc_mtrcs_pca_plot_png
    doc: |
      PC2 and PC3 from the QC metrics PCA (not filtered).
      PNG format

  raw_cells_count_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_cells_count_plot_png
    doc: |
      Number of cells per dataset (not filtered).
      PNG format

  raw_umi_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_umi_dnst_plot_png
    doc: |
      UMI per cell density (not filtered).
      PNG format

  raw_gene_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_gene_dnst_plot_png
    doc: |
      Genes per cell density (not filtered).
      PNG format

  raw_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_gene_umi_corr_plot_png
    doc: |
      Genes vs UMI per cell correlation (not filtered).
      PNG format

  raw_mito_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_mito_dnst_plot_png
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (not filtered).
      PNG format

  raw_nvlt_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_nvlt_dnst_plot_png
    doc: |
      Novelty score per cell density (not filtered).
      PNG format

  raw_qc_mtrcs_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_qc_mtrcs_dnst_plot_png
    doc: |
      QC metrics per cell density (not filtered).
      PNG format

  raw_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_umi_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition UMI per cell density (not filtered).
      PNG format

  raw_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_gene_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition genes per cell density (not filtered).
      PNG format

  raw_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_mito_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (not filtered).
      PNG format

  raw_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_nvlt_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the novelty score per cell density (not filtered).
      PNG format

  fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_1_2_qc_mtrcs_pca_plot_png
    doc: |
      PC1 and PC2 from the QC metrics PCA (filtered).
      PNG format

  fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_2_3_qc_mtrcs_pca_plot_png
    doc: |
      PC2 and PC3 from the QC metrics PCA (filtered).
      PNG format

  fltr_cells_count_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_cells_count_plot_png
    doc: |
      Number of cells per dataset (filtered).
      PNG format

  fltr_umi_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_umi_dnst_plot_png
    doc: |
      UMI per cell density (filtered).
      PNG format

  fltr_gene_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_gene_dnst_plot_png
    doc: |
      Genes per cell density (filtered).
      PNG format

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_gene_umi_corr_plot_png
    doc: |
      Genes vs UMI per cell correlation (filtered).
      PNG format

  fltr_mito_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_mito_dnst_plot_png
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (filtered).
      PNG format

  fltr_nvlt_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_nvlt_dnst_plot_png
    doc: |
      Novelty score per cell density (filtered).
      PNG format

  fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_qc_mtrcs_dnst_plot_png
    doc: |
      QC metrics per cell density (filtered).
      PNG format

  fltr_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_umi_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition UMI per cell density (filtered).
      PNG format

  fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_gene_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition genes per cell density (filtered).
      PNG format

  fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_mito_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (filtered).
      PNG format

  fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_nvlt_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the novelty score per cell density (filtered).
      PNG format

  sc_rna_filter_stdout_log:
    type: File
    outputSource: sc_rna_filter/stdout_log
    doc: |
      stdout log generated by sc_rna_filter step

  sc_rna_filter_stderr_log:
    type: File
    outputSource: sc_rna_filter/stderr_log
    doc: |
      stderr log generated by sc_rna_filter step

  elbow_plot_png:
    type: File?
    outputSource: sc_rna_reduce/elbow_plot_png
    doc: |
      Elbow plot (from cells PCA).
      PNG format

  qc_dim_corr_plot_png:
    type: File?
    outputSource: sc_rna_reduce/qc_dim_corr_plot_png
    doc: |
      Correlation plots between QC metrics and cells PCA components.
      PNG format

  umap_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_qc_mtrcs_plot_png
    doc: |
      QC metrics on cells UMAP.
      PNG format

  umap_spl_mito_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_mito_plot_png
    doc: |
      Split by the percentage of transcripts mapped to mitochondrial genes cells UMAP.
      PNG format

  umap_spl_umi_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_umi_plot_png
    doc: |
      Split by the UMI per cell counts cells UMAP.
      PNG format

  umap_spl_gene_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_gene_plot_png
    doc: |
      Split by the genes per cell counts cells UMAP.
      PNG format

  umap_gr_cnd_spl_ph_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_ph_plot_png
    doc: |
      Grouped by condition split by cell cycle cells UMAP.
      PNG format

  umap_gr_cnd_spl_mito_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_mito_plot_png
    doc: |
      Grouped by condition split by the percentage of transcripts mapped to mitochondrial genes cells UMAP.
      PNG format

  umap_gr_cnd_spl_umi_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_umi_plot_png
    doc: |
      Grouped by condition split by the UMI per cell counts cells UMAP.
      PNG format

  umap_gr_cnd_spl_gene_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_gene_plot_png
    doc: |
      Grouped by condition split by the genes per cell counts cells UMAP.
      PNG format

  sc_rna_reduce_stdout_log:
    type: File
    outputSource: sc_rna_reduce/stdout_log
    doc: |
      stdout log generated by sc_rna_reduce step

  sc_rna_reduce_stderr_log:
    type: File
    outputSource: sc_rna_reduce/stderr_log
    doc: |
      stderr log generated by sc_rna_reduce step

  umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/umap_res_plot_png
    doc: |
      Clustered cells UMAP.
      PNG format

  slh_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/slh_res_plot_png
    doc: |
      Silhouette scores. Downsampled to max 500 cells per cluster.
      PNG format

  umap_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/umap_spl_idnt_res_plot_png
    doc: |
      Split by dataset clustered cells UMAP.
      PNG format

  cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_clst_spl_idnt_res_plot_png
    doc: |
      Grouped by cluster split by dataset cells composition plot. Downsampled.
      PNG format

  cmp_gr_idnt_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_idnt_spl_clst_res_plot_png
    doc: |
      Grouped by dataset split by cluster cells composition plot. Downsampled.
      PNG format

  umap_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/umap_spl_cnd_res_plot_png
    doc: |
      Split by grouping condition clustered cells UMAP.
      PNG format

  cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_clst_spl_cnd_res_plot_png
    doc: |
      Grouped by cluster split by condition cells composition plot. Downsampled.
      PNG format

  cmp_gr_cnd_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_cnd_spl_clst_res_plot_png
    doc: |
      Grouped by condition split by cluster cells composition plot. Downsampled.
      PNG format

  umap_spl_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/umap_spl_ph_res_plot_png
    doc: |
      Split by cell cycle phase clustered cells UMAP.
      PNG format

  cmp_gr_ph_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_ph_spl_idnt_res_plot_png
    doc: |
      Grouped by cell cycle phase split by dataset cells composition plot. Downsampled.
      PNG format

  cmp_gr_ph_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/cmp_gr_ph_spl_clst_res_plot_png
    doc: |
      Grouped by cell cycle phase split by cluster cells composition plot. Downsampled.
      PNG format

  xpr_avg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/xpr_avg_res_plot_png
    doc: |
      Log normalized scaled average gene expression per cluster.
      PNG format

  xpr_per_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/xpr_per_cell_res_plot_png
    doc: |
      Log normalized gene expression on cells UMAP.
      PNG format

  xpr_dnst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_rna_cluster/xpr_dnst_res_plot_png
    doc: |
      Log normalized gene expression density per cluster.
      PNG format

  gene_markers_tsv:
    type: File?
    outputSource: sc_rna_cluster/gene_markers_tsv
    doc: |
      Differentially expressed genes between each pair of clusters for all resolutions.
      TSV format

  ucsc_cb_html_data:
    type: Directory?
    outputSource: sc_rna_cluster/ucsc_cb_html_data
    doc: |
      Directory with UCSC Cellbrowser html data.

  seurat_data_rds:
    type: File
    outputSource: sc_rna_cluster/seurat_data_rds
    doc: |
      Processed Seurat data in RDS format

  sc_rna_cluster_stdout_log:
    type: File
    outputSource: sc_rna_cluster/stdout_log
    doc: |
      stdout log generated by sc_rna_cluster step

  sc_rna_cluster_stderr_log:
    type: File
    outputSource: sc_rna_cluster/stderr_log
    doc: |
      stderr log generated by sc_rna_cluster step


steps:

  uncompress_feature_bc_matrices:
    doc: |
      Extracts the content of TAR file into a folder
    run: ../tools/tar-extract.cwl
    in:
      file_to_extract: feature_bc_matrices_folder
    out:
    - extracted_folder

  sc_rna_filter:
    doc: |
      Filters single-cell RNA-Seq datasets based on the common QC metrics
    run: ../tools/sc-rna-filter.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/extracted_folder
      aggregation_metadata: aggregation_metadata
      grouping_data: grouping_data
      barcodes_data: barcodes_data
      rna_minimum_cells: rna_minimum_cells
      minimum_genes: minimum_genes
      maximum_genes: maximum_genes
      rna_minimum_umi: rna_minimum_umi
      minimum_novelty_score: minimum_novelty_score
      mito_pattern: mito_pattern
      maximum_mito_perc: maximum_mito_perc
      output_prefix:
        default: "step_1"
      parallel_memory_limit: parallel_memory_limit
      vector_memory_limit: vector_memory_limit
      threads: threads
    out:
    - raw_1_2_qc_mtrcs_pca_plot_png
    - raw_2_3_qc_mtrcs_pca_plot_png
    - raw_cells_count_plot_png
    - raw_umi_dnst_plot_png
    - raw_gene_dnst_plot_png
    - raw_gene_umi_corr_plot_png
    - raw_mito_dnst_plot_png
    - raw_nvlt_dnst_plot_png
    - raw_qc_mtrcs_dnst_plot_png
    - raw_umi_dnst_spl_cnd_plot_png
    - raw_gene_dnst_spl_cnd_plot_png
    - raw_mito_dnst_spl_cnd_plot_png
    - raw_nvlt_dnst_spl_cnd_plot_png
    - fltr_1_2_qc_mtrcs_pca_plot_png
    - fltr_2_3_qc_mtrcs_pca_plot_png
    - fltr_cells_count_plot_png
    - fltr_umi_dnst_plot_png
    - fltr_gene_dnst_plot_png
    - fltr_gene_umi_corr_plot_png
    - fltr_mito_dnst_plot_png
    - fltr_nvlt_dnst_plot_png
    - fltr_qc_mtrcs_dnst_plot_png
    - fltr_umi_dnst_spl_cnd_plot_png
    - fltr_gene_dnst_spl_cnd_plot_png
    - fltr_mito_dnst_spl_cnd_plot_png
    - fltr_nvlt_dnst_spl_cnd_plot_png
    - seurat_data_rds
    - stdout_log
    - stderr_log

  sc_rna_reduce:
    doc: |
      Integrates multiple single-cell RNA-Seq datasets,
      reduces dimensionality using PCA
    run: ../tools/sc-rna-reduce.cwl
    in:
      query_data_rds: sc_rna_filter/seurat_data_rds
      cell_cycle_data: cell_cycle_data
      normalization_method:
        default: "sctglm"
      highly_var_genes_count: highly_var_genes_count
      regress_mito_perc: regress_mito_perc
      regress_rna_umi: regress_rna_umi
      regress_genes: regress_genes
      regress_cellcycle: regress_cellcycle
      dimensions: dimensions
      low_memory:
        default: true
      output_prefix:
        default: "step_2"
      parallel_memory_limit: parallel_memory_limit
      vector_memory_limit: vector_memory_limit
      threads: threads
    out:
    - elbow_plot_png
    - qc_dim_corr_plot_png
    - umap_qc_mtrcs_plot_png
    - umap_spl_mito_plot_png
    - umap_spl_umi_plot_png
    - umap_spl_gene_plot_png
    - umap_gr_cnd_spl_ph_plot_png
    - umap_gr_cnd_spl_mito_plot_png
    - umap_gr_cnd_spl_umi_plot_png
    - umap_gr_cnd_spl_gene_plot_png
    - seurat_data_rds
    - stdout_log
    - stderr_log

  sc_rna_cluster:
    doc: |
      Clusters single-cell RNA-Seq datasets, identifies gene markers
    run: ../tools/sc-rna-cluster.cwl
    in:
      query_data_rds: sc_rna_reduce/seurat_data_rds
      dimensions: dimensions
      resolution: resolution
      genes_of_interest: genes_of_interest
      identify_diff_genes:
        default: true
      minimum_logfc: minimum_logfc
      minimum_pct: minimum_pct
      only_positive_diff_genes:
        default: true
      export_ucsc_cb:
        default: true
      output_prefix:
        default: "step_3"
      parallel_memory_limit: parallel_memory_limit
      vector_memory_limit: vector_memory_limit
      threads: threads
    out:
    - umap_res_plot_png
    - slh_res_plot_png
    - umap_spl_idnt_res_plot_png
    - cmp_gr_clst_spl_idnt_res_plot_png
    - cmp_gr_idnt_spl_clst_res_plot_png
    - umap_spl_cnd_res_plot_png
    - cmp_gr_clst_spl_cnd_res_plot_png
    - cmp_gr_cnd_spl_clst_res_plot_png
    - umap_spl_ph_res_plot_png
    - cmp_gr_ph_spl_idnt_res_plot_png
    - cmp_gr_ph_spl_clst_res_plot_png
    - xpr_avg_res_plot_png
    - xpr_per_cell_res_plot_png
    - xpr_dnst_res_plot_png
    - gene_markers_tsv
    - ucsc_cb_html_data
    - seurat_data_rds
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell RNA-Seq Analyze"
s:name: "Single-cell RNA-Seq Analyze"
s:alternateName: |
  Runs filtering, normalization, scaling, integration (optionally) and
  clustering for a single or aggregated single-cell RNA-Seq datasets

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/scRNA-Seq-Analysis/main/workflows/sc-rna-analyze.cwl
s:codeRepository: https://github.com/Barski-lab/scRNA-Seq-Analysis
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Single-cell RNA-Seq Analyze

  Runs filtering, normalization, scaling, integration (optionally) and
  clustering for a single or aggregated single-cell RNA-Seq datasets.