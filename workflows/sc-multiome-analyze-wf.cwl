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
      Path to the compressed folder with feature-barcode matrix from Cell Ranger ARC Count/Aggregate
      experiment in MEX format. The rows consist of all the genes and peaks concatenated
      together and the columns are restricted to those barcodes that are identified as cells.

  aggregation_metadata:
    type: File?
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities. If '--mex' points to
      the Cell Ranger ARC Aggregate outputs, the aggr.csv file can be used. If input is not
      provided, the default dummy_metadata.csv will be used instead.

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    doc: |
      Count and barcode information for every ATAC fragment observed in the experiment in TSV
      format. Tbi-index file is required.

  annotation_gtf_file:
    type: File
    doc: |
      Path to the genome annotation file in GTF format.

  grouping_data:
    type: File?
    doc: |
      Path to the TSV/CSV file to define datasets grouping.
      First column - 'library_id' with the values and order
      that correspond to the 'library_id' column from the '
      --identity' file, second column 'condition'.
      Default: each dataset is assigned to its own group.

  blacklist_regions_file:
    type: File?
    doc: |
      Path to the optional BED file with the genomic blacklist regions.

  barcodes_data:
    type: File?
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata be selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present.
      Default: all cells used, no extra metadata is added

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
      Include cells where at least this many UMI (RNA transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file.
      Default: 500 (applied to all datasets)

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

  minimum_novelty_score:
    type:
    - "null"
    - float
    - float[]
    doc: |
      Include cells with the novelty score not lower than this value, calculated for
      as log10(genes)/log10(UMI) for RNA assay. If multiple values provided, each of them will
      be applied to the correspondent dataset from the '--mex' input based on the
      '--identity' file.
      Default: 0.8 (applied to all datasets)

  atac_minimum_cells:
    type: int?
    doc: |
      Include only peaks detected in at least this many cells.
      Default: 5 (applied to all datasets)

  atac_minimum_umi:
    type:
    - "null"
    - int
    - int[]
    doc: |
      Include cells where at least this many UMI (ATAC transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file.
      Default: 1000 (applied to all datasets)

  maximum_nucl_signal:
    type:
    - "null"
    - float
    - float[]
    doc: |
      Include cells with the nucleosome signal not bigger than this value.
      Nucleosome signal quantifies the approximate ratio of mononucleosomal
      to nucleosome-free fragments. If multiple values provided, each of
      them will be applied to the correspondent dataset from the '--mex' input
      based on the '--identity' file.
      Default: 4 (applied to all datasets)

  minimum_tss_enrich:
    type:
    - "null"
    - float
    - float[]
    doc: |
      Include cells with the TSS enrichment score not lower than this value.
      Score is calculated based on the ratio of fragments centered at the TSS
      to fragments in TSS-flanking regions. If multiple values provided, each
      of them will be applied to the correspondent dataset from the '--mex' input
      based on the '--identity' file.
      Default: 2 (applied to all datasets)

  minimum_frip:
    type:
    - "null"
    - float
    - float[]
    doc: |
      Include cells with the FRiP not lower than this value. If multiple values
      provided, each of them will be applied to the correspondent dataset from the
      '--mex' input based on the '--identity' file. FRiP is calculated for fragments.
      Default: 0.15 (applied to all datasets)

  maximum_blacklist_fraction:
    type:
    - "null"
    - float
    - float[]
    doc: |
      Include cells with the fraction of fragments in genomic blacklist regions
      not bigger than this value. If multiple values provided, each of them
      will be applied to the correspondent dataset from the '--mex' input based
      on the '--identity' file.
      Default: 0.05 (applied to all datasets)

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

  regress_cellcycle:
    type: boolean?
    doc: |
      Regress cell cycle scores as a confounding source of variation.
      Ignored if --cellcycle is not provided.
      Default: false

  rna_dimensions:
    type:
    - "null"
    - int
    - int[]
    doc: |
      Dimensionality from the 'pca' reduction to use when constructing weighted
      nearest-neighbor graph before clustering (from 1 to 50). If single value N
      is provided, use from 1 to N dimensions. If multiple values are provided,
      subset to only selected dimensions.
      Default: from 1 to 10

  minimum_var_peaks_perc:
    type: int?
    doc: |
      Minimum percentile for identifying the top most common peaks as highly variable.
      For example, setting to 5 will use the the top 95 percent most common among all cells
      peaks as highly variable. These peaks are used for datasets integration, scaling
      and dimensionality reduction.
      Default: 0 (use all available peaks)

  atac_dimensions:
    type:
    - "null"
    - int
    - int[]
    doc: |
      Dimensionality from the 'atac_lsi' reduction to use when constructing weighted
      nearest-neighbor graph before clustering (from 1 to 50). If single value N
      is provided, use from 2 to N dimensions. If multiple values are provided,
      subset to only selected dimensions.
      Default: from 2 to 10

  resolution:
    type:
    - "null"
    - float
    - float[]
    doc: |
      Clustering resolution applied to the constructed weighted nearest-neighbor
      graph. Can be set as an array.
      Default: 0.3, 0.5, 1.0

  genes_of_interest:
    type:
    - "null"
    - string
    - string[]
    doc: |
      Genes of interest to build gene expression and Tn5 insertion frequency plots
      for the nearest peaks. If '--fragments' is not provided only gene expression
      plots will be built.
      Default: None

  rna_minimum_logfc:
    type: float?
    doc: |
      For putative gene markers identification include only those genes that
      on average have log fold change difference in expression between every
      tested pair of clusters not lower than this value. Ignored if '--diffgenes'
      is not set.
      Default: 0.25

  rna_minimum_pct:
    type: float?
    doc: |
      For putative gene markers identification include only those genes that
      are detected in not lower than this fraction of cells in either of the
      two tested clusters. Ignored if '--diffgenes' is not set.
      Default: 0.1

  atac_minimum_logfc:
    type: float?
    doc: |
      For differentially accessible peaks identification include only those peaks that
      on average have log fold change difference in the chromatin accessibility between
      every tested pair of clusters not lower than this value. Ignored if '--diffpeaks'
      is not set.
      Default: 0.25

  atac_minimum_pct:
    type: float?
    doc: |
      For differentially accessible peaks identification include only those peaks that
      are detected in not lower than this fraction of cells in either of the two tested
      clusters. Ignored if '--diffpeaks' is not set.
      Default: 0.05

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
    outputSource: sc_multiome_filter/raw_1_2_qc_mtrcs_pca_plot_png
    doc: |
      PC1 and PC2 from the QC metrics PCA (not filtered).
      PNG format

  raw_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_2_3_qc_mtrcs_pca_plot_png
    doc: |
      PC2 and PC3 from the QC metrics PCA (not filtered).
      PNG format

  raw_cells_count_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_cells_count_plot_png
    doc: |
      Number of cells per dataset (not filtered).
      PNG format

  raw_rna_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_rna_umi_dnst_plot_png
    doc: |
      UMI per cell density for RNA assay (not filtered).
      PNG format

  raw_gene_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gene_dnst_plot_png
    doc: |
      Genes per cell density (not filtered).
      PNG format

  raw_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gene_umi_corr_plot_png
    doc: |
      Genes vs UMI per cell correlation for RNA assay (not filtered).
      PNG format

  raw_mito_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_mito_dnst_plot_png
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (not filtered).
      PNG format

  raw_nvlt_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_nvlt_dnst_plot_png
    doc: |
      Novelty score per cell density for RNA assay (not filtered).
      PNG format

  raw_atac_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_atac_umi_dnst_plot_png
    doc: |
      UMI per cell density for ATAC assay (not filtered).
      PNG format

  raw_peak_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_peak_dnst_plot_png
    doc: |
      Peaks per cell density (not filtered).
      PNG format

  raw_blck_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_blck_dnst_plot_png
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (not filtered).
      PNG format

  raw_rna_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_rna_atac_umi_corr_plot_png
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (not filtered).
      PNG format

  raw_tss_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_tss_atac_umi_corr_plot_png
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (not filtered).
      PNG format

  raw_qc_mtrcs_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_qc_mtrcs_dnst_plot_png
    doc: |
      QC metrics per cell density (not filtered).
      PNG format

  raw_tss_nrch_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_tss_nrch_plot_png
    doc: |
      TSS enrichment score (not filtered).
      PNG format

  raw_frgm_hist_png:
    type: File?
    outputSource: sc_multiome_filter/raw_frgm_hist_png
    doc: |
      Fragments length histogram (not filtered).
      PNG format

  raw_rna_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_rna_umi_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (not filtered).
      PNG format

  raw_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gene_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition genes per cell density (not filtered).
      PNG format

  raw_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_mito_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (not filtered).
      PNG format

  raw_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_nvlt_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (not filtered).
      PNG format

  raw_atac_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_atac_umi_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (not filtered).
      PNG format

  raw_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_peak_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition peaks per cell density (not filtered).
      PNG format

  raw_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_blck_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (not filtered).
      PNG format

  fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_1_2_qc_mtrcs_pca_plot_png
    doc: |
      PC1 and PC2 from the QC metrics PCA (filtered).
      PNG format

  fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_2_3_qc_mtrcs_pca_plot_png
    doc: |
      PC2 and PC3 from the QC metrics PCA (filtered).
      PNG format

  fltr_cells_count_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_cells_count_plot_png
    doc: |
      Number of cells per dataset (filtered).
      PNG format

  fltr_rna_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_rna_umi_dnst_plot_png
    doc: |
      UMI per cell density for RNA assay (filtered).
      PNG format

  fltr_gene_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gene_dnst_plot_png
    doc: |
      Genes per cell density (filtered).
      PNG format

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gene_umi_corr_plot_png
    doc: |
      Genes vs UMI per cell correlation for RNA assay (filtered).
      PNG format

  fltr_mito_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_mito_dnst_plot_png
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (filtered).
      PNG format

  fltr_nvlt_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_nvlt_dnst_plot_png
    doc: |
      Novelty score per cell density for RNA assay (filtered).
      PNG format

  fltr_atac_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_atac_umi_dnst_plot_png
    doc: |
      UMI per cell density for ATAC assay (filtered).
      PNG format

  fltr_peak_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_peak_dnst_plot_png
    doc: |
      Peaks per cell density (filtered).
      PNG format

  fltr_blck_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_blck_dnst_plot_png
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (filtered).
      PNG format

  fltr_rna_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_rna_atac_umi_corr_plot_png
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (filtered).
      PNG format

  fltr_tss_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_tss_atac_umi_corr_plot_png
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (filtered).
      PNG format

  fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_qc_mtrcs_dnst_plot_png
    doc: |
      QC metrics per cell density (filtered).
      PNG format

  fltr_tss_nrch_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_tss_nrch_plot_png
    doc: |
      TSS enrichment score (filtered).
      PNG format

  fltr_frgm_hist_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_frgm_hist_png
    doc: |
      Fragments length histogram (filtered).
      PNG format

  fltr_rna_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_rna_umi_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (filtered).
      PNG format

  fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gene_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition genes per cell density (filtered).
      PNG format

  fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_mito_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (filtered).
      PNG format

  fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_nvlt_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (filtered).
      PNG format

  fltr_atac_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_atac_umi_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (filtered).
      PNG format

  fltr_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_peak_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition peaks per cell density (filtered).
      PNG format

  fltr_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_blck_dnst_spl_cnd_plot_png
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (filtered).
      PNG format

  sc_multiome_filter_stdout_log:
    type: File
    outputSource: sc_multiome_filter/stdout_log
    doc: |
      stdout log generated by sc_multiome_filter step

  sc_multiome_filter_stderr_log:
    type: File
    outputSource: sc_multiome_filter/stderr_log
    doc: |
      stderr log generated by sc_multiome_filter step

  rna_elbow_plot_png:
    type: File?
    outputSource: sc_rna_reduce/elbow_plot_png
    doc: |
      Elbow plot from cells PCA of RNA assay.
      PNG format

  rna_qc_dim_corr_plot_png:
    type: File?
    outputSource: sc_rna_reduce/qc_dim_corr_plot_png
    doc: |
      Correlation plots between QC metrics and cells PCA components
      of RNA assay.
      PNG format

  rna_umap_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_qc_mtrcs_plot_png
    doc: |
      QC metrics on cells UMAP for RNA assay.
      PNG format

  rna_umap_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_plot_png
    doc: |
      Cells UMAP for RNA assay.
      PNG format

  rna_umap_spl_ph_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_ph_plot_png
    doc: |
      Split by cell cycle phase cells UMAP for RNA assay.
      PNG format

  rna_umap_spl_mito_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_mito_plot_png
    doc: |
      Split by the percentage of transcripts mapped to mitochondrial
      genes cells UMAP for RNA assay.
      PNG format

  rna_umap_spl_umi_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_umi_plot_png
    doc: |
      Split by the UMI per cell counts cells UMAP for RNA assay.
      PNG format

  rna_umap_spl_gene_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_gene_plot_png
    doc: |
      Split by the genes per cell counts cells UMAP for RNA assay.
      PNG format

  rna_umap_spl_idnt_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_idnt_plot_png
    doc: |
      Split by dataset cells UMAP for RNA assay.
      PNG format

  rna_umap_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_cnd_plot_png
    doc: |
      Split by grouping condition cells UMAP for RNA assay.
      PNG format

  rna_umap_gr_cnd_spl_ph_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_ph_plot_png
    doc: |
      Grouped by condition split by cell cycle cells UMAP for RNA assay.
      PNG format

  rna_umap_gr_cnd_spl_mito_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_mito_plot_png
    doc: |
      Grouped by condition split by the percentage of transcripts
      mapped to mitochondrial genes cells UMAP for RNA assay.
      PNG format

  rna_umap_gr_cnd_spl_umi_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_umi_plot_png
    doc: |
      Grouped by condition split by the UMI per cell counts cells UMAP
      for RNA assay.
      PNG format

  rna_umap_gr_cnd_spl_gene_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_gene_plot_png
    doc: |
      Grouped by condition split by the genes per cell counts cells UMAP
      for RNA assay.
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

  atac_qc_dim_corr_plot_png:
    type: File?
    outputSource: sc_atac_reduce/qc_dim_corr_plot_png
    doc: |
      Correlation plots between QC metrics and cells LSI dimensions
      of ATAC assay.
      PNG format

  atac_umap_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_qc_mtrcs_plot_png
    doc: |
      QC metrics on cells UMAP for ATAC assay.
      PNG format

  atac_umap_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_plot_png
    doc: |
      Cells UMAP for ATAC assay.
      PNG format

  atac_umap_spl_idnt_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_idnt_plot_png
    doc: |
      Split by dataset cells UMAP for ATAC assay.
      PNG format

  atac_umap_spl_cnd_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_cnd_plot_png
    doc: |
      Split by grouping condition cells UMAP for ATAC assay.
      PNG format

  sc_atac_reduce_stdout_log:
    type: File
    outputSource: sc_atac_reduce/stdout_log
    doc: |
      stdout log generated by sc_atac_reduce step

  sc_atac_reduce_stderr_log:
    type: File
    outputSource: sc_atac_reduce/stderr_log
    doc: |
      stderr log generated by sc_atac_reduce step

  wnn_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/umap_res_plot_png
    doc: |
      Clustered cells UMAP for integrated ATAC and RNA assays.
      PNG format

  wnn_umap_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/umap_spl_idnt_res_plot_png
    doc: |
      Split by dataset clustered cells UMAP for integrated
      ATAC and RNA assays.
      PNG format

  wnn_cmp_gr_clst_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cmp_gr_clst_spl_idnt_res_plot_png
    doc: |
      Grouped by cluster split by dataset cells composition plot for
      integrated ATAC and RNA assays. Downsampled.
      PNG format

  wnn_cmp_gr_idnt_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cmp_gr_idnt_spl_clst_res_plot_png
    doc: |
      Grouped by dataset split by cluster cells composition plot for
      integrated ATAC and RNA assays. Downsampled.
      PNG format

  wnn_umap_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/umap_spl_cnd_res_plot_png
    doc: |
      Split by grouping condition clustered cells UMAP for
      integrated ATAC and RNA assays.
      PNG format

  wnn_cmp_gr_clst_spl_cnd_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cmp_gr_clst_spl_cnd_res_plot_png
    doc: |
      Grouped by cluster split by condition cells composition plot for
      integrated ATAC and RNA assays. Downsampled.
      PNG format

  wnn_cmp_gr_cnd_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cmp_gr_cnd_spl_clst_res_plot_png
    doc: |
      Grouped by condition split by cluster cells composition plot for
      integrated ATAC and RNA assays. Downsampled.
      PNG format

  wnn_umap_spl_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/umap_spl_ph_res_plot_png
    doc: |
      Split by cell cycle phase clustered cells UMAP for
      integrated ATAC and RNA assays.
      PNG format

  wnn_cmp_gr_ph_spl_idnt_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cmp_gr_ph_spl_idnt_res_plot_png
    doc: |
      Grouped by cell cycle phase split by dataset cells composition plot
      for integrated ATAC and RNA assays. Downsampled.
      PNG format

  wnn_cmp_gr_ph_spl_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cmp_gr_ph_spl_clst_res_plot_png
    doc: |
      Grouped by cell cycle phase split by cluster cells composition plot
      for integrated ATAC and RNA assays. Downsampled.
      PNG format

  wnn_xpr_avg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/xpr_avg_res_plot_png
    doc: |
      Log normalized scaled average gene expression per cluster for
      integrated ATAC and RNA assays.
      PNG format

  wnn_xpr_per_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/xpr_per_cell_res_plot_png
    doc: |
      Log normalized gene expression on cells UMAP for integrated
      ATAC and RNA assays.
      PNG format

  wnn_xpr_per_cell_sgnl_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/xpr_per_cell_sgnl_res_plot_png
    doc: |
      Log normalized gene expression density on cells UMAP.
      PNG format

  wnn_xpr_dnst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/xpr_dnst_res_plot_png
    doc: |
      Log normalized gene expression density per cluster for
      integrated ATAC and RNA assays.
      PNG format

  wnn_cvrg_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_wnn_cluster/cvrg_res_plot_png
    doc: |
      Tn5 insertion frequency plot around gene for integrated
      ATAC and RNA assays.
      PNG format

  gene_markers_tsv:
    type: File?
    outputSource: sc_wnn_cluster/gene_markers_tsv
    doc: |
      Differentially expressed genes between each pair of clusters for
      all resolutions of clustered integrated ATAC and RNA assays.
      TSV format

  peak_markers_tsv:
    type: File?
    outputSource: sc_wnn_cluster/peak_markers_tsv
    doc: |
      Differentially accessible peaks between each pair of clusters for
      all resolutions of clustered integrated ATAC and RNA assays.
      TSV format

  ucsc_cb_html_data:
    type: Directory?
    outputSource: sc_wnn_cluster/ucsc_cb_html_data
    doc: |
      Directory with UCSC Cellbrowser html data.

  seurat_data_rds:
    type: File
    outputSource: sc_wnn_cluster/seurat_data_rds
    doc: |
      Processed Seurat data in RDS format

  seurat_data_h5ad:
    type: File?
    outputSource: sc_wnn_cluster/seurat_data_h5ad
    doc: |
      Reduced Seurat data in h5ad format

  sc_wnn_cluster_stdout_log:
    type: File
    outputSource: sc_wnn_cluster/stdout_log
    doc: |
      stdout log generated by sc_wnn_cluster step

  sc_wnn_cluster_stderr_log:
    type: File
    outputSource: sc_wnn_cluster/stderr_log
    doc: |
      stderr log generated by sc_wnn_cluster step


steps:

  uncompress_feature_bc_matrices:
    doc: |
      Extracts the content of TAR file into a folder
    run: ../tools/tar-extract.cwl
    in:
      file_to_extract: feature_bc_matrices_folder
    out:
    - extracted_folder

  sc_multiome_filter:
    doc: |
      Filters single-cell multiome ATAC and RNA-Seq datasets
      based on the common QC metrics
    run: ../tools/sc-multiome-filter.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/extracted_folder
      aggregation_metadata: aggregation_metadata
      atac_fragments_file: atac_fragments_file
      annotation_gtf_file: annotation_gtf_file
      grouping_data: grouping_data
      blacklist_regions_file: blacklist_regions_file
      barcodes_data: barcodes_data
      rna_minimum_cells: rna_minimum_cells
      minimum_genes: minimum_genes
      maximum_genes: maximum_genes
      rna_minimum_umi: rna_minimum_umi
      mito_pattern: mito_pattern
      maximum_mito_perc: maximum_mito_perc
      minimum_novelty_score: minimum_novelty_score
      atac_minimum_cells: atac_minimum_cells
      atac_minimum_umi: atac_minimum_umi
      maximum_nucl_signal: maximum_nucl_signal
      minimum_tss_enrich: minimum_tss_enrich
      minimum_frip: minimum_frip
      maximum_blacklist_fraction: maximum_blacklist_fraction
      output_prefix:
        default: "s_1"
      parallel_memory_limit: parallel_memory_limit
      vector_memory_limit: vector_memory_limit
      threads: threads
    out:
    - raw_1_2_qc_mtrcs_pca_plot_png
    - raw_2_3_qc_mtrcs_pca_plot_png
    - raw_cells_count_plot_png
    - raw_rna_umi_dnst_plot_png
    - raw_gene_dnst_plot_png
    - raw_gene_umi_corr_plot_png
    - raw_mito_dnst_plot_png
    - raw_nvlt_dnst_plot_png
    - raw_atac_umi_dnst_plot_png
    - raw_peak_dnst_plot_png
    - raw_blck_dnst_plot_png
    - raw_rna_atac_umi_corr_plot_png
    - raw_tss_atac_umi_corr_plot_png
    - raw_qc_mtrcs_dnst_plot_png
    - raw_tss_nrch_plot_png
    - raw_frgm_hist_png
    - raw_rna_umi_dnst_spl_cnd_plot_png
    - raw_gene_dnst_spl_cnd_plot_png
    - raw_mito_dnst_spl_cnd_plot_png
    - raw_nvlt_dnst_spl_cnd_plot_png
    - raw_atac_umi_dnst_spl_cnd_plot_png
    - raw_peak_dnst_spl_cnd_plot_png
    - raw_blck_dnst_spl_cnd_plot_png
    - fltr_1_2_qc_mtrcs_pca_plot_png
    - fltr_2_3_qc_mtrcs_pca_plot_png
    - fltr_cells_count_plot_png
    - fltr_rna_umi_dnst_plot_png
    - fltr_gene_dnst_plot_png
    - fltr_gene_umi_corr_plot_png
    - fltr_mito_dnst_plot_png
    - fltr_nvlt_dnst_plot_png
    - fltr_atac_umi_dnst_plot_png
    - fltr_peak_dnst_plot_png
    - fltr_blck_dnst_plot_png
    - fltr_rna_atac_umi_corr_plot_png
    - fltr_tss_atac_umi_corr_plot_png
    - fltr_qc_mtrcs_dnst_plot_png
    - fltr_tss_nrch_plot_png
    - fltr_frgm_hist_png
    - fltr_rna_umi_dnst_spl_cnd_plot_png
    - fltr_gene_dnst_spl_cnd_plot_png
    - fltr_mito_dnst_spl_cnd_plot_png
    - fltr_nvlt_dnst_spl_cnd_plot_png
    - fltr_atac_umi_dnst_spl_cnd_plot_png
    - fltr_peak_dnst_spl_cnd_plot_png
    - fltr_blck_dnst_spl_cnd_plot_png
    - seurat_data_rds
    - stdout_log
    - stderr_log

  sc_rna_reduce:
    doc: |
      Integrates multiple single-cell RNA-Seq datasets,
      reduces dimensionality using PCA
    run: ../tools/sc-rna-reduce.cwl
    in:
      query_data_rds: sc_multiome_filter/seurat_data_rds
      cell_cycle_data: cell_cycle_data
      normalization_method:
        default: "sctglm"
      integration_method:
        default: "seurat"
      highly_var_genes_count: highly_var_genes_count
      regress_mito_perc: regress_mito_perc
      regress_cellcycle: regress_cellcycle
      dimensions: rna_dimensions
      low_memory:
        default: true
      output_prefix:
        default: "s_2"
      parallel_memory_limit: parallel_memory_limit
      vector_memory_limit: vector_memory_limit
      threads: threads
    out:
    - elbow_plot_png
    - qc_dim_corr_plot_png
    - umap_qc_mtrcs_plot_png
    - umap_plot_png
    - umap_spl_ph_plot_png
    - umap_spl_mito_plot_png
    - umap_spl_umi_plot_png
    - umap_spl_gene_plot_png
    - umap_spl_idnt_plot_png
    - umap_spl_cnd_plot_png
    - umap_gr_cnd_spl_ph_plot_png
    - umap_gr_cnd_spl_mito_plot_png
    - umap_gr_cnd_spl_umi_plot_png
    - umap_gr_cnd_spl_gene_plot_png
    - seurat_data_rds
    - stdout_log
    - stderr_log

  sc_atac_reduce:
    doc: |
      Integrates multiple single-cell ATAC-Seq datasets,
      reduces dimensionality using LSI
    run: ../tools/sc-atac-reduce.cwl
    in:
      query_data_rds: sc_rna_reduce/seurat_data_rds
      normalization_method:
        default: "log-tfidf"
      integration_method:
        default: "signac"
      minimum_var_peaks_perc: minimum_var_peaks_perc
      dimensions: atac_dimensions
      output_prefix:
        default: "s_3"
      parallel_memory_limit: parallel_memory_limit
      vector_memory_limit: vector_memory_limit
      threads: threads
    out:
    - qc_dim_corr_plot_png
    - umap_qc_mtrcs_plot_png
    - umap_plot_png
    - umap_spl_idnt_plot_png
    - umap_spl_cnd_plot_png
    - seurat_data_rds
    - stdout_log
    - stderr_log

  sc_wnn_cluster:
    doc: |
      Clusters multiome ATAC and RNA-Seq datasets, identifies
      gene markers and differentially accessible peaks
    run: ../tools/sc-wnn-cluster.cwl
    in:
      query_data_rds: sc_atac_reduce/seurat_data_rds
      rna_dimensions: rna_dimensions
      atac_dimensions: atac_dimensions
      resolution: resolution
      atac_fragments_file: atac_fragments_file
      genes_of_interest: genes_of_interest
      identify_diff_genes:
        default: true
      identify_diff_peaks:
        default: true
      rna_minimum_logfc: rna_minimum_logfc
      rna_minimum_pct: rna_minimum_pct
      only_positive_diff_genes:
        default: true
      atac_minimum_logfc: atac_minimum_logfc
      atac_minimum_pct: atac_minimum_pct
      export_ucsc_cb:
        default: true
      export_h5ad_data:
        default: true
      output_prefix:
        default: "s_4"
      parallel_memory_limit: parallel_memory_limit
      vector_memory_limit: vector_memory_limit
      threads: threads
    out:
    - umap_res_plot_png
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
    - xpr_per_cell_sgnl_res_plot_png
    - xpr_dnst_res_plot_png
    - cvrg_res_plot_png
    - gene_markers_tsv
    - peak_markers_tsv
    - ucsc_cb_html_data
    - seurat_data_rds
    - seurat_data_h5ad
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell Multiome ATAC and RNA-Seq Analyze"
s:name: "Single-cell Multiome ATAC and RNA-Seq Analyze"
s:alternateName: |
  Runs filtering, normalization, scaling, integration (optionally) and
  clustering for a single or aggregated single-cell Multiome ATAC and
  RNA-Seq datasets

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/sc-seq-analysis/main/workflows/sc-multiome-analyze-wf.cwl
s:codeRepository: https://github.com/Barski-lab/sc-seq-analysis
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
  Single-cell Multiome ATAC and RNA-Seq Analyze

  Runs filtering, normalization, scaling, integration (optionally) and
  clustering for a single or aggregated single-cell Multiome ATAC and
  RNA-Seq datasets.