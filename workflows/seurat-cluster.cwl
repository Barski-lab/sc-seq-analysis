cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_features = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };
    - var split_numbers = function(line) {
          let splitted_line = line?line.split(/[\s,]+/).map(parseFloat):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


inputs:

  filtered_feature_bc_matrix_folder:
    type: File
    doc: "Compressed folder with aggregated filtered feature-barcode matrices in MEX format"

  aggregation_metadata:
    type: File
    doc: "Aggregation metadata in CSV format"

  minimum_cells:
    type: int?
    default: 5
    doc: "Include genes detected in at least this many cells. Applied to all datasets together"

  minimum_features:
    type: string?
    default: "250"
    doc: |
      Include cells where at least this many genes are detected.
      If multiple values provided each of them will be applied to
      the correspondent dataset.

  maximum_features:
    type: string?
    default: "5000"
    doc: |
      Include cells with the number of genes not bigger than this value.
      If multiple values provided each of them will be applied to the
      correspondent dataset.

  minimum_umis:
    type: string?
    default: "500"
    doc: |
      Include cells where at least this many UMIs are detected.
      If multiple values provided each of them will be applied
      to the correspondent dataset.

  minimum_novelty_score:
    type: string?
    default: "0.8"
    doc: |
      Include cells with the novelty score (the ratio of genes per cell over UMIs per cell)
      not lower than this value (calculated as log10(genes)/log10(UMIs)). If multiple values
      provided each of them will be applied to the correspondent dataset.

  maximum_mito_perc:
    type: float?
    default: 5
    doc: |
      Include cells with the percentage of transcripts mapped to mitochondrial genes
      not bigger than this value.

  mito_pattern:
    type: string?
    default: "^Mt-"
    doc: "Pattern to identify mitochondrial genes"

  high_var_features_count:
    type: int?
    default: 3000
    doc: "Number of highly variable genes to detect (used for dataset integration and dimensional reduction)"

  dimensionality:
    type: int?
    default: 10
    doc: "Number of principal components to use in UMAP projection and clustering (from 1 to 50)"

  umap_spread:
    type: float?
    default: 1
    doc: |
      The effective scale of embedded points on UMAP. In combination with mindist
      this determines how clustered/clumped the embedded points are.

  umap_mindist:
    type: float?
    default: 0.3
    doc: |
      Controls how tightly the embedding is allowed compress points together on UMAP.
      Larger values ensure embedded points are moreevenly distributed, while smaller
      values allow the algorithm to optimise more accurately with regard to local structure.
      Sensible values are in the range 0.001 to 0.5.

  umap_nneighbors:
    type: int?
    default: 30
    doc: |
      Determines the number of neighboring points used in UMAP. Larger values will result
      in more global structure being preserved at the loss of detailed local structure.
      In general this parameter should often be in the range 5 to 50.

  resolution:
    type: string?
    default: "0.1"
    doc: "Comma or space separated list of clustering resolutions"

  minimum_logfc:
    type: float?
    default: 0.25
    doc: |
      Include only those genes that on average have log fold change difference in
      expression between every tested pair of clusters not lower than this value.

  minimum_pct:
    type: float?
    default: 0.1
    doc: |
      Include only those genes that are detected in not lower than
      this fraction of cells in either of the two tested clusters.

  test_use:
    type:
    - "null"
    - type: enum
      symbols:
      - "wilcox"
      - "bimod"
      - "roc"
      - "t"
      - "negbinom"
      - "poisson"
      - "LR"
      - "MAST"
      - "DESeq2"
    default: "wilcox"
    doc: "Statistical test to use for gene markers identification"

  threads:
    type: int?
    default: 6
    doc: "Threads number to use"

  species:
    type:
    - "null"
    - type: enum
      symbols:
      - "hs"
      - "mm"
      - "none"
    default: "none"
    doc: "Species for gene name conversion when running cell type prediction"

  regress_cellcycle:
    type: boolean?
    default: false
    doc: "Regress cell cycle as a confounding source of variation"

  regress_mito_perc:
    type: boolean?
    default: false
    doc: "Regress mitochondrial gene expression as a confounding source of variation"

  only_positive_markers:
    type: boolean?
    default: false
    doc: "Report only positive gene markers"

  selected_features:
    type: string?
    default: null
    doc: "Comma or space separated list of genes of interest"

  conditions_data:
    type: File?
    doc: |
      Path to the TSV/CSV file to define datasets conditions
      for grouping. First column - 'library_id' with values
      from the --identity file, second column 'condition'.
      If not provided, each dataset is assigned to its own
      biological condition

  barcodes_data:
    type: File?
    doc: |
      Path to the headerless TSV/CSV file with selected barcodes
      (one per line) to prefilter input feature-barcode matrices.
      If not provided, use all cells

  cell_cycle_data:
    type: File?
    doc: |
      TSV/CSV file with cell cycle data. First column - 'phase', second column 'gene_id'.
      If not provided, skip cell cycle score assignment

  classifier_rds:
    type: File?
    doc: |
      Path to the Garnett classifier rds file for cell type prediction.
      If not provided, skip cell type prediction


outputs:

  raw_cell_count_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_cell_count_plot_png
    doc: |
      Number of cells per dataset (not filtered).
      PNG format

  raw_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_umi_dnst_spl_by_cond_plot_png
    doc: |
      Split by condition UMI density per cell (not filtered).
      PNG format

  raw_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_gene_dnst_spl_by_cond_plot_png
    doc: |
      Split by condition gene density per cell (not filtered).
      PNG format

  raw_gene_umi_corr_spl_by_ident_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_gene_umi_corr_spl_by_ident_plot_png
    doc: |
      Split by identity genes vs UMIs per cell correlation (not filtered).
      PNG format

  raw_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_mito_perc_dnst_spl_by_cond_plot_png
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PNG format

  raw_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_nvlt_score_dnst_spl_by_cond_plot_png
    doc: |
      Split by condition novelty score density per cell (not filtered).
      PNG format

  raw_qc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_qc_mtrcs_plot_png
    doc: |
      QC metrics densities per cell (not filtered).
      PNG format

  raw_qc_mtrcs_gr_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/raw_qc_mtrcs_gr_by_cond_plot_png
    doc: |
      Grouped by condition QC metrics densities per cell (not filtered).
      PNG format

  fltr_cell_count_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_cell_count_plot_png
    doc: |
      Number of cells per dataset (filtered).
      PNG format

  fltr_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_umi_dnst_spl_by_cond_plot_png
    doc: |
      Split by condition UMI density per cell (filtered).
      PNG format

  fltr_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_gene_dnst_spl_by_cond_plot_png
    doc: |
      Split by condition gene density per cell (filtered).
      PNG format

  fltr_gene_umi_corr_spl_by_ident_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_gene_umi_corr_spl_by_ident_plot_png
    doc: |
      Split by identity genes vs UMIs per cell correlation (filtered).
      PNG format

  fltr_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_mito_perc_dnst_spl_by_cond_plot_png
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered).
      PNG format

  fltr_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_nvlt_score_dnst_spl_by_cond_plot_png
    doc: |
      Split by condition novelty score density per cell (filtered).
      PNG format

  fltr_qc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_qc_mtrcs_plot_png
    doc: |
      QC metrics densities per cell (filtered).
      PNG format

  fltr_qc_mtrcs_gr_by_cond_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_qc_mtrcs_gr_by_cond_plot_png
    doc: |
      Grouped by condition QC metrics densities per cell (filtered).
      PNG format

  fltr_pca_spl_by_ph_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_pca_spl_by_ph_plot_png
    doc: |
      Split by cell cycle phase PCA of filtered unintegrated/scaled datasets.
      PNG format

  fltr_pca_spl_by_mito_perc_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_pca_spl_by_mito_perc_plot_png
    doc: |
      Split by level of mitochondrial gene expression PCA of filtered unintegrated/scaled datasets.
      PNG format

  fltr_umap_spl_by_idnt_plot_png:
    type: File?
    outputSource: seurat_cluster/fltr_umap_spl_by_idnt_plot_png
    doc: |
      Split by identity UMAP projected PCA of filtered unintegrated/scaled datasets.
      PNG format

  ntgr_elbow_plot_png:
    type: File?
    outputSource: seurat_cluster/ntgr_elbow_plot_png
    doc: |
      Elbow plot from PCA of filtered integrated/scaled datasets.
      PNG format

  ntgr_pca_plot_png:
    type: File?
    outputSource: seurat_cluster/ntgr_pca_plot_png
    doc: |
      PCA of filtered integrated/scaled datasets.
      PNG format

  ntgr_pca_heatmap_png:
    type: File?
    outputSource: seurat_cluster/ntgr_pca_heatmap_png
    doc: |
      Genes per cells expression heatmap sorted by their PC scores from PCA of filtered integrated/scaled datasets.
      PNG format

  ntgr_pca_loadings_plot_png:
    type: File?
    outputSource: seurat_cluster/ntgr_pca_loadings_plot_png
    doc: |
      PC scores of the most variant genes from PCA of filtered integrated/scaled datasets.
      PNG format

  ntgr_umap_spl_by_idnt_plot_png:
    type: File?
    outputSource: seurat_cluster/ntgr_umap_spl_by_idnt_plot_png
    doc: |
      Split by identity UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_res_plot_png
    doc: |
      Clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_umap_spl_by_cond_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_spl_by_cond_res_plot_png
    doc: |
      Split by condition clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_umap_ctype_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_ctype_res_plot_png
    doc: |
      Grouped by predicted cell types UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_umap_spl_by_ph_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_umap_spl_by_ph_res_plot_png
    doc: |
      Split by cell cycle phase clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  clst_qc_mtrcs_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/clst_qc_mtrcs_res_plot_png
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  expr_avg_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_avg_per_clst_res_plot_png
    doc: |
      Scaled average log normalized gene expression per cluster of filtered integrated/scaled datasets.
      PNG format

  expr_per_clst_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_per_clst_cell_res_plot_png
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets.
      PNG format

  expr_clst_heatmap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_clst_heatmap_res_plot_png
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets.
      PNG format

  expr_dnst_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_dnst_per_clst_res_plot_png
    doc: |
      Log normalized gene expression densities per cluster of filtered integrated/scaled datasets.
      PNG format

  expr_avg_per_ctype_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_avg_per_ctype_res_plot_png
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets.
      PNG format

  expr_per_ctype_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_per_ctype_cell_res_plot_png
    doc: |
      Log normalized gene expression per cell of clustered filtered/scaled integrated datasets with predicted cell types.
      PNG format

  expr_ctype_heatmap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_ctype_heatmap_res_plot_png
    doc: |
      Log normalized gene expression heatmap of clustered filtered integrated/scaled datasets with predicted cell types.
      PNG format

  expr_dnst_per_ctype_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_cluster/expr_dnst_per_ctype_res_plot_png
    doc: |
      Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets.
      PNG format

  seurat_clst_data_rds:
    type: File
    outputSource: seurat_cluster/seurat_clst_data_rds
    doc: |
      Clustered filtered integrated Seurat data.
      RDS format

  clst_pttv_gene_markers:
    type: File
    outputSource: seurat_cluster/clst_pttv_gene_markers
    doc: |
      Putative gene markers file for all clusters and all resolutions.
      TSV format

  clst_csrvd_gene_markers:
    type: File
    outputSource: seurat_cluster/clst_csrvd_gene_markers
    doc: |
      Conserved gene markers file for all clusters and all resolutions.
      TSV format

  compressed_cellbrowser_config_data:
    type: File
    outputSource: compress_cellbrowser_config_data/compressed_folder
    doc: |
      Compressed directory with UCSC Cellbrowser configuration data

  cellbrowser_html_data:
    type: Directory
    outputSource: seurat_cluster/cellbrowser_html_data
    doc: |
      Directory with UCSC Cellbrowser formatted html data

  cellbrowser_html_file:
    type: File
    outputSource: seurat_cluster/cellbrowser_html_file
    doc: |
      HTML index file from the directory with UCSC Cellbrowser formatted html data

  seurat_cluster_stdout_log:
    type: File
    outputSource: seurat_cluster/stdout_log
    doc: "Stdout log generated by Seurat"

  seurat_cluster_stderr_log:
    type: File
    outputSource: seurat_cluster/stderr_log
    doc: "Stderr log generated by Seurat"


steps:

  uncompress_feature_bc_matrices:
    in:
      compressed: filtered_feature_bc_matrix_folder
    out:
    - uncompressed
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        compressed:
          type: File
          inputBinding:
            position: 1
      outputs:
        uncompressed:
          type: Directory
          outputBinding:
            glob: "*"
      baseCommand: ["tar", "xzf"]

  seurat_cluster:
    run: ../tools/seurat-cluster.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/uncompressed
      aggregation_metadata: aggregation_metadata
      cell_cycle_data: cell_cycle_data
      conditions_data: conditions_data
      classifier_rds: classifier_rds
      species: species
      barcodes_data: barcodes_data
      minimum_cells: minimum_cells
      minimum_features:
        source: minimum_features
        valueFrom: $(split_numbers(self))
      maximum_features:
        source: maximum_features
        valueFrom: $(split_numbers(self))
      selected_features:
        source: selected_features
        valueFrom: $(split_features(self))
      minimum_umis:
        source: minimum_umis
        valueFrom: $(split_numbers(self))
      minimum_novelty_score:
        source: minimum_novelty_score
        valueFrom: $(split_numbers(self))
      maximum_mito_perc: maximum_mito_perc
      mito_pattern: mito_pattern
      regress_cellcycle: regress_cellcycle
      regress_mito_perc: regress_mito_perc
      high_var_features_count: high_var_features_count
      dimensionality: dimensionality
      umap_spread: umap_spread
      umap_mindist: umap_mindist
      umap_nneighbors: umap_nneighbors
      resolution:
        source: resolution
        valueFrom: $(split_numbers(self))
      minimum_logfc: minimum_logfc
      minimum_pct: minimum_pct
      only_positive_markers: only_positive_markers
      test_use: test_use
      export_pdf_plots:
        default: false
      export_rds_data:
        default: true
      threads: threads
    out:
    - raw_cell_count_plot_png
    - raw_umi_dnst_spl_by_cond_plot_png
    - raw_gene_dnst_spl_by_cond_plot_png
    - raw_gene_umi_corr_spl_by_ident_plot_png
    - raw_mito_perc_dnst_spl_by_cond_plot_png
    - raw_nvlt_score_dnst_spl_by_cond_plot_png
    - raw_qc_mtrcs_plot_png
    - raw_qc_mtrcs_gr_by_cond_plot_png
    - fltr_cell_count_plot_png
    - fltr_umi_dnst_spl_by_cond_plot_png
    - fltr_gene_dnst_spl_by_cond_plot_png
    - fltr_gene_umi_corr_spl_by_ident_plot_png
    - fltr_mito_perc_dnst_spl_by_cond_plot_png
    - fltr_nvlt_score_dnst_spl_by_cond_plot_png
    - fltr_qc_mtrcs_plot_png
    - fltr_qc_mtrcs_gr_by_cond_plot_png
    - fltr_pca_spl_by_ph_plot_png
    - fltr_pca_spl_by_mito_perc_plot_png
    - fltr_umap_spl_by_idnt_plot_png
    - ntgr_elbow_plot_png
    - ntgr_pca_plot_png
    - ntgr_pca_heatmap_png
    - ntgr_pca_loadings_plot_png
    - ntgr_umap_spl_by_idnt_plot_png
    - clst_umap_res_plot_png
    - clst_umap_spl_by_cond_res_plot_png
    - clst_umap_ctype_res_plot_png
    - clst_umap_spl_by_ph_res_plot_png
    - clst_qc_mtrcs_res_plot_png
    - clst_pttv_gene_markers   
    - clst_csrvd_gene_markers
    - expr_avg_per_clst_res_plot_png
    - expr_per_clst_cell_res_plot_png
    - expr_clst_heatmap_res_plot_png
    - expr_dnst_per_clst_res_plot_png
    - expr_avg_per_ctype_res_plot_png
    - expr_per_ctype_cell_res_plot_png
    - expr_ctype_heatmap_res_plot_png
    - expr_dnst_per_ctype_res_plot_png
    - seurat_clst_data_rds
    - cellbrowser_config_data
    - cellbrowser_html_data
    - cellbrowser_html_file
    - stdout_log
    - stderr_log

  compress_cellbrowser_config_data:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: seurat_cluster/cellbrowser_config_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Seurat Cluster"
s:name: "Seurat Cluster"
s:alternateName: "Runs filtering, integration, and clustering analyses for Cell Ranger Count Gene Expression or Cell Ranger Aggregate experiments"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/scRNA-Seq-Analysis/master/workflows/seurat-cluster.cwl
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
  Seurat Cluster
  ==============

  Runs filtering, integration, and clustering analyses for Cell Ranger
  Count Gene Expression or Cell Ranger Aggregate experiments.