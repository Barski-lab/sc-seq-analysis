cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_label = function(i) {
        var rootname = inputs.molecule_info_h5[i].basename.split('.').slice(0,-1).join('.');
        rootname = (rootname=="")?inputs.molecule_info_h5[i].basename:rootname;
        return inputs.gem_well_labels?inputs.gem_well_labels[i].replace(/,/g, "_"):rootname;
    };
- class: InitialWorkDirRequirement
  listing: |
    ${
      var entry = "library_id,molecule_h5\n"
      for (var i=0; i < inputs.molecule_info_h5.length; i++){
        entry += get_label(i) + "," + inputs.molecule_info_h5[i].path + "\n"
      }
      return [{
        "entry": entry,
        "entryname": "metadata.csv"
      }];
    }


hints:
- class: DockerRequirement
  dockerPull: cumulusprod/cellranger:4.0.0


inputs:
  
  molecule_info_h5:
    type: File[]
    doc: |
      Array of molecule-level information files in HDF5 format.
      Outputs from "cellranger count" command
  
  gem_well_labels:
    type:
    - "null"
    - string[]
    doc: |
      Array of GEM well identifiers to be used for labeling purposes only.
      If not provided use rootnames of files from the molecule_info_h5 input

  normalization_mode:
    type:
    - "null"
    - type: enum
      name: "normalization"
      symbols: ["none", "mapped"]
    inputBinding:
      position: 5
      prefix: "--normalize"
    doc: |
      Library depth normalization mode: mapped, none.
      Default: mapped

  threads:
    type: int?
    inputBinding:
      position: 6
      prefix: "--localcores"
    doc: |
      Set max cores the pipeline may request at one time.
      Default: all available

  memory_limit:
    type: int?
    inputBinding:
      position: 7
      prefix: "--localmem"
    doc: |
      Set max GB the pipeline may request at one time
      Default: all available

  virt_memory_limit:
    type: int?
    inputBinding:
      position: 8
      prefix: "--localvmem"
    doc: |
      Set max virtual address space in GB for the pipeline
      Default: all available


outputs:

  web_summary_report:
    type: File
    outputBinding:
      glob: "aggregated/outs/web_summary.html"
    doc: |
      Aggregated run summary metrics and charts in HTML format

  metrics_summary_report_json:
    type: File
    outputBinding:
      glob: "aggregated/outs/summary.json"
    doc: |
      Aggregated run summary metrics in JSON format
  
  secondary_analysis_report_folder:
    type: Directory
    outputBinding:
      glob: "aggregated/outs/analysis"
    doc: |
      Folder with secondary analysis results including dimensionality reduction,
      cell clustering, and differential expression for aggregated results

  filtered_feature_bc_matrix_folder:
    type: Directory
    outputBinding:
      glob: "aggregated/outs/filtered_feature_bc_matrix"
    doc: |
      Folder with aggregated filtered feature-barcode matrices containing only cellular barcodes in MEX format

  filtered_feature_bc_matrix_h5:
    type: File
    outputBinding:
      glob: "aggregated/outs/filtered_feature_bc_matrix.h5"
    doc: |
      Aggregated filtered feature-barcode matrices containing only cellular barcodes in HDF5 format

  raw_feature_bc_matrices_folder:
    type: Directory
    outputBinding:
      glob: "aggregated/outs/raw_feature_bc_matrix"
    doc: |
      Folder with aggregated unfiltered feature-barcode matrices containing all barcodes in MEX format

  raw_feature_bc_matrices_h5:
    type: File
    outputBinding:
      glob: "aggregated/outs/raw_feature_bc_matrix.h5"
    doc: |
      Aggregated unfiltered feature-barcode matrices containing all barcodes in HDF5 format

  aggregation_metadata:
    type: File
    outputBinding:
      glob: "aggregated/outs/aggregation.csv"
    doc: |
      Copy of the input aggregation CSV file

  loupe_browser_track:
    type: File
    outputBinding:
      glob: "aggregated/outs/cloupe.cloupe"
    doc: |
      Loupe Browser visualization and analysis file for aggregated results

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["cellranger", "aggr", "--disable-ui", "--id", "aggregated", "--csv", "metadata.csv"]


stdout: cellranger_aggr_stdout.log
stderr: cellranger_aggr_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger Aggregate"
s:name: "Cell Ranger Aggregate"
s:alternateName: "Aggregates outputs from multiple runs of Cell Ranger Count Gene Expression"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/sc-seq-analysis/main/tools/cellranger-aggr.cwl
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
  Cell Ranger Aggregate

  Aggregates outputs from multiple runs of Cell Ranger Count Gene Expression.

  Parameters set by default:
  --disable-ui - no need in any UI when running in Docker container
  --id - hardcoded to `aggregated` as we want to return the content
         of the outputs folder as separate outputs

  Skipped parameters:
  --nosecondary
  --dry
  --noexit
  --nopreflight
  --description
  --jobmode
  --mempercore
  --maxjobs
  --jobinterval
  --overrides
  --uiport

  Not supported features:
  - Batch correction caused by different versions of the Single Cell Gene
    Expression chemistry is not supported as the generated metadata file
    doesn't include "batch" field.