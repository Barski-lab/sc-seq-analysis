cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      shopt -s nocaseglob
      set -- "$0" "$@"

      function extract {
        FILE=$1
        COMBINED=$2
        T=`file -b "${FILE}" | awk '{print $1}'`
        case "${T}" in
          "bzip2"|"gzip"|"Zip")
            7z e -so "${FILE}" >> "${COMBINED}"
            ;;
          "ASCII")
            cat "${FILE}" >> "${COMBINED}" || true
            ;;
          *)
            echo "Error: file type unknown"
            rm -f "${COMBINED}"
            exit 1
        esac
      } 
      
      COMBINED=""
      for FILE in "$@"; do
          BASENAME=$(basename "$FILE")
          ROOTNAME="${BASENAME%.*}"
          EXT_LIST=( ".fastq" ".fq" )
          for ITEM in $EXT_LIST; do
            if [[ $ROOTNAME == *$ITEM ]]; then
              ROOTNAME="${ROOTNAME%.*}"
            fi
          done
          if [[ $COMBINED == "" ]]; then
            COMBINED="${ROOTNAME}"
          else
            COMBINED="${COMBINED}"_"${ROOTNAME}"
          fi
      done
      COMBINED="${COMBINED}".fastq

      for FILE in "$@"; do
          echo "Extracting:" $FILE;
          extract "${FILE}" "${COMBINED}"
      done;

    inputBinding:
      position: 5
    doc: |
      Bash script to extract compressed FASTQ file

  compressed_file:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      position: 6
    doc: |
      Compressed or uncompressed FASTQ file(s)


outputs:

  fastq_file:
    type: File
    outputBinding:
      glob: "*"


baseCommand: [bash, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Extract FASTQ"
s:name: "Extract FASTQ"
s:alternateName: "Extracts bzip2-/gzip-/zip-compressed FASTQ file"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/sc-seq-analysis/main/tools/extract-fastq.cwl
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
  Extract FASTQ

  Extracts bzip2-/gzip-/zip-compressed FASTQ file.
  If several FASTQ files are provided, they will be concatenated
  in the order that corresponds to files in input.
  Bash script's logic:
  - disable case sensitive glob check
  - check if root name of input file already include '.fastq' or
    '.fq' extension.
    If yes, set DEFAULT_EXT to "", otherwise use '.fastq'
  - check file type, decompress if needed
  - return 1, if file type is not recognized
  This script also works of input file doesn't have any extension
  at all.