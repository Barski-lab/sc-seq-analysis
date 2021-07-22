cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastqc:v0.11.5


inputs:
  reads_file:
    type:
      - File
    inputBinding:
      position: 50
    doc: |
      Input bam,sam,bam_mapped,sam_mapped or fastq file

  format_enum:
    type:
      - "null"
      - type: enum
        name: "format"
        symbols: ['bam','sam','bam_mapped','sam_mapped','fastq']
    inputBinding:
      position: 6
      prefix: '--format'
    doc: |
      Bypasses the normal sequence file format detection and
      forces the program to use the specified format.  Valid
      formats are bam,sam,bam_mapped,sam_mapped and fastq

  threads:
    type:
      - "null"
      - int
    inputBinding:
      position: 7
      prefix: '--threads'
    doc: |
      Specifies the number of files which can be processed
      simultaneously.  Each thread will be allocated 250MB of
      memory so you shouldn't run more threads than your
      available memory will cope with, and not more than
      6 threads on a 32 bit machine

  contaminants:
    type:
      - "null"
      - File
    inputBinding:
      position: 8
      prefix: '--contaminants'
    doc: |
      Specifies a non-default file which contains the list of
      contaminants to screen overrepresented sequences against.
      The file must contain sets of named contaminants in the
      form name[tab]sequence.  Lines prefixed with a hash will
      be ignored.

  adapters:
    type:
      - "null"
      - File
    inputBinding:
      position: 9
      prefix: '--adapters'
    doc: |
      Specifies a non-default file which contains the list of
      adapter sequences which will be explicity searched against
      the library. The file must contain sets of named adapters
      in the form name[tab]sequence.  Lines prefixed with a hash
      will be ignored.

  limits:
    type:
      - "null"
      - File
    inputBinding:
      position: 10
      prefix: '--limits'
    doc: |
      Specifies a non-default file which contains a set of criteria
      which will be used to determine the warn/error limits for the
      various modules.  This file can also be used to selectively
      remove some modules from the output all together.  The format
      needs to mirror the default limits.txt file found in the
      Configuration folder.

  kmers:
    type:
      - "null"
      - int
    inputBinding:
      position: 11
      prefix: '--kmers'
    doc: |
      Specifies the length of Kmer to look for in the Kmer content
      module. Specified Kmer length must be between 2 and 10. Default
      length is 7 if not specified.

  casava:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 13
      prefix: '--casava'
    doc: |
      Files come from raw casava output. Files in the same sample
      group (differing only by the group number) will be analysed
      as a set rather than individually. Sequences with the filter
      flag set in the header will be excluded from the analysis.
      Files must have the same names given to them by casava
      (including being gzipped and ending with .gz) otherwise they
      won't be grouped together correctly.

  nofilter:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 14
      prefix: '--nofilter'
    doc: |
      If running with --casava then don't remove read flagged by
      casava as poor quality when performing the QC analysis.

  hide_group:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 15
      prefix: '--nogroup'
    doc: |
      Disable grouping of bases for reads >50bp. All reports will
      show data for every base in the read.  WARNING: Using this
      option will cause fastqc to crash and burn if you use it on
      really long reads, and your plots may end up a ridiculous size.
      You have been warned!


outputs:
  zipped_file:
    type:
      - File
    outputBinding:
      glob: '*.zip'
  html_file:
    type:
      - File
    outputBinding:
      glob: '*.html'
  summary_file:
    type:
      - File
    outputBinding:
      glob: |
        ${
          return "*/summary.txt";
        }


baseCommand: [fastqc, --extract, --outdir, .]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "fastqc"
s:name: "fastqc"
s:alternateName: "fastqc"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/scRNA-Seq-Analysis/master/tools/fastqc.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
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
  Tool runs FastQC from Babraham Bioinformatics

s:about: |
  SYNOPSIS

      fastqc seqfile1 seqfile2 .. seqfileN

      fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
             [-c contaminant file] seqfile1 .. seqfileN

  DESCRIPTION

      FastQC reads a set of sequence files and produces from each one a quality
      control report consisting of a number of different modules, each one of
      which will help to identify a different potential type of problem in your
      data.

      If no files to process are specified on the command line then the program
      will start as an interactive graphical application.  If files are provided
      on the command line then the program will run with no user interaction
      required.  In this mode it is suitable for inclusion into a standardised
      analysis pipeline.

      The options for the program as as follows:

      -h --help       Print this help file and exit

      -v --version    Print the version of the program and exit

      -o --outdir     Create all output files in the specified output directory.
                      Please note that this directory must exist as the program
                      will not create it.  If this option is not set then the
                      output file for each sequence file is created in the same
                      directory as the sequence file which was processed.

      --casava        Files come from raw casava output. Files in the same sample
                      group (differing only by the group number) will be analysed
                      as a set rather than individually. Sequences with the filter
                      flag set in the header will be excluded from the analysis.
                      Files must have the same names given to them by casava
                      (including being gzipped and ending with .gz) otherwise they
                      won't be grouped together correctly.

      --nano          Files come from naopore sequences and are in fast5 format. In
                      this mode you can pass in directories to process and the program
                      will take in all fast5 files within those directories and produce
                      a single output file from the sequences found in all files.

      --nofilter      If running with --casava then don't remove read flagged by
                      casava as poor quality when performing the QC analysis.

      --extract       If set then the zipped output file will be uncompressed in
                      the same directory after it has been created.  By default
                      this option will be set if fastqc is run in non-interactive
                      mode.

      -j --java       Provides the full path to the java binary you want to use to
                      launch fastqc. If not supplied then java is assumed to be in
                      your path.

      --noextract     Do not uncompress the output file after creating it.  You
                      should set this option if you do not wish to uncompress
                      the output when running in non-interactive mode.

      --nogroup       Disable grouping of bases for reads >50bp. All reports will
                      show data for every base in the read.  WARNING: Using this
                      option will cause fastqc to crash and burn if you use it on
                      really long reads, and your plots may end up a ridiculous size.
                      You have been warned!

      -f --format     Bypasses the normal sequence file format detection and
                      forces the program to use the specified format.  Valid
                      formats are bam,sam,bam_mapped,sam_mapped and fastq

      -t --threads    Specifies the number of files which can be processed
                      simultaneously.  Each thread will be allocated 250MB of
                      memory so you shouldn't run more threads than your
                      available memory will cope with, and not more than
                      6 threads on a 32 bit machine

      -c              Specifies a non-default file which contains the list of
      --contaminants  contaminants to screen overrepresented sequences against.
                      The file must contain sets of named contaminants in the
                      form name[tab]sequence.  Lines prefixed with a hash will
                      be ignored.

      -a              Specifies a non-default file which contains the list of
      --adapters      adapter sequences which will be explicity searched against
                      the library. The file must contain sets of named adapters
                      in the form name[tab]sequence.  Lines prefixed with a hash
                      will be ignored.

      -l              Specifies a non-default file which contains a set of criteria
      --limits        which will be used to determine the warn/error limits for the
                      various modules.  This file can also be used to selectively
                      remove some modules from the output all together.  The format
                      needs to mirror the default limits.txt file found in the
                      Configuration folder.

     -k --kmers       Specifies the length of Kmer to look for in the Kmer content
                      module. Specified Kmer length must be between 2 and 10. Default
                      length is 7 if not specified.

     -q --quiet       Supress all progress messages on stdout and only report errors.

     -d --dir         Selects a directory to be used for temporary files written when
                      generating report images. Defaults to system temp directory if
                      not specified.
