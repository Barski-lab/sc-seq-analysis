#!/bin/bash
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
INPUTS_DIR="${BASE_DIR}/data/multiome_inputs"
TOOLS_DIR="${BASE_DIR}/../tools"
JOBS_DIR=$BASE_DIR

# mkdir -p $INPUTS_DIR && cd $INPUTS_DIR
# echo "Fetching data for tests into ${INPUTS_DIR} directory. Previously downloaded files will be overwritten."
# wget -O pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz
# wget -O pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
# wget -O pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
# wget -O gencode.v40.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz
# wget -O hg38-blacklist.v2.bed.gz https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
# gzip -d hg38-blacklist.v2.bed.gz
# gzip -d gencode.v40.annotation.gtf.gz
# tar xzf pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz && rm -f pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz

WORKING_DIR="${BASE_DIR}/data/multiome_outputs/sc_multiome_filter"
mkdir -p $WORKING_DIR && cd $WORKING_DIR
echo "Running sc-multiome-filter.cwl in ${WORKING_DIR} directory."
cwltool --debug "$TOOLS_DIR/sc-multiome-filter.cwl" "$JOBS_DIR/sc-multiome-filter.yaml"

