#!/bin/bash
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
INPUTS_DIR="${BASE_DIR}/data/rna_inputs"
TOOLS_DIR="${BASE_DIR}/../tools"
JOBS_DIR=$BASE_DIR

mkdir -p $INPUTS_DIR && cd $INPUTS_DIR
echo "Fetching data for tests into ${INPUTS_DIR} directory. Previously downloaded files will be overwritten."
wget -q --show-progress -O filtered_feature_bc_matrix.tar.gz https://figshare.com/ndownloader/files/34819513
tar xzf filtered_feature_bc_matrix.tar.gz && rm -f filtered_feature_bc_matrix.tar.gz
wget -q --show-progress -O aggregation.csv https://figshare.com/ndownloader/files/34819516
wget -q --show-progress -O condition.csv https://figshare.com/ndownloader/files/34819519
wget -q --show-progress -O mouse_cell_cycle_genes.csv https://figshare.com/ndownloader/files/34822054

WORKING_DIR="${BASE_DIR}/data/rna_outputs/sc_rna_filter"
mkdir -p $WORKING_DIR && cd $WORKING_DIR
echo "Running sc-rna-filter.cwl in ${WORKING_DIR} directory."
cwltool --debug "$TOOLS_DIR/sc-rna-filter.cwl" "$JOBS_DIR/sc-rna-filter.yaml"

WORKING_DIR="${BASE_DIR}/data/rna_outputs/sc_rna_reduce"
mkdir -p $WORKING_DIR && cd $WORKING_DIR
echo "Running sc-rna-reduce.cwl in ${WORKING_DIR} directory."
cwltool --debug "$TOOLS_DIR/sc-rna-reduce.cwl" "$JOBS_DIR/sc-rna-reduce.yaml"

WORKING_DIR="${BASE_DIR}/data/rna_outputs/sc_rna_cluster"
mkdir -p $WORKING_DIR && cd $WORKING_DIR
echo "Running sc-rna-cluster.cwl in ${WORKING_DIR} directory."
cwltool --debug "$TOOLS_DIR/sc-rna-cluster.cwl" "$JOBS_DIR/sc-rna-cluster.yaml"

WORKING_DIR="${BASE_DIR}/data/rna_outputs/sc_ctype_assign"
mkdir -p $WORKING_DIR && cd $WORKING_DIR
echo -e "cluster\tctype" > cell_types.tsv
for i in {0..18}
do
    echo -e "$i\tctype_$i" >> cell_types.tsv
done
echo "Running sc-ctype-assign.cwl in ${WORKING_DIR} directory."
cwltool --debug "$TOOLS_DIR/sc-ctype-assign.cwl" "$JOBS_DIR/sc-ctype-assign.yaml"

WORKING_DIR="${BASE_DIR}/data/rna_outputs/sc_rna_de_pseudobulk"
mkdir -p $WORKING_DIR && cd $WORKING_DIR
echo "Running sc-rna-de-pseudobulk.cwl in ${WORKING_DIR} directory."
cwltool --debug "$TOOLS_DIR/sc-rna-de-pseudobulk.cwl" "$JOBS_DIR/sc-rna-de-pseudobulk.yaml"

WORKING_DIR="${BASE_DIR}/data/rna_outputs/sc_rna_da_cells"
mkdir -p $WORKING_DIR && cd $WORKING_DIR
echo "Running sc-rna-da-cells.cwl in ${WORKING_DIR} directory."
cwltool --debug "$TOOLS_DIR/sc-rna-da-cells.cwl" "$JOBS_DIR/sc-rna-da-cells.yaml"