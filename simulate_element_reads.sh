set -ex

GRAPH_DIR=/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs
READ_DIR=/private/groups/patenlab/xhchang/reads/
mkdir -p "${READ_DIR}/sim/element/HG002"

# Prepare a separete GBWT
vg gbwt -o "${GRAPH_DIR}/hprc-chm-hg002-2024-03-25-mc.full.gbwt" -g "${GRAPH_DIR}/hprc-chm-hg002-2024-03-25-mc.full.gg" -Z "${GRAPH_DIR}/hprc-chm-hg002-2024-03-25-mc.full.gbz"
# Prepare an xg graph
vg convert --drop-haplotypes --xg-out "${GRAPH_DIR}//hprc-chm-hg002-2024-03-25-mc.full.gbz" >"${GRAPH_DIR}//hprc-chm-hg002-2024-03-25-mc.full.xg"

# Simulate reads (assuming you already put real reads in ${READ_DIR}/real/element/HG002/element.GAT-LI-C044.fq.gz)
# Error rates taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC11316309/
vg sim -r -n 2500000 -a -s 12345 -p 550 -v 100 -e 0.0003 -i 0.00000001 -x "${GRAPH_DIR}//hprc-chm-hg002-2024-03-25-mc.full.xg" -g "${GRAPH_DIR}//hprc-chm-hg002-2024-03-25-mc.full.gbwt" --sample-name HG002 -F "${READ_DIR}/real/element/HG002/element.GAT-LI-C044.fq.gz" --multi-position > "${READ_DIR}/sim/element/HG002/HG002-sim-element.gam"
# Subset reads
for READ_COUNT in 100 1000 10000 100000 1000000 ; do
    vg filter --interleaved -t1 --max-reads "${READ_COUNT}" "${READ_DIR}/sim/element/HG002/HG002-sim-element.gam" >"${READ_DIR}/sim/element/HG002/HG002-sim-element-${READ_COUNT}.gam"
done
