#!/usr/bin/env bash
# simulate_illumina_reads.sh: Simulate reads to look like Illumina sequencing data.

set -ex

: "${GRAPH_DIR:=/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs}"
: "${READ_DIR:=/private/groups/patenlab/xhchang/reads/}"
: "${GRAPH_NAME:=hprc-chm-hg002-v2.0.full}"
: "${SAMPLE_NAME:=HG002}"
: "${SAMPLE_FASTQ:=${READ_DIR}/real/illumina/HG002/HG002.novaseq.pcr-free.40x.3m.fq}"
: "${VG:=vg}"

mkdir -p "${READ_DIR}/sim/illumina/${SAMPLE_NAME}"

if [[ ! -e "${GRAPH_DIR}/${GRAPH_NAME}.gbwt" || ! -e "${GRAPH_DIR}/${GRAPH_NAME}.gg" ]] ; then
    # Prepare a separete GBWT
    vg gbwt -o "${GRAPH_DIR}/${GRAPH_NAME}.gbwt.tmp" -g "${GRAPH_DIR}/${GRAPH_NAME}.gg.tmp" -Z "${GRAPH_DIR}/${GRAPH_NAME}.gbz"
    mv "${GRAPH_DIR}/${GRAPH_NAME}.gbwt.tmp" "${GRAPH_DIR}/${GRAPH_NAME}.gbwt"
    mv "${GRAPH_DIR}/${GRAPH_NAME}.gg.tmp" "${GRAPH_DIR}/${GRAPH_NAME}.gg"
fi
if [[ ! -e "${GRAPH_DIR}/${GRAPH_NAME}.xg" ]] ; then
    # Prepare an xg graph
    vg convert --drop-haplotypes --xg-out "${GRAPH_DIR}/${GRAPH_NAME}.gbz" >"${GRAPH_DIR}/${GRAPH_NAME}.xg.tmp"
    mv "${GRAPH_DIR}/${GRAPH_NAME}.xg.tmp" "${GRAPH_DIR}/${GRAPH_NAME}.xg"
fi

# Simulate reads
"${VG}" sim -r -n 2500000 -a -s 12345 -p 570 -v 165 -i 0.00029 -x "${GRAPH_DIR}/${GRAPH_NAME}.xg" -g "${GRAPH_DIR}/${GRAPH_NAME}.gbwt" --sample-name ${SAMPLE_NAME} -F "${SAMPLE_FASTQ}" --multi-position > "${READ_DIR}/sim/illumina/${SAMPLE_NAME}/${SAMPLE_NAME}-sim-illumina.gam"
# Subset reads
for READ_COUNT in 100 1000 10000 100000 1000000 ; do
    "${VG}" filter --interleaved -t1 --max-reads "${READ_COUNT}" "${READ_DIR}/sim/illumina/${SAMPLE_NAME}/${SAMPLE_NAME}-sim-illumina.gam" >"${READ_DIR}/sim/illumina/${SAMPLE_NAME}/${SAMPLE_NAME}-sim-illumina-${READ_COUNT}.gam"
done
