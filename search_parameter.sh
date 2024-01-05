#!/usr/bin/env bash
# search_parameter.sh: run a Giraffe grid search on one parameter and plot the results.

set -ex

# Here we use : and := to set variables to default values if not present in the environment.
# You can set these in the environment to override them and I don't have to write a CLI option parser.
# See https://stackoverflow.com/a/28085062

: "${WORK_DIR:=output/gridsearch}"
: "${GRAPH_BASE:=/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v1.1-mc-chm13.d9}"
: "${MINPARAMS:=k29.w11.W}"
: "${INPUT_READS:=/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads/sim/illumina/HG002/HG002-sim-illumina-100000.gam}"

: "${OPTION_NAME:=downsample-min}"
: "${OPTION_RANGE:=30:200:20}"
: "${PRESET:=sr}"


SEARCH_DIR="${WORK_DIR}/${OPTION_NAME}"

mkdir -p "${SEARCH_DIR}"

srun -c20 --time 01:00:00 --partition=medium --mem 300G --job-name giraffe-run vg giraffe -t16 --parameter-preset ${PRESET} --progress --track-provenance -Z ${GRAPH_BASE}.gbz -d ${GRAPH_BASE}.dist -m ${GRAPH_BASE}.${MINPARAMS}.withzip.min -z ${GRAPH_BASE}.${MINPARAMS}.zipcodes -G ${INPUT_READS} --${OPTION_NAME} ${OPTION_RANGE} --output-basename ${SEARCH_DIR}/search

SLURM_JOBS=()
for GAM_FILE in ${SEARCH_DIR}/*.gam ; do
    if [[ ! -e "${GAM_FILE%.gam}.compare.txt" ]] ; then
        SLURM_JOBS+=($(sbatch --parsable -c4 --time 01:00:00 --partition=medium --mem 25G --job-name giraffe-eval --wrap "vg annotate -a ${GAM_FILE} -x ${GRAPH_BASE}.gbz -m | vg gamcompare --range 200 - ${INPUT_READS} -T -a "${GAM_FILE}" > ${GAM_FILE%.gam}.compared.tsv 2>${GAM_FILE%.gam}.compare.txt"))
    fi
    if [[ ! -e "${GAM_FILE%.gam}.gamstats.txt" ]] ; then
        SLURM_JOBS+=($(sbatch --parsable -c2 --time 01:00:00 --partition=medium --mem 10G --job-name giraffe-stats --wrap "vg stats -a "${GAM_FILE}" >${GAM_FILE%.gam}.gamstats.txt"))
    fi
done

# If there are any jobs in the array, wait for them to finish
if [[ "${#SLURM_JOBS[@]}" != "0" ]] ; then
    # Wait for these particular jobs to finish 
    QUEUE_LINES=0
    while [[ "${QUEUE_LINES}" != "1" ]] ; do
        sleep 2
        QUEUE_LINES="$(squeue -u $USER -j $(IFS=,; echo "${SLURM_JOBS[*]}") | wc -l)"
    done
fi

COMPARISON_SCRATCH="${WORK_DIR}/combined.tsv"
printf "correct\tmq\taligner\tread\teligible\n" >"${COMPARISON_SCRATCH}"
cat "${SEARCH_DIR}"/*.compared.tsv | sed "s_${SEARCH_DIR}[-/a-zA-Z]*_v_g" | sed 's_.gam__g' | grep -v "^correct" >>"${COMPARISON_SCRATCH}"

# Plot PR and QQ
srun -c2 --mem 10G Rscript ./plot-pr.R "${COMPARISON_SCRATCH}" ${WORK_DIR}/${OPTION_NAME}.compared.png
srun -c2 --mem 10G Rscript ./plot-qq.R "${COMPARISON_SCRATCH}" ${WORK_DIR}/${OPTION_NAME}.qq.png

# Dump QC stats and speed
rm -f ${WORK_DIR}/${OPTION_NAME}.stats.txt
for GAM_FILE in ${SEARCH_DIR}/*.gam ; do
    echo "${GAM_FILE}" >>${WORK_DIR}/${OPTION_NAME}.stats.txt
    cat ${GAM_FILE%.gam}.compare.txt >>${WORK_DIR}/${OPTION_NAME}.stats.txt
    cat ${GAM_FILE%.gam}.gamstats.txt | grep "Speed:" >>${WORK_DIR}/${OPTION_NAME}.stats.txt
    # Add a blank line at the end
    echo "" >>${WORK_DIR}/${OPTION_NAME}.stats.txt
done

