#!/usr/bin/env bash
#SBATCH --array 1-14%7
#SBATCH -c8
#SBATCH --mem 128G
#SBATCH --partition long
#SBATCH --time 40:00:00
set -e
set -x
set -u
export SINGULARITY_DOCKER_HUB_MIRROR=http://10.50.100.88
CONDITION_NUMBER=0
EXPERIMENT_NAME="r10y2025-training"

for BAM_URL in \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG001/PAW79146/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG001/PAW81754/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG002/PAW71238/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG003/PAY87794/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG003/PAY87954/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG004/PAY87778/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG004/PAY88428/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG005/PAW87816/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG005/PAW88001/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG006/PAY77227/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG006/PBA16846/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG007/PAY12990/calls.sorted.bam \
https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG007/PBA20413/calls.sorted.bam \
; do
    ((CONDITION_NUMBER+=1))
    if [[ "${CONDITION_NUMBER}" != "${SLURM_ARRAY_TASK_ID}" ]] ; then
        # Not our condition
        continue
    fi
    
    BASE_NAME="$(echo "${BAM_URL}" | rev | cut -f2,3 -d'/' | rev | tr '/' '_')"
    SAMPLE_NAME="$(echo "${BAM_URL}" | rev | cut -f3 -d'/' | rev)"
    
    mkdir -p /private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/r10y2025
    INPUT_JSON="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/r10y2025/${BASE_NAME}.input.json"
    JOB_STORE="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/r10y2025/${BASE_NAME}.tree"
    BUS_PATH="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/r10y2025/${BASE_NAME}.bus"
    OUT_DIR="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/r10y2025/${BASE_NAME}.dir"
    OUT_FILE="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/r10y2025/${BASE_NAME}.output.json"
    LOG_FILE="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/r10y2025/${BASE_NAME}.log"
    LOG_DIR="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/r10y2025/${BASE_NAME}.batchlog"
    mkdir -p "${LOG_DIR}"
    
    cat >"${INPUT_JSON}" <<EOF
{
    "Giraffe.INPUT_BAM_FILE": "${BAM_URL}",
    "Giraffe.SAMPLE_NAME": "${SAMPLE_NAME}",
    "Giraffe.PAIRED_READS": false,
    "Giraffe.HAPLOTYPE_SAMPLING": true,
    "Giraffe.DIPLOID": false,
    "Giraffe.SET_REFERENCE": "GRCh38",
    "Giraffe.HAPLOTYPE_NUMBER": 16,
    "Giraffe.VG_DOCKER": "quay.io/vgteam/vg:ci-1977-791547fbb3bf130aaf5b6c3a61fbda9a63b4623b",
    "Giraffe.GIRAFFE_PRESET": "r10",
    "Giraffe.READS_PER_CHUNK": 150000,
    "Giraffe.GBZ_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v2.0-mc-grch38-eval.ec1M.gbz",
    "Giraffe.HAPL_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v2.0-mc-grch38-eval.ec1M.fragment.hapl",
    "Giraffe.REFERENCE_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    "Giraffe.REFERENCE_INDEX_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
    "Giraffe.REFERENCE_DICT_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict",
    "Giraffe.PATH_LIST_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/grch38-pansn.fa.paths.callable.txt",
    "Giraffe.REFERENCE_PREFIX": "GRCh38#0#",
    "Giraffe.REALIGN_INDELS": false,
    "Giraffe.LEFTALIGN_BAM": true,
    "Giraffe.INDEX_MINIMIZER_MEM": 600,
    "Giraffe.SPLIT_READ_MEM": 120 
}
EOF

    . /private/home/anovak1/workspace/toil/venv/bin/activate 

    if [[ -e "${JOB_STORE}" ]] ; then
        toil-wdl-runner --restart --jobStore "${JOB_STORE}" --wdlOutputDialect miniwdl --wdlOutputDirectory "${OUT_DIR}" --wdlOutputFile "${OUT_FILE}" --clean=onSuccess --batchSystem slurm --slurmTime 11:59:59 --statePollingWait 30 --disableProgress --caching=False --retryCount 0 --stopOnFirstFailure=True --logFile "${LOG_FILE}" --batchLogsDir "${LOG_DIR}" --writeMessages "${BUS_PATH}" --runImportsOnWorkers '#workflow/github.com/vgteam/vg_wdl/Giraffe:lr-giraffe' "${INPUT_JSON}"
    else
        toil-wdl-runner --jobStore "${JOB_STORE}" --wdlOutputDialect miniwdl --wdlOutputDirectory "${OUT_DIR}" --wdlOutputFile "${OUT_FILE}" --clean=onSuccess --batchSystem slurm --slurmTime 11:59:59 --statePollingWait 30 --disableProgress --caching=False --retryCount 0 --stopOnFirstFailure=True --logFile "${LOG_FILE}" --batchLogsDir "${LOG_DIR}" --writeMessages "${BUS_PATH}" --runImportsOnWorkers '#workflow/github.com/vgteam/vg_wdl/Giraffe:lr-giraffe' "${INPUT_JSON}"
    fi

done
