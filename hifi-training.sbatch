#!/usr/bin/env bash
#SBATCH --array 1
#SBATCH -c8
#SBATCH --mem 64G
#SBATCH --partition long
#SBATCH --time 40:00:00
set -e
set -x
set -u
export SINGULARITY_DOCKER_HUB_MIRROR=http://10.50.100.88
CONDITION_NUMBER=0
EXPERIMENT_NAME="hifi-training"

for FASTQ_URL in \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio/HG002.m84005_220827_014912_s1.GRCh38.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio/HG002.m84005_220919_232112_s2.GRCh38.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio/HG002.m84011_220902_175841_s1.GRCh38.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio/HG003.m84010_220919_235306_s2.GRCh38.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio/HG004.m84010_220919_232145_s1.GRCh38.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio/NA12878.GRCh38.haplotagged.35x.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio/NA12878.GRCh38.haplotagged.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio_amp/HG001.Twist.2_cells.32x.grch38.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio_amp/HG002.Twist.2_cells.33x.grch38.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/revio_amp/HG002.ULI.2_cells.60x.grch38.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/sequel2/HG002.pfda_challenge.grch38.phased.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/sequel2/HG003.pfda_challenge.grch38.phased.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/sequel2/HG004.pfda_challenge.grch38.phased.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/sprq/aligned.sorted.hg002.m84029_240730_192004_s1.hifi_reads.ccs.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/sprq/aligned.sorted.hg002.m84034_240730_194815_s1.hifi_reads.ccs.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/sprq/aligned.sorted.hg003.m84028_240802_214905_s1.hifi_reads.ccs.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/sprq/aligned.sorted.hg003.m84029_240730_212142_s2.hifi_reads.ccs.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/sprq/aligned.sorted.hg004.m84030_240730_192221_s1.hifi_reads.ccs.fastq.gz \
https://storage.googleapis.com/brain-genomics/awcarroll/pacbio_training_2024/fastq/sprq/aligned.sorted.hg004.m84030_240802_210346_s1.hifi_reads.ccs.fastq.gz \
; do
    
    ((CONDITION_NUMBER+=1))
    if [[ "${CONDITION_NUMBER}" != "${SLURM_ARRAY_TASK_ID}" ]] ; then
        # Not our condition
        continue
    fi
    
    FASTQ_NAME="$(basename "${FASTQ_URL}")"
    BASE_NAME="${FASTQ_NAME%.fastq.gz}"
    SAMPLE_NAME="$(echo "${BASE_NAME}" | tr '.' '\n' | tr 'a-z' 'A-Z' | grep "^\(HG\|NA\)")"
    
    mkdir -p /private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/hifi
    INPUT_JSON="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/hifi/${BASE_NAME}.input.json"
    JOB_STORE="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/hifi/${BASE_NAME}.tree"
    OUT_DIR="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/hifi/${BASE_NAME}.dir"
    OUT_FILE="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/hifi/${BASE_NAME}.output.json"
    LOG_FILE="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/hifi/${BASE_NAME}.log"
    LOG_DIR="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/training-data/hifi/${BASE_NAME}.batchlog"
    mkdir -p "${LOG_DIR}"

    cat >"${INPUT_JSON}" <<EOF
{
    "Giraffe.INPUT_READ_FILE_FIRST": "${FASTQ_URL}",
    "Giraffe.SAMPLE_NAME": "${SAMPLE_NAME}",
    "Giraffe.PAIRED_READS": false,
    "Giraffe.HAPLOTYPE_SAMPLING": true,
    "Giraffe.DIPLOID": false,
    "Giraffe.SET_REFERENCE": "GRCh38",
    "Giraffe.HAPLOTYPE_NUMBER": 16,
    "Giraffe.VG_DOCKER": "quay.io/vgteam/vg:ci-1792-344434952cfede3c8a34cb7584f1732db9832a4a",
    "Giraffe.GIRAFFE_PRESET": "hifi",
    "Giraffe.READS_PER_CHUNK": 150000,
    "Giraffe.GBZ_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v2.prerelease.3-mc-chm13-eval-jan29.gbz",
    "Giraffe.HAPL_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v2.prerelease.3-mc-chm13-eval-jan29.fragment.hapl",
    "Giraffe.REFERENCE_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/grch38-pansn.fa",
    "Giraffe.REFERENCE_INDEX_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/grch38-pansn.fa.fai",
    "Giraffe.REFERENCE_DICT_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/grch38-pansn.fa.paths.callable.dict",
    "Giraffe.PATH_LIST_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/grch38-pansn.fa.paths.callable.txt",
    "Giraffe.REFERENCE_PREFIX": "GRCh38#0#",
    "Giraffe.REALIGN_INDELS": false,
    "Giraffe.LEFTALIGN_BAM": true
}
EOF

    /private/home/anovak/workspace/toil/venv/bin/toil-wdl-runner --jobStore "${JOB_STORE}" --wdlOutputDialect miniwdl --wdlOutputDirectory "${OUT_DIR}" --wdlOutputFile "${OUT_FILE}" --clean=onSuccess --batchSystem slurm --slurmTime 11:59:59 --statePollingWait 30 --disableProgress --caching=False --retryCount 2 --logFile "${LOG_FILE}" --batchLogsDir "${LOG_DIR}" '#workflow/github.com/vgteam/vg_wdl/Giraffe:lr-giraffe' "${INPUT_JSON}"

done
