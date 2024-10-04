#!/usr/bin/env bash
#SBATCH --array 1
#SBATCH -c16 
#SBATCH --mem 32G 
#SBATCH --partition long
#SBATCH --time 7-00:00:00
set -e
set -x
set -u

export SINGULARITY_DOCKER_HUB_MIRROR=http://10.50.100.88

CONDITION_NUMBER=0
EXPERIMENT_NAME=dvtest6
for VG_VERSION in new; do
    for DV_VERSION in 1.5 ; do
        for DV_MODEL in default ; do
            if [[ "${DV_MODEL}" == "trained" && "${DV_VERSION}" == "1.5" ]] ; then
                continue
            fi
            if [[ "${DV_MODEL}" == "1.7_default" && ("${DV_VERSION}" == "1.5" || "${DV_VERSION}" == "1.6.1") ]] ; then
                continue
            fi
            if [[ "${DV_MODEL}" == "1.6.1_default" && "${DV_VERSION}" == "1.5" ]] ; then
                continue
            fi
            if [[ "${DV_MODEL}" == "default" && ("${DV_VERSION}" != "1.5" && "${DV_VERSION}" != "1.6.1") ]] ; then
                continue
            fi
            for MIN_MAPQ in 0; do
                ((CONDITION_NUMBER+=1))
                if [[ "${CONDITION_NUMBER}" != "${SLURM_ARRAY_TASK_ID}" ]] ; then
                    # Not our condition
                    continue
                fi
                CONDITION_NAME="${VG_VERSION}_vg-${DV_VERSION}_dv-${DV_MODEL}_model-${MIN_MAPQ}_min_mapq"
                INPUTS_FILE="inputs-${EXPERIMENT_NAME}-${CONDITION_NAME}.json"
                cat >"${INPUTS_FILE}" <<EOF
{
"GiraffeDeepVariant.CONTIGS": ["CHM13#0#chr1", "CHM13#0#chr2", "CHM13#0#chr3",  "CHM13#0#chr4", "CHM13#0#chr5", "CHM13#0#chr6", "CHM13#0#chr7", "CHM13#0#chr8", "CHM13#0#chr9", "CHM13#0#chr10", "CHM13#0#chr11", "CHM13#0#chr12", "CHM13#0#chr13", "CHM13#0#chr14", "CHM13#0#chr15", "CHM13#0#chr16", "CHM13#0#chr17", "CHM13#0#chr18", "CHM13#0#chr19", "CHM13#0#chr20", "CHM13#0#chr21", "CHM13#0#chr22", "CHM13#0#chrX"],
"GiraffeDeepVariant.DIST_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v1.1-mc-chm13.d9.dist",
"GiraffeDeepVariant.DV_KEEP_LEGACY_AC": false,
"GiraffeDeepVariant.DV_MODEL_TYPE": "PACBIO",
"GiraffeDeepVariant.DV_NORM_READS": true,
"GiraffeDeepVariant.EVALUATION_REGIONS_BED": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
"GiraffeDeepVariant.GBZ_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v1.1-mc-chm13.d9.gbz",
"GiraffeDeepVariant.GIRAFFE_PRESET": "hifi",
"GiraffeDeepVariant.INPUT_READ_FILE_1": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads/real/hifi/HG002/HiFi_DC_v1.2_HG002_combined_unshuffled.fq.gz",
"GiraffeDeepVariant.LEFTALIGN_BAM": true,
"GiraffeDeepVariant.MIN_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v1.1-mc-chm13.d9.k31.w50.W.withzip.min",
"GiraffeDeepVariant.OUTPUT_SINGLE_BAM": true,
"GiraffeDeepVariant.OUTPUT_CALLING_BAMS": false,
"GiraffeDeepVariant.PAIRED_READS": false,
"GiraffeDeepVariant.PRUNE_LOW_COMPLEXITY": true,
"GiraffeDeepVariant.READS_PER_CHUNK": 150000,
"GiraffeDeepVariant.REALIGN_INDELS": false,
"GiraffeDeepVariant.REFERENCE_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/chm13v2.0.fa",
"GiraffeDeepVariant.REFERENCE_PREFIX": "CHM13#0#",
"GiraffeDeepVariant.RESTRICT_REGIONS_BED": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/chm13v2.0.fa.callable.from.chm13.bed",
"GiraffeDeepVariant.RUN_STANDALONE_VCFEVAL": false,
"GiraffeDeepVariant.SAMPLE_NAME": "HG002",
"GiraffeDeepVariant.TRUTH_VCF": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz",
"GiraffeDeepVariant.TRUTH_VCF_INDEX": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz.tbi",
"GiraffeDeepVariant.EVAL_MEM": 100,
"GiraffeDeepVariant.ZIPCODES_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v1.1-mc-chm13.d9.k31.w50.W.zipcodes",
EOF
                
                if [[ "${VG_VERSION}" == "old" ]] ; then
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.VG_DOCKER": "quay.io/adamnovak/vg:bd8ded6",'
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.VG_SURJECT_DOCKER": "quay.io/adamnovak/vg:bd8ded6",'
                elif [[ "${VG_VERSION}" == "new" ]] ; then
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.VG_DOCKER": "quay.io/adamnovak/vg:beec239",'
                fi

                if [[ "${DV_VERSION}" == "668585654" || "${DV_VERSION}" == "680683261" ]] ; then
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_GPU_DOCKER": "gcr.io/deepvariant-docker/deepvariant:head'${DV_VERSION}'",'
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_NO_GPU_DOCKER": "gcr.io/deepvariant-docker/deepvariant:head'${DV_VERSION}'",'
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_IS_1_7_OR_NEWER": true,'
                elif [[ "${DV_VERSION}" == "1.6.1" ]] ; then
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_GPU_DOCKER": "google/deepvariant:'${DV_VERSION}'-gpu",'
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_NO_GPU_DOCKER": "google/deepvariant:'${DV_VERSION}'",'
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_IS_1_7_OR_NEWER": false,'
                fi

                if [[ "${DV_MODEL}" == "trained" ]] ; then
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_MODEL_FILES": ["/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/models/hifi/2024-09-18/example_info.json", "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/models/hifi/2024-09-18/weights-67-0.994181.ckpt.index", "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/models/hifi/2024-09-18/weights-67-0.994181.ckpt.data-00000-of-00001"],'
                elif [[ "${DV_MODEL}" == "1.7_default" ]] ; then
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_MODEL_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/fingerprint.pb", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/saved_model.pb", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/example_info.json"],'
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_MODEL_VARIABLES_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.data-00000-of-00001", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.index"],'
                elif [[ "${DV_MODEL}" == "1.6.1_default" ]] ; then
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_MODEL_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/fingerprint.pb", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/saved_model.pb", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/example_info.json"],'
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_MODEL_VARIABLES_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.data-00000-of-00001", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.index"],'
                elif [[ "${DV_MODEL}" == "1.5_default" ]] ; then
                    echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.DV_MODEL_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.5.0/DeepVariant-inception_v3-1.5.0+data-pacbio_standard/model.ckpt.data-00000-of-00001","https://storage.googleapis.com/deepvariant/models/DeepVariant/1.5.0/DeepVariant-inception_v3-1.5.0+data-pacbio_standard/model.ckpt.index","https://storage.googleapis.com/deepvariant/models/DeepVariant/1.5.0/DeepVariant-inception_v3-1.5.0+data-pacbio_standard/model.ckpt.meta","https://storage.googleapis.com/deepvariant/models/DeepVariant/1.5.0/DeepVariant-inception_v3-1.5.0+data-pacbio_standard/model.ckpt.example_info.json"],'
                fi
                
                echo >>"${INPUTS_FILE}" '"GiraffeDeepVariant.MIN_MAPQ": '${MIN_MAPQ}

                echo >>"${INPUTS_FILE}" '}'
                
                OUTPUTS_FILE="toil-result-${EXPERIMENT_NAME}-${CONDITION_NAME}.json"
                
                echo "Running ${CONDITION_NAME}"
                mkdir -p ./output/${EXPERIMENT_NAME}
                MINIWDL__CALL_CACHE__GET=true MINIWDL__CALL_CACHE__PUT=true MINIWDL__CALL_CACHE__DIR=$(pwd)/output/miniwdl-cache /private/home/anovak/workspace/toil/venv/bin/toil-wdl-runner https://raw.githubusercontent.com/vgteam/vg_wdl/fde5f8f17e5f406fce75daa7a5d1b52f08710de0/workflows/giraffe_and_deepvariant.wdl ${INPUTS_FILE} --wdlOutputDialect miniwdl --wdlOutputDirectory ./output/${EXPERIMENT_NAME}/${CONDITION_NAME} --wdlOutputFile ${OUTPUTS_FILE} --jobStore ./output/${EXPERIMENT_NAME}/${CONDITION_NAME}-tree --clean=always --batchSystem slurm --slurmTime 11:59:59 --statePollingWait 30 --disableProgress --caching=False --retryCount 0 --logDebug --logFile ./output/${EXPERIMENT_NAME}/${CONDITION_NAME}.log --writeMessages ./output/${EXPERIMENT_NAME}/${CONDITION_NAME}.bus
            done
        done
    done
done

