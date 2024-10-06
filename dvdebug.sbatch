#!/usr/bin/env bash


#SBATCH --array 1-6
#SBATCH -c16 
#SBATCH --mem 32G 
#SBATCH --partition long
#SBATCH --time 7-00:00:00
set -e
set -x
set -u

export SINGULARITY_DOCKER_HUB_MIRROR=http://10.50.100.88

CONDITION_NUMBER=0
EXPERIMENT_NAME=dvdebug
for MAPPER in giraffe minimap2 ; do
    for REFERENCE in chm13 GRCh38 ; do
        for BENCHMARK in beta v4.2.1 ; do

            VG_VERSION=new
            DV_VERSION=1.6.1
            DV_MODEL=default
            MIN_MAPQ=0

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
            if [[ "${REFERENCE}" == "chm13" && "${BENCHMARK}" == "v4.2.1" ]] ; then
                continue
            fi

            ((CONDITION_NUMBER+=1))
            if [[ "${CONDITION_NUMBER}" != "${SLURM_ARRAY_TASK_ID}" ]] ; then
                # Not our condition
                continue
            fi
            CONDITION_NAME="${MAPPER}-${REFERENCE}-${BENCHMARK}_benchmark"
            INPUTS_FILE="inputs-${EXPERIMENT_NAME}-${CONDITION_NAME}.json"

             if [[ "${MAPPER}" == "giraffe" ]] ; then
                WF_NAME="GiraffeDeepVariant"
                WF_FILE="giraffe_and_deepvariant.wdl"
            elif [[ "${MAPPER}" == "minimap2" ]] ; then
                WF_NAME="DeepVariant"
                WF_FILE="deepvariant.wdl"
            else
                exit 1
            fi

            cat >"${INPUTS_FILE}" <<EOF
{
"${WF_NAME}.DV_KEEP_LEGACY_AC": false,
"${WF_NAME}.DV_MODEL_TYPE": "PACBIO",
"${WF_NAME}.DV_NORM_READS": true,
"${WF_NAME}.LEFTALIGN_BAM": true,
"${WF_NAME}.OUTPUT_SINGLE_BAM": false,
"${WF_NAME}.OUTPUT_CALLING_BAMS": false,
"${WF_NAME}.REALIGN_INDELS": false,
"${WF_NAME}.RUN_STANDALONE_VCFEVAL": false,
"${WF_NAME}.SAMPLE_NAME": "HG002",
EOF

            if [[ "${DV_VERSION}" == "668585654" || "${DV_VERSION}" == "680683261" ]] ; then
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_GPU_DOCKER": "gcr.io/deepvariant-docker/deepvariant:head'${DV_VERSION}'",'
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_NO_GPU_DOCKER": "gcr.io/deepvariant-docker/deepvariant:head'${DV_VERSION}'",'
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_IS_1_7_OR_NEWER": true,'
            elif [[ "${DV_VERSION}" == "1.6.1" ]] ; then
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_GPU_DOCKER": "google/deepvariant:'${DV_VERSION}'-gpu",'
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_NO_GPU_DOCKER": "google/deepvariant:'${DV_VERSION}'",'
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_IS_1_7_OR_NEWER": false,'
            fi

            if [[ "${MAPPER}" == "giraffe" ]] ; then
                if [[ "${VG_VERSION}" == "old" ]] ; then
                    echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.VG_DOCKER": "quay.io/adamnovak/vg:bd8ded6",'
                    echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.VG_SURJECT_DOCKER": "quay.io/adamnovak/vg:bd8ded6",'
                elif [[ "${VG_VERSION}" == "new" ]] ; then
                    echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.VG_DOCKER": "quay.io/adamnovak/vg:beec239",'
                else
                    exit 1
                fi

                cat >>"${INPUTS_FILE}" <<EOF
"${WF_NAME}.INPUT_READ_FILE_1": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads/real/hifi/HG002/HiFi_DC_v1.2_HG002_combined_unshuffled.fq.gz",
"${WF_NAME}.PAIRED_READS": false,
"${WF_NAME}.READS_PER_CHUNK": 150000,
"${WF_NAME}.DIST_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v1.1-mc-chm13.d9.dist",
"${WF_NAME}.ZIPCODES_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v1.1-mc-chm13.d9.k31.w50.W.zipcodes",
"${WF_NAME}.MIN_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v1.1-mc-chm13.d9.k31.w50.W.withzip.min",
"${WF_NAME}.GBZ_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs/hprc-v1.1-mc-chm13.d9.gbz",
"${WF_NAME}.GIRAFFE_PRESET": "hifi",
"${WF_NAME}.PRUNE_LOW_COMPLEXITY": true,
EOF
            elif [[ "${MAPPER}" == "minimap2" ]] ; then
                if [[ "${REFERENCE}" == "chm13" ]] ; then
                    BAM_FILE="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/control-alignments/chm13v2.0/minimap2-lr:hq/real/hifi/HiFi_DC_v1.2_HG002_combined_unshuffled.sorted.nosecondaries.bam"
                    
                elif [[ "${REFERENCE}" == "GRCh38" ]] ; then
                    BAM_FILE="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/control-alignments/GCA_000001405.15_GRCh38_no_alt_analysis_set/minimap2-lr:hq/real/hifi/HiFi_DC_v1.2_HG002_combined_unshuffled.sorted.nosecondaries.bam"
                else
                    exit 1
                fi
                cat >>"${INPUTS_FILE}" <<EOF
"${WF_NAME}.MERGED_BAM_FILE": "${BAM_FILE}",
"${WF_NAME}.MERGED_BAM_FILE_INDEX": "${BAM_FILE}.bai",
EOF
            else
                exit 1
            fi
            

            if [[ "${BENCHMARK}" == "v4.2.1" ]] ; then
                if [[ "${REFERENCE}" == "GRCh38" ]] ; then
                    TRUTH_VCF="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
                    TRUTH_BED="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
                else
                    exit 1
                fi
            elif [[ "${BENCHMARK}" == "beta" ]] ; then
                if [[ "${REFERENCE}" == "chm13" ]] ; then
                    TRUTH_VCF="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz"
                    TRUTH_BED="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed"
                elif [[ "${REFERENCE}" == "GRCh38" ]] ; then
                    TRUTH_VCF="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/GRCh38_HG2-T2TQ100-V1.1.vcf.gz"
                    TRUTH_BED="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed"
                else
                    exit 1
                fi
            fi
            cat >>"${INPUTS_FILE}" <<EOF
"${WF_NAME}.TRUTH_VCF": "${TRUTH_VCF}",
"${WF_NAME}.TRUTH_VCF_INDEX": "${TRUTH_VCF}.tbi",
"${WF_NAME}.EVALUATION_REGIONS_BED": "${TRUTH_BED}",
EOF

            if [[ "${REFERENCE}" == "chm13" ]] ; then
                cat >>"${INPUTS_FILE}" <<EOF
"${WF_NAME}.REFERENCE_PREFIX": "CHM13#0#",
"${WF_NAME}.CONTIGS": ["CHM13#0#chr1", "CHM13#0#chr2", "CHM13#0#chr3",  "CHM13#0#chr4", "CHM13#0#chr5", "CHM13#0#chr6", "CHM13#0#chr7", "CHM13#0#chr8", "CHM13#0#chr9", "CHM13#0#chr10", "CHM13#0#chr11", "CHM13#0#chr12", "CHM13#0#chr13", "CHM13#0#chr14", "CHM13#0#chr15", "CHM13#0#chr16", "CHM13#0#chr17", "CHM13#0#chr18", "CHM13#0#chr19", "CHM13#0#chr20", "CHM13#0#chr21", "CHM13#0#chr22", "CHM13#0#chrX"],
"${WF_NAME}.REFERENCE_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/chm13v2.0.fa",
"${WF_NAME}.RESTRICT_REGIONS_BED": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/chm13v2.0.fa.callable.from.chm13.bed",
EOF
            elif [[ "${REFERENCE}" == "GRCh38" ]] ; then
                 cat >>"${INPUTS_FILE}" <<EOF
"${WF_NAME}.REFERENCE_PREFIX": "GRCh38#0#",
"${WF_NAME}.CONTIGS": ["GRCh38#0#chr1", "GRCh38#0#chr2", "GRCh38#0#chr3",  "GRCh38#0#chr4", "GRCh38#0#chr5", "GRCh38#0#chr6", "GRCh38#0#chr7", "GRCh38#0#chr8", "GRCh38#0#chr9", "GRCh38#0#chr10", "GRCh38#0#chr11", "GRCh38#0#chr12", "GRCh38#0#chr13", "GRCh38#0#chr14", "GRCh38#0#chr15", "GRCh38#0#chr16", "GRCh38#0#chr17", "GRCh38#0#chr18", "GRCh38#0#chr19", "GRCh38#0#chr20", "GRCh38#0#chr21", "GRCh38#0#chr22", "GRCh38#0#chrX", "GRCh38#0#chrY"],
"${WF_NAME}.REFERENCE_FILE": "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
EOF
            else
                exit 1
            fi
                
            if [[ "${DV_MODEL}" == "trained" ]] ; then
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_MODEL_FILES": ["/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/models/hifi/2024-09-18/example_info.json", "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/models/hifi/2024-09-18/weights-67-0.994181.ckpt.index", "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/models/hifi/2024-09-18/weights-67-0.994181.ckpt.data-00000-of-00001"],'
            elif [[ "${DV_MODEL}" == "1.7_default" ]] ; then
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_MODEL_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/fingerprint.pb", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/saved_model.pb", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/example_info.json"],'
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_MODEL_VARIABLES_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.data-00000-of-00001", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.7.0/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.index"],'
            elif [[ "${DV_MODEL}" == "1.6.1_default" ]] ; then
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_MODEL_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/fingerprint.pb", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/saved_model.pb", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/example_info.json"],'
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_MODEL_VARIABLES_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.data-00000-of-00001", "https://storage.googleapis.com/deepvariant/models/DeepVariant/1.6.1/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.index"],'
            elif [[ "${DV_MODEL}" == "1.5_default" ]] ; then
                echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.DV_MODEL_FILES": ["https://storage.googleapis.com/deepvariant/models/DeepVariant/1.5.0/DeepVariant-inception_v3-1.5.0+data-pacbio_standard/model.ckpt.data-00000-of-00001","https://storage.googleapis.com/deepvariant/models/DeepVariant/1.5.0/DeepVariant-inception_v3-1.5.0+data-pacbio_standard/model.ckpt.index","https://storage.googleapis.com/deepvariant/models/DeepVariant/1.5.0/DeepVariant-inception_v3-1.5.0+data-pacbio_standard/model.ckpt.meta","https://storage.googleapis.com/deepvariant/models/DeepVariant/1.5.0/DeepVariant-inception_v3-1.5.0+data-pacbio_standard/model.ckpt.example_info.json"],'
            fi
                
            echo >>"${INPUTS_FILE}" '"'${WF_NAME}'.MIN_MAPQ": '${MIN_MAPQ}

            echo >>"${INPUTS_FILE}" '}'
                
            OUTPUTS_FILE="toil-result-${EXPERIMENT_NAME}-${CONDITION_NAME}.json"
            
            echo "Running ${CONDITION_NAME}"
            mkdir -p ./output/${EXPERIMENT_NAME}
            MINIWDL__CALL_CACHE__GET=true MINIWDL__CALL_CACHE__PUT=true MINIWDL__CALL_CACHE__DIR=$(pwd)/output/miniwdl-cache /private/home/anovak/workspace/toil/venv/bin/toil-wdl-runner "https://raw.githubusercontent.com/vgteam/vg_wdl/90940196ef371652f2dc3aa2c0ead0d20234dbf7/workflows/${WF_FILE}" ${INPUTS_FILE} --wdlOutputDialect miniwdl --wdlOutputDirectory ./output/${EXPERIMENT_NAME}/${CONDITION_NAME} --wdlOutputFile ${OUTPUTS_FILE} --jobStore ./output/${EXPERIMENT_NAME}/${CONDITION_NAME}-tree --clean=always --batchSystem slurm --slurmTime 11:59:59 --statePollingWait 30 --disableProgress --caching=False --retryCount 0 --logDebug --logFile ./output/${EXPERIMENT_NAME}/${CONDITION_NAME}.log --writeMessages ./output/${EXPERIMENT_NAME}/${CONDITION_NAME}.bus
        done
    done
done

