#!/usr/bin/env bash
# make_pbsim_reads.sh: script to simulate reads with pbsim2.
# Intended to run on UCSC behind-the-firewall systems
# You may also need to CFLAGS=-fPIC pip3 install --user bioconvert

set -ex

# Here we use : and := to set variables to default values if not present in the environment.
# You can set these in the environment to override them and I don't have to write a CLI option parser.
# See https://stackoverflow.com/a/28085062

# Graph to simulate from. Can be S3 URLs or local file paths. If GRAPH_GBZ_URL
# is set, GRAPH_XG_URL and GRAPH_GBWT_URL are not used.
: "${GRAPH_XG_URL:=s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.xg}"
: "${GRAPH_GBWT_URL:=s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.gbwt}"
: "${GRAPH_GBZ_URL:=""}"
# Name to use for graph when downloaded
: "${GRAPH_NAME:=hprc-v1.0-mc-grch38}"
# Sample to simulate from
: "${SAMPLE_NAME:=HG00741}"
# Sample to name output as coming from (in case you need multiple replicates of a sample)
: "${SAMPLE_NAME_OUT:=${SAMPLE_NAME}}"
# Technology name to use in output filenames
: "${TECH_NAME:=hifi}"
# FASTQ to use as a template, or "/dev/null"
# MUST NOT be a .fastq.gz; no error will be raised but it will not work!
: "${SAMPLE_FASTQ:=/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads/real/hifi/HiFi_reads_100k.fq}"
# HMM model to use instead of a FASTQ, or "/dev/null"
: "${PBSIM_HMM:=/dev/null}"
# This needs to be the pbsim2 binary, which might not be in $PATH.
# It can be installed with
# git clone https://github.com/yukiteruono/pbsim2.git
# cd pbsim2
# git checkout eeb5a19420534a0f672c81db2670117e62a9ee38
# autoupdate
# automake --add-missing
# autoreconf 
# ./configure --prefix=$HOME/.local && make
# The binary will be in src/pbsim
: "${PBSIM:=pbsim}"
# Parameters to use with pbsim for simulating reads for each contig. Parameters are space-separated and internal spaces must be escaped.
: "${PBSIM_PARAMS:=--depth 4 --accuracy-min 0.00 --length-min 10000 --difference-ratio 6:50:54}"
# This needs to be the vg binary. Since this script may run for a long time and
# on one node in a cluster with a shared filesystem, it needs to be a vg binary
# that will not be replaced while the script is running, or vg may crash with a
# bus error.
: "${VG:=vg}"
# Directory to save results in
: "${OUT_DIR:=./reads/sim/${TECH_NAME}/${SAMPLE_NAME_OUT}}"
# Number of MAFs to convert at once
: "${MAX_JOBS:=10}"
# Group to make output directory writable to
: "${OUT_DIR_GROUP=""}"

if [[ "${WORK_DIR}" == "" ]] ; then
    # Make a work directory
    WORK_DIR="$(mktemp -d)"
    CLEAN_WORK_DIR=1
else
    # Let the user send one in in the environment.
    CLEAN_WORK_DIR=0
fi


# Make sure scratch directory exists
mkdir -p "${WORK_DIR}"

if [[ -z "${GRAPH_GBZ_URL}" ]] ; then

    # Fetch graph
    if [[ ! -e "${WORK_DIR}/${GRAPH_NAME}.xg" ]] ; then
        # This comparison require Bash 3 or later. See <https://stackoverflow.com/a/2172365>
        if [[ ${GRAPH_XG_URL} =~ ^s3:.* ]]; then
            # Download from S3
            aws s3 cp "${GRAPH_XG_URL}" "${WORK_DIR}/${GRAPH_NAME}.xg.tmp"
            mv "${WORK_DIR}/${GRAPH_NAME}.xg.tmp" "${WORK_DIR}/${GRAPH_NAME}.xg"
        else
            # Use local symlink
            ln -s "$(realpath "${GRAPH_XG_URL}")" "${WORK_DIR}/${GRAPH_NAME}.xg"
        fi
    fi
    if [[ ! -e "${WORK_DIR}/${GRAPH_NAME}.gbwt" ]] ; then
        if [[ ${GRAPH_GBWT_URL} =~ ^s3:.* ]]; then
            # Download from S3
            aws s3 cp "${GRAPH_GBWT_URL}" "${WORK_DIR}/${GRAPH_NAME}.gbwt.tmp"
            mv "${WORK_DIR}/${GRAPH_NAME}.gbwt.tmp" "${WORK_DIR}/${GRAPH_NAME}.gbwt"
        else
            # Use local symlink
            ln -s "$(realpath "${GRAPH_GBWT_URL}")" "${WORK_DIR}/${GRAPH_NAME}.gbwt"
        fi
    fi

    if [[ ! -e "${WORK_DIR}/${GRAPH_NAME}.gbz" ]] ; then
        # Make it one file
        time "${VG}" gbwt -x "${WORK_DIR}/${GRAPH_NAME}.xg" "${WORK_DIR}/${GRAPH_NAME}.gbwt" --gbz-format -g "${WORK_DIR}/${GRAPH_NAME}.gbz.tmp"
        mv "${WORK_DIR}/${GRAPH_NAME}.gbz.tmp" "${WORK_DIR}/${GRAPH_NAME}.gbz"
    fi

elif [[ ! -e "${WORK_DIR}/${GRAPH_NAME}.gbz" ]] ; then
    # Fetch the GBZ
    if [[ ${GRAPH_GBZ_URL} =~ ^s3:.* ]]; then
        # Download from S3
        aws s3 cp "${GRAPH_GBZ_URL}" "${WORK_DIR}/${GRAPH_NAME}.gbz.tmp"
        mv "${WORK_DIR}/${GRAPH_NAME}.gbz.tmp" "${WORK_DIR}/${GRAPH_NAME}.gbz"
    else
        # Use local symlink
        ln -s "$(realpath "${GRAPH_GBZ_URL}")" "${WORK_DIR}/${GRAPH_NAME}.gbz"
    fi

fi

if [[ ! -e "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz" ]] ; then
    # Make it have our sample as the reference
    "${VG}" gbwt -Z ${WORK_DIR}/${GRAPH_NAME}.gbz --set-tag "reference_samples=${SAMPLE_NAME}" --gbz-format -g "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz.tmp"
    mv "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz.tmp" "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz"
fi

if [[ ! -e "${WORK_DIR}/${SAMPLE_NAME}.fa" ]] ; then
    # Extract sample assembly FASTA from graph where sample is the *reference*. If
    # we do it from the one where the sample is haplotypes, we get different path
    # name strings and we can't inject without hacking them up. We leave the code
    # to hack them up anyway though, for reference later.
    "${VG}" paths -x "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME}-as-ref.gbz" \
        --sample "${SAMPLE_NAME}" \
        --extract-fasta \
      > "${WORK_DIR}/${SAMPLE_NAME}.fa.tmp"
    mv "${WORK_DIR}/${SAMPLE_NAME}.fa.tmp" "${WORK_DIR}/${SAMPLE_NAME}.fa"
fi

# Make a directory for our sample and tech so multiple jobs can share a reference.
SAMPLE_WORK_DIR="${WORK_DIR}/${SAMPLE_NAME_OUT}-${TECH_NAME}-reads"

if [[ -d "${SAMPLE_WORK_DIR}" && "$(ls "${SAMPLE_WORK_DIR}/"sim_*.maf | wc -l)" == "0" ]] ; then
    # Sim directory exists but has no MAFs. Shouldn't have any files at all.
    rmdir "${SAMPLE_WORK_DIR}"
fi

if [[ ! -d "${SAMPLE_WORK_DIR}" ]] ; then
    rm -Rf "${SAMPLE_WORK_DIR}.tmp"
    mkdir "${SAMPLE_WORK_DIR}.tmp"
    
    if [[ "${PBSIM_HMM}" != "/dev/null" ]] ; then
        if [[ "${SAMPLE_FASTQ}" != "/dev/null" ]] ; then
            echo "Can't use both a PBSIM_HMM and a SAMPLE_FASTQ"
            exit 1
        fi
        # Using an HMM to make qualities.
        QUAL_SOURCE_ARGS=(--hmm_model "${SAMPLE_FASTQ}")
    else
        # Using a FASTQ to make qualities.
        # No read may be over 1 megabase or pbsim2 will crash.
        QUAL_SOURCE_ARGS=(--sample-fastq "${SAMPLE_FASTQ}")
    fi
    
    # Simulate reads
    time "${PBSIM}" \
        ${PBSIM_PARAMS} \
       "${QUAL_SOURCE_ARGS[@]}" \
       --prefix "${SAMPLE_WORK_DIR}.tmp/sim" \
       "${WORK_DIR}/${SAMPLE_NAME}.fa"
    
    mv "${SAMPLE_WORK_DIR}.tmp" "${SAMPLE_WORK_DIR}"
fi

function do_job() {
    # Run this file in a job
    set -e
    
    SAM_NAME="${MAF_NAME%.maf}.sam"
    FASTQ_NAME="${MAF_NAME%.maf}.fastq"
    REF_NAME="${MAF_NAME%.maf}.ref"
    RENAMED_BAM_NAME="${MAF_NAME%.maf}.renamed.bam"
    # Get the contig name in the format it would be as a reference sense path.
    # It may already be a reference sense path.
    # Can't run under pipefail because some of these may not match.
    CONTIG_NAME="$(cat "${REF_NAME}" | head -n1 | sed 's/^>//' | sed 's/ .*//' | sed 's/#\([0-9]*\)$/[\1]/')"
    # Haplotype paths can end in a 0 offset/fragment but reference paths don't include that in the name.
    CONTIG_NAME="${CONTIG_NAME%\[0\]}"
    if [[ ! -e "${RENAMED_BAM_NAME}" ]] ; then
        echo "Making ${RENAMED_BAM_NAME}..."
        if [[ ! -e "${SAM_NAME}" ]] ; then
            echo "Making SAM ${SAM_NAME}..."
            /usr/bin/time -v bioconvert maf2sam --force "${MAF_NAME}" "${SAM_NAME}.tmp" 
            mv "${SAM_NAME}.tmp" "${SAM_NAME}"
        fi
        set -o pipefail
        python3 "$(dirname -- "${BASH_SOURCE[0]}")/reinsert_qualities.py" -s "${SAM_NAME}" -f "${FASTQ_NAME}" | sed "s/ref/${CONTIG_NAME}/g" | samtools view -b - > "${RENAMED_BAM_NAME}.tmp"
        set +o pipefail
        mv "${RENAMED_BAM_NAME}.tmp" "${RENAMED_BAM_NAME}"
    else
        echo "Already have ${RENAMED_BAM_NAME}..."
    fi
}

# Convert all the reads to BAM in the space of the sample as a primary reference
for MAF_NAME in "${SAMPLE_WORK_DIR}/"sim_*.maf ; do
    if [[ "${MAX_JOBS}" == "1" ]] ; then
        # Serial mode
        do_job
    else
        # Parallel mode
        while [[ "$(jobs -p | wc -l)" -ge "${MAX_JOBS}" ]] ; do
            # Don't do too much in parallel
            # Fake wait on any job without wait -n
            sleep 0.5
        done
        (
            do_job
        ) &
        ((RUNNING_JOBS += 1))
    fi
done
# Wait on all jobs
wait

if [[ "$(ls "${SAMPLE_WORK_DIR}"/sim_*.tmp 2>/dev/null | wc -l)" != "0" ]] ; then
    # Make sure all the per-file temp files got moved 
    echo "Loose temp files; failure detected."
    exit 1
fi

if [[ ! -e "${SAMPLE_WORK_DIR}/merged.bam" ]] ; then
    # Combine all the BAM files
    time samtools merge -n "${SAMPLE_WORK_DIR}"/sim_*.renamed.bam -o "${SAMPLE_WORK_DIR}/merged.bam.tmp" --threads 14
    mv "${SAMPLE_WORK_DIR}/merged.bam.tmp" "${SAMPLE_WORK_DIR}/merged.bam"
fi

if [[ ! -e "${SAMPLE_WORK_DIR}/injected.gam" ]] ; then
    # Move reads into graph space
    time "${VG}" inject -x "${WORK_DIR}/${GRAPH_NAME}-${SAMPLE_NAME_OUT}-as-ref.gbz" "${SAMPLE_WORK_DIR}/merged.bam" -t 16 >"${SAMPLE_WORK_DIR}/injected.gam.tmp"
    mv "${SAMPLE_WORK_DIR}/injected.gam.tmp" "${SAMPLE_WORK_DIR}/injected.gam"
fi

if [[ ! -e "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.gam" ]] ; then
    # Annotate reads with linear reference positions
    time "${VG}" annotate -x "${WORK_DIR}/${GRAPH_NAME}.gbz" -a "${SAMPLE_WORK_DIR}/injected.gam" --multi-position --search-limit=-1 -t 16 >"${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.gam.tmp"
    mv "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.gam.tmp" "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.gam"
fi

# Work out howe many reads there are
TOTAL_READS="$("${VG}" stats -a "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.gam" | grep "^Total alignments:" | cut -f2 -d':' | tr -d ' ')"

if [[ "${TOTAL_READS}" -lt 1000500 ]] ; then
    echo "Only ${TOTAL_READS} reads were simulated. Cannot subset to 1000000 reads with buffer!"
    exit 1
fi
echo "Simulated ${TOTAL_READS} reads overall"

SUBSAMPLE_SEED=1
for READ_COUNT in 100 1000 10000 100000 1000000 ; do
    # Subset to manageable sizes (always)
    # Get the fraction of reads to keep, overestimated, with no leading 0, to paste onto subsample seed.
    FRACTION="$(echo "(${READ_COUNT} + 500)/${TOTAL_READS}" | bc -l | sed 's/^[0-9]*//g')"
    
    if [[ ! -e "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}-${READ_COUNT}.gam" ]] ; then
        "${VG}" filter -d "${SUBSAMPLE_SEED}${FRACTION}" "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.gam" >"${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.coarse.gam"
        "${VG}" gamsort --shuffle "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.coarse.gam" >"${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.coarse.shuffled.gam"
        "${VG}" filter -t1 --max-reads "${READ_COUNT}" "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.coarse.shuffled.gam" >"${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}-${READ_COUNT}.gam.tmp"
        mv "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}-${READ_COUNT}.gam.tmp" "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}-${READ_COUNT}.gam"
    fi
   
    ((SUBSAMPLE_SEED+=1))    
done

# Output them
mkdir -p "${OUT_DIR}"
cp "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}.gam" "${SAMPLE_WORK_DIR}/${SAMPLE_NAME_OUT}-sim-${TECH_NAME}-"*".gam" "${OUT_DIR}/"

if [[ ! -z "${OUT_DIR_GROUP}" ]] ; then
    # Make output directory owned by and writable to this group, so Snakemake
    # can make subset files.
    chgrp -R "${OUT_DIR_GROUP}" "${OUT_DIR}"
    chmod g+w "${OUT_DIR}"
fi

if [[ "${CLEAN_WORK_DIR}" == "1" ]] ; then
    # Clean up the work directory
    rm -Rf "${WORK_DIR}"
fi
