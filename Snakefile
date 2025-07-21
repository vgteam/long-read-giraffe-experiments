###########################################################
# Experimental procedure for evaluating Long Read Giraffe #
###########################################################

import parameter_search

import functools
import tempfile
import os
import numpy as np

# Set a default config file. This can be overridden with --configfile.
# See the config file for how to define experiments.
configfile: "lr-config.yaml"

# Where are the input graphs?
#
# For each reference (here "chm13"), this directory must contain:
#
# hprc-v1.1-mc-chm13.d9.gbz
# hprc-v1.1-mc-chm13.d9.dist
#
# Also, it must either be writable, or already contain zipcode and minimizer
# indexes for each set of minimizer indexing parameters (here "k31.w50.W"),
# named like:
#
# hprc-v1.1-mc-chm13.d9.k31.w50.W.withzip.min
# hprc-v1.1-mc-chm13.d9.k31.w50.W.zipcodes
#
# For comparing against old vg versions, it must also ber writable or have a
# non-zipcode-based minimizer index for each graph:
#
# hprc-v1.1-mc-chm13.d9.k31.w50.W.nozip.min
#
# It must also either be writable or contain an hg and gfa file
# for running GraphAligner, named like:
#
# hprc-v1.1-mc-chm13.hg
# hprc-v1.1-mc-chm13.gfa
#
# And for haplotype sampling, the full graph (without any .dXX) must be available:
#
# hprc-v1.1-mc-chm13.gbz
#
# A top-level chains only "distance" index will be made if not present:
#
# hprc-v1.1-mc-chm13.tcdist
#
# As will the haplotype sampling indexes:
#
# hprc-v1.1-mc-chm13.ri
# hprc-v1.1-mc-chm13.fragment.hapl
#
GRAPHS_DIR = config.get("graphs_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs"

#
# For SV calling, we need a truth vcf and a BED file with confident regions
#
# giab6_chm13.vcf.gz
# giab6_chm13.confreg.bed
#
# And a tandem repeat bed for sniffles
# human_chm13v2.0_maskedY_rCRS.trf.bed (from https://raw.githubusercontent.com/PacificBiosciences/pbsv/refs/heads/master/annotations/human_chm13v2.0_maskedY_rCRS.trf.bed)
#
# TODO: For HG001 we might need https://platinum-pedigree-data.s3.amazonaws.com/variants/merged_sv_truthset/GRCh38/merged_hg38.svs.sort.oa.vcf.gz

SV_DATA_DIR = config.get("sv_data_dir", None) or "/private/home/jmonlong/workspace/lreval/data"

# Where are the reads to use?
#
# This directory must have "real" and "sim" subdirectories. Within each, there
# must be a subdirectory for the sequencing technology, and within each of
# those, a subdirectory for the sample.
#
# For real reads, each sample directory must have a ".fq.gz" or ".fastq.gz" file.
# The name of the file must contain the sample name. If the directory is not
# writable, and you want to trim adapters off nanopore reads, there must also
# be a ".trimmed.fq.gz" or ".trimmed.fastq.gz" version of this file, with the
# first 100 and last 10 bases trimmed off. The workflow will generate
# "{basename}.{subset}.fq.gz" files for each subset size in reads ("1k", "1m:,
# etc.) that you want to work with.
#
# For simulated reads, each sample directory must have files
# "{sample}-sim-{tech}-{subset}.gam" for each subset size as a number (100, 1000,
# 1000000, etc.) that you want to work with. If the directory is not writable,
# it must already have abbreviated versions ("1k" or "1m" instead of the full
# number) of the GAM files, and the corresponding extracted ".fq" files.
#
# There can also be a "{sample}-sim-{tech}.{category}.txt" file with the names
# of reads in a gategory (like "centromeric") for analysis of different subsets
# of reads.
#
# Simulated reads should be made with the "make_pbsim_reads.sh" script in this
# repository, or, for paired-end read simulation for Illumina reads:
#
# vg sim -r -n 2500000 -a -s 12345 -p 570 -v 165 -i 0.00029 --multi-position
#
# Using an XG and a GBWT for the target sample-and-reference graph,
# --sample-name with the target sample's name, and the full real read set for
# the training FASTQ.
#
# A fully filled out reads directory might look like:
#.
#├── real
#│   ├── hifi
#│   │   └── HG002
#│   │       ├── HiFi_DC_v1.2_HG002_combined_unshuffled.1k.fq
#│   │       └── HiFi_DC_v1.2_HG002_combined_unshuffled.fq.gz
#│   └── r10
#│       └── HG002
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.fastq.gz
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.fastq.gz
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.10k.fq.gz
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.1k.fq.gz
#│           └── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.1m.fq.gz
#└── sim
#    ├── hifi
#    │   └── HG002
#    │       ├── HG002-sim-hifi.centromeric.txt
#    │       ├── HG002-sim-hifi-1000.gam
#    │       ├── HG002-sim-hifi-10000.gam
#    │       ├── HG002-sim-hifi-1000000.gam
#    │       ├── HG002-sim-hifi-10k.fq
#    │       ├── HG002-sim-hifi-10k.gam
#    │       ├── HG002-sim-hifi-1k.fq
#    │       ├── HG002-sim-hifi-1k.gam
#    │       ├── HG002-sim-hifi-1m.fq
#    │       └── HG002-sim-hifi-1m.gam
#    └── r10
#        └── HG002
#            ├── HG002-sim-r10.centromeric.txt
#            ├── HG002-sim-r10-1000.gam
#            ├── HG002-sim-r10-10000.gam
#            ├── HG002-sim-r10-1000000.gam
#            ├── HG002-sim-r10-10k.fq
#            ├── HG002-sim-r10-10k.gam
#            ├── HG002-sim-r10-1k.fq
#            ├── HG002-sim-r10-1k.gam
#            ├── HG002-sim-r10-1m.fq
#            └── HG002-sim-r10-1m.gam
#
READS_DIR = config.get("reads_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads"

# Where are the linear reference files?
#
# For each reference name (here "chm13") this directory must contain:
#
# A FASTA file with PanSN-style (CHM13#0#chr1) contig names: 
# chm13-pansn-newY.fa
#
# (For grch38 we don't use the -newY part.)
#
# For the calling references (chm13v2.0 and grch38) we also need a plain .fa
# and .fa.fai without pansn names, and _PAR.bed files with the pseudo-autosomal
# regions.
#
# We also use, but can generate:
#
# Index files for Minimap2 for each preset (here "hifi", can also be "ont" or "sr", and can be generated from the FASTA):
# chm13-pansn-newY.hifi.mmi
# 
# A Winnowmap repetitive kmers file:
# chm13-pansn-newY.repetitive_k15.txt
#
# Minimap2 and BWA-MEM indexes
#
REFS_DIR = config.get("refs_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references"

# Where are variant call truth set files kept (for when there isn't a handy hosted URL somewhere).
# These are organized by reference and then sample
TRUTH_DIR = config.get("truth_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/truth-sets"

# Where are custom trained DeepVariant models for Giraffe?
# Models should be in directories by tech and then by date in YYYY-MM-DD format.
MODELS_DIR = config.get("models_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/models/"

# When we "snakemake all_paper_figures", where should the results go?
ALL_OUT_DIR = config.get("all_out_dir", None) or "/private/groups/patenlab/project-lrg"

# What stages does the Giraffe mapper report times for?
STAGES = ["minimizer", "seed", "tree", "fragment", "chain", "align", "winner"]

# What stages does the Giraffe mapper report times for on the non-chainign codepath?
NON_CHAINING_STAGES = ["minimizer", "seed", "cluster", "extend", "align", "pairing", "winner"]

# What aligner and read part combinations does Giraffe report statistics for?
ALIGNER_PARTS = ["wfa_tail", "dozeu_tail", "wfa_middle", "bga_middle"]

# To allow for splitting and variable numbers of output files, we need to know
# the available subset values to generate rules.
KNOWN_SUBSETS = ["100", "1k", "10k", "100k", "1m"]
CHUNK_SIZE = 10000

# For each Slurm partition name, what is its max wall time in minutes?
# TODO: Put this in the config
SLURM_PARTITIONS = [
    ("short", 60),
    ("medium", 12 * 60),
    ("long", 7 * 24 * 60)
]

PARAM_SEARCH = parameter_search.ParameterSearch()

# Where is a large temp directory?
LARGE_TEMP_DIR = config.get("large_temp_dir", "/data/tmp")

#Different phoenix nodes seem to run at different speeds, so we can specify which node to run
#This gets added as a slurm_extra for all the real read runs
REAL_SLURM_EXTRA = config.get("real_slurm_extra", None) or ""

# If set to True, jobs where we care about speed will demand entire nodes.
# If False, they will just use one thread per core.
EXCLUSIVE_TIMING = config.get("exclusive_timing", True)

# How many threads do we want mapping to use?
# TODO: If we're reserving whole nodes we want to use all of the node actually, but only on exclusive runs.
MAPPER_THREADS = 64 if EXCLUSIVE_TIMING else 32

# We may not want to populate the MiniWDL task cache because it makes us take
# more shared disk space.
FILL_WDL_CACHE = "true" if config.get("fill_wdl_cache", True) else "false"

# Figure out what columns to put in a table comparing all the conditions in an experiment.
# TODO: Make this be per-experiment and let multiple tables be defined
IMPORTANT_STATS_TABLE_COLUMNS=config.get("important_stats_table_columns", ["speed_from_log", "softclipped_or_unmapped", "accuracy", "indel_f1", "snp_f1"])

# What versions of Giraffe should be run with the old, non-zipcode-aware indexes
NON_ZIPCODE_GIRAFFE_VERSIONS = set(config.get("non_zipcode_giraffe_versions")) if "non_zipcode_giraffe_versions" in config else set()

# What version of vg should be used to make fragment-aware haplotype indexes?
VG_FRAGMENT_HAPLOTYPE_INDEXING_VERSION="v1.64.1"
# What version of vg should be used to haplotype-sample graphs?
VG_HAPLOTYPE_SAMPLING_VERSION="v1.64.1"
# What version of vg should be used to haplotype-sample graphs when we want to keep just one reference?
VG_HAPLOTYPE_SAMPLING_ONEREF_VERSION="v1.64.1"

wildcard_constraints:
    expname="[^/]+",
    refgraphbase="[^/]+?",
    reference="chm13|grch38|chm13v1",
    # We can have multiple versions of graphs with different modifications and clipping regimes
    modifications="(-[^.-]+(\\.trimmed)?(\\.clip\\.[0-9]*\\.[0-9]*)?(\\.ec[0-9kmKM]+)?)*",
    clipping="\\.d[0-9]+|",
    full="\\.full|",
    chopping="\\.unchopped|",
    # After a full graph, we might have sampling stuff in a graph name
    sampling="(-[^-]+)*",
    trimmedness="\\.trimmed|",
    sample=".+(?<!\\.trimmed)",
    basename=".+(?<!\\.trimmed)",
    subset="([0-9]+[km]?|full)",
    category="((not_)?(centromeric))?|",
    # We can restrict calling to a small region for testing
    region="(|chr21)",
    # We use this for an optional separating dot, so we can leave it out if we also leave the field empty
    dot="\\.?",
    callparams="(|(\\.model[0-9-]+|\\.nomodel)(\\.legacy)?)",
    tech="[a-zA-Z0-9]+",
    statname="[a-zA-Z0-9_]+(?<!compared)(?<!sv_summary)(?<!mapping_stats_real)(?<!mapping_stats_sim)(.mean|.total)?",
    statnamex="[a-zA-Z0-9_]+(?<!compared)(?<!sv_summary)(?<!mapping_stats_real)(?<!mapping_stats_sim)(.mean|.total)?",
    statnamey="[a-zA-Z0-9_]+(?<!compared)(?<!sv_summary)(?<!mapping_stats_real)(?<!mapping_stats_sim)(.mean|.total)?",
    realness="(real|sim)",
    realnessx="(real|sim)",
    realnessy="(real|sim)",

def auto_mapping_threads(wildcards):
    """
    Choose the number of threads to use map reads, from subset.
    """
    number = subset_to_number(wildcards["subset"])
    mapping_threads = 0
    if number >= 100000:
        mapping_threads= MAPPER_THREADS
    elif number >= 10000:
        mapping_threads= 16
    else:
        mapping_threads= 8

    if wildcards.get("mapper", "").startswith("graphaligner") and wildcards.get("realness", "") == "sim":
        #Graphaligner is really slow so for simulated reads where we don't care about time
        #double the number of threads
        #At most 128 because it errors with too many threads sometimes
        return min(mapping_threads * 2, 128) 
    else:
        return mapping_threads

def auto_mapping_slurm_extra(wildcards):
    """
    Determine Slurm extra arguments for a timed, real-read mapping job.
    """
    if EXCLUSIVE_TIMING:
        return "--exclusive " + REAL_SLURM_EXTRA
    else:
        return REAL_SLURM_EXTRA

def auto_mapping_full_cluster_nodes(wildcards):
    """
    Determine number of full cluster nodes for a timed, real-read mapping job.

    TODO: Is this really used by Slurm?
    """
    number = subset_to_number(wildcards["subset"])
    if EXCLUSIVE_TIMING:
        return 1
    else:
        return 0

def auto_mapping_memory(wildcards):
    """
    Determine the memory to use for Giraffe mapping, in MB, from subset and tech.
    """
    thread_count = auto_mapping_threads(wildcards)

    base_mb = 60000

    if wildcards["tech"] == "illumina":
        scale_mb = 200000
    elif wildcards["tech"] == "hifi":
        scale_mb = 240000
    else:
        scale_mb = 600000

    # Scale down memory with threads
    return scale_mb / 64 * thread_count + base_mb

def choose_partition(minutes):
    """
    Get a Slurm partition that can fit a job running for the given number of
    minutes, or raise an error.
    """
    for name, limit in SLURM_PARTITIONS:
        if minutes <= limit:
            return name
    raise ValueError(f"No Slurm partition accepts jobs that run for {minutes} minutes")

import snakemake
if int(snakemake.__version__.split(".")[0]) >= 8:
    # Remote providers have been replaced with storage plugins.

    # TODO: test this on Snakemake 8
    # TODO: Really depend on snakemake-storage-plugin-http

    def remote_or_local(url):
        """
        Wrap a URL as a Snakemake "remote file", but pass a local path through.
        """

        if url.startswith("https://"):
            # Looks like a remote.
            return storage.http(url)
        else:
            return url

    def to_local(possibly_remote_file):
        """
        Given a result of remote_or_local, turn it into a local filesystem path.

        Snakemake must have already downloaded it for us if needed.
        """

        # TODO: Do we still have the list problem in Python code with storage
        # plugins?

        if isinstance(possibly_remote_file, list):
            return possibly_remote_file[0]
        else:
            return possibly_remote_file

else:
    # This is for Snakemake 7 and below. Snakemake 8 replaces Snakemake 7
    # remote providers with storage plugins, which aren't available in
    # Snakemake 7.

    # Also, Snakemake remote providers are terrible, and can fail a successful job with something like:
    # Error recording metadata for finished job (HTTPSConnectionPool(host='ftp-trace.ncbi.nlm.nih.gov', port=443): Max retries exceeded with url: /ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed (Caused by SSLError(SSLEOFError(8, '[SSL: UNEXPECTED_EOF_WHILE_READING] EOF occurred in violation of protocol (_ssl.c:1007)')))). Please ensure write permissions for the directory /private/home/anovak/workspace/lr-giraffe/.snakemake
    # Because for some reason they want to run web requests at the end of a job.

    # So we just use a magic local "remote_cache_https" directory

    rule remote_cache_https:
        output:
            "remote_cache_https/{host}/{path}"
        wildcard_constraints:
            host="[^/]*"
        threads: 1
        resources:
            mem_mb=2000,
            runtime=30,
            slurm_partition=choose_partition(30)
        shell:
            "wget https://{wildcards.host}/{wildcards.path} -O {output}"

    def remote_or_local(url):
        """
        Wrap a URL as a usable Snakemake input, but pass a local path through.
        """

        if url.startswith("https://"):
            # Looks like a remote.
            return  "remote_cache_https/" + url[8:]
        else:
            return url

    def to_local(possibly_remote_file):
        """
        Given a result of remote_or_local, turn it into a local filesystem path.

        Snakemake must have already downloaded it for us if needed.
        """

        return possibly_remote_file

def subset_to_number(subset):
    """
    Take a subset like 1m or full and turn it into a number.
    """
    if subset == "full":
        return float("inf")
    elif subset.endswith("m"):
        multiplier = 1000000
        subset = subset[:-1]
    elif subset.endswith("k"):
        multiplier = 1000
        subset = subset[:-1]
    else:
        multiplier = 1

    return int(subset) * multiplier

def repetitive_kmers(wildcards):
    """
    Find the Winnowmap repetitive kmers file from a reference.
    """
    return reference_basename(wildcards) + ".repetitive_k15.txt"

def minimap_derivative_mode(wildcards):
    """
    Determine the right Minimap2/Winnowmap preset (map-pb, etc.) from minimapmode or tech.
    """
    explicit_mode = wildcards.get("minimapmode", None)
    if explicit_mode is not None:
        return explicit_mode

    MODE_BY_TECH = {
        "r9": "map-ont",
        "r10": "map-ont",
        "r10y2025": "map-ont",
        "hifi": "map-pb",
        "illumina": "sr" # Only Minimap2 has this one, Winnowmap doesn't.
    }

    return MODE_BY_TECH[wildcards["tech"]]

def minimap2_index(wildcards):
    """
    Find the minimap2 index from reference and tech.
    """
    
    mode_part = minimap_derivative_mode(wildcards)
    return reference_basename(wildcards) + "." + mode_part + ".mmi"

def pbmm2_hifi_index(wildcards):
    """
    Find the index for pbmm2's HiFI preset, from reference and tech.

    This will always be the same as minimap2's map-hifi preset index, because
    they use the same k, w, and holopolymer compression. It will NOT be the
    same as minimap2's map-pb preset index, which uses different settings.
    """
    return reference_basename(wildcards) + ".map-hifi.mmi"

def reference_basename(wildcards):
    """
    Find the linear reference base name without extension from a reference.

    This reference is used for mapping. Calling may be against a different
    calling reference (with possibly a different Y) or renamed contigs.

    This reference will use PanSN contig names.
    """
    parts = [wildcards["reference"], "pansn"]
    if wildcards["reference"] == "chm13":
        # We want to use a version of the reference FASTA with the "new"
        # non-HG002, non-GRCh38 Y contig.
        parts.append("newY")
    elif wildcards["reference"] == "chm13v1":
        # Use this for the older version, so just chm13 without the v1 and without the newY
        parts[0] = "chm13"
    return os.path.join(REFS_DIR, "-".join(parts))

def reference_fasta(wildcards):
    """
    Find the linear reference FASTA from a reference.
    """
    return reference_basename(wildcards) + ".fa"

def reference_dict(wildcards):
    """
    Find the linear reference FASTA dictionary from a reference.
    """
    return reference_fasta(wildcards) + ".dict"

def reference_path_list_callable(wildcards):
    """
    Find the path list file for a linear reference that we can actually call on, from reference and region.
   
    We "can't" call on chrY for CHM13 because the one we use in the graphs
    isn't the same as the one in CHM13v2.0 where the calling happens.
    """
    return reference_fasta(wildcards) + ".paths" + wildcards.get("region", "") + ".callable.txt"

def reference_path_dict_callable(wildcards):
    """
    Find the path dict file for a linear reference that we can actually call on, from reference and region.
   
    We need this because when surjecting to a non-base reference we can't infer path lengths from the graph.
    """
    return reference_fasta(wildcards) + ".paths" + wildcards.get("region", "") + ".callable.dict"

def bwa_index_set(wildcards):
    """
    Find a reference indexed for BWA from reference.

    Doesn't actually include the FASTA itself.
    """
    fasta_name = reference_fasta(wildcards)
    return {ext: fasta_name + "." + ext for ext in ("amb", "ann", "bwt", "pac", "sa")}

def reference_prefix(wildcards):
    """
    Find the PanSN prefix we need to remove to convert form PanSN names to
    non-PanSN names, from reference.
    """
    return {
        "chm13": "CHM13#0#",
        "chm13v1": "CHM13#0#",
        "grch38": "GRCh38#0#"
    }[wildcards["reference"]]

def reference_sample(wildcards):
    """
    Get the reference as a sample name
    """
    return {
        "chm13": "CHM13",
        "chm13v1": "CHM13",
        "grch38": "GRCh38"
    }[wildcards["truthref"]]

def calling_reference_fasta(wildcards):
    """
    Find the linear reference FASTA with non-PanSN names from a reference (for
    interpreting VCFs).

    For CHM13, we always use CHM13v2.0 as the calling reference since that's
    the one we can get a truth on.
    """
    match wildcards["reference"]:
        case "chm13":
            return os.path.join(REFS_DIR, "chm13v2.0.fa")
        case "chm13v1":
            return os.path.join(REFS_DIR, "chm13v2.0.fa")
        case reference:
            return os.path.join(REFS_DIR, reference + ".fa")

def calling_reference_fasta_index(wildcards):
    """
    Find the index for the linear calling (non-PanSN) reference, from reference.
    """
    return calling_reference_fasta(wildcards) + ".fai"

def calling_reference_restrict_bed(wildcards):
    """
    Find the BED for the linear calling (non-PanSN) reference region we think we can call on, from reference.
    """
    return calling_reference_fasta(wildcards) + ".callable.from." + wildcards["reference"] + ".bed"

def calling_reference_par_bed(wildcards):
    """
    Find the BED for the psuedo-autosomal regions of a reference, from reference.
    """
    return os.path.splitext(calling_reference_fasta(wildcards))[0] + "_PAR.bed"
def sv_calling_vntr_bed(wildcards):
    """
    Find the vntr bed file for sniffles
    """

    match wildcards["reference"]:
        case "chm13":
            return os.path.join(SV_DATA_DIR, "human_chm13v2.0_maskedY_rCRS.trf.bed")
        case "chm13v1":
            return os.path.join(SV_DATA_DIR, "human_chm13v2.0_maskedY_rCRS.trf.bed")
        case reference:
            return os.path.join(SV_DATA_DIR, "human_GRCh38_no_alt_analysis_set.trf.bed")


def uncallable_contig_regex(wildcards):
    """
    Get a grep regex matching a substring in all uncallable contigs in the
    calling or PanSN reference, from reference or basename.

    We expect basename to start with a reference identifier like "chm13", maybe
    with other stuff after a "-".

    This will be a regex compatible with non-E grep.
    """

    wc_dict = dict(wildcards)
    if "reference" in wc_dict:
        reference = wc_dict["reference"]
    else:
        reference = wc_dict["basename"].split("-")[0]
    match reference:
        case "chm13":
            # TODO: We don't want to try and call on Y on CHM13 because it's
            # not the same Y as CHM13v2.0, where the truth set is and where we
            # will do the calling.
            return "chr[YM]"
        case "chm13v1":
            # TODO: We don't want to try and call on Y on CHM13 because it's
            # not the same Y as CHM13v2.0, where the truth set is and where we
            # will do the calling.
            return "chr[YM]"
        case "grch38":
            # We don't want to try and call on the _random/Un/EBV contigs that
            # probably aren't in the graph at all.
            return "chrM\\|.*_random\\|chrUn\\|chrEBV"
        case _:
            return "chrM"

def haploid_contigs(wildcards):
    """
    Get a list of all haploid contigs in a sample, with no prefixes, from sample.

    Will be either [] or ["chrX", "chrY"].
    """
    XX_SAMPLES = {"HG001"}
    XY_SAMPLES = {"HG002"}
    if wildcards["sample"] in XX_SAMPLES:
        return []
    elif wildcards["sample"] in XY_SAMPLES:
        return ["chrX", "chrY"]
    else:
        raise RuntimeError(f"Unknown karyotype for {wildcards['sample']}")

def wdl_cache(wildcards):
    """
    Get a WDL workflow step cache directory path from root.
    """
    return os.path.join(os.path.abspath("."), wildcards["root"], "miniwdl-cache")

def model_files(wildcards):
    """
    Get the DV model files to use, or None, from tech, mapper, and callparams.
    """

    if wildcards.tech in ("hifi", "r10y2025") and "giraffe" in wildcards.mapper:
        if wildcards.callparams.startswith(".model"):
            # Grab the date; we assume we always use dates as model names.
            model_name = wildcards.callparams[len(".model"):len(".model") + 10]
        elif wildcards.callparams.startswith(".nomodel"):
            # Just use the model as shipped in DV
            return None
        else:
            if wildcards.tech == "hifi":
                # Use a particular trained model as the default for HiFi.
                model_name = "2025-03-26"
            else:
                # For r10, keep defaulting to no trained model for now.
                return None
        model_base_path = os.path.join(MODELS_DIR, wildcards.tech, model_name)
        # Get all the files as absolute paths in sorted order, minus any README.
        # We want to be the same as the manual order so we don't get WDL cache misses.
        model_files = sorted([os.path.join(model_base_path, f) for f in os.listdir(model_base_path) if f != "README.txt"])
        return model_files
    return None

def truth_vcf_url(wildcards):
    """
    Find the URL or local file for the variant calling truth VCF, from reference and sample.
    """

    if wildcards["sample"] == "HG002":
        # These are available online directly
        return  {
            "chm13": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz",
            "chm13v1": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz",
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
        }[wildcards["reference"]]
    elif wildcards["sample"] == "HG001":
        return  {
            # On CHM13 we don't have a real benchmark set, so we have to use the raw Platinum Pedigree dipcall calls.
            "chm13": os.path.join(TRUTH_DIR, wildcards["reference"], wildcards["sample"], wildcards["sample"] + ".dip.vcf.gz"),
            # On CHM13 we don't have a real benchmark set, so we have to use the raw Platinum Pedigree dipcall calls.
            "chm13v1": os.path.join(TRUTH_DIR, wildcards["reference"], wildcards["sample"], wildcards["sample"] + ".dip.vcf.gz"),
            # On GRCh38 there's a GIAB truth set
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
        }[wildcards["reference"]]
    else:
        raise RuntimeError("Unsupported sample: " + wildcards["sample"])

def truth_vcf_index_url(wildcards):
    """
    Find the URL or local file for the variant calling truth VCF index, from reference and sample.
    """
    return truth_vcf_url(wildcards) + ".tbi"

def truth_bed_url(wildcards):
    """
    Find the URL or local file for the variant calling truth high confidence BED, from reference.

    If compressed, must end in ".gz".
    """

    if wildcards["sample"] == "HG002":
        # For HG002, these are available online directly.
        return {
            "chm13": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            "chm13v1": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
        }[wildcards["reference"]]
    elif wildcards["sample"] == "HG001":
        return  {
            # TODO: On CHM13 we don't have Platinum Pedigree high-confidence regions, so
            # we need to just use the HG002 ones for other samples and hope they're close enough.
            "chm13": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            # TODO: On CHM13 we don't have Platinum Pedigree high-confidence regions, so
            # we need to just use the HG002 ones for other samples and hope they're close enough.
            "chm13v1": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed",
            # On GRCh38 there's a GIAB truth set
            "grch38": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
        }[wildcards["reference"]]
    else:
        raise RuntimeError("Unsupported sample: " + wildcards["sample"])

def surjectable_gam(wildcards):
    """
    Find a GAM mapped to a graph built on a reference that we can use to
    surject to the reference we're interested in.
    """
    # TODO: Make similar redirects for the benchmarks!

    format_data = dict(wildcards)

    return "{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam".format(**format_data)

def graph_base(wildcards):
    """
    Find the base name for a collection of graph files from reference and either refgraph or all of refgraphbase, modifications, and clipping.

    For graphs ending in "-sampled" and sampling parameters, autodetects the right haplotype-sampled full graph name to use from tech and sample.

    For GraphAligner, selects an unchopped version of the graph based on mapper.
    """

    # For membership testing, we need a set of wildcard keys
    wc_keys = set(wildcards.keys())
    reference = wildcards["reference"]
    # TODO: This sucks
    if reference == "chm13v1":
        reference = "chm13"
    modifications = []

    if "refgraphbase" in wc_keys and "modifications" in wc_keys:
        # We already have the reference graph base (hprc-v1.1-mc) and clipping
        # (.d9) and chopping (.unchopped) and full (.full) if allowed cut
        # apart. Sampling info may also be there (we assume it already
        # describes what to sample for specifically).
        refgraphbase = wildcards["refgraphbase"]
        modifications.append(wildcards["modifications"])
        if "clipping" in wc_keys:
            modifications.append(wildcards["clipping"])
        if "full" in wc_keys:
            modifications.append(wildcards["full"])
        if "chopping" in wc_keys:
            modifications.append(wildcards["chopping"])
        if "sampling" in wc_keys:
            modifications.append(wildcards["sampling"])
    else:
        assert "refgraph" in wc_keys, f"No refgraph wildcard in: {wc_keys}"
        # We need to handle hprc-v1.1-mc and hprc-v1.1-mc.full and hprc-v1.1-mc-d9 and hprc-v2.prereease-mc-R2-d32 and hprc-v2.prereease-mc-R2-sampled10d.
        # Also probaby primary.
        # They need the reference inserted after the -mc. but before the other stuff, and -d32 needs to become .d32.
        # And -sampled10d needs to be expanded to say what sample and read set we are haplotype sampling from.

        # TODO: If it's not -mc, where would the reference go?
        
        refgraph = wildcards["refgraph"]

        parts = refgraph.split("-")
        last = parts[-1]
        if re.fullmatch("unchopped", last):
            # We have a choppedness modifier, which gets a dot.
            modifications.append("." + last)
            parts.pop()

        last = parts[-1]
        if re.fullmatch("sampled[0-9]+d?o?", last):
            # We have a generic haplotype sampling flag.
            # Autodetect the right haplotype-sampled graph to use.

            # We need to convert this into the graph sampled for the right sample.
            # Which means we need some info about the sample we are working on.
            assert "tech" in wc_keys, "No tech known for haplotype sampling graph " + wildcards["refgraph"]
            assert "sample" in wc_keys, "No sample known for haplotype sampling graph " + wildcards["refgraph"]
           
            # TODO: We always sample for the full real trimmed-if-R10 version
            # of whatever reads we're going to map, so we can consistently use
            # one graph. Note r10y2025 doesn't get trimmed.
            sampling_trimmedness = ".trimmed" if wildcards["tech"] == "r10" else ""
            modifications.append(f"-{last}_fragmentlinked-for-real-{wildcards['tech']}-{wildcards['sample']}{sampling_trimmedness}-full")
            parts.pop()
        elif re.fullmatch("d[0-9]+", last):
            # We have a clipping modifier, which gets a dot.
            modifications.append("." + last)
            parts.pop()
        
        while len(parts) > 3:
            # We have more than just the 3-tuple of name, version, algorithm. Take the last thing into modifications.
            modifications.append("-" + parts[-1])
            parts.pop()
        
        last = parts[-1]
        if last.endswith(".full"):
            # Move a .full over the reference
            modifications.append(".full")
            last = last[:-5]
            parts[-1] = last
        

        # Now we have all the modifications. Flip them around the right way.
        modifications.reverse()
        
        # The first 3 or fewer parts are the graph base name.
        refgraphbase = "-".join(parts)

    result = os.path.join(GRAPHS_DIR, refgraphbase + "-" + reference + "".join(modifications))
    return result

def gbz(wildcards):
    """
    Find a graph GBZ file from reference.
    """
    return graph_base(wildcards) + ".gbz"

def hg(wildcards):
    """
    Find a graph hg file from reference.
    """
    return graph_base(wildcards) + ".hg"

def gfa(wildcards):
    """
    Find a graph GFA file from reference.
    """
    return graph_base(wildcards) + ".gfa"

def snarls(wildcards):
    """
    Find a graph snarls file from reference.
    """
    return graph_base(wildcards) + ".snarls"

def minimizer_k(wildcards):
    """
    Find the minimizer kmer size from mapper.
    """
    if wildcards["mapper"].startswith("giraffe"):
        # Looks like "giraffe-k31.w50.W-lr-default-noflags".
        # So get second part on - and first part of that on . and number-ify it after the k.
        return int(wildcards["mapper"].split("-")[1].split(".")[0][1:])
    else:
        mode = minimap_derivative_mode(wildcards)
        match mode:
            # See minimap2 man page
            case "map-ont":
                return 15
            case "map-pb":
                return 19
            case "sr":
                return 21
        raise RuntimeError("Unimplemented mode: " + mode)

def dist_indexed_graph(wildcards):
    """
    Find a GBZ and its dist index from reference.
    """
    base = graph_base(wildcards)
    return {
        "gbz": gbz(wildcards),
        "dist": base + ".dist"
    }

def indexed_graph(wildcards):
    """
    Find an indexed graph and all its indexes from reference and minparams.

    Also checks vgversion to see if we need a no-zipcodes index for old vg.
    """
    base = graph_base(wildcards)
    indexes = dist_indexed_graph(wildcards)
    # Some versions of Giraffe can't use zipcodes. But all the Giraffe versions
    # we test can use the same distance indexes.
    if "vgversion" in dict(wildcards).keys() and wildcards["vgversion"] in NON_ZIPCODE_GIRAFFE_VERSIONS:
        # Don't use zipcodes
        new_indexes = {
            "minfile": base + "." + wildcards["minparams"] + ".nozip.min",
        }
    else:
        # Use zipcodes
        new_indexes = {
            "minfile": base + "." + wildcards["minparams"] + ".withzip.min",
            "zipfile": base + "." + wildcards["minparams"] + ".zipcodes"
        }
    new_indexes.update(indexes)
    return new_indexes

def r_and_snarl_indexed_graph(wildcards):
    """
    Find a GBZ and its snarl/distance and ri indexes.

    Uses the full distance index if present or the snarls-only one otherwise.
    """
    base = graph_base(wildcards)
    return {
        "gbz": gbz(wildcards),
        "ri": base + ".ri",
        "snarls": base + ".dist" if os.path.exists(base + ".dist") else base + ".tcdist"
    }

def haplotype_indexed_graph(wildcards):
    """
    Find a GBZ and its hapl index.

    Distance and ri indexes are not needed for haplotype sampling.
    """
    base = graph_base(wildcards)
    return {
        "gbz": gbz(wildcards),
        "hapl": base + ".fragment.hapl"
    }
    return result

def base_fastq_gz(wildcards):
    """
    Find a full compressed FASTQ for a real sample, based on realness, sample,
    tech, and trimmedness.

    If an untrimmed version exists and the trimmed version does not, returns
    the name of the trimmed version to make.
    """
    import glob
    full_gz_pattern = os.path.join(READS_DIR, "{realness}/{tech}/{sample}/*{sample}*{trimmedness}.f*q.gz".format(**wildcards))
    results = glob.glob(full_gz_pattern)
    if wildcards["trimmedness"] != ".trimmed":
        # Don't match trimmed files when not trimmed.
        results = [r for r in results if ".trimmed" not in r]
    # Skip any subset files
    subset_regex = re.compile("([0-9]+[km]?|full)")
    results = [r for r in results if not subset_regex.fullmatch(r.split(".")[-3])]
    # Skip any chunk files
    chunk_regex = re.compile("chunk[0-9]*")
    results = [r for r in results if not chunk_regex.fullmatch(r.split(".")[-3])]
    if len(results) == 0:
        # Can't find it
        if wildcards["trimmedness"] == ".trimmed":
            # Look for an untrimmed version
            untrimmed_pattern = os.path.join(READS_DIR, "{realness}/{tech}/{sample}/*{sample}*.f*q.gz".format(**wildcards))
            results = glob.glob(untrimmed_pattern)
            # Skip any subset files
            results = [r for r in results if not subset_regex.fullmatch(r.split(".")[-3])]
            if len(results) == 1:
                # We can use this untrimmed one to make a trimmed one
                without_gz = os.path.splitext(results[0])[0]
                without_fq, fq_ext = os.path.splitext(without_gz)
                trimmed_base = without_fq + ".trimmed" + fq_ext + ".gz"
                return trimmed_base
        raise FileNotFoundError(f"No  non-chunk, non-subset files found matching {full_gz_pattern}")
    elif len(results) > 1:
        raise RuntimeError("Multiple non-chunk, non-subset files matched " + full_gz_pattern)
    return results[0]

def fastq_finder(wildcards, compressed = False):
    """
    Find a FASTQ or compressed FASTQ from realness, tech, sample, trimmedness, and subset.

    Works even if there is extra stuff in the name besides sample. Accounts for
    being able to make a FASTQ from a GAM.
    """

    extra_ext = ".gz" if compressed else ""

    #If we chunked the reads, we need to add it. Otherwise there isn't a chunk
    readchunk = wildcards.get("readchunk", "")

    import glob
    fastq_by_sample_pattern = os.path.join(READS_DIR, ("{realness}/{tech}/{sample}/*{sample}*{trimmedness}[._-]{subset}" + readchunk+ ".f*q" + extra_ext).format(**wildcards))
    results = glob.glob(fastq_by_sample_pattern)
    if wildcards["trimmedness"] != ".trimmed":
        # Don't match trimmed files when not trimmed.
        results = [r for r in results if ".trimmed" not in r]
    if len(results) == 0:
        if wildcards["realness"] == "real":
            # Make sure there's a full .fq.gz to extract from (i.e. this doesn't raise)
            full_file = base_fastq_gz(wildcards)
            # And compute the subset name
            without_gz = os.path.splitext(full_file)[0]
            without_fq = os.path.splitext(without_gz)[0]
            return without_fq + (".{subset}" + readchunk+ ".fq" + extra_ext).format(**wildcards)
        elif wildcards["realness"] == "sim":
            # Assume we can get this FASTQ.
            # For simulated reads we assume the right subset GAM is there. We
            # don't want to deal with the 1k/1000 difference here.
            return os.path.join(READS_DIR, ("{realness}/{tech}/{sample}/{sample}-{realness}-{tech}{trimmedness}-{subset}" + readchunk+ ".fq" + extra_ext).format(**wildcards))
        else:
            raise FileNotFoundError(f"No files found matching {fastq_by_sample_pattern}")
    elif len(results) > 1:
        raise AmbiguousRuleException("Multiple files matched " + fastq_by_sample_pattern)
    return results[0]

def fastq(wildcards):
    """
    Find a FASTQ from realness, tech, sample, trimmedness, and subset.

    Works even if there is extra stuff in the name besides sample. Accounts for
    being able to make a FASTQ from a GAM.
    """
    return fastq_finder(wildcards, compressed=False)

def fastq_gz(wildcards):
    """
    Find a compressed FASTQ from realness, tech, sample, trimmedness, and subset.

    Works even if there is extra stuff in the name besides sample. Accounts for
    being able to make a FASTQ from a GAM.
    """
    return fastq_finder(wildcards, compressed=True)

def kmer_counts(wildcards):
    """
    Find a compressed FASTQ from realness, tech, sample, trimmedness, and subset, and get the corresponding kmer file name.
    """
    result = re.sub("\\.f.*?q.gz$", ".kff", fastq_gz(wildcards))
    assert result.endswith(".kff")
    return result

def kmer_counts_benchmark(wildcards):
    return re.sub(".kff", ".kmer_counting.benchmark", kmer_counts(wildcards))

def mapper_stages(wildcards):
    """
    Find the list of mapping stages from mapper.
    """

    if wildcards["mapper"].startswith("giraffe"):
        parts = wildcards["mapper"].split("-")
        assert len(parts) >= 3
        # Should be giraffe, then the minimizer parameters, then the parameter preset.
        if parts[3] == "default":
            # Default mapping preset is non-chaining
            return NON_CHAINING_STAGES
        else:
            return STAGES
    else:
        return []


def all_experiment_conditions(expname, filter_function=None, debug=False):
    """
    Yield dictionaries of all conditions for the given experiment.
    
    The config file should have a dict in "experiments", of which the given
    expname should be a key. The value is the experiment dict.

    The experiment dict should have a "control" dict, listing names and values
    of variables to keep constant.

    The experiment dict should have a "vary" dict, listing names and values
    lists of variables to vary. All combinations will be generated.

    The experiment dict should have a "constrain" list. Each item in the list
    is a "pass", which is a list of constraints. Each item in the pass is a
    dict of variable names and values (or lists of values). A condition must
    match *at least* one of these dicts on *all* values in the dict in order to
    survive the pass. And it must survive all passes in order to be run.

    If optional variables (defined internally) don't appear in the experiment
    "vary" or "control" sections, they will be given an empty-string value.

    If filter_function is provided, only yields conditions that the filter
    function is true for.

    Yields variable name to value dicts for all passing conditions for the
    given experiment.
    """

    if "experiments" not in config:
        raise RuntimeError(f"No experiments section in configuration; cannot run experiment {expname}")
    all_experiments = config["experiments"]
    
    if expname not in all_experiments:
        raise RuntimeError(f"Experiment {expname} not in configuration")
    exp_dict = all_experiments[expname]

    # Make a base dict of all controlled variables.
    base_condition = exp_dict.get("control", {})

    to_vary = exp_dict.get("vary", {})

    constraint_passes = exp_dict.get("constrain", [])

    total_conditions = 0
    for condition in augmented_with_all(base_condition, to_vary):
        # For each combination of independent variables on top of the base condition

        # Fill in optional experiment variables with empty strings
        for optional_key in ["trimmedness", "callparams"]:
            if optional_key not in condition:
                condition[optional_key] = ""

        # We need to see if this is a combination we want to do
        if matches_all_constraint_passes(condition, constraint_passes):
            if not filter_function or filter_function(condition):
                total_conditions += 1
                yield condition
            else:
                if debug:
                    print(f"Condition {condition} does not match requested filter function")
        else:
            if debug:
                print(f"Condition {condition} does not match a constraint in some pass")
    print(f"Experiment {expname} has {total_conditions} eligible conditions")
    

def augmented_with_each(base_dict, new_key, possible_values):
    """
    Yield copies of base_dict with each value from possible_values under new_key.
    """

    for value in sorted(possible_values):
        clone = dict(base_dict)
        clone[new_key] = value
        yield clone

def augmented_with_all(base_dict, keys_and_values):
    """
    Yield copies of base_dict augmented with all combinations of values from
    keys_and_values, under the corresponding keys.
    """

    if len(keys_and_values) == 0:
        # Base case: nothing to add
        yield base_dict
    else:
        # Break off one facet
        first_key = next(iter(keys_and_values.keys()))
        first_values = keys_and_values[first_key]
        rest = dict(keys_and_values)
        del rest[first_key]

        for with_first in augmented_with_each(base_dict, first_key, first_values):
            # Augment with this key
            for with_rest in augmented_with_all(with_first, rest):
                # And augment with the rest
                yield with_rest

def matches_constraint_value(query, value):
    """
    Returns True if query equals value, except if value is a list, query has to
    be in the list instead.
    """

    if isinstance(value, list):
        return query in value
    else:
        return query == value

def matches_constraint(condition, constraint, debug=False):
    """
    Returns True if all keys in constraint are in condition with the same
    values, or with values in the list in constraint.
    """
    for key, match in constraint.items():
        if key not in condition or not matches_constraint_value(condition[key], match):
            if debug:
                print(f"Condition {condition} mismatched constraint {constraint} on {key}")
            return False
    return True

def matches_any_constraint(condition, constraints):
    """
    Return True if, for some constraint dict, the condition dict matches all
    values in the constraint dict.
    """

    for constraint in constraints:
        if matches_constraint(condition, constraint):
            return True
    return False

def matches_all_constraint_passes(condition, passes):
    """
    Return True if the condfition matches some constraint in each pass in passes.
    """
    
    if len(passes) > 0 and not isinstance(passes[0], list) and isinstance(passes[0], dict):
        # Old style config where there's just one pass of constraints. Fix it up.
        passes = [passes]

    for constraints in passes:
        if not matches_any_constraint(condition, constraints):
            return False
    return True

def wildcards_to_condition(all_wildcards):
    """
    Filter down wildcards to just the condition parameters for the experiment in expname.
    
    Raises an error if any varied variable in the experiment cannot be determined.
    """

    exp_dict = config.get("experiments", {}).get(all_wildcards["expname"], {})
    base_condition = exp_dict.get("control", {})
    to_vary = exp_dict.get("vary", {})
    
    # For membership testing, we need a set of wildcard keys
    wc_keys = set(all_wildcards.keys())

    condition = {}

    for var in base_condition.keys():
        if var in wc_keys:
            # Constant variables across the whole experiment don't need to
            # actually be in the condition name unless we actually have them
            # available.
            condition[var] = all_wildcards[var]

    for var in to_vary.keys():
        if len(all_wildcards.get(var, "")) != 0:
            condition[var] = all_wildcards[var]
        elif var == "mapper":
            #If we haven't defined a mapper, but we have the components of a giraffe mapper
            minimizer_params = all_wildcards.get("minparams", "k"+all_wildcards["k"]+".w"+all_wildcards["w"]+all_wildcards["weightedness"])
            condition[var] = "giraffe-" + minimizer_params + "-" + all_wildcards["preset"] + "-" + all_wildcards["vgversion"] + "-" + all_wildcards["vgflag"]
        elif var == "refgraph" and len(all_wildcards.get("refgraphbase", "")) != 0:
            #If we haven't defined refgraph, but we have the components
            #TODO: This isn't great but it works in the one case I need it to
            graph = all_wildcards.get("refgraphbase")
            if len(all_wildcards.get("sampling", "")) != 0:
                graph += "-"
                graph += all_wildcards.get("sampling")
            graph += all_wildcards.get("clipping","")
            graph += all_wildcards.get("full","")
            graph += all_wildcards.get("chopping","")

            condition[var] = graph
        else:
            #Catch any case where it fails
            condition[var] = all_wildcards[var]



    return condition

def condition_name(wildcards):
    """
    Determine a human-readable condition name from expname and the experiment's variable values.
    """
 
    # Get what changes in the experiment
    exp_dict = config.get("experiments", {}).get(wildcards["expname"], {})
    to_vary = exp_dict.get("vary", {})

    # Get the condition dict in use here
    condition = wildcards_to_condition(wildcards)

    name_parts = []

    for varied_key in to_vary:

        #Don't include realness
        if varied_key == "realness":
            continue

        # Look at the value we have for this varied variable
        condition_value = condition[varied_key]
        # And the other possible values that are used
        alternatives = set(to_vary.get(varied_key, []))

        if varied_key == "mapper" and "giraffe" in condition_value:
            # If we're working on a Giraffe mapper name, only compare against other Giraffe mapper names
            alternatives = {a for a in alternatives if "giraffe" in a}

        # Find all the name parts used in alternatives.
        parts_in_alternatives = [set(alternative.split("-")) for alternative in alternatives]
        # Find those used in all alternatives
        universal_alternative_parts = functools.reduce(lambda a, b: a & b, parts_in_alternatives)

        # Find all the name parts used in us
        condition_value_parts = condition_value.split("-")
       
        # Drop parts that are universal among applicable alternatives, keeping
        # only the parts that represent differences.
        interesting_parts = [p for p in condition_value_parts if p not in universal_alternative_parts]
        
        if varied_key == "mapper" and "giraffe" in condition_value:
            # Make sure to mark Giraffe conditions even when all of them have "giraffe" in them.
            interesting_parts = ["giraffe"] + interesting_parts

        # And add that to the name
        name_parts.append("-".join(interesting_parts))

    return ",".join(name_parts)

def all_experiment(wildcard_values, pattern, filter_function=None, empty_ok=False, debug=False):
    """
    Produce all values of pattern substituted with the wildcards and the experiment conditions' values, from expname.

    If provided, restricts to conditions passing the filter function.

    Throws an error if nothing is produced and empty_ok is not set.

    Needs to be used like:
        lambda w: all_experiment(w, "your pattern")
    """

    empty = True
    for condition in all_experiment_conditions(wildcard_values["expname"], filter_function=filter_function):
        merged = dict(wildcard_values)
        merged.update(condition)

        # TODO: Hackily fill in the optional {dot} which should be set if {category} is set.
        # There's no way in a Snakemake expandion to tack on a leader sequence.
        if "category" in merged and len(merged["category"]) > 0:
            merged["dot"] = "."
        else:
            merged["dot"] = ""
        if "category" not in merged:
            # Category might appear in templates but in experiments it is meant to be optional.
            merged["category"] = ""

        if debug:
            print(f"Evaluate {pattern} in {merged} from {wildcard_values} and {condition}")
        filename = pattern.format(**merged)
        yield filename
        empty = False
    if empty:
        if debug:
            print("Produced no values for " + pattern + " in experiment!")
        if not empty_ok:
            if debug:
                print("Failing!")
            raise RuntimeError("Produced no values for " + pattern + " in experiment!")

def has_stat_filter(stat_name):
    """
    Produce a filter function for conditions that might have the stat stat_name.

    Applies to stat files like:
    {root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv

    Use with all_experiment() when aggregating stats to compare, to avoid
    trying to agregate from conditions for which the stat cannot be measured.
    """

    def filter_function(condition):
        """
        Return True if the given condition dict should have the stat named stat_name.
        """

        if stat_name in {"correct", "accuracy", "wrong"}:
            # These stats only exist for conditions with a truth set (i.e. simulated ones)
            if condition["realness"] != "sim":
                return False

        if stat_name.startswith("indel_") or stat_name.startswith("snp_"):
            # These are calling stats, and we shoudl only do calling for real reads.
            if condition["realness"] != "real":
                return False

        if stat_name.startswith("time_used") or stat_name in ("mapping_speed", "chain_coverage"):
            # This is a Giraffe time used stat or mean thereof. We need to be a
            # Giraffe condition.
            if not condition["mapper"].startswith("giraffe"):
                return False

        return True

    return filter_function

def get_vg_flags(wildcard_flag):
    match wildcard_flag:
        case "gapExt":
            return "--do-gapless-extension"
        case "mqCap":
            return "--explored-cap"
        case downsample_number if downsample_number[0:10] == "downsample":
            return "--downsample-min " + downsample_number[10:]
        case "mapqscale":
            return "--mapq-score-scale 0.01"
        case scorescale_number if scorescale_number[0:10] == "scorescale":
            return "--mapq-score-scale " + scorescale_number[10:]
        case "moreseeds":
            return "--downsample-window-length 400"
        case "mqWindow":
            return "--mapq-score-scale 1 --mapq-score-window 150"
        case "noflags":
            return ""
        case unknown:
            #otherwise this is a hash and we get the flags from ParameterSearch
            return PARAM_SEARCH.hash_to_parameter_string(wildcard_flag)

def get_graphaligner_flags(wildcard_flag):
    match wildcard_flag:
        case "default":
            return "-x vg --multimap-score-fraction 1.0"
        case "fast":
            return "--seeds-mxm-length 30 --seeds-mem-count 10000 --bandwidth 15 --multimap-score-fraction 0.99 --precise-clipping 0.85 --min-alignment-score 100 --clip-ambiguous-ends 100 --overlap-incompatible-cutoff 0.15 --max-trace-count 5 --mem-index-no-wavelet-tree"

def get_vg_version(wildcard_vgversion):
    if wildcard_vgversion == "default":
        return "vg"
    else:
        return "./vg_"+wildcard_vgversion

def cap_reference(wildcards):
    """
    Get the partially-capitalized reference sample name from reference.

    We expect reference to be all lower-case.
    """

    return {"chm13": "CHM13", "chm13v1": "CHM13", "grch38": "GRCh38"}[wildcards["reference"]]


def haplotype_sampling_flags(wildcards):
    """
    Return a string of command line flags for haplotype-sampling a GBZ, from
    hapcount, diploidtag, onereftag, and reference.
    """
    
    parts = []
    
    # Start with the number of haplotypes to sample
    parts.append("--num-haplotypes")
    parts.append(wildcards["hapcount"])

    match wildcards["diploidtag"]:
        case "d":
            # Use diploid sampling
            parts.append("--diploid-sampling")
        case "":
            # Don't use diploid sampling. Nothing to do.
            pass
        case unknown:
            raise RuntimeError("Unimplemented diploid tag " + unknown)

    match wildcards["onereftag"]:
        case "o":
            # Use single-reference sampling
            parts.append("--set-reference")
            parts.append(cap_reference(wildcards))
        case "":
            # Don't use single-reference sampling. Nothing to do.
            pass
        case unknown:
            raise RuntimeError("Unimplemented oneref tag " + unknown)

    return " ".join(parts)

    


    


def param_search_tsvs(wildcards, statname="time_used.mean", realness="real"):
    """
    Get the combined (i.e. mean) TSVs for the conditions in the parameter search.

    TSVs are in the same order as PARAM_SEARCH.get_hashes().

    Needs to be used like:
        lambda w: param_search_tsv(w, "time_used.mean")
    """
    
    # Note r10y2025 doesn't get trimmed
    trimmedness = ".trimmed" if wildcards["tech"] in ("r9", "r10", "q27") and realness == "real" else wildcards["trimmedness"]
    values = dict(wildcards)
    values["trimmedness"] = trimmedness
    values["param_hash"] = PARAM_SEARCH.get_hashes()
    values["realness"] = realness
    values["statname"] = statname
    
    return expand("{root}/stats/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{param_hash}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv", **values)

rule hg_index_graph:
    input:
        gbz="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.gbz"
    output:
        hgfile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.hg"
    threads: 26
    resources:
        mem_mb=120000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "vg convert -a -H -t {threads} {input.gbz} > {output.hgfile}"

rule unchop_hg_graph:
    input:
        hgfile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}.hg"
    output:
        unchopped_hg="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}.unchopped.hg"
    threads: 26
    resources:
        mem_mb=120000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "vg mod -u {input.hgfile} -t {threads} > {output.unchopped_hg}"

rule hg_to_gfa:
    input:
        hgfile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}.hg"
    output:
        gfa="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}.gfa"
    wildcard_constraints:
        chopping="\.unchopped|"
    threads: 26
    resources:
        mem_mb=120000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "vg convert -f -t {threads} {input.hgfile} > {output.gfa}"

rule make_primary_graph:
    input:
        reference_fasta=reference_fasta
    output:
        vg="{graphs_dir}/primary-{reference}.vg"
    threads: 16
    resources:
        mem_mb=60000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg construct -t {threads} -m 1024 -p -r {input.reference_fasta} >{output.vg}"

rule gbz_index_primary_graph:
    input:
        vg="{graphs_dir}/primary-{reference}.vg"
    output:
        gbz="{graphs_dir}/primary-{reference}.gbz"
    threads: 16
    resources:
        mem_mb=60000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg gbwt --progress --num-threads {threads} -x {input.vg} -E --gbz-format -g {output.gbz}"

rule distance_index_graph:
    input:
        gbz="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.gbz"
    output:
        distfile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.dist"
    benchmark: "{graphs_dir}/indexing_benchmarks/distance_indexing_{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.benchmark"
    # TODO: Distance indexing only really uses 1 thread
    threads: 1
    resources:
        mem_mb=lambda w: 120000 if w["full"] == "" else 2000000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "vg index -t {threads} -j {output.distfile} {input.gbz}"

rule precompute_snarls:
    input:
        gbz="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.gbz"
    output:
        snarls="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.snarls"
    benchmark: "{graphs_dir}/indexing_benchmarks/snarls_{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.benchmark"
    threads: 4
    resources:
        mem_mb=120000,
        runtime=180,
        slurm_partition=choose_partition(180)
    shell: "vg snarls -t {threads} -T {input.gbz} > {output.snarls}"

rule tcdist_index_graph:
    input:
        gbz="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}.gbz"
    output:
        tcdistfile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}.tcdist"
    benchmark: "{graphs_dir}/indexing_benchmarks/tcdistance_{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}.benchmark"

    # TODO: Distance indexing only really uses 1 thread
    threads: 1
    resources:
        # TODO: If we could do no-Dijkstra indexing we'd only need 240 minutes
        # and 120G memory, but we need top-level chain distances.
        mem_mb=500000,
        runtime=2880,
        slurm_partition=choose_partition(2880)
    shell:
        "vg index -t {threads} --no-nested-distance -j {output.tcdistfile} {input.gbz}"

rule r_index_graph:
    input:
        # We don't operate on clipped graphs
        gbz="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}.gbz"
    output:
        rifile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}.ri"
    threads: 16
    resources:
        mem_mb=lambda w: 120000 if w["full"] == "" else 240000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "vg gbwt -p --num-threads {threads} -r {output.rifile} -Z {input.gbz}"

rule haplotype_index_graph:
    input:
        unpack(r_and_snarl_indexed_graph),
    output:
        haplfile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}.fragment.hapl"
    params:
        vg_binary=VG_FRAGMENT_HAPLOTYPE_INDEXING_VERSION
    threads: 16
    resources:
        mem_mb=1000000,
        runtime=5760,
        slurm_partition=choose_partition(5760)
    shell:
        "{params.vg_binary} haplotypes -v 2 -t 16 -d {input.snarls} -r {input.ri} {input.gbz} -H {output.haplfile}"

rule xg_index_graph:
    input:
        gbz="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.gbz"
    output:
        xg="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.xg"
    threads: 8
    resources:
        mem_mb=100000,
        runtime=1440,
        slurm_partition=choose_partition(1440)
    shell:
        "vg convert -t {threads} -x --drop-haplotypes {input.gbz} >{output.xg}"

rule kmer_count_full_sample:
    input:
        base_fastq_gz=base_fastq_gz
    output:
        kmer_counts="{reads_dir}/{realness}/{tech}/{sample}/{basename}{trimmedness}.kff"
    params:
        output_basename="{reads_dir}/{realness}/{tech}/{sample}/{basename}{trimmedness}"
    benchmark: "{reads_dir}/{realness}/{tech}/{sample}/{basename}{trimmedness}.kmer_counting.benchmark"
    threads: 8
    resources:
        mem_mb=160000,
        runtime=240,
        slurm_partition=choose_partition(240),
        tmpdir=LARGE_TEMP_DIR
    shell:
        # Need to cut the extension off the output file because KMC will supply it.
        "kmc -k29 -m128 -okff -t{threads} -hp {input.base_fastq_gz} {params.output_basename} \"$TMPDIR\""

rule haplotype_sample_graph:
    input:
        unpack(haplotype_indexed_graph),
        kmer_counts=kmer_counts
    output:
        # Need to sample back into the graphs directory so e.g. minimizer indexing and mapping can work.
        sampled_gbz="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}-sampled{hapcount}{diploidtag}{onereftag}_fragmentlinked-for-{realness}-{tech}-{sample}{trimmedness}-{subset}.gbz"
    benchmark: "{graphs_dir}/indexing_benchmarks/haplotype_sampling_{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}-sampled{hapcount}{diploidtag}{onereftag}_fragmentlinked-for-{realness}-{tech}-{sample}{trimmedness}-{subset}.benchmark"
    wildcard_constraints:
        hapcount="[0-9]+",
        diploidtag="d?",
        onereftag="o?"
    params:
        haplotype_sampling_flags=haplotype_sampling_flags,
        vg_binary=lambda w: get_vg_version(VG_HAPLOTYPE_SAMPLING_VERSION if w["onereftag"] == "" else VG_HAPLOTYPE_SAMPLING_ONEREF_VERSION)
    threads: 8
    resources:
        mem_mb=lambda w: 120000 if w["full"] == "" else 240000,
        runtime=lambda w: 60 if int(w["hapcount"]) <= 32 else 240,
        slurm_partition=lambda w: choose_partition(60) if int(w["hapcount"]) <= 32 else choose_partition(240)
    shell:
        "{params.vg_binary} haplotypes -v 2 -t {threads} --include-reference {params.haplotype_sampling_flags} -i {input.hapl} -k {input.kmer_counts} {input.gbz} -g {output.sampled_gbz}"


rule minimizer_index_graph:
    input:
        unpack(dist_indexed_graph)
    output:
        minfile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.k{k}.w{w}{weightedness}.withzip.min",
        zipfile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.k{k}.w{w}{weightedness}.zipcodes"
    benchmark: "{graphs_dir}/indexing_benchmarks/minimizer_indexing_{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.k{k}.w{w}{weightedness}.benchmark"
    wildcard_constraints:
        weightedness="\\.W|",
        k="[0-9]+",
        w="[0-9]+"
    params:
        weighting_option=lambda w: "--weighted" if w["weightedness"] == ".W" else ""
    threads: 16
    resources:
        mem_mb=lambda w: 600000 if ("hprc-v2" in w["refgraphbase"]  or "hprc-v1.1-nov.11.2024" in w["refgraphbase"]) else 320000 if w["weightedness"] == ".W" else 80000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "vg minimizer --progress -k {wildcards.k} -w {wildcards.w} {params.weighting_option} -t {threads} -p -d {input.dist} -z {output.zipfile} -o {output.minfile} {input.gbz}"

rule non_zipcode_minimizer_index_graph:
    input:
        unpack(dist_indexed_graph),
        non_zipcode_vg="vg_v1.62.0"
    output:
        minfile="{graphs_dir}/{refgraphbase}-{reference}{modifications}{clipping}{full}{chopping}{sampling}.k{k}.w{w}{weightedness}.nozip.min",
    wildcard_constraints:
        weightedness="\\.W|",
        k="[0-9]+",
        w="[0-9]+"
    params:
        weighting_option=lambda w: "--weighted" if w["weightedness"] == ".W" else ""
    threads: 16
    resources:
        mem_mb=lambda w: 600000 if ("hprc-v2" in w["refgraphbase"]  or "hprc-v1.1-nov.11.2024" in w["refgraphbase"]) else 320000 if w["weightedness"] == ".W" else 80000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "./{input.non_zipcode_vg} minimizer --progress -k {wildcards.k} -w {wildcards.w} {params.weighting_option} -t {threads} -p -d {input.dist} -o {output.minfile} {input.gbz}"

rule get_release_vg:
    output:
        "vg_{releaseversion}"
    wildcard_constraints:
        releaseversion=r"v[0-9]+\.[0-9]+\.[0-9]+"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "wget https://github.com/vgteam/vg/releases/download/{wildcards.releaseversion}/vg -O {output} && chmod +x {output}"

# See https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets and https://www.biorxiv.org/content/10.1101/2024.10.02.616333v1.full.pdf
rule get_dipcall_truth_vcf:
    output:
        vcf=TRUTH_DIR + "/{reference}/{sample}/{sample}.dip.vcf.gz",
        index=TRUTH_DIR + "/{reference}/{sample}/{sample}.dip.vcf.gz.tbi"
    wildcard_constraints:
        sample="HG001",
        reference="(chm13|grch38|chm13v1)"
    params:
        cap_reference=cap_reference,
        na_sample=lambda w: {"HG001": "NA12878"}[w["sample"]]
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "curl https://platinum-pedigree-data.s3.amazonaws.com/variants/assembly-based/dipcall/{params.cap_reference}/{params.na_sample}.dip.vcf.gz | bcftools reheader --samples <(echo {wildcards.sample}) | bcftools sort -O z >{output.vcf} && tabix -p vcf {output.vcf}"

rule get_pedigree_truth_vcf:
    output:
        vcf=TRUTH_DIR + "/{reference}/{sample}/{sample}.family-truthset.ov.vcf.gz",
        index=TRUTH_DIR + "/{reference}/{sample}/{sample}.family-truthset.ov.vcf.gz.tbi"
    wildcard_constraints:
        sample="HG001",
        reference="grch38"
    params:
        cap_reference=cap_reference,
        na_sample=lambda w: {"HG001": "NA12878"}[w["sample"]]
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "curl https://platinum-pedigree-data.s3.amazonaws.com/variants/small_variant_truthset/{params.cap_reference}/CEPH1463.GRCh38.family-truthset.ov.vcf.gz | bcftools view --samples {params.na_sample} | bcftools reheader --samples <(echo {wildcards.sample}) | bcftools sort -O z >{output.vcf} && tabix -p vcf {output.vcf}"

rule alias_gam_k:
    input:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{part_subset}000.gam"
    output:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{part_subset}k.gam"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "ln {input.gam} {output.gam}"

rule alias_gam_m:
    input:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{part_subset}000000.gam"
    output:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{part_subset}m.gam" 
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "ln {input.gam} {output.gam}"

rule trim_base_fastq_gz:
    input:
        fq_gz="{reads_dir}/real/{tech}/{sample}/{basename}.{fq_ext}.gz"
    output:
        fq_gz="{reads_dir}/real/{tech}/{sample}/{basename}.trimmed.{fq_ext}.gz"
    wildcard_constraints:
        fq_ext="fq|fastq"
    threads: 8
    resources:
        mem_mb=10000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "seqkit subseq -j {threads} -r 100:-10 {input.fq_gz} -o {output.fq_gz}"
    
rule subset_base_fastq_gz:
    input:
        base_fastq_gz=base_fastq_gz
    output:
        fastq_gz="{reads_dir}/{realness}/{tech}/{sample}/{basename}{trimmedness}.{subset}.fq.gz"
    wildcard_constraints:
        realness="real",
        subset="[0-9]+[km]?"
    params:
        lines=lambda w: str(subset_to_number(w["subset"]) * 4)
    threads: 8
    resources:
        mem_mb=10000,
        runtime=120,
        slurm_partition=choose_partition(120)
    shell:
        # We need to account for bgzip getting upset that we close the pipe before it is done writing.
        "(bgzip -d <{input.base_fastq_gz} || true) | head -n {params.lines} | pigz -p8 -c >{output.fastq_gz}"

rule decompress_fastq:
    input:
        fastq_gz="{reads_dir}/{realness}/{tech}/{sample}/{basename}{trimmedness}.{subset}.fq.gz"
    output:
        fastq="{reads_dir}/{realness}/{tech}/{sample}/{basename}{trimmedness}.{subset}.fq"
    threads: 8
    resources:
        mem_mb=1000,
        runtime=120,
        slurm_partition=choose_partition(120)
    shell:
        "pigz -d -c -p 8 {input.fastq_gz} > {output.fastq}"

rule subset_alias_base_fastq_gz:
    input:
        base_fastq_gz=base_fastq_gz
    output:
        fastq_gz="{reads_dir}/{realness}/{tech}/{sample}/{basename}{trimmedness}.{subset}.fq.gz"
    wildcard_constraints:
        realness="real",
        subset="full"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "ln {input.base_fastq_gz} {output.fastq_gz}"

rule extract_fastq_from_gam:
    input:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"
    output:
        fastq="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.fq"
    threads: 16
    resources:
        mem_mb=10000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg view --fastq-out --threads {threads} {input.gam} >{output.fastq}"

rule extract_fastq_gz_from_gam:
    input:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"
    output:
        fastq_gz="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.fq.gz"
    threads: 24
    resources:
        mem_mb=10000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg view --fastq-out --threads 16 {input.gam} | pigz -p 8 -c >{output.fastq_gz}"

rule extract_fastq_from_full_gam:
    input:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}.gam"
    output:
        fastq="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}.fq"
    threads: 16
    resources:
        mem_mb=10000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg view --fastq-out --threads {threads} {input.gam} >{output.fastq}"

rule extract_fastq_gz_from_full_gam:
    input:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}.gam"
    output:
        fastq_gz="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}.fq.gz"
    threads: 24
    resources:
        mem_mb=10000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg view --fastq-out --threads 16 {input.gam} | pigz -p 8 -c >{output.fastq_gz}"


# From Jean, chunk and unchunk reads, because graphaligner takes too long and errors when given too many threads
# number of read chunks (for GraphAligner)
if 'n_read_chunk' not in config:
    config['n_read_chunk'] = 20
READ_CHUNKS = list(range(config['n_read_chunk']))

#Split read fq.gz file into chunks
rule split_reads:
    input: "{reads_dir}/{realness}/{tech}/{sample}/{basename}{trimmedness}.{subset}.fq.gz"
    output: temp(expand('{{reads_dir}}/{{realness}}/{{tech}}/{{sample}}/{{basename}}{{trimmedness}}.{{subset}}.chunk{chunk}.fq.gz', chunk=READ_CHUNKS))
    threads: 4
    resources:
        mem_mb=12000,
        runtime=600,
        slurm_partition=choose_partition(600)
    container: 'docker://quay.io/biocontainers/fastqsplitter:1.2.0--py312hf67a6ed_6'
    params:
        outfiles=lambda wildcards, output: '-o ' + ' -o '.join(output)
    shell:
        """
        fastqsplitter -i {input} {params.outfiles}
        """

# Put mapped chunks of reads back together
rule merge_graphaligner_gams:
    input: expand('{{root}}/aligned-secsup/{{reference}}/{{refgraph}}/{{mapper}}-{{graphalignerflag}}/{{realness}}/{{tech}}/{{sample}}{{trimmedness}}.{{subset}}.chunk{chunk}.gam', chunk=READ_CHUNKS)
    output:
        gam=temp("{root}/aligned-secsup/{reference}/{refgraph}/{mapper}-{graphalignerflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam")
    wildcard_constraints:
        realness="real",
        mapper="graphaligner.*",
        subset="full",
        refgraph="hprc-v1.1-mc"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell: "cat {input} > {output}"

# Prefer to chunk gams for graphaligner, but only works for real reads on the full hprc graph
ruleorder: merge_graphaligner_gams > graphaligner_real_reads

rule dict_index_reference:
    input:
        reference_fasta=REFS_DIR + "/{basename}.fa"
    output:
        index=REFS_DIR + "/{basename}.fa.dict"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        # TODO: Needs picard.jar sitting in this directory
        "java -jar ./picard.jar CreateSequenceDictionary R={input.reference_fasta} O={output.index}"

rule paths_index_reference:
    input:
        reference_dict=REFS_DIR + "/{basename}.fa.dict"
    output:
        index=REFS_DIR + "/{basename}.fa.paths{region}.txt"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.reference_dict} | grep '^@SQ' | sed 's/^@SQ.*[[:blank:]]SN:\\([^[:blank:]]*\\).*/\\1/g' | grep '{wildcards.region}$' > {output.index}"

rule callable_paths_index_reference:
    input:
        paths=REFS_DIR + "/{basename}.fa.paths{region}.txt"
    output:
        paths=REFS_DIR + "/{basename}.fa.paths{region}.callable.txt"
    params:
        uncallable_contig_regex=uncallable_contig_regex
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell: "cat {input.paths} | grep -v '{params.uncallable_contig_regex}' > {output.paths}"

rule callable_paths_dict_reference:
    input:
        paths=REFS_DIR + "/{basename}.fa.paths{region}.callable.txt",
        reference_dict=REFS_DIR + "/{basename}.fa.dict"
    output:
        callable_dict=REFS_DIR + "/{basename}.fa.paths{region}.callable.dict",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    run:
        wanted_paths = set((line.strip() for line in open(input.paths) if line.strip()))
        with open(output.callable_dict, "w") as out:
            for line in open(input.reference_dict):
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if parts[0] == "@HD":
                    out.write(line)
                    out.write("\n")
                elif parts[0] == "@SQ":
                    # Retain only @SQ lines from the dict for callable paths
                    keep = False
                    for part in parts:
                        if part.startswith("SN:"):
                            name = part[3:]
                            if name in wanted_paths:
                                keep = True
                                break
                    if keep:
                        out.write(line)
                        out.write("\n")

rule callable_bed_index_calling_reference:
    input:
        # First two columns of FAI are name and length
        fai=REFS_DIR + "/{calling_reference}.fa.fai"
    output:
        bed=REFS_DIR + "/{calling_reference}.fa.callable.from.{reference}.bed"
    params:
        uncallable_contig_regex=uncallable_contig_regex
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell: "cat {input.fai} | cut -f1,2 | grep -v '{params.uncallable_contig_regex}' | sed 's/\\t/\\t0\\t/g' > {output.bed}"

rule giraffe_real_reads:
    input:
        unpack(indexed_graph),
        fastq_gz=fastq_gz,
    output:
        # Giraffe can dump out pre-annotated reads at annotation range -1.
        gam="{root}/aligned/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    log:"{root}/aligned/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    benchmark: "{root}/aligned/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    wildcard_constraints:
        realness="real"
    threads: auto_mapping_threads
    params:
        exclusive_timing=EXCLUSIVE_TIMING 
    resources:
        mem_mb=auto_mapping_memory,
        runtime=1200,
        slurm_partition=choose_partition(1200),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    run:
        vg_binary = get_vg_version(wildcards.vgversion)
        flags=get_vg_flags(wildcards.vgflag)
        pairing_flag="-i" if wildcards.preset == "default" else ""
        zipcodes_flag=f"-z {input.zipfile}" if "zipfile" in dict(input).keys() else ""

        shell(vg_binary + " giraffe -t{threads} --parameter-preset {wildcards.preset} --progress -Z {input.gbz} -d {input.dist} -m {input.minfile} -f {input.fastq_gz} " + zipcodes_flag + " " + flags + " " + pairing_flag + " >{output.gam} 2>{log}")

rule giraffe_sim_reads:
    input:
        unpack(indexed_graph),
        gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/annotated-1/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    log:"{root}/annotated-1/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="sim"
    threads: auto_mapping_threads
    resources:
        mem_mb=auto_mapping_memory,
        runtime=600,
        slurm_partition=choose_partition(600)
    run:
        vg_binary = get_vg_version(wildcards.vgversion)
        flags=get_vg_flags(wildcards.vgflag)
        pairing_flag="-i" if wildcards.preset == "default" else ""
        zipcodes_flag=f"-z {input.zipfile}" if "zipfile" in dict(input).keys() else ""

        shell(vg_binary + " giraffe -t{threads} --parameter-preset {wildcards.preset} --progress --track-provenance --set-refpos -Z {input.gbz} -d {input.dist} -m {input.minfile} -G {input.gam} " + zipcodes_flag + " " + flags + " " + pairing_flag + " >{output.gam} 2>{log}")

rule giraffe_sim_reads_with_correctness:
    input:
        unpack(indexed_graph),
        gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/correctness/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    log:"{root}/correctness/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="sim"
    threads: auto_mapping_threads
    resources:
        mem_mb=auto_mapping_memory,
        runtime=600,
        slurm_partition=choose_partition(600)
    run:
        vg_binary = get_vg_version(wildcards.vgversion)
        flags=get_vg_flags(wildcards.vgflag)
        pairing_flag="-i" if wildcards.preset == "default" else ""
        zipcodes_flag=f"-z {input.zipfile}" if "zipfile" in dict(input).keys() else ""

        shell(vg_binary + " giraffe -t{threads} --parameter-preset {wildcards.preset} --progress --track-provenance --track-correctness --set-refpos -Z {input.gbz} -d {input.dist} -m {input.minfile} -G {input.gam} " + zipcodes_flag + " " + flags + " " + pairing_flag + " >{output.gam} 2>{log}")

rule winnowmap_repetitive_kmers:
    input:
        fasta=REFS_DIR + "/{basename}.fa"
    output:
        kmers=REFS_DIR + "/{basename}.repetitive_k15.txt",
        db=temp(directory(REFS_DIR + "/{basename}.db"))
    shell:
        "meryl count k=15 output {output.db} {input.fasta} && meryl print greater-than distinct=0.9998 {output.db} > {output.kmers}"

rule winnowmap_sim_reads:
    input:
        reference_fasta=reference_fasta,
        repetitive_kmers=repetitive_kmers,
        fastq=fastq
    params:
        mode=minimap_derivative_mode,
    output:
        sam=temp("{root}/aligned-secsup/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam")
    log:"{root}/aligned-secsup/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="sim"
    wildcard_constraints:
        # Winnowmap doesn't have a short read preset, so we can't do Illumina reads.
        # So match any string but that. See https://stackoverflow.com/a/14683066
        tech="(?!illumina).+"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "winnowmap -t {threads} -W {input.repetitive_kmers} -ax {params.mode} {input.reference_fasta} {input.fastq} >{output.sam} 2>{log}"

rule winnowmap_real_reads:
    input:
        reference_fasta=reference_fasta,
        repetitive_kmers=repetitive_kmers,
        fastq_gz=fastq_gz
    params:
        mode=minimap_derivative_mode,
    output:
        sam=temp("{root}/aligned-secsup/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam")
    log:"{root}/aligned-secsup/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    benchmark: "{root}/aligned-secsup/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    wildcard_constraints:
        realness="real"
    wildcard_constraints:
        # Winnowmap doesn't have a short read preset, so we can't do Illumina reads.
        # So match any string but that. See https://stackoverflow.com/a/14683066
        tech="(?!illumina).+"
    threads: auto_mapping_threads
    params:
        exclusive_timing=EXCLUSIVE_TIMING
    resources:
        mem_mb=300000,
        runtime=1600,
        slurm_partition=choose_partition(1600),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "winnowmap -t {threads} -W {input.repetitive_kmers} -ax {params.mode} {input.reference_fasta} {input.fastq_gz} >{output.sam} 2>{log}"

rule minimap2_index_reference:
    input:
        reference_fasta=REFS_DIR + "/{basename}.fa"
    output:
        index=REFS_DIR + "/{basename}.{preset}.mmi"
    threads: 16
    resources:
        mem_mb=16000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
         "minimap2 -t {threads} -x {wildcards.preset} -d {output.index} {input.reference_fasta}"


rule minimap2_sim_reads:
    input:
        minimap2_index=minimap2_index,
        fastq=fastq
    output:
        sam=temp("{root}/aligned-secsup/{reference}/minimap2-{minimapmode}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam")
    log:"{root}/aligned-secsup/{reference}/minimap2-{minimapmode}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="sim"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "minimap2 -t {threads} -ax {wildcards.minimapmode} --secondary=no {input.minimap2_index} {input.fastq} >{output.sam} 2> {log}"

rule minimap2_real_reads:
    input:
        minimap2_index=minimap2_index,
        fastq_gz=fastq_gz
    output:
        sam=temp("{root}/aligned-secsup/{reference}/minimap2-{minimapmode}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam")
    benchmark: "{root}/aligned-secsup/{reference}/minimap2-{minimapmode}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    log: "{root}/aligned-secsup/{reference}/minimap2-{minimapmode}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="real"
    threads: auto_mapping_threads
    params:
        exclusive_timing=EXCLUSIVE_TIMING
    resources:
        mem_mb=300000,
        runtime=1200,
        slurm_partition=choose_partition(1200),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "minimap2 -t {threads} -ax {wildcards.minimapmode} --secondary=no {input.minimap2_index} {input.fastq_gz} >{output.sam} 2> {log}"

# pbmm2 uses the same indexes as minimap2
# TODO: pbmm2 makes BAM output but we're calling it SAM here.
rule pbmm2_sim_reads:
    input:
        pbmm2_index=pbmm2_hifi_index,
        fastq=fastq
    output:
        sam=temp("{root}/aligned-secsup/{reference}/pbmm2/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam")
    log:"{root}/aligned-secsup/{reference}/pbmm2/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="sim",
        tech="hifi"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "pbmm2 align --num-threads {threads} --preset HiFi --best-n 1 --unmapped {input.pbmm2_index} {input.fastq} >{output.sam} 2> {log}"

rule pbmm2_real_reads:
    input:
        pbmm2_index=pbmm2_hifi_index,
        fastq_gz=fastq_gz
    output:
        sam=temp("{root}/aligned-secsup/{reference}/pbmm2/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam")
    benchmark: "{root}/aligned-secsup/{reference}/pbmm2/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    log: "{root}/aligned-secsup/{reference}/pbmm2/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="real",
        tech="hifi"
    threads: auto_mapping_threads
    params:
        exclusive_timing=EXCLUSIVE_TIMING
    resources:
        mem_mb=300000,
        runtime=1200,
        slurm_partition=choose_partition(1200),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "pbmm2 align --num-threads {threads} --preset HiFi --best-n 1 --unmapped {input.pbmm2_index} {input.fastq_gz} >{output.sam} 2> {log}"

# The BWA index file names are implicit but can still be dependencies.
rule bwa_index_reference:
    input:
        reference_fasta=REFS_DIR + "/{basename}.fa"
    output:
        amb=REFS_DIR + "/{basename}.fa.amb",
        ann=REFS_DIR + "/{basename}.fa.ann",
        bwt=REFS_DIR + "/{basename}.fa.bwt",
        pac=REFS_DIR + "/{basename}.fa.pac",
        sa=REFS_DIR + "/{basename}.fa.sa"
    threads: 2
    resources:
        mem_mb=16000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
         "bwa index {input.reference_fasta}"

rule bwa_sim_reads:
    input:
        unpack(bwa_index_set),
        reference_fasta=reference_fasta,
        fastq_gz=fastq_gz,
    output:
        sam=temp("{root}/aligned-secsup/{reference}/bwa{pairing}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam")
    log:"{root}/aligned-secsup/{reference}/bwa{pairing}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="sim",
        tech="illumina",
        pairing="(-pe|)"
    params:
        pairing_flag=lambda w: "-p" if w["pairing"] == "-pe" else ""
    threads: auto_mapping_threads
    resources:
        mem_mb=30000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "bwa mem -t {threads} -R '@RG\\tID:1\\tLB:lib1\\tSM:{wildcards.sample}\\tPL:{wildcards.tech}\\tPU:unit1' {params.pairing_flag} {input.reference_fasta} {input.fastq_gz} >{output.sam} 2> {log}"

rule bwa_real_reads:
    input:
        unpack(bwa_index_set),
        reference_fasta=reference_fasta,
        fastq_gz=fastq_gz
    output:
        sam=temp("{root}/aligned-secsup/{reference}/bwa{pairing}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam")
    benchmark: "{root}/aligned-secsup/{reference}/bwa{pairing}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    log: "{root}/aligned-secsup/{reference}/bwa{pairing}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="real",
        pairing="(-pe|)"
    params:
        pairing_flag=lambda w: "-p" if w["pairing"] == "-pe" else "",
        exclusive_timing=EXCLUSIVE_TIMING
    threads: auto_mapping_threads
    resources:
        mem_mb=30000,
        runtime=1200,
        slurm_partition=choose_partition(1200),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "bwa mem -t {threads} -R '@RG\\tID:1\\tLB:lib1\\tSM:{wildcards.sample}\\tPL:{wildcards.tech}\\tPU:unit1' {params.pairing_flag} {input.reference_fasta} {input.fastq_gz} >{output.sam} 2> {log}"

# Minimap2 and Winnowmap include secondary alignments in the output by default, and Winnowmap doesn't quite have a way to limit them (minimap2 has -N)
# Also they only speak SAM and we don't want to benchmark the BAM-ification time.
# So drop secondary and supplementary alignments and convert to BAM.
# TODO: Get the downstream stats tools to be able to do things like measure average softclips when there are supplementary alignments.
rule drop_secondary_and_supplementary:
    input:
        sam="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam"
    output:
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    wildcard_constraints:
        mapper="(minimap2-.*|winnowmap|bwa(-pe)?|pbmm2)"
    threads: 16
    resources:
        mem_mb=30000,
        runtime=600,
        slurm_partition=choose_partition(600),
    shell:
        "samtools view --threads 16 -h -F 2048 -F 256 --bam {input.sam} >{output.bam}"

rule graphaligner_sim_reads:
    input:
        gfa=gfa,
        fastq=fastq
    output:
        gam="{root}/aligned-secsup/{reference}/{refgraph}/{mapper}-{graphalignerflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        realness="sim",
        # Need to hint mapper to auto_mapping_threads
        mapper="graphaligner.*"
    threads: auto_mapping_threads
    params:
        mapping_threads=lambda wildcards, threads: threads if threads <= 2 else threads-2,
        flags=lambda wildcards: get_graphaligner_flags(wildcards.graphalignerflag)
    container: "docker://quay.io/biocontainers/graphaligner:1.0.20--h06902ac_0"
    resources:
        mem_mb=300000,
        runtime=3000,
        slurm_partition=choose_partition(3000)
    shell:
        "GraphAligner -t {params.mapping_threads} -g {input.gfa} -f {input.fastq} {params.flags} -a {output.gam}"


rule graphaligner_real_reads:
    input:
        gfa=gfa,
        fastq_gz=fastq_gz
    output:
        gam=temp("{root}/aligned-secsup/{reference}/{refgraph}/{mapper}-{graphalignerflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}{readchunk}.gam")
    benchmark: "{root}/aligned/{reference}/{refgraph}/{mapper}-{graphalignerflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}{readchunk}.benchmark"
    log: "{root}/aligned/{reference}/{refgraph}/{mapper}-{graphalignerflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}{readchunk}.log"
    wildcard_constraints:
        realness="real",
        mapper="graphaligner.*",
        readchunk=".chunk[0-9]*|"
    threads: auto_mapping_threads
    params:
        mapping_threads=lambda wildcards, threads: threads if threads <= 2 else threads-2,
        flags=lambda wildcards: get_graphaligner_flags(wildcards.graphalignerflag),
        exclusive_timing=EXCLUSIVE_TIMING
    container: "docker://quay.io/biocontainers/graphaligner:1.0.20--h06902ac_0"
    resources:
        mem_mb=900000,
        runtime=9000,
        slurm_partition=choose_partition(9000),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "GraphAligner -t {params.mapping_threads} -g {input.gfa} -f {input.fastq_gz} {params.flags} -a {output.gam} 2> {log}"

rule minigraph_sim_reads:
    input:
        gfa=gfa,
        fastq=fastq
    output:
        gaf=temp("{root}/aligned-secsup/{reference}/{refgraph}/minigraph/{realness}/{tech}/{sample}{trimmedness}.{subset}.gaf")
    wildcard_constraints:
        realness="sim"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "minigraph --vc --secondary=no -cx lr -t {threads} {input.gfa} {input.fastq} >{output.gaf}"

rule minigraph_real_reads:
    input:
        gfa=gfa,
        fastq_gz=fastq_gz
    output:
        gaf=temp("{root}/aligned-secsup/{reference}/{refgraph}/minigraph/{realness}/{tech}/{sample}{trimmedness}.{subset}.gaf")
    benchmark: "{root}/aligned/{reference}/{refgraph}/minigraph/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    log: "{root}/aligned/{reference}/{refgraph}/minigraph/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="real"
    threads: auto_mapping_threads
    params:
        exclusive_timing=EXCLUSIVE_TIMING
    resources:
        mem_mb=300000,
        runtime=1200,
        slurm_partition=choose_partition(1200),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "minigraph --vc --secondary=no -cx lr -t {threads} {input.gfa} {input.fastq_gz} >{output.gaf} 2>{log}"


rule panaligner_sim_reads:
    input:
        gfa=gfa,
        fastq=fastq
    output:
        gaf="{root}/aligned/{reference}/{refgraph}/panaligner/{realness}/{tech}/{sample}{trimmedness}.{subset}.gaf"
    wildcard_constraints:
        realness="sim"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "PanAligner --vc -cx lr -t {threads} {input.gfa} {input.fastq} >{output.gaf}"

rule panaligner_real_reads:
    input:
        gfa=gfa,
        fastq=fastq
    output:
        gaf="{root}/aligned/{reference}/{refgraph}/panaligner/{realness}/{tech}/{sample}{trimmedness}.{subset}.gaf"
    benchmark: "{root}/aligned/{reference}/{refgraph}/panaligner/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    log: "{root}/aligned/{reference}/{refgraph}/panaligner/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    wildcard_constraints:
        realness="real"
    threads: auto_mapping_threads
    params:
        exclusive_timing=EXCLUSIVE_TIMING
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "PanAligner --vc -cx lr -t {threads} {input.gfa} {input.fastq} >{output.gaf} 2>{log}"

rule gaf_to_gam:
    input:
        hg=hg,
        gaf="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gaf"
    output:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        mapper="(minigraph|panaligner)"
    threads: 16
    resources:
        mem_mb=100000,
        runtime=300,
        slurm_partition=choose_partition(300)
    shell:
        "vg convert -t {threads} --gaf-to-gam {input.gaf} {input.hg} >{output.gam}"

rule gam_to_gaf:
    input:
        hg=hg,
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    output:
        gaf="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gaf"
    wildcard_constraints:
        mapper="(graphaligner-.*)"
    threads: 16
    resources:
        mem_mb=100000,
        runtime=300,
        slurm_partition=choose_partition(300)
    shell:
        "vg convert -t {threads} --gam-to-gaf {input.gam} {input.hg} >{output.gaf}"

rule select_first_duplicate_read_gam:
    input:
        gam="{root}/aligned-secsup/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    output:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        mapper="(graphaligner-.*)"
    threads: 2
    resources:
        mem_mb=30000,
        runtime=1600,
        slurm_partition=choose_partition(1600)
    shell:
        "vg filter -l {input.gam} > {output.gam}"

rule select_best_duplicate_read_gaf:
    input:
        gaf="{root}/aligned-secsup/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gaf"
    output:
        gaf="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gaf"
    wildcard_constraints:
        mapper="(minigraph)"
    threads: 2
    resources:
        mem_mb=30000,
        runtime=1600,
        slurm_partition=choose_partition(1600)
    shell:
        "cat {input.gaf} | python3 select_best_gaf.py > {output.gaf}"

rule inject_bam:
    input:
        gbz=gbz,
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        mapper="(minimap2.+|winnowmap|bwa(-pe)?|pbmm2)"
    threads: 64
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg inject --threads {threads} -x {input.gbz} {input.bam} >{output.gam}"

rule surject_gam:
    input:
        gbz=gbz,
        # We surject onto all target reference paths, even those we can't call
        # on, like Y in CHM13 (due to different Ys being used in different
        # graphs), and the _random contigs in GRCh38. Otherwise mappers mapping
        # against GRCh38 with its _random contigs have an advantage when
        # calling there, because Giraffe will stick all the _random reads in
        # the main genome.
        reference_dict=reference_dict,
        gam=surjectable_gam
    output:
        # Don't keep this around because we're going to keep the sorted version and use that for most things.
        bam=temp("{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam")
    wildcard_constraints:
        mapper="(giraffe.*|graphaligner-.*)"
    params:
        # The default preset is paired
        paired_flag=lambda w: "-i" if re.match("giraffe-[^-]*-default-[^-]*-[^-]*", w["mapper"]) else ""
    threads: 40
    resources:
        mem_mb=lambda w: (600000 / 64 * 40) if w["tech"] in ("r10", "r10y2025") else (150000 / 64 * 40),
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg surject -F {input.reference_dict} -x {input.gbz} -t {threads} --bam-output --sample {wildcards.sample} --read-group \"ID:1 LB:lib1 SM:{wildcards.sample} PL:{wildcards.tech} PU:unit1\" --prune-low-cplx {params.paired_flag} {input.gam} > {output.bam}"

rule sort_gam:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    output: 
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.gam",
        gam_index="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.gam.gai"
    wildcard_constraints:
        mapper="(giraffe.*|graphaligner-.*)"
    threads: 8
    resources:
        mem_mb=60000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg gamsort -t {threads} -i {output.gam_index} {input.gam} >{output.gam}"

rule alias_bam_graph:
    # For BAM-generating mappers we can view their BAMs as if they mapped to any reference graph for a reference
    input:
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        bam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    wildcard_constraints:
        mapper="(minimap2.+|winnowmap|bwa(-pe)?|pbmm2)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "ln {input.bam} {output.bam}"

rule sort_bam:
    input:
        bam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        bam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.bam",
        bai="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.bam.bai"
    threads: 16
    resources:
        mem_mb=16000,
        runtime=600,
        slurm_partition=choose_partition(600)
    run:
        with tempfile.TemporaryDirectory() as sort_scratch:
            shell("samtools sort -T " + os.path.join(sort_scratch, "scratch") + " --threads {threads} {input.bam} -O BAM > {output.bam} && samtools index -b {output.bam} {output.bai}")

rule reheader_bam:
    input:
        bam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.bam"
    output:
        bam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.reheadered.bam",
        bai="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.reheadered.bam.bai"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "samtools reheader -c 'sed -e \"s/SN:[^\\t#]*#0#/SN:/g\"' {input.bam} >{output.bam} && samtools index -b {output.bam} {output.bai}"

rule call_variants_dv:
    input:
        sorted_bam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.bam",
        sorted_bam_index="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sorted.bam.bai",
        reference_path_list_callable=reference_path_list_callable,
        calling_reference_fasta=calling_reference_fasta,
        calling_reference_fasta_index=calling_reference_fasta_index,
        calling_reference_restrict_bed=calling_reference_restrict_bed,
        calling_reference_par_bed=calling_reference_par_bed,
        truth_vcf=lambda w: remote_or_local(truth_vcf_url(w)),
        truth_vcf_index=lambda w: remote_or_local(truth_vcf_index_url(w)),
        truth_bed=lambda w: remote_or_local(truth_bed_url(w))
    output:
        wdl_output_file="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.json",
        # TODO: make this temp so we can delete it?
        wdl_output_directory=directory("{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.wdlrun"),
        vcf="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.vcf.gz",
        vcf_index="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.vcf.gz.tbi",
        happy_evaluation_archive="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.happy_results.tar.gz"
    log:
        logfile="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.log"
    params:
        wdl_input_file="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.input.json",
        job_store="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.jobstore",
        batch_logs="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.batchlogs",
        reference_prefix=reference_prefix,
        haploid_contigs=haploid_contigs,
        wdl_cache=wdl_cache,
        model_files=model_files
    threads: 2
    resources:
        mem_mb=4000,
        # Since Toil tasks need to schedule through the cluster again, this
        # might need to wait a long time.
        # TODO: Run single machine so Slurm can't plan to run our jobs *after*
        # the leader finishes?
        runtime=8640,
        slurm_partition=choose_partition(8640)
    run:
        import json
        wf_source = "#workflow/github.com/vgteam/vg_wdl/DeepVariant:lr-giraffe"
        wf_inputs = {
            "DeepVariant.MERGED_BAM_FILE": input.sorted_bam,
            "DeepVariant.MERGED_BAM_FILE_INDEX": input.sorted_bam_index,
            "DeepVariant.SAMPLE_NAME": wildcards.sample,
            "DeepVariant.PATH_LIST_FILE": input.reference_path_list_callable,
            "DeepVariant.REFERENCE_PREFIX": params.reference_prefix,
            "DeepVariant.REFERENCE_PREFIX_ON_BAM": True,
            "DeepVariant.REFERENCE_FILE": input.calling_reference_fasta,
            # TODO: Should we not left align for minimap2 or GraphAligner?
            "DeepVariant.LEFTALIGN_BAM": True,
            # Indel realignment tools aren't available for long reads.
            # Also, we can't realign indels on our CHM13 non-2.0 mapping
            # reference when using our CHM13v2.0 calling reference.
            # TODO: Turn indel realignment back on if we ever we manage to use
            # CHM13v2.0. Which we probably won't because it uses HG002 Y and
            # we hold that out of the eval graphs.
            "DeepVariant.REALIGN_INDELS": False,
            "DeepVariant.DV_MODEL_TYPE": {"hifi": "PACBIO", "r10": "ONT_R104", "r10y2025": "ONT_R104", "illumina": "WGS"}[wildcards.tech],
            # We don't send MIN_MAPQ and let the model example info define it.
            # TODO: Should we use legacy AC like in the paper?
            # We can't actually turn off legacy AC if it's on in the model example_info. So don't pass DV_KEEP_LEGACY_AC either.
            # Read normalization should work for long reads as of gcr.io/deepvariant-docker/deepvariant:head756846963 (which is post-1.9)
            # But we also don't send DV_NORM_READS because the model defines it as of google/deepvariant:CL782981885
            "DeepVariant.DV_GPU_DOCKER": "google/deepvariant:CL782981885-gpu",
            "DeepVariant.DV_NO_GPU_DOCKER": "google/deepvariant:CL782981885",
            "DeepVariant.TRUTH_VCF": to_local(input.truth_vcf),
            "DeepVariant.TRUTH_VCF_INDEX": to_local(input.truth_vcf_index),
            "DeepVariant.EVALUATION_REGIONS_BED": to_local(input.truth_bed),
            "DeepVariant.RESTRICT_REGIONS_BED": input.calling_reference_restrict_bed,
            "DeepVariant.PAR_REGIONS_BED_FILE": input.calling_reference_par_bed,
            "DeepVariant.HAPLOID_CONTIGS": params.haploid_contigs,
            # We just need hap.py; we don't need a stand-alone vcfeval run.
            "DeepVariant.RUN_STANDALONE_VCFEVAL": False,
            # We also don't need BAMs.
            "DeepVariant.OUTPUT_SINGLE_BAM": False,
            "DeepVariant.OUTPUT_CALLING_BAMS": False,
            "DeepVariant.CALL_CORES": 8 * 4,
            "DeepVariant.CALL_MEM": 50 * 4,
            "DeepVariant.MAKE_EXAMPLES_MEM": 50
        }
        if params.model_files:
            # Use a model that's not built in.
            wf_inputs["DeepVariant.DV_MODEL_FILES"] = params.model_files
        json.dump(wf_inputs, open(params["wdl_input_file"], "w"))
        shell("rm -Rf {params.job_store}")
        # Run and keep the first manageable amount of logs not sent to the log
        # file in case we can't start. Don't stop when we hit the log limit.
        # See https://superuser.com/a/1531706
        toil_command = "MINIWDL__CALL_CACHE__GET=true MINIWDL__CALL_CACHE__PUT={FILL_WDL_CACHE} MINIWDL__CALL_CACHE__DIR={params.wdl_cache} toil-wdl-runner '" + wf_source + "' {params.wdl_input_file} --clean=onSuccess --jobStore file:{params.job_store} --batchLogsDir {params.batch_logs} --wdlOutputDirectory {output.wdl_output_directory} --wdlOutputFile {output.wdl_output_file} --batchSystem slurm --slurmTime 11:59:59 --disableProgress --caching=False --logFile={log.logfile} 2>&1 | (head -c1000000; cat >/dev/null)"
        print("Running Toil: " + toil_command)
        shell(toil_command)

        wdl_result=json.load(open(output.wdl_output_file))
        shell("cp " + wdl_result["DeepVariant.output_vcf"] + " {output.vcf}")
        shell("cp " + wdl_result["DeepVariant.output_vcf_index"] + " {output.vcf_index}")
        shell("cp " + wdl_result["DeepVariant.output_happy_evaluation_archive"] + " {output.happy_evaluation_archive}")
        
rule extract_happy_summary:
    input:
        happy_evaluation_archive="{root}/called/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.happy_results.tar.gz"
    output:
        happy_evaluation_summary="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.eval.summary.csv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "tar -xOf {input.happy_evaluation_archive} happy_results/eval.summary.csv >{output.happy_evaluation_summary}"

rule stat_from_happy_summary:
    input:
        happy_evaluation_summary="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}.eval.summary.csv"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.{vartype}_{colname}.tsv"
    wildcard_constraints:
        vartype="(snp|indel)",
        colname="(f1|precision|recall|fn|fp|tp)"
    params:
        colnum=lambda w: {"f1": 14, "precision": 12, "recall": 11, "fn": 5, "fp": 7, "tp":4}[w["colname"]]
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    run:
        shell("cat {input.happy_evaluation_summary} | grep '^" + wildcards["vartype"].upper() + ",PASS' | cut -f{params.colnum} -d',' >{output.tsv}")

rule total_errors_from_fp_and_fn:
    input:
        snp_fn="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.snp_fn.tsv",
        snp_fp="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.snp_fp.tsv",
        indel_fn="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.indel_fn.tsv",
        indel_fp="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.indel_fp.tsv"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.total_errors.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.snp_fn} {input.snp_fp} {input.indel_fn} {input.indel_fp} | awk '{{sum += $1 }} END {{ print sum }}' >{output.tsv}"

rule vartype_errors_from_fp_and_fn:
    input:
        vartype_fn="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.{vartype}_fn.tsv",
        vartype_fp="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.{vartype}_fp.tsv",
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{region}{callparams}{dot}{category}.{vartype}_errors.tsv"
    wildcard_constraints:
        vartype="(snp|indel)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.vartype_fn} {input.vartype_fp} | awk '{{sum += $1 }} END {{ print sum }}' >{output.tsv}"

#This outputs a tsv of: tp, fn, fp, recall, precision, f1
rule dv_summary_by_condition:
    input:
        happy_evaluation_summary="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{callparams}.eval.summary.csv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/stats/deepvariant/{mapper}/{realness}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}{callparams}.dv_{vartype}_summary.tsv"
    wildcard_constraints:
        refgraph="[^/_]+",
        vartype="(snp|indel)"
    threads: 1
    resources:
        mem_mb=200,
        runtime=10,
        slurm_partition=choose_partition(10)
    run:
        shell("echo \"{params.condition_name}\t$(cat {input.happy_evaluation_summary} | grep '^" + wildcards["vartype"].upper() + ",PASS' | cut -f 4,5,7,11,12,14 -d',' | sed 's/,/\\t/g')\" >{output.tsv}")

rule dv_summary_table:
    input:
        tsv=lambda w: all_experiment(w, "{root}/experiments/{expname}/stats/deepvariant/{mapper}/{realness}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}{callparams}.dv_{vartype}_summary.tsv")
    output:
        tsv="{root}/experiments/{expname}/results/dv_{vartype}_summary.tsv"
    wildcard_constraints:
        vartype="(snp|indel)"
    threads: 1
    resources:
        mem_mb=200,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        """
        printf "condition\tTP\tFN\tFP\trecall\tprecision\tF1\n" >> {output.tsv}
        cat {input} >>{output.tsv}
        """



rule compare_alignments:
    input:
        gam="{root}/annotated-1/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        truth_gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        tsv="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv",
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    params:
        # In the v2.0 CHM13-based graphs we now use a non-HG002 Y, and that's where our
        # Y truth positions are. For other graphs (like v2 prerelease 3 or v1.1) we won't
        # have the right CHM13 Y.
        # So we just ignore Y in CHM13. 
        ignore_flag="--ignore 'CHM13#0#chrY'"
    threads: 16
    resources:
        mem_mb=200000,
        runtime=800,
        slurm_partition=choose_partition(800)
    shell:
        "vg gamcompare --threads 16 --range 200 {params.ignore_flag} {input.gam} {input.truth_gam} --output-gam {output.gam} -T -a {wildcards.mapper} > {output.tsv} 2>{output.compare}"

rule compare_alignments_category:
    input:
        gam="{root}/annotated-1/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        truth_gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
        category_list=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}.{category}.txt"),
    output:
        gam="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.{category}.gam",
        tsv="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.{category}.compared.tsv",
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.{category}.compare.txt"
    params:
        ignore_flag="--ignore 'CHM13#0#chrY'"
    threads: 24
    resources:
        mem_mb=200000,
        runtime=800,
        slurm_partition=choose_partition(800)
    shell:
        "vg gamcompare --threads 16 --range 200 {params.ignore_flag} <(vg filter --threads 4 --exact-name -N {input.category_list} {input.gam}) <(vg filter --threads 4 --exact-name -N {input.category_list} {input.truth_gam}) --output-gam {output.gam} -T -a {wildcards.mapper} > {output.tsv} 2>{output.compare}"

rule compare_alignments_not_category:
    input:
        gam="{root}/annotated-1/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        truth_gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
        category_list=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}.{category}.txt"),
    output:
        gam="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.not_{category}.gam",
        tsv="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.not_{category}.compared.tsv",
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.not_{category}.compare.txt"
    params:
        ignore_flag="--ignore 'CHM13#0#chrY'"
    threads: 24
    resources:
        mem_mb=200000,
        runtime=800,
        slurm_partition=choose_partition(800)
    shell:
        "vg gamcompare --threads 16 --range 200 {params.ignore_flag} <(vg filter --complement --threads 4 --exact-name -N {input.category_list} {input.gam}) <(vg filter --complement --threads 4 --exact-name -N {input.category_list} {input.truth_gam}) --output-gam {output.gam} -T -a {wildcards.mapper} > {output.tsv} 2>{output.compare}"

rule annotate_alignments:
    input:
        gbz=gbz,
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    output:
        gam="{root}/annotated-1/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    threads: 16
    resources:
        mem_mb=100000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg annotate -t16 -a {input.gam} -x {input.gbz} -m --search-limit=-1 >{output.gam}"
ruleorder: giraffe_sim_reads > annotate_alignments

rule annotate_alignments_against_hg:
    input:
        hg=hg,
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    output:
        gam="{root}/annotated-1/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    threads: 16
    resources:
        mem_mb=100000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg annotate -t16 -a {input.gam} -x {input.hg} -m --search-limit=-1 >{output.gam}"
ruleorder: annotate_alignments > annotate_alignments_against_hg
ruleorder: giraffe_sim_reads > annotate_alignments_against_hg

rule de_annotate_sim_alignments:
    input:
        gam="{root}/annotated-1/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    output:
        gam="{root}/aligned/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        realness="sim"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "ln {input.gam} {output.gam}"

rule correct_from_comparison:
    input:
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.compare.txt"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.correct.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.compare} | grep -o '[0-9]* reads correct' | cut -f1 -d' ' >{output.tsv}"

rule accuracy_from_comparison:
    input:
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.compare.txt"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.accuracy.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.compare} | grep -o '[0-9%.]* accuracy' | cut -f1 -d' ' | tr -d '%' >{output.tsv}"

rule eligible_from_comparison:
    input:
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.compare.txt"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.eligible.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.compare} | grep -o '[0-9]* reads eligible' | cut -f1 -d' ' >{output.tsv}"

rule wrong_from_correct_and_eligible:
    input:
        correct="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.correct.tsv",
        eligible="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.eligible.tsv"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.wrong.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"$(cat {input.eligible}) - $(cat {input.correct})\" | bc -l >{output.tsv}"


rule overall_fraction:
    input:
        number="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.{state}.tsv",
        all_comparison="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.compared.tsv",
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.overall_fraction_{state}.tsv"
    wildcard_constraints:
        state="(correct|eligible|wrong)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"$(cat {input.number}) / ($(wc -l {input.all_comparison} | cut -f1 -d' ') - 1)\" | bc -l >{output.tsv}"

rule speed_from_log_giraffe_stats:
    input:
        giraffe_log="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.speed_from_log.tsv"
    wildcard_constraints:
        realness="real",
        mapper="giraffe.*"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"$(cat {input.giraffe_log} | grep \"reads per CPU-second\" | sed \'s/Achieved \([0-9]*\.[0-9]*\) reads per CPU-second.*/\\1/g\')\" >{output.tsv}"

rule speed_from_log_giraffe:
    input:
        giraffe_log="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.speed_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="giraffe.*"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"{params.condition_name}\t$(cat {input.giraffe_log} | grep \"reads per CPU-second\" | sed \'s/Achieved \([0-9]*\.[0-9]*\) reads per CPU-second.*/\\1/g\')\" >{output.tsv}"

#Put the mapper name and memory into a tsv in experiments directory for the experiment
#This makes it easier to find for different mappers that may or may not have a refgraph in the path
rule memory_from_log_giraffe_experiment:
    input:
        giraffe_log="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="giraffe.*"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"{params.condition_name}\t$(cat {input.giraffe_log} | grep \"Memory footprint\" | sed \'s/Memory footprint: \([0-9]*\.[0-9]*\) GB.*/\\1/g\')\" >{output.tsv}"

#Put just the memory use in the stats folder for parameter search
rule memory_from_log_giraffe_stat:
    input:
        giraffe_log="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_from_log.tsv"
    wildcard_constraints:
        realness="real",
        mapper="giraffe.*"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"$(cat {input.giraffe_log} | grep \"Memory footprint\" | sed \'s/Memory footprint: \([0-9]*\.[0-9]*\) GB.*/\\1/g\')\" >{output.tsv}"

rule speed_from_log_bwa:
    input:
        bwa_log="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.speed_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="bwa(-pe)?"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        mapped_count=$(cat {input.bwa_log} | grep "Processed" | awk '{{sum+=$3}} END {{print sum}}')
        total_time=$(cat {input.bwa_log} | grep "Processed" | sed 's/.*in \([0-9]*\.[0-9]*\) CPU sec.*/\\1/g' | awk '{{sum+=$1}} END {{print sum}}')
        echo "{params.condition_name}\t$(echo "$mapped_count / $total_time" | bc -l)" >{output.tsv}
        """


rule speed_from_log_bam:
    input:
        minimap2_log="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.speed_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(minimap2-.*|winnowmap|pbmm2)"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        mapped_count=$(cat {input.minimap2_log} | grep "mapped" | awk '{{sum+=$3}} END {{print sum}}')
        total_cpu_time=$(cat {input.minimap2_log} | grep "\\[M::main\\] Real time" | sed 's/.*CPU: \([0-9]*\.[0-9]*\) sec.*/\\1/g')
        startup_cpu_time_expr=$(cat {input.minimap2_log} | grep "loaded/built the index" | sed 's/.M::main::\([0-9]*\.[0-9]*\*[0-9]*\.[0-9]*\).*/\\1 /g')
        echo "{params.condition_name}\t$(echo "$mapped_count / ($total_cpu_time - $startup_cpu_time_expr)" | bc -l)" >{output.tsv}
        """


rule speed_from_log_minigraph:
    input:
        log="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.speed_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(minigraph|panaligner)"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        mapped_count=$(cat {input.log} | grep "mapped" | awk '{{sum+=$3}} END {{print sum}}')
        total_cpu_time=$(cat {input.log} | grep "\\[M::main\\] Real time" | sed 's/.*CPU: \([0-9]*\.[0-9]*\) sec.*/\\1/g')
        startup_cpu_time_expr=$(cat {input.log} | grep "indexed the graph" | sed 's/.M::mg_index::\([0-9]*\.[0-9]*\*[0-9]*\.[0-9]*\).*/\\1 /g')
        echo "{params.condition_name}\t$(echo "$mapped_count / ($total_cpu_time - $startup_cpu_time_expr)" | bc -l)" >{output.tsv}
        """

#We need a speed_from_log.tsv file but graphaligner doesn't have a log so just make a dummy file
rule speed_from_log_graphaligner:
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.speed_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="graphaligner-.*"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        echo "{params.condition_name}\tNA" >{output.tsv}
        """

rule memory_from_log_bam:
    input:
        minimap2_log="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(minimap2-.*|winnowmap|bwa(-pe)?|pbmm2)"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"{params.condition_name}\t$(cat {input.minimap2_log} | grep \"Peak RSS\" | sed \'s/.*Peak RSS: \([0-9]*\.[0-9]*\) GB.*/\\1/g\')\" >{output.tsv}"

rule memory_from_log_minigraph:
    input:
        log="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(minigraph|panaligner)"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"{params.condition_name}\t$(cat {input.log} | grep \"Peak RSS\" | sed \'s/.*Peak RSS: \([0-9]*\.[0-9]*\) GB.*/\\1/g\')\" >{output.tsv}"


#We need memory_from_log for all mappers but graphaligner doesn't have a log so make a dummy file
rule memory_from_log_graphaligner:
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="graphaligner-.*"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"{params.condition_name}\tNA\" >{output.tsv}"


#output tsv of:
#condition name, index load time in minutes
rule index_load_time_from_log_giraffe:
    input:
        giraffe_log="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.index_load_time_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="giraffe.*"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        index_sec=$(cat {input.giraffe_log} | grep "Loading and initialization" | sed 's/Loading and initialization: \([0-9]*\.[0-9]*\) second.*/\\1/g')
        echo \"{params.condition_name}\t$(echo "$index_sec / 60" | bc -l)\" >{output.tsv}
        """

#Index load time in minutes
rule index_load_time_from_log_bam:
    input:
        minimap2_log="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.index_load_time_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(minimap2-.*|winnowmap)"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        startup_cpu_time_expr=$(cat {input.minimap2_log} | grep "loaded/built the index" | sed 's/.M::main::\([0-9]*\.[0-9]*\*[0-9]*\.[0-9]*\).*/\\1 /g')
        echo "{params.condition_name}\t$(echo "($startup_cpu_time_expr) / 60" | bc -l)" >{output.tsv}
        """
rule index_load_time_from_log_minigraph:
    input:
        log="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.log"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.index_load_time_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(minigraph|panaligner)"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        startup_cpu_time_expr=$(cat {input.log} | grep "indexed the graph" | sed 's/.M::mg_index::\([0-9]*\.[0-9]*\*[0-9]*\.[0-9]*\).*/\\1 /g')
        echo "{params.condition_name}\t$(echo "($startup_cpu_time_expr) / 60" | bc -l)" >{output.tsv}
        """

#a dummy file for the index load time since graphaligner doesn't have a log
rule index_load_time_from_log_graphaligner:
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.index_load_time_from_log.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="graphaligner-.*|bwa(-pe)?|pbmm2"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        echo "{params.condition_name}\t0" >{output.tsv}
        """

#empty time for anything not haplotype sampled
rule haplotype_sampling_time_empty:
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmednes}.{subset}.haplotype_sampling_time.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    run:
        shell("echo \"{params.condition_name}\t0\" >{output.tsv}")

rule haplotype_sampling_time_giraffe:
    input:
        kmer_counting=kmer_counts_benchmark,
        haplotype_sampling=os.path.join(GRAPHS_DIR, "indexing_benchmarks/haplotype_sampling_{refgraphbase}-{reference}{full}-{sampling}_fragmentlinked-for-{realness}-{tech}-{sample}{trimmedness}-{subset}.benchmark"),
        distance_indexing=os.path.join(GRAPHS_DIR, "indexing_benchmarks/distance_indexing_{refgraphbase}-{reference}{full}-{sampling}_fragmentlinked-for-{realness}-{tech}-{sample}{trimmedness}-{subset}.benchmark"),
        minimizer_indexing=os.path.join(GRAPHS_DIR, "indexing_benchmarks/minimizer_indexing_{refgraphbase}-{reference}{full}-{sampling}_fragmentlinked-for-{realness}-{tech}-{sample}{trimmedness}-{subset}.k{k}.w{w}{weightedness}.benchmark")
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraphbase}-{sampling}{clipping}{full}{chopping}/giraffe-k{k}.w{w}{weightedness}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.haplotype_sampling_time.tsv"
    wildcard_constraints:
        sampling="sampled[0-9]+d?o?",
        weightedness="\\.W|",
        k="[0-9]+",
        w="[0-9]+",
        # We're thinking about haplotype sampling only and these don't go into the inputs actually
        clipping="",
        chopping=""
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    run:
        runtime = 0
        for infile in [input.kmer_counting, input.haplotype_sampling, input.distance_indexing, input.minimizer_indexing]:
            f = open(infile)
            assert(f.readline().split("\t")[1] == "h:m:s")
            runtime_list = f.readline().split("\t")[1].split()
            hms_list = runtime_list[-1].split(":")
            days = 0 if len(runtime_list) == 1 else int(runtime_list[0])
            runtime = (int(hms_list[0]) * 60) + int(hms_list[1]) + (int(hms_list[2]) / 60)
            runtime += 24 * 60 * days
            f.close()
        shell("echo \"{params.condition_name}\t{runtime}\" >{output.tsv}")

ruleorder: haplotype_sampling_time_giraffe > haplotype_sampling_time_empty

#output tsv of:
#condition name, index load time plus haplotype sampling time in minutes
#input tsv have condition name and index load time / sampling time
rule index_load_and_sampling_time:
    input:
        startup_time="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.index_load_time_from_log.tsv",
        sampling_time="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.haplotype_sampling_time.tsv",
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.index_load_and_sampling_time.tsv"
    wildcard_constraints:
        realness="real"
    threads: 1
    resources:
        mem_mb=200,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        paste {input.startup_time} {input.sampling_time} | awk -v OFS='\t' '{{print $1,($2+$4)}}' >{output.tsv}
        """

# Some experiment stats can come straight from stats for the individual conditions
rule condition_experiment_stat:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{callparams}{dot}{category}.{conditionstat}.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{callparams}{dot}{category}.{conditionstat}.tsv"
    wildcard_constraints:
        refgraph="[^/_]+",
        conditionstat="((overall_fraction_)?(wrong|correct|eligible)|accuracy|(snp|indel)_(f1|precision|recall|fn|fp)|(snp|indel|total)_errors)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && cat {input.tsv} >>{output.tsv}"

rule experiment_stat_table:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{callparams}{dot}{category}.{statname}.tsv", filter_function=has_stat_filter(w["statname"]))
    output:
        table="{root}/experiments/{expname}/results/{statname}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "cat {input} >{output.table}"

ruleorder: dv_summary_table > experiment_stat_table

rule experiment_important_stats_table:
    input:
        lambda w: expand("{{root}}/experiments/{{expname}}/results/{stat}.tsv", stat=IMPORTANT_STATS_TABLE_COLUMNS)
    output:
        table="{root}/experiments/{expname}/tables/important_stats.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
        slurm_partition=choose_partition(10)
    run:
        with open(output.table, "w") as fp:
            # Write a header
            fp.write("\t".join(["#condition"] + IMPORTANT_STATS_TABLE_COLUMNS) + "\n")

        # Join all the parts
        command_parts = ["cat ", input[0]]
        for input_file in input[1:]:
            command_parts.append("| join -a 1 -a 2 -e 'N/A' -o auto - " + input_file)
        command_parts.append(" >> {output.table}")
        shell("".join(command_parts))

rule experiment_correctness_plot:
    input:
        tsv="{root}/experiments/{expname}/results/correct.tsv"
    output:
        "{root}/experiments/{expname}/plots/correct.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Correctness' --y_label 'Correct Reads' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_wrongness_plot:
    input:
        tsv="{root}/experiments/{expname}/results/wrong.tsv"
    output:
        "{root}/experiments/{expname}/plots/wrong.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Wrongness' --y_label 'Wrong Reads' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_accuracy_plot:
    input:
        tsv="{root}/experiments/{expname}/results/accuracy.tsv"
    output:
        "{root}/experiments/{expname}/plots/accuracy.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Accuracy' --y_label 'Percentage Correct of Eligible' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_overall_fraction_wrong_plot:
    input:
        tsv="{root}/experiments/{expname}/results/overall_fraction_wrong.tsv"
    output:
        "{root}/experiments/{expname}/plots/overall_fraction_wrong.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Overall Fraction Wrong' --y_label 'Fraction Eligible and Wrong' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_overall_fraction_correct_plot:
    input:
        tsv="{root}/experiments/{expname}/results/overall_fraction_correct.tsv"
    output:
        "{root}/experiments/{expname}/plots/overall_fraction_correct.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Overall Fraction Correct' --y_label 'Fraction Eligible and Correct' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_overall_fraction_eligible_plot:
    input:
        tsv="{root}/experiments/{expname}/results/overall_fraction_eligible.tsv"
    output:
        "{root}/experiments/{expname}/plots/overall_fraction_eligible.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Overall Fraction Eligible' --y_label 'Fraction Eligible' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule compared_named_from_compared:
    input:
        tsv="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv",
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compared.tsv"
    wildcard_constraints:
        refgraph="[^/_]+"
    threads: 3
    resources:
        mem_mb=1000,
        runtime=120,
        slurm_partition=choose_partition(120)
    shell:
        "printf 'correct\\tmq\\taligner\\tread\\teligible\\n' >{output.tsv} && cat {input.tsv} | grep -v '^correct' | awk -F '\\t' -v OFS='\\t' '{{ $3 = \"{params.condition_name}\"; print }}' >>{output.tsv}"


rule experiment_compared_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{callparams}{dot}{category}.compared.tsv", lambda condition: condition["realness"] == "sim")
    output:
        tsv="{root}/experiments/{expname}/results/compared.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "printf 'correct\\tmq\\taligner\\tread\\teligible\\n' >{output.tsv} && cat {input} | grep -v '^correct' >>{output.tsv}"

rule experiment_qq_plot_from_compared:
    input:
        tsv="{root}/experiments/{expname}/results/compared.tsv"
    output:
        "{root}/experiments/{expname}/plots/qq.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "Rscript plot-qq.R {input.tsv} {output}"

rule experiment_pr_plot_from_compared:
    input:
        tsv="{root}/experiments/{expname}/results/compared.tsv"
    output:
        "{root}/experiments/{expname}/plots/pr.pdf"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "Rscript plot-pr.R {input.tsv} {output}"

rule experiment_roc_plot_from_compared:
    input:
        tsv="{root}/experiments/{expname}/results/compared.tsv"
    output:
        "{root}/experiments/{expname}/plots/roc.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "Rscript plot-roc.R {input.tsv} {output}"

rule experiment_calling_plot:
    input:
        tsv="{root}/experiments/{expname}/results/{vartype}_{colname}.tsv"
    output:
        "{root}/experiments/{expname}/plots/{vartype}_{colname}.{ext}"
    wildcard_constraints:
        vartype="(snp|indel)",
        colname="(f1|precision|recall|fn|fp)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} {wildcards.vartype} {wildcards.colname}' --y_label '{wildcards.colname}' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_total_errors_plot:
    input:
        tsv="{root}/experiments/{expname}/results/total_errors.tsv"
    output:
        "{root}/experiments/{expname}/plots/total_errors.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Total Point Variant Errors' --y_label 'Total Errors' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_snp_errors_plot:
    input:
        tsv="{root}/experiments/{expname}/results/snp_errors.tsv"
    output:
        "{root}/experiments/{expname}/plots/snp_errors.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} SNP Errors' --y_label 'Errors' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_indel_errors_plot:
    input:
        tsv="{root}/experiments/{expname}/results/indel_errors.tsv"
    output:
        "{root}/experiments/{expname}/plots/indel_errors.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Indel Errors' --y_label 'Errors' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule experiment_calling_summary_plot:
    input:
        tsv="{root}/experiments/{expname}/results/{vartype}_summary.tsv"
    output:
        "{root}/experiments/{expname}/plots/{vartype}_summary.{ext}"
    wildcard_constraints:
        vartype="(dv_snp|dv_indel|sv)",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "Rscript plot-calling-results.R {input.tsv} {output}"


rule experiment_speed_from_log_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.speed_from_log.tsv", lambda condition: condition["realness"] == "real" and any(m in condition["mapper"] for m in ("giraffe", "minimap2", "winnowmap", "bwa", "pbmm2", "minigraph", "panaligner")))
    output:
        tsv="{root}/experiments/{expname}/results/speed_from_log.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cat {input} >>{output.tsv}"

ruleorder: experiment_speed_from_log_tsv > experiment_stat_table

rule experiment_speed_from_log_plot:
    input:
        tsv=rules.experiment_speed_from_log_tsv.output.tsv
    output:
        "{root}/experiments/{expname}/plots/speed_from_log.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Speed From Log' --y_label 'Speed (reads/second/thread)' --x_label 'Mapper' --x_sideways --no_n --save {output}"

rule experiment_memory_from_log_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.memory_from_log.tsv", lambda condition: condition["realness"] == "real" and any(m in condition["mapper"] for m in ("giraffe", "minimap2", "winnowmap", "bwa", "pbmm2", "minigraph", "panaligner")))
    output:
        tsv="{root}/experiments/{expname}/results/memory_from_log.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cat {input} >>{output.tsv}"

rule experiment_memory_from_log_plot:
    input:
        tsv=rules.experiment_memory_from_log_tsv.output.tsv
    output:
        "{root}/experiments/{expname}/plots/memory_from_log.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Memory From Log' --y_label 'Memory use (GB)' --x_label 'Mapper' --x_sideways --no_n --save {output}"

rule experiment_runtime_from_benchmark_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.runtime_from_benchmark.tsv", lambda condition: condition["realness"] == "real")
    output:
        tsv="{root}/experiments/{expname}/results/runtime_from_benchmark.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cat {input} >>{output.tsv}"

ruleorder: experiment_runtime_from_benchmark_tsv > experiment_stat_table

rule experiment_run_and_sampling_time_from_benchmark_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.run_and_sampling_time_from_benchmark.tsv", lambda condition: condition["realness"] == "real"),
    output:
        tsv="{root}/experiments/{expname}/results/run_and_sampling_time_from_benchmark.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cat {input} >>{output.tsv}"

ruleorder: experiment_run_and_sampling_time_from_benchmark_tsv > experiment_stat_table



rule experiment_runtime_from_benchmark_plot:
    input:
        tsv=rules.experiment_runtime_from_benchmark_tsv.output.tsv
    output:
        "{root}/experiments/{expname}/plots/runtime_from_benchmark.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Runtime From Benchmark' --y_label 'Runtime (minutes)' --x_label 'Mapper' --x_sideways --no_n --save {output}"


#This includes the time for haplotype sampling
rule experiment_index_load_and_sampling_time_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.index_load_and_sampling_time_from_log.tsv", lambda condition: condition["realness"] == "real")
    output:
        tsv="{root}/experiments/{expname}/results/index_load_time.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cat {input} >>{output.tsv}"

rule experiment_index_load_time_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.index_load_time_from_log.tsv", lambda condition: condition["realness"] == "real")
    output:
        tsv="{root}/experiments/{expname}/results/index_load_time.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cat {input} >>{output.tsv}"

#This includes the time for sampling in both the runtime and index time
#Runtime is really the total time
rule experiment_run_and_index_time_plot:
    input:
        runtime=rules.experiment_run_and_sampling_time_from_benchmark_tsv.output.tsv,
        index_time=rules.experiment_index_load_and_sampling_time_tsv.output.tsv
    output:
        "{root}/experiments/{expname}/plots/run_and_index_time.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "python3 barchart.py {input.runtime} --width 8 --height 8 --divisions {input.index_time} --title '{wildcards.expname} Runtime and Index Load Time' --y_label 'Time (minutes)' --x_label 'Mapper' --x_sideways --no_n --save {output}"

rule experiment_run_and_index_time_split_plot:
    input:
        runtime=rules.experiment_run_and_sampling_time_from_benchmark_tsv.output.tsv,
        index_time=rules.experiment_index_load_and_sampling_time_tsv.output.tsv
    output:
        slow="{root}/experiments/{expname}/plots/run_and_index_slow_time.{ext}",
        fast="{root}/experiments/{expname}/plots/run_and_index_fast_time.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    run:
        runtimes = []
        with open(input.runtime) as in_file:
            for line in in_file:
                runtimes.append(float(line.split()[1]))


        iqr = np.percentile(runtimes, 75) - np.percentile(runtimes, 25)
        cutoff = np.percentile(runtimes, 75) + (1.5 * iqr)
        bigs = list(filter(lambda x : x > cutoff, runtimes))

        if not len(bigs) == 0:

            lower_limit = min(bigs) * 0.80

            shell("python3 barchart.py {input.runtime} --width 8 --height 8 --min " + str(lower_limit) + " --divisions {input.index_time} --title '{wildcards.expname} Runtime and Index Load Time' --y_label 'Time (minutes)' --x_label 'Mapper' --x_sideways --no_n --save {output.slow}")

        smalls = list(filter(lambda x : x <= cutoff, runtimes))

        if not len(smalls) == 0:

            upper_limit = max(smalls) * 1.1

            shell("python3 barchart.py {input.runtime} --width 8 --height 8 --max " + str(upper_limit) + " --divisions {input.index_time} --title '{wildcards.expname} Runtime and Index Load Time' --y_label 'Time (minutes)' --x_label 'Mapper' --x_sideways --no_n --save {output.fast}")

rule experiment_run_and_index_time_hours_plot:
    input:
        runtime=rules.experiment_run_and_sampling_time_from_benchmark_tsv.output.tsv,
        index_time=rules.experiment_index_load_and_sampling_time_tsv.output.tsv
    output:
        "{root}/experiments/{expname}/plots/run_and_index_time_hours.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "python3 barchart.py <(awk -v OFS='\\t' '{{print $1,$2/60}}' {input.runtime}) --divisions <(awk -v OFS='\\t' '{{print $1,$2/60}}' {input.index_time}) --width 8 --height 8 --title '{wildcards.expname} Runtime and Index Load Time' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --save {output}"

#Plot the slower and faster times separately
rule experiment_run_and_index_time_split_hours_plot:
    input:
        runtime=rules.experiment_run_and_sampling_time_from_benchmark_tsv.output.tsv,
        index_time=rules.experiment_index_load_and_sampling_time_tsv.output.tsv
    output:
        slow="{root}/experiments/{expname}/plots/run_and_index_slow_time_hours.{ext}",
        fast="{root}/experiments/{expname}/plots/run_and_index_fast_time_hours.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    run:
        runtimes = []
        with open(input.runtime) as in_file:
            for line in in_file:
                runtimes.append(float(line.split()[1]))


        iqr = np.percentile(runtimes, 75) - np.percentile(runtimes, 25)
        cutoff = np.percentile(runtimes, 75) + (1.5 * iqr)
        bigs = list(filter(lambda x : x > cutoff, runtimes))

        if not len(bigs) == 0:

            lower_limit = (min(bigs) / 60.0 * 0.80)

            shell("python3 barchart.py <(awk -v OFS='\\t' '{{print $1,$2/60}}' {input.runtime}) --min " + str(lower_limit) + " --divisions <(awk -v OFS='\\t' '{{print $1,$2/60}}'      {input.index_time}) --width 8 --height 8 --title '{wildcards.expname} Runtime and Index Load Time' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --save {output.slow}")

        smalls = list(filter(lambda x : x <= cutoff, runtimes))


        if not len(smalls) == 0:

            upper_limit = (max(smalls) * 1.1) / 60

            shell("python3 barchart.py <(awk -v OFS='\\t' '{{print $1,$2/60}}' {input.runtime}) --max " + str(upper_limit) + " --divisions <(awk -v OFS='\\t' '{{print $1,$2/60}}'      {input.index_time}) --width 8 --height 8 --title '{wildcards.expname} Runtime and Index Load Time' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --save {output.fast}")


rule experiment_memory_from_benchmark_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.memory_from_benchmark.tsv", lambda condition: condition["realness"] == "real")
    output:
        tsv="{root}/experiments/{expname}/results/memory_from_benchmark.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cat {input} >>{output.tsv}"
ruleorder: experiment_memory_from_benchmark_tsv > experiment_stat_table


rule experiment_memory_from_benchmark_plot:
    input:
        tsv=rules.experiment_memory_from_benchmark_tsv.output.tsv
    output:
        "{root}/experiments/{expname}/plots/memory_from_benchmark.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Memory From Benchmark' --y_label 'Memory (GB)' --x_label 'Mapper' --x_sideways --no_n --save {output}"

#Get the accuracy from simulated reads for one condition
rule experiment_mapping_stats_sim_tsv_from_stats:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_accuracy.tsv"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_accuracy.tsv"
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="sim"
    params:
        condition_name=condition_name
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "echo \"{params.condition_name}\t$(cat {input})\" >>{output.tsv}"

#Get the accuracy from simulated reads for all conditions in the experiment
rule experiment_mapping_stats_sim_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.mapping_accuracy.tsv", lambda condition: condition["realness"] == "sim", empty_ok=True)
    output:
        tsv="{root}/experiments/{expname}/results/mapping_stats_sim.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        """
        printf "condition\tcorrect\tmapq60\twrong_mapq60\twrong_eligible\n" >> {output.tsv}
        cat {input} /dev/null >>{output.tsv}
        """

#Get the speed, memory use, and softclips from real reads for each condition
rule experiment_mapping_stats_real_tsv_from_stats:
    input:
        startup_time="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.index_load_time_from_log.tsv",
        sampling_time="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.haplotype_sampling_time.tsv",
        runtime_from_benchmark="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.runtime_from_benchmark.tsv",
        memory_from_benchmark="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_from_benchmark.tsv",
        softclipped_or_unmapped="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclipped_or_unmapped.tsv"

    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_stats_real.tsv"
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real"
    params:
        condition_name=condition_name
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        """
        echo "{params.condition_name}\t$(cat {input.runtime_from_benchmark} | cut -f 2)\t$(cat {input.startup_time} | cut -f 2)\t$(cat {input.sampling_time} | cut -f 2)\t$(cat {input.memory_from_benchmark} | cut -f 2)\t$(cat {input.softclipped_or_unmapped} | cut -f 2)" >>{output.tsv}
        """


#Get the speed, memory use, and softclips from real reads
rule experiment_mapping_stats_real_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}{dot}{category}.mapping_stats_real.tsv", lambda condition: condition["realness"] == "real", empty_ok=True)
    output:
        tsv="{root}/experiments/{expname}/results/mapping_stats_real.tsv"
    wildcard_constraints:
        realness="real"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        """
        printf "condition\truntime(min)\tstartup_time(min)\tsampling_time(min)\tmemory(GB)\tsoftclipped_or_unmapped\n" >> {output.tsv} 
        cat {input} /dev/null >>{output.tsv}
        """
ruleorder: experiment_mapping_stats_real_tsv > experiment_stat_table

#Get mapping stats from real and sim reads
rule experiment_mapping_stats_tsv:
    input:
        real="{root}/experiments/{expname}/results/mapping_stats_real.tsv",
        sim="{root}/experiments/{expname}/results/mapping_stats_sim.tsv"
    output:
        tsv="{root}/experiments/{expname}/results/mapping_stats.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        # Do a full outer join
        "join -t '\t' -a 1 -a 2 -o auto {input.sim} {input.real} >> {output.tsv}"
ruleorder: experiment_mapping_stats_tsv > experiment_stat_table

rule stats_from_alignments:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        stats="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gamstats.txt"
    threads: 16
    resources:
        mem_mb=10000,
        runtime=90,
        slurm_partition=choose_partition(90)
    shell:
        "vg stats -p {threads} -a {input.gam} >{output.stats}"

rule facts_from_alignments_with_correctness:
    input:
        gam="{root}/correctness/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        facts="{root}/stats/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.facts.txt",
        facts_dir=directory("{root}/stats/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.facts")
    threads: 64
    resources:
        mem_mb=10000,
        runtime=90,
        slurm_partition=choose_partition(90)
    run:
        vg_binary = get_vg_version(wildcards.vgversion)
        shell("python3 giraffe-facts.py --threads {threads} {input.gam} {output.facts_dir} --stage --filter-help --vg " + vg_binary + " >{output.facts}")


rule mapping_rate_from_stats:
    input:
        stats="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gamstats.txt"
    params:
        condition_name=condition_name
    output:
        rate="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_rate.tsv"
    wildcard_constraints:
        refgraph="[^/_]+"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.rate} && cat {input.stats} | grep 'Total aligned:' | cut -f2 -d':' | tr -d ' ' >>{output.rate}"

rule experiment_mapping_rate_plot:
    input:
        tsv="{root}/experiments/{expname}/results/mapping_rate.tsv"
    output:
        "{root}/experiments/{expname}/plots/mapping_rate.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Mapping Rate' --y_label 'Mapped Reads' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule unmapped_from_stats:
    input:
        stats="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gamstats.txt",
        mapped="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_rate.tsv"
    params:
        condition_name=condition_name
    output:
        unmapped="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.unmapped.tsv"
    wildcard_constraints:
        refgraph="[^/_]+"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.unmapped} && echo $(( $(cat {input.stats} | grep 'Total alignments:' | cut -f2 -d':' | tr -d ' ') - $(cut -f2 {input.mapped}) )) >>{output.unmapped}"

rule experiment_unmapped_plot:
    input:
        tsv="{root}/experiments/{expname}/results/unmapped.tsv"
    output:
        "{root}/experiments/{expname}/plots/unmapped.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Unmapped Reads' --y_label 'Unmapped Reads' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule mapping_speed_from_mean_time_used:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.time_used.mean.tsv"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_speed.tsv"
    wildcard_constraints:
        mapper="giraffe-.+"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"1 / $(cat {input.tsv})\" | bc -l >>{output.tsv}"

rule mapping_speed_from_stats:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_speed.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_speed.tsv"
    wildcard_constraints:
        refgraph="[^/_]+",
        mapper="giraffe-.+"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && cat {input.tsv} >>{output.tsv}"

rule experiment_mapping_speed_plot:
    input:
        tsv="{root}/experiments/{expname}/results/mapping_speed.tsv"
    output:
        "{root}/experiments/{expname}/plots/mapping_speed.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Speed' --y_label 'Reads per Second' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule softclips_from_mean_softclips:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips.mean.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips.tsv"
    wildcard_constraints:
        refgraph="[^/_]+"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && cat {input.tsv} >>{output.tsv}"

rule experiment_softclips_plot:
    input:
        tsv="{root}/experiments/{expname}/results/softclips.tsv"
    output:
        "{root}/experiments/{expname}/plots/softclips.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Read End Softclips' --y_label 'Mean Softclip (bp/end)' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule softclipped_or_unmapped_from_softclipped_or_unmapped:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclipped_or_unmapped.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclipped_or_unmapped.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && cat {input.tsv} >>{output.tsv}"

rule experiment_softclipped_or_unmapped_plot:
    input:
        tsv="{root}/experiments/{expname}/results/softclipped_or_unmapped.tsv"
    output:
        "{root}/experiments/{expname}/plots/softclipped_or_unmapped.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Softclipped or Unmapped Bases' --y_label 'Total (bp)' --x_label 'Condition' --x_sideways --no_n --save {output}"


rule softclipped_or_unmapped_percent:
    input:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclipped_or_unmapped.tsv",
        unmapped="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclipped_or_unmapped_percent.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        """
        printf '{params.condition_name}\\t' >{output.tsv} && echo "$(cut -f 2 {input.tsv}) * 100 / $(awk '{{SUM += $2}} END {{print SUM}}' {input.unmapped})" | bc -l  >>{output.tsv}
        """

#Plot the unmapped bases but split off the high outliers 
rule experiment_softclipped_or_unmapped_percent_split_plot:
    input:
        tsv="{root}/experiments/{expname}/results/softclipped_or_unmapped_percent.tsv"
    output:
        low="{root}/experiments/{expname}/plots/softclipped_or_unmapped_percent_low.{ext}",
        high="{root}/experiments/{expname}/plots/softclipped_or_unmapped_percent_high.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    run:
        values = []
        with open(input.tsv) as in_file:
            for line in in_file:
                values.append(float(line.split()[1]))


        iqr = np.percentile(values, 75) - np.percentile(values, 25)
        cutoff = np.percentile(values, 75) + (1.5 * iqr)
        bigs = list(filter(lambda x : x > cutoff, values))
        smalls = list(filter(lambda x : x <= cutoff, values))

        if not len(bigs) == 0:

            lower_limit = min(bigs) * 0.80

            shell("python3 barchart.py {input.tsv} --min " + str(lower_limit) + " --title '{wildcards.expname} Softclipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --save {output.low}")


        if not len(smalls) == 0:

            upper_limit = max(smalls) * 1.1

            shell("python3 barchart.py {input.tsv} --max " + str(upper_limit) + " --title '{wildcards.expname} Softclipped or unmapped bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --save {output.high}")


rule chain_coverage_from_mean_best_chain_coverage:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.mean.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.chain_coverage.tsv"
    wildcard_constraints:
        refgraph="[^/_]+"
    wildcard_constraints:
        mapper="giraffe-.+"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && cat {input.tsv} >>{output.tsv}"

rule experiment_chain_coverage_plot:
    input:
        tsv="{root}/experiments/{expname}/results/chain_coverage.tsv"
    output:
        "{root}/experiments/{expname}/plots/chain_coverage.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "python3 barchart.py {input.tsv} --width 8 --height 8 --title '{wildcards.expname} Chain Coverage' --y_label 'Best Chain Coverage (fraction)' --x_label 'Condition' --x_sideways --no_n --save {output}"


rule best_chain_coverage:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"annotation.best_chain.coverage\" {input.gam} | grep -v \"#\" >{output}"

rule chain_anchors_by_name_giraffe:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_anchors_by_name.tsv"
    wildcard_constraints:
        mapper="giraffe.*"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"name;annotation.best_chain.anchors\" {input.gam} | grep -v \"#\" >{output}"

rule chain_anchors_by_name_other:
    input:
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_anchors_by_name.tsv"
    wildcard_constraints:
        mapper="(?!giraffe).+"
    threads: 7
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        r"samtools view {input.bam} | sed 's/^\([^\t]*\).*\tcm:i:\([0-9]*\).*$/\1\t\2/' > {output}"

rule chain_anchor_bases_by_name:
    input:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_anchors_by_name.tsv"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_anchor_bases_by_name.tsv"
    params:
        minimizer_k=minimizer_k
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "awk '{{ print $1 \"\\t\" $2 * {params.minimizer_k} }}' {input} >{output}"

rule remove_names:
    input:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{allowedstat}_by_name.tsv"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{allowedstat}.tsv"
    wildcard_constraints:
        allowedstat="(best_chain_anchors|best_chain_anchor_bases)"  
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cut -f2 {input} > {output}"

rule mapq_giraffe:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapq.tsv"
    wildcard_constraints:
        mapper="giraffe.*"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"mapping_quality\" {input.gam} | grep -v \"#\" >{output}"

rule mapq_other:
    input:
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapq.tsv"
    wildcard_constraints:
        mapper="(?!giraffe).+"
    threads: 7
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        r"samtools view {input.bam} | cut -f5 > {output}"

rule time_used:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.time_used.tsv"
    threads: 10
    resources:
        mem_mb=4000,
        runtime=120,
        slurm_partition=choose_partition(120)
    shell:
        "vg filter -t {threads} -T \"time_used\" {input.gam} | grep -v \"#\" >{output}"

rule stage_time_sim:
    input:
        # Simulated reads will have provenance tracked and stage time recorded
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.stage_{stage}_time.tsv"
    wildcard_constraints:
        realness="sim"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"annotation.stage.{wildcards.stage}.time\" {input.gam} | grep -v \"#\" >{output}"

rule aligner_part_stat:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.aligner_{aligner}_{part}_{stat}.tsv"
    wildcard_constraints:
        stat="(time|bases|invocations)"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"annotation.aligner_stats.per_read.{wildcards.part}.{wildcards.stat}.{wildcards.aligner}\" {input.gam} | grep -v \"#\" >{output}"

# We also want to get the fraction of the read processed by each aligner
rule aligner_part_fraction:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.aligner_{aligner}_{part}_fraction.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"annotation.aligner_stats.per_read.{wildcards.part}.bases.{wildcards.aligner};length\" {input.gam} | grep -v \"#\" | awk '{{ print $1 / $2 }}' >{output}"

# And the speed in bases per second
rule aligner_part_speed:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.aligner_{aligner}_{part}_speed.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"annotation.aligner_stats.per_read.{wildcards.part}.bases.{wildcards.aligner};annotation.aligner_stats.per_read.{wildcards.part}.time.{wildcards.aligner}\" {input.gam} | grep -v \"#\" | grep -v \"^0.0$\" | awk '{{ print $1 / $2 }}' >{output}"

# And the problem sizes
rule aligner_part_probsize:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.aligner_{aligner}_{part}_probsize.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"annotation.aligner_stats.per_read.{wildcards.part}.bases.{wildcards.aligner};annotation.aligner_stats.per_read.{wildcards.part}.invocations.{wildcards.aligner}\" {input.gam} | grep -v \"#\" | grep -v \"^0.0$\" | awk '{{ print $1 / $2 }}' >{output}"

rule length:
    input:
        gam="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"length\" {input.gam} >{output}"

# tsv of "mapped"/"unmapped" and the length of the original read
rule length_by_mapping:
    input:
        fastq=fastq,
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv",
        length_by_name=temp("{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_name.tsv"),
        mapped_names=temp("{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapped_names.tsv")
    threads: 2
    resources:
        mem_mb=28000,
        runtime=720,
        slurm_partition=choose_partition(720)
    shell:
        """
        vg filter -t {threads} -T "name" {input.gam} | grep -v "#" | sort -k 1b,1 | sed 's/\/[1-2]\$//g' > {output.mapped_names}
        seqkit fx2tab -n -l {input.fastq} | sort -k 1b,1 | awk -v OFS='\t' '{{print $1,$NF}}' > {output.length_by_name}
        join -j 1 -a 2 {output.length_by_name} {output.mapped_names} | awk -v OFS='\t' '{{if (NF == 1) {{print "unmapped",$1}} else {{print "mapped",$2}}}}' > {output.tsv}
        """


rule unmapped_length:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.unmapped_length.tsv"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cat {input.tsv} | {{ grep '^unmapped' || test $? = 1; }} | cut -f2 >{output}"

rule length_by_correctness:
    input:
        gam="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_correctness.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=120,
        slurm_partition=choose_partition(120)
    shell:
        "vg filter -t {threads} -T \"correctness;length\" {input.gam} | grep -v \"#\" > {output}"

rule effective_mismatch_rate_by_correctness:
    input:
        gam="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.effective_mismatch_rate_by_correctness.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=120,
        slurm_partition=choose_partition(120)
    shell:
        "vg filter -t {threads} -T \"correctness;annotation.effective_mismatch_rate\" {input.gam} | grep -v \"#\" > {output}"

rule mapq_by_correctness:
    input:
        gam="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapq_by_correctness.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=120,
        slurm_partition=choose_partition(120)
    shell:
        "vg filter -t {threads} -T \"correctness;mapping_quality\" {input.gam} | grep -v \"#\" > {output}"

rule softclips_by_name_gam:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips_by_name.tsv"
    wildcard_constraints:
        mapper="(minigraph|graphaligner-.*|giraffe.*)"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"name;softclip_start;softclip_end\" {input.gam} | grep -v \"#\" > {output}"

#Graphaligner doesn't always map the entire read. This outputs a tsv of the read name and the number of bases that didn't get put in the final alignment
rule unmapped_ends_by_name:
    input:
        fastq=fastq,
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.unmapped_ends_by_name.tsv",
        read_length_by_name=temp("{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.read_length_by_name.tsv"),
        mapped_length_by_name=temp("{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapped_length_by_name.tsv")
    threads: 2
    resources:
        mem_mb=28000,
        runtime=720,
        slurm_partition=choose_partition(720)
    shell:
        """
        vg filter -t {threads} -T "name;length" {input.gam} | grep -v "#" | sort -k 1b,1 | sed 's/\/[1-2]\$//g' > {output.mapped_length_by_name}
        seqkit fx2tab -n -l {input.fastq} | sort -k 1b,1 | awk -v OFS='\t' '{{print $1,$NF}}' > {output.read_length_by_name}
        join -a 2 {output.read_length_by_name} {output.mapped_length_by_name} | awk -v OFS='\t' '{{print $1,$2-$3}}' > {output.tsv}
        """


rule softclips_by_name_other:
    input:
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips_by_name.tsv"
    wildcard_constraints:
        mapper="(?!(giraffe|graphaligner)).+"
    threads: 7
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        r"samtools view {input.bam} | cut -f1,2,6 | sed 's/\t\(\([0-9]*\)S\)\?\([0-9]*[IDMH=X]\|\*\)*\(\([0-9]*\)S\)\?$/\t\2\t\5/g' | sed 's/\t\t/\t0\t/g' | sed 's/\t$/\t0/g' | sed 's/16\t\([0-9]*\)\t\([0-9]*\)/\2\t\1/g' | sed 's/\t[0-9]\+\t\([0-9]*\t[0-9]*\)$/\t\1/g' > {output}"

ruleorder: softclips_by_name_gam > softclips_by_name_other

rule softclips:
    input:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips_by_name.tsv"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        r"sed 's/^.*\t\([0-9]*\)\t\([0-9]*\)$/\1\n\2/' {input} > {output}"

rule unmapped_ends:
    input:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.unmapped_ends_by_name.tsv"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.unmapped_ends.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        r"sed 's/^.*\t\([0-9]*\)$/\1/' {input} > {output}"

rule softclipped_or_unmapped:
    input:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips.total.tsv",
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.unmapped_length.total.tsv",
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.unmapped_ends.total.tsv"
    output:
         "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclipped_or_unmapped.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    run:
        # Sum the one-column TSVs
        total = 0
        for file in input:
            for line in open(file):
                line = line.strip()
                if line:
                    total += float(line)
        with open(output[0], "w") as f:
            f.write(f"{total}\n")
        

#Getting this for short reads is super slow and we don't really care so skip it
rule softclipped_or_unmapped_empty:
    output:
         "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclipped_or_unmapped.tsv"
    wildcard_constraints:
        tech="illumina"
    threads: 1
    resources:
        mem_mb=400,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        """
        echo "0" >{output}
        """
ruleorder: softclipped_or_unmapped_empty > softclipped_or_unmapped



rule memory_usage_gam:
    input:
        bench="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_usage.tsv"
    wildcard_constraints:
        realness="real",
        mapper="(giraffe.+|graphaligner-.*|minigraph|panaligner)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        # max_rss happens to be column 3, but try to check
        "cat {input.bench} | cut -d \"\\t\" -f3 | grep -A1 max_rss | tail -n1 >{output.tsv}"

rule memory_usage_sam:
    input:
        bench="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_usage.tsv"
    wildcard_constraints:
        realness="real",
        mapper="(minimap2.+|winnowmap|bwa(-pe)?|pbmm2)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        # max_rss happens to be column 3, but try to check
        "cat {input.bench} | cut -d \"\\t\" -f3 | grep -A1 max_rss | tail -n1 >{output.tsv}"

rule runtime_from_benchmark_bam:
    input:
        bench="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.runtime_from_benchmark.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(minimap2.+|winnowmap|bwa(-pe)?|pbmm2)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    run:
        f = open(input.bench)
        assert(f.readline().split("\t")[1] == "h:m:s")
        runtime_list = f.readline().split("\t")[1].split()
        hms_list = runtime_list[-1].split(":")
        days = 0 if len(runtime_list) == 1 else int(runtime_list[0])
        runtime = (int(hms_list[0]) * 60) + int(hms_list[1]) + (int(hms_list[2]) / 60)
        runtime += 24 * 60 * days
        f.close()

        shell("echo \"{params.condition_name}\t{runtime}\" >{output.tsv}")

rule runtime_from_benchmark_gam:
    input:
        bench="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.runtime_from_benchmark.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(giraffe.*|graphaligner-.*|minigraph|panaligner)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    run:
        f = open(input.bench)
        assert(f.readline().split("\t")[1] == "h:m:s")
        runtime_list = f.readline().split("\t")[1].split()
        hms_list = runtime_list[-1].split(":")
        days = 0 if len(runtime_list) == 1 else int(runtime_list[0])
        runtime = (int(hms_list[0]) * 60) + int(hms_list[1]) + (int(hms_list[2]) / 60)
        runtime += 24 * 60 * days
        f.close()

        shell("echo \"{params.condition_name}\t{runtime}\" >{output.tsv}")

rule run_and_sampling_time_from_benchmark:
    input:
        runtime="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.runtime_from_benchmark.tsv",
        sampling="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.haplotype_sampling_time.tsv"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.run_and_sampling_time_from_benchmark.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        """
        paste {input.runtime} {input.sampling} | awk -v OFS='\t' '{{print $1,($2+$4)}}' >{output.tsv}
        """

rule memory_from_benchmark_sam:
    input:
        bench="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_from_benchmark.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(minimap2.+|winnowmap|bwa(-pe)?|pbmm2)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    run:
        f = open(input.bench)
        assert(f.readline().split("\t")[2] == "max_rss")
        memory = float(f.readline().split("\t")[2]) / 1000
        f.close()

        shell("echo \"{params.condition_name}\t{memory}\" >{output.tsv}")

rule memory_from_benchmark_gam:
    input:
        bench="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_from_benchmark.tsv"
    params:
        condition_name=condition_name
    wildcard_constraints:
        refgraph="[^/_]+",
        realness="real",
        mapper="(giraffe.*|graphaligner-.*|minigraph|panaligner)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    run:
        f = open(input.bench)
        assert(f.readline().split("\t")[2] == "max_rss")
        memory = float(f.readline().split("\t")[2]) / 1000
        f.close()

        shell("echo \"{params.condition_name}\t{memory}\" >{output.tsv}")

rule mean_stat:
    input:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.mean.tsv"
    threads: 1
    resources:
        mem_mb=512,
        runtime=20,
        slurm_partition=choose_partition(20)
    run:
        # Average the one-column TSV
        total = 0
        count = 0
        for line in open(input[0]):
            line = line.strip()
            if line:
                if line != "null":
                    total += float(line)
                count += 1
        with open(output[0], "w") as f:
            f.write(f"{total/count}\n")

rule total_stat:
    input:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.total.tsv"
    threads: 1
    resources:
        mem_mb=512,
        runtime=20,
        slurm_partition=choose_partition(20)
    run:
        # Sum the one-column TSV
        total = 0
        for line in open(input[0]):
            line = line.strip()
            if line:
                total += float(line)
        with open(output[0], "w") as f:
            f.write(f"{total}\n")

rule wfa_portion:
    input:
        lambda w: expand("{{root}}/stats/{{reference}}/{{refgraph}}/{{mapper}}/{{realness}}/{{tech}}/{{sample}}{{trimmedness}}.{{subset}}.aligner_{alignerpart}_{{stat}}.total.tsv", alignerpart=[p for p in ALIGNER_PARTS if p.endswith(w["tailmiddle"])])
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.wfa_{tailmiddle}_portion_of_{stat}.tsv"
    threads: 1
    resources:
        mem_mb=512,
        runtime=20,
        slurm_partition=choose_partition(20)
    run:
        # Make a TSV with just the fraction of this stat's total that belongs
        # to WFA for tail or middle, whichever this is.
        aligner_parts = [p for p in ALIGNER_PARTS if p.endswith(wildcards["tailmiddle"])]
        total = 0
        wfa = 0
        for (aligner_part, filename) in zip(aligner_parts, input):
            value = float(open(filename).read().strip())
            total += value
            if aligner_part.startswith("wfa"):
                wfa += value
        
        with open(output[0], "w") as out_stream:
            out_stream.write(f"{wfa/total}\n")


rule average_stage_time_table:
    input:
        # Input files must be in the same order as stages, which is dynamic. But we can't use params here.
        lambda w: expand("{{root}}/stats/{{reference}}/{{refgraph}}/{{mapper}}/{{realness}}/{{tech}}/{{sample}}{{trimmedness}}.{{subset}}.stage_{stage}_time.mean.tsv", stage=mapper_stages(w))
    output:
        "{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_stage_time.tsv"
    params:
        mapper_stages=mapper_stages
    threads: 1
    resources:
        mem_mb=512,
        runtime=20,
        slurm_partition=choose_partition(20)
    run:
        # Make a TSV of stage name and its average value
        with open(output[0], "w") as out_stream:
            for (stage, filename) in zip(params.mapper_stages, input):
                out_stream.write(f"{stage}\t{open(filename).read().strip()}\n")

rule combine_aligner_stat_table:
    input:
        expand("{{root}}/stats/{{reference}}/{{refgraph}}/{{mapper}}/{{realness}}/{{tech}}/{{sample}}{{trimmedness}}.{{subset}}.aligner_{alignerpart}_{{stat}}.mean.tsv", alignerpart=ALIGNER_PARTS)
    output:
        "{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_aligner_{stat}.tsv"
    threads: 1
    resources:
        mem_mb=512,
        runtime=20,
        slurm_partition=choose_partition(20)
    run:
        # Make a TSV of aligner and part name and its average value
        with open(output[0], "w") as out_stream:
            for (alignerpart, filename) in zip(ALIGNER_PARTS, input):
                out_stream.write(f"{alignerpart}\t{open(filename).read().strip()}\n")

rule total_aligner_stat_table:
    input:
        expand("{{root}}/stats/{{reference}}/{{refgraph}}/{{mapper}}/{{realness}}/{{tech}}/{{sample}}{{trimmedness}}.{{subset}}.aligner_{alignerpart}_{{stat}}.total.tsv", alignerpart=ALIGNER_PARTS)
    output:
        "{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.total_aligner_{stat}.tsv"
    threads: 1
    resources:
        mem_mb=512,
        runtime=20,
        slurm_partition=choose_partition(20)
    run:
        # Make a TSV of aligner and part name and its total value
        with open(output[0], "w") as out_stream:
            for (alignerpart, filename) in zip(ALIGNER_PARTS, input):
                out_stream.write(f"{alignerpart}\t{open(filename).read().strip()}\n")

rule best_chain_coverage_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/best_chain_coverage-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title '{wildcards.tech} {wildcards.realness} Fraction Covered' --y_label 'Items' --x_label 'Coverage' --no_n --save {output}"

rule time_used_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.time_used.tsv",
        mean="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.time_used.mean.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/time_used-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --log_counts --bins 100 --title \"{wildcards.tech} {wildcards.realness} Time Used, Mean=$(cat {input.mean})\" --y_label 'Items' --x_label 'Time (s)' --no_n --save {output}"

rule stage_time_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.stage_{stage}_time.tsv",
        mean="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.stage_{stage}_time.mean.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/stage_{stage}_time-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title \"{wildcards.tech} {wildcards.realness} Stage {wildcards.stage} Time, Mean=$(cat {input.mean})\" --y_label 'Items' --x_label 'Time (s)' --no_n --save {output}"

rule average_stage_time_barchart:
    input:
        tsv="{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_stage_time.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/average_stage_time-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    params:
        mapper_stages=mapper_stages
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {params.mapper_stages} --title '{wildcards.tech} {wildcards.realness} Mean Stage Times' --y_label 'Time (s)' --x_label 'Stage' --no_n --save {output}"

rule average_aligner_time_barchart:
    input:
        tsv="{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_aligner_time.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_time-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {ALIGNER_PARTS} --title '{wildcards.tech} {wildcards.realness} Mean Aligner Times' --y_label 'Time per Read (s)' --x_label 'Aligner and Part' --no_n --save {output}"

rule total_aligner_time_barchart:
    input:
        tsv="{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.total_aligner_time.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/total_aligner_time-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {ALIGNER_PARTS} --title '{wildcards.tech} {wildcards.realness} Total Aligner Times' --y_label 'Time (s)' --x_label 'Aligner and Part' --no_n --save {output}"


rule average_aligner_bases_barchart:
    input:
        tsv="{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_aligner_bases.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_bases-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {ALIGNER_PARTS} --title '{wildcards.tech} {wildcards.realness} Mean Aligner Bases' --y_label 'Bases Processed per Read (bp)' --x_label 'Aligner and Part' --no_n --save {output}"

rule average_aligner_invocations_barchart:
    input:
        tsv="{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_aligner_invocations.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_invocations-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {ALIGNER_PARTS} --title '{wildcards.tech} {wildcards.realness} Mean Aligner Invocations' --y_label 'Invocations per Read (count)' --x_label 'Aligner and Part' --no_n --save {output}"

rule average_aligner_fraction_barchart:
    input:
        tsv="{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_aligner_fraction.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_fraction-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {ALIGNER_PARTS} --title '{wildcards.tech} {wildcards.realness} Mean Aligner Fraction' --y_label 'Read Processed (fraction)' --x_label 'Aligner and Part' --no_n --save {output}"

rule average_aligner_speed_barchart:
    input:
        tsv="{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_aligner_speed.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_speed-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {ALIGNER_PARTS} --title '{wildcards.tech} {wildcards.realness} Mean Aligner Speed' --y_label 'Average of Read Speeds (bp/s)' --x_label 'Aligner and Part' --no_n --save {output}"

rule average_aligner_probsize_barchart:
    input:
        tsv="{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_aligner_probsize.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_probsize-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {ALIGNER_PARTS} --title '{wildcards.tech} {wildcards.realness} Mean Aligner Problem Size' --y_label 'Average of Read Average Problem Sizes (bp)' --x_label 'Aligner and Part' --no_n --save {output}"

rule length_by_mapping_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv",
        mean="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length.mean.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/length_by_mapping-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title \"{wildcards.tech} {wildcards.realness} Read Length, Mean=$(cat {input.mean})\" --y_label 'Items' --x_label 'Length (bp)' --no_n --categories mapped unmapped --category_labels Mapped Unmapped --legend_overlay 'best' --save {output}"


rule length_by_correctness_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_correctness.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/length_by_correctness-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --log_counts --bins 100 --title '{wildcards.tech} {wildcards.realness} Read Length for {wildcards.mapper}' --y_label 'Items' --x_label 'Length (bp)' --no_n --categories correct incorrect off-reference --category_labels Correct Incorrect 'Off Reference' --legend_overlay 'best' --stack --save {output}"

rule effective_mismatch_rate_by_correctness_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.effective_mismatch_rate_by_correctness.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/effective_mismatch_rate_by_correctness-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --normalize --line --bins 100 --title '{wildcards.tech} {wildcards.realness} Effective Mismatch Rate for {wildcards.mapper}' --y_label 'Density' --x_label 'Effective Mismatch Rate' --no_n --categories correct incorrect off-reference --category_labels Correct Incorrect 'Off Reference' --legend_overlay 'best' --save {output}"


rule softclips_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips.tsv",
        mean="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips.mean.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/softclips-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title \"{wildcards.tech} {wildcards.realness} Softclip Length, Mean=$(cat {input.mean})\" --y_label 'Ends' --x_label 'Softclip Length (bp)' --no_n --log_counts --save {output}"

rule mapq_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapq.tsv",
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/mapq-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --log_counts --bins 61 --title \"{wildcards.tech} {wildcards.realness} Mapping Quality)\" --y_label 'Items' --x_label 'MAPQ (Phred)' --no_n --save {output}"

rule mapq_by_correctness_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapq_by_correctness.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/mapq_by_correctness-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --log_counts --bins 61 --title '{wildcards.tech} {wildcards.realness} MAPQ for {wildcards.mapper}' --y_label 'Items' --x_label 'Length (bp)' --no_n --categories correct incorrect off-reference --category_labels Correct Incorrect 'Off Reference' --legend_overlay 'best' --line --save {output}"


rule mapping_accuracy:
    input:
        compared_tsv="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.mapping_accuracy.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=100,
        slurm_partition=choose_partition(100)
    run:

        correct_count = 0
        mapq60_count = 0
        wrong_mapq60_count = 0
        wrong_eligible_count = 0

        f = open(input.compared_tsv)
        f.readline()
        for line in f:
            l = line.split()
            if l[0] == "1":
                correct_count+=1
            if int(l[1]) == 60:
                mapq60_count += 1
                if int(l[0]) == 0 and int(l[4]) == 1:
                    wrong_mapq60_count+=1
            if int(l[0]) == 0 and int(l[4]) == 1:
                wrong_eligible_count+=1
        f.close()
        shell("printf \"" + str(correct_count) + "\t" + str(mapq60_count) + "\t" + str(wrong_mapq60_count) + "\t" + str(wrong_eligible_count) + "\" > {output}")
        

rule parameter_search_mapping_stats:
    input:
        times = lambda w: param_search_tsvs(w, "time_used.mean"),
        speed_log = lambda w: param_search_tsvs(w, "speed_from_log"),
        memory = lambda w: param_search_tsvs(w, "memory_from_log"),
        softclips = lambda w : param_search_tsvs(w, "softclips.mean"),
        unmapped = lambda w : param_search_tsvs(w, "softclipped_or_unmapped"),
        mapping_accuracy = expand("{{root}}/stats/{{reference}}/{{refgraph}}/giraffe-{{minparams}}-{{preset}}-{{vgversion}}-{param_hash}/sim/{{tech}}/{{sample}}{{trimmedness}}.{{subset}}.mapping_accuracy.tsv",param_hash=PARAM_SEARCH.get_hashes())
    output:
        outfile="{root}/parameter_search/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.parameter_mapping_stats.tsv"
    log: "{root}/parameter_search/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.param_search_mapping_stats.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=100,
        slurm_partition=choose_partition(100)
    run:
        f = open(output.outfile, "w")
        f.write("#correct\tmapq60\twrong_mapq60\tsoftclips\tsoftclipped_or_unmapped\tspeed(r/s/t)\tspeed_from_log\tmemory(GB)\t" + '\t'.join([param.name for param in PARAM_SEARCH.parameters]))
        for param_hash, stats_file, times_file, speed_file, memory_file, softclips_file, unmapped_file in zip(PARAM_SEARCH.get_hashes(), input.mapping_accuracy, input.times, input.speed_log, input.memory, input.softclips, input.unmapped):


            param_f = open(stats_file)
            l = param_f.readline().split()
            correct_count = l[0]
            mapq60_count = l[1]
            wrong_mapq60_count = l[2]
            param_f.close()

            softclips_f = open(softclips_file)
            l = softclips_f.readline().split()
            softclips = l[0]
            softclips_f.close()

            unmapped_f = open(unmapped_file)
            l = unmapped_f.readline().split()
            unmapped = l[0]
            unmapped_f.close()

            time_f = open(times_file)
            l = time_f.readline().split()
            speed = str(1/float(l[0]))
            time_f.close()

            speed_f = open(speed_file)
            l = speed_f.readline().split()
            log_speed=l[0]
            speed_f.close()

            memory_f = open(memory_file)
            l = memory_f.readline().split()
            memory=l[0]
            memory_f.close()

            parameters = PARAM_SEARCH.hash_to_parameters[param_hash]
            f.write("\n" + correct_count + "\t" + mapq60_count + "\t" + wrong_mapq60_count + "\t" + softclips + "\t" + unmapped + "\t" + speed + "\t" + log_speed + "\t" + memory + "\t" + '\t'.join([str(x) for x in parameters])) 
        f.close()

rule plot_correct_speed_vs_parameter:
    input:
        tsv = "{root}/parameter_search/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.parameter_mapping_stats.tsv"
    output:
        plot = "{root}/parameter_search/plots/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.correct_speed_vs_{parameter}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    run:
        infile = open(input.tsv)
        header = infile.readline().split()
        parameter_col = str(header.index(wildcards.parameter)+1) 
        infile.close()
        shell("cat <(cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $6}}' | sed 's/^/RPS /g') <(cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $1}}' | sed 's/^/Correct /g') | ./scatter.py --title 'Speed vs Correctness vs " + wildcards.parameter + "' --x_label " + wildcards.parameter + "  --y_per_category --categories 'RPS' 'Correct' --y_label 'Reads per Second' 'Reads Correct' --legend_overlay 'best' --save {output.plot} /dev/stdin")


rule parameter_search_stat:
    input:
        stat_files=lambda w: param_search_tsvs(w, w["statname"], w["realness"]),
    output:
        outfile="{root}/parameter_search/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.{realness}.parameter_{statname}.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=100,
        slurm_partition=choose_partition(100)
    run:
        f = open(output.outfile, "w")
        f.write("#" + wildcards["statname"] + "\t" + '\t'.join([param.name for param in PARAM_SEARCH.parameters]))
        for param_hash, stat_file in zip(PARAM_SEARCH.get_hashes(), input.stat_files):

            param_f = open(stat_file)
            l = param_f.readline().split()
            stat_value = l[0]
            param_f.close()

            parameters = PARAM_SEARCH.hash_to_parameters[param_hash]
            f.write("\n" + stat_value + "\t" + '\t'.join([str(x) for x in parameters])) 
        f.close()

rule plot_stat_vs_parameter:
    input:
        tsv = "{root}/parameter_search/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.{realness}.parameter_{statname}.tsv"
    output:
        plot = "{root}/parameter_search/plots/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.{realness}.{statname}_vs_{parameter}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    run:
        infile = open(input.tsv)
        header = infile.readline().split()
        parameter_col = str(header.index(wildcards.parameter)+1) 
        infile.close()
        # TODO: Aren't wildcards available here with {}?
        shell("cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $1}}' | ./scatter.py --title '" + wildcards.statname + " vs. " + wildcards.parameter + "' --x_label " + wildcards.parameter + " --y_label '" + wildcards.statname + "' --legend_overlay 'best' --save {output.plot} /dev/stdin")

rule parameter_search_parametric_stats:
    input:
        stat_files_y=lambda w: param_search_tsvs(w, w["statnamey"], w["realnessy"]),
        stat_files_x=lambda w: param_search_tsvs(w, w["statnamex"], w["realnessx"]),
    output:
        outfile="{root}/parameter_search/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.parametric.{realnessy}.{statnamey}_vs_{realnessx}.{statnamex}.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=100,
        slurm_partition=choose_partition(100)
    run:
        f = open(output.outfile, "w")
        f.write("#" + wildcards["statnamex"] + "\t" + wildcards["statnamey"] + "\t" + '\t'.join([param.name for param in PARAM_SEARCH.parameters]))
        for param_hash, stat_file_x, stat_file_y in zip(PARAM_SEARCH.get_hashes(), input.stat_files_x, input.stat_files_y):

            param_f = open(stat_file_x)
            l = param_f.readline().split()
            stat_value_x = l[0]
            param_f.close()

            param_f = open(stat_file_y)
            l = param_f.readline().split()
            stat_value_y = l[0]
            param_f.close()

            parameters = PARAM_SEARCH.hash_to_parameters[param_hash]
            f.write("\n" + stat_value_x + "\t" + stat_value_y + "\t" + '\t'.join([str(x) for x in parameters])) 
        f.close()

rule plot_stat_vs_stat:
    input:
        tsv = "{root}/parameter_search/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.parametric.{realnessy}.{statnamey}_vs_{realnessx}.{statnamex}.tsv"
    output:
        plot = "{root}/parameter_search/plots/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}/{sample}{trimmedness}.{subset}/{tech}.parametric.{realnessy}.{statnamey}_vs_{realnessx}.{statnamex}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "cat {input.tsv} | grep -v '#' | cut -f1-2 | ./scatter.py --title '{wildcards.statnamey} vs. {wildcards.statnamex}' --x_label '{wildcards.statnamex}' --y_label '{wildcards.statnamey}' --save {output.plot} /dev/stdin"

rule chain_anchors_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_anchors.tsv",
        mean="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_anchors.mean.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/chain_anchors-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title \"{wildcards.tech} {wildcards.realness} Chained Anchors, Mean=$(cat {input.mean})\" --y_label 'Reads' --x_label 'Chained Anchors' --no_n --save {output}"

rule chain_anchor_bases_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_anchor_bases.tsv",
        mean="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_anchor_bases.mean.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/chain_anchor_length-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title \"{wildcards.tech} {wildcards.realness} Chain Anchor Length, Mean=$(cat {input.mean})\" --y_label 'Reads' --x_label 'Chained Anchor Length (bp)' --no_n --save {output}"

rule chain_coverage_histogram:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.tsv",
        mean="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.mean.tsv"
    output:
        "{root}/plots/{reference}/{refgraph}/{mapper}/chain_coverage-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    wildcard_constraints:
        mapper="giraffe.*"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title \"{wildcards.tech} {wildcards.realness} Best Chain Coverage, Mean=$(cat {input.mean})\" --y_label 'Reads' --x_label 'Best Chain Coverage' --no_n --save {output}"


rule add_mapper_to_plot:
    input:
        "{root}/plots/{reference}/{refgraph}/{mapper}/{plotname}-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    output:
        "{root}/plots/{reference}/{refgraph}/{plotname}-for-{mapper}-on-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "cp {input} {output}"

rule latex_table_from_tsv:
    input:
        tsv="{root}/experiments/{expname}/results/{table_name}.tsv",
    output:
        tsv="{root}/experiments/{expname}/results/{table_name}.latex.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "cat {input.tsv} | sed 's/\t/\t\& /g' > {output.tsv}"

##########################################################################
## Genotyping with vg call
rule vgpack:
    input:
        gam='{root}/aligned/{reference}/{refgraph}/{mapper}/real/{tech}/{sample}{trimmedness}.{subset}.gam',
        gbz=gbz
    output:
        pack='{root}/svcall/vgcall/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.pack'
    benchmark: '{root}/svcall/vgcall/{mapper}/benchmark.call.vgcall_pack.{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.tsv'
    resources:
        mem_mb=300000,
        runtime=1400,
        slurm_partition=choose_partition(1400)
    threads: 8
    run:
        gam_path = os.path.abspath(input.gam)
        gbz_path = os.path.abspath(input.gbz)
        out_path = os.path.abspath(output.pack)
        shell("cd {LARGE_TEMP_DIR} && vg pack -e -x " + gbz_path + " -o " + out_path + " -g " + gam_path + " -Q 5 -t {threads}")

rule vgcall:
    input:
        pack='{root}/svcall/vgcall/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.pack',
        graph=gbz,
        snarls=snarls
    output:
        vcf='{root}/svcall/vgcall/unpopped/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.called_on_{truthref}.vcf.gz'
    benchmark: '{root}/svcall/vgcall/unpopped/{mapper}/benchmark.call.vgcall_call.{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.called_on_{truthref}.tsv'
    log:
        logfile='{root}/svcall/vgcall/unpopped/{mapper}/log.call.vgcall_call.{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.called_on_{truthref}.log'
    wildcard_constraints:
        mapper="giraffe.*|graphaligner.*"
    threads: 4
    params:
        reference_sample=reference_sample
    resources:
        mem_mb=160000,
        runtime=360,
        slurm_partition=choose_partition(360)
    run:
        pack_path = os.path.abspath(input.pack)
        graph_path = os.path.abspath(input.graph)
        snarls_path = os.path.abspath(input.snarls)
        out_path = os.path.abspath(output.vcf)
        log_path = os.path.abspath(log.logfile)
        shell("cd {LARGE_TEMP_DIR} && vg call -Az -s {wildcards.sample} -S {params.reference_sample} -c 30 -k " + pack_path + " -r " + snarls_path + " -t {threads} " + graph_path + " | gzip > " + out_path + " 2> " + log_path)

rule decompose_nested_variants:
    input:
        vcf='{root}/svcall/vgcall/unpopped/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.called_on_{truthref}.vcf.gz'
    output:
        vcf='{root}/svcall/vgcall/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.called_on_{truthref}.vcf.gz',
        temp_vcf=temp('{root}/svcall/vgcall/unpopped/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.called_on_{truthref}.temp.vcf.gz'),
        temp_tbi=temp('{root}/svcall/vgcall/unpopped/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.called_on_{truthref}.temp.vcf.gz.tbi'),

    wildcard_constraints:
        mapper="giraffe.*|graphaligner.*"
    threads: 2
    resources:
        mem_mb=100000,
        runtime=600,
        slurm_partition=choose_partition(600)
    container: 'docker://quay.io/jmonlong/vcfwave-vcfbub:1.0.12_0.1.1'
    shell:
        """
        export TMPDIR={wildcards.root}/temp
        zcat {input.vcf} | bgzip > {output.temp_vcf}
        tabix -p vcf {output.temp_vcf}
        vcfbub -l 0 -a 100000 --input {output.temp_vcf} | vcfwave -I 1000 | bcftools sort -O z -o {output.vcf}
        """

# Call SVs with linear caller sniffles
# Copied from Jean
rule call_svs_sniffles:
    input:
        bam="{root}/aligned/{reference}/{refgraph}/{mapper}/real/{tech}/{sample}{trimmedness}.{subset}.sorted.bam",
        bai="{root}/aligned/{reference}/{refgraph}/{mapper}/real/{tech}/{sample}{trimmedness}.{subset}.sorted.bam.bai",
        vntr=sv_calling_vntr_bed,
        ref_fasta=reference_fasta
    output: "{root}/svcall/sniffles/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.called_on_{reference}.vcf.gz"
    threads: 16
    log: "{root}/svcall/sniffles/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.log"
    benchmark: "{root}/svcall/sniffles/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.benchmark"
    resources:
        mem_mb=32000,
        runtime=120,
        slurm_partition=choose_partition(120)
    container: 'docker://quay.io/biocontainers/sniffles:2.5.3--pyhdfd78af_0'
    shell:
        """
        sniffles -i {input.bam} -v {output}.temp.vcf.gz -t {threads} --tandem-repeats {input.vntr} --reference {input.ref_fasta} 2>&1 > {log}
        zcat {output}.temp.vcf.gz | sed 's/CHM13#0#//g' | sed 's/GRCh38#0#//g' | gzip >{output}
        rm {output}.temp.vcf.gz
        """

rule truvari:
    input:
        truth_vcf=SV_DATA_DIR+'/{truthset}_{truthref}.vcf.gz',
        # TODO: Rename caller to svcaller to avoid confusing it with the point variant caller (always DV)
        sample_vcf='{root}/svcall/{caller}/{mapper}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.called_on_{truthref}.vcf.gz',
        confreg=SV_DATA_DIR+'/{truthset}_{truthref}.confreg.bed',
        ref=calling_reference_fasta
    output:
        summary="{root}/svcall/{caller}/{mapper}/eval/{truthset}_{truthref}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.truvari.summary.json",
        refine_var="{root}/svcall/{caller}/{mapper}/eval/{truthset}_{truthref}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.truvari.refine.variant_summary.json",
        refine_reg="{root}/svcall/{caller}/{mapper}/eval/{truthset}_{truthref}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}.truvari.refine.region_summary.json",
        # These aren't actual outputs but are scratch files which we rely on Snakemake to manage cleanup for
        odir=temp(directory('{root}/temp/{caller}.{mapper}.truvari_{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}')),
        bvcf=temp('{root}/temp/{caller}.{mapper}.truvari_{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.base.vcf.gz'),
        bvcf_index=temp('{root}/temp/{caller}.{mapper}.truvari_{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.base.vcf.gz.tbi'),
        cvcf_temp=temp('{root}/temp/{caller}.{mapper}.temp.truvari_{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.call.vcf.gz'),
        cvcf_temp_index=temp('{root}/temp/{caller}.{mapper}.temp.truvari_{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.call.vcf.gz.tbi'),
        cvcf=temp('{root}/temp/{caller}.{mapper}.truvari_{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.call.nomulti.vcf.gz'),
        cvcf_index=temp('{root}/temp/{caller}.{mapper}.truvari_{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.call.nomulti.vcf.gz.tbi')
    resources:
        mem_mb=12000,
        runtime=360,
        slurm_partition=choose_partition(360)
    container: 'docker://quay.io/jmonlong/truvari:v4.3.1'
    threads: 4
    shell:
        """
        # TODO: Use LARGE_TEMP_DIR?
        export TMPDIR={wildcards.root}/temp

        bcftools norm -m -both {input.truth_vcf} -O z -o {output.bvcf}
        tabix -p vcf {output.bvcf}
        zcat {input.sample_vcf} | bgzip > {output.cvcf_temp}
        tabix -p vcf {output.cvcf_temp}
        bcftools norm -m -both {output.cvcf_temp} -O z -o {output.cvcf}
        tabix -p vcf {output.cvcf}
        
        truvari bench -b {output.bvcf} -c {output.cvcf} -o {output.odir} --pick ac --passonly -r 2000 -C 5000 --includebed {input.confreg}

        truvari refine --recount --use-region-coords --use-original-vcfs --align mafft --reference {input.ref} --regions {output.odir}/candidate.refine.bed {output.odir}

        cp {output.odir}/summary.json {output.summary}
        cp {output.odir}/refine.variant_summary.json {output.refine_var}
        cp {output.odir}/refine.region_summary.json {output.refine_reg}
        """

#Print summary statistics from sv calling with the format:
#This outputs a tsv of: tp, fn, fp, recall, precision, f1
#The input json has "TP-base", "TP-comp", "FP", "FN", "precision", "recall", "f1", "base cnt", "comp cnt"
# This uses TP-base, the number of matching variants as counted by the truth set
rule sv_summary_by_condition:
    input:
        json="{root}/svcall/{caller}/{mapper}/eval/{truthset}_{truthref}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.truvari.refine.variant_summary.json"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/stats/{caller}/{mapper}/{truthset}_{truthref}/{realness}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.sv_summary.tsv"
    wildcard_constraints:
        refgraph="[^/_]+"
    threads: 1
    resources:
        mem_mb=200,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "echo \"{params.condition_name}\t$(jq -r '[.[\"TP-base\"],.FN,.FP,.recall,.precision,.f1] | @tsv' {input.json})\" >{output.tsv}"
rule sv_summary_table:
    input:
        tsv=lambda w: all_experiment(w, "{root}/experiments/{expname}/stats/{caller}/{mapper}/{truthset}_{truthref}/{realness}/{sample}{trimmedness}.{subset}.{tech}.{reference}.{refgraph}.{truthset}_{truthref}.sv_summary.tsv")
    output:
        tsv="{root}/experiments/{expname}/results/sv_summary.tsv"
    threads: 1
    resources:
        mem_mb=200,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        """
        printf "condition\tTP\tFN\tFP\trecall\tprecision\tF1\n" >> {output.tsv}
        cat {input} >>{output.tsv}
        """
ruleorder: sv_summary_table > experiment_stat_table





#Make the figures we want for the paper
#Experiments are specified in the config file
#TODO: idk what we want for variant calling
rule all_paper_figures:
    input:
        mapping_stats_real=expand(ALL_OUT_DIR + "/experiments/{expname}/results/mapping_stats_real.latex.tsv", expname=config["real_exps"]),
        softclippeds=expand(ALL_OUT_DIR + "/experiments/{expname}/results/softclipped_or_unmapped_percent.tsv", expname=config["headline_real_exps"]),
        runtime=expand(ALL_OUT_DIR + "/experiments/{expname}/results/run_and_sampling_time_from_benchmark.tsv", expname=config["headline_real_exps"]),
        index_time=expand(ALL_OUT_DIR + "/experiments/{expname}/results/index_load_time.tsv", expname=config["headline_real_exps"]),
        memory=expand(ALL_OUT_DIR + "/experiments/{expname}/results/memory_from_benchmark.tsv", expname=config["headline_real_exps"]),
        mapping_stats_sim=expand(ALL_OUT_DIR + "/experiments/{expname}/results/mapping_stats_sim.tsv", expname=config["headline_sim_exps"]),
        compared_sim=expand(ALL_OUT_DIR + "/experiments/{expname}/results/compared.tsv", expname=config["headline_sim_exps"]),
        dv_indel=expand(ALL_OUT_DIR + "/experiments/{expname}/plots/dv_indel_summary.pdf", expname=config["dv_exps"]),
        dv_snp=expand(ALL_OUT_DIR + "/experiments/{expname}/plots/dv_snp_summary.pdf", expname=config["dv_exps"]),
        svs=expand(ALL_OUT_DIR + "/experiments/{expname}/plots/sv_summary.pdf", expname=config["sv_exps"])

