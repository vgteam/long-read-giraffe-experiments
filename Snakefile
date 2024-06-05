###########################################################
# Experimental procedure for evaluating Long Read Giraffe #
###########################################################

import parameter_search

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
GRAPHS_DIR = config.get("graphs_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs"

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
# first 100 and last 10 bases trimmed off. Also, there must be
# "{basename}-{subset}.fq" files for each subset size in reads ("1k", "1m:,
# etc.) that you want to work with. "_" and "." are also accepted for setting
# off the subset, and that ".fastq" is alos accepted as the extension.
#
# For simulated reads, each sample directory must have files
# "{sample}-sim-{tech}-{subset}.gam" for each subset size as a number (100, 1000,
# 1000000, etc.) that you want to work with. If the directory is not writable,
# it must already have abbreviated versions ("1k" or "1m" instead of the full
# number) of the GAM files, and the corresponding extracted ".fq" files.
#
# Simulated reads should be made with the "make_pbsim_reads.sh" script in this
# repository.
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
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.10k.fastq
#│           ├── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.1k.fastq
#│           └── HG002_1_R1041_UL_Guppy_6.3.7_5mc_cg_sup_prom_pass.trimmed.1m.fastq
#└── sim
#    ├── hifi
#    │   └── HG002
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
# chm13-pansn.fa
#
# Index files for Minimap2 for each preset (here "hifi", can also be "ont" or "sr", and can be generated from the FASTA):
# chm13-pansn.hifi.mmi
# 
# A Winnowmap repetitive kmers file:
# chm13-pansn.repetitive_k15.txt
#
# TODO: Right now these indexes must be manually generated.
#
REFS_DIR = config.get("refs_dir", None) or "/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references"

# What stages does the Giraffe mapper report times for?
STAGES = ["minimizer", "seed", "tree", "fragment", "chain", "align", "winner"]

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

# How many threads do we want mapping to use?
MAPPER_THREADS=64

PARAM_SEARCH = parameter_search.ParameterSearch()

#Different phoenix nodes seem to run at different speeds, so we can specify which node to run
#This gets added as a slurm_extra for all the real read runs
REAL_SLURM_EXTRA = config.get("real_slurm_extra", None) or ""

# If set to True, jobs where we care about speed will demand entire nodes.
# If False, they will just use one thread per core.
EXCLUSIVE_TIMING = config.get("exclusive_timing", True)

wildcard_constraints:
    trimmedness="\\.trimmed|",
    sample=".+(?<!\\.trimmed)",
    basename=".+(?<!\\.trimmed)",
    subset="[0-9]+[km]?",
    tech="[a-zA-Z0-9]+",
    statname="[a-zA-Z0-9_]+(?<!compared)(.mean|.total)?",
    statnamex="[a-zA-Z0-9_]+(?<!compared)(.mean|.total)?",
    statnamey="[a-zA-Z0-9_]+(?<!compared)(.mean|.total)?",
    realness="(real|sim)",
    realnessx="(real|sim)",
    realnessy="(real|sim)",

def auto_mapping_threads(wildcards):
    """
    Choose the number of threads to use map reads, from subset.
    """
    number = subset_to_number(wildcards["subset"])
    if number >= 100000:
        return MAPPER_THREADS
    elif number >= 10000:
        return 16
    else:
        return 8

def auto_mapping_slurm_extra(wildcards):
    """
    Determine Slurm extra arguments for a timed, real-read mapping job.
    """
    if EXCLUSIVE_TIMING:
        return "--exclusive " + REAL_SLURM_EXTRA
    else:
        return "--threads-per-core 1 " + REAL_SLURM_EXTRA

def auto_mapping_full_cluster_nodes(wildcards):
    """
    Determine number of full cluster nodes for a timed, real-read mapping job.

    TODO: Is this really used by Slurm?
    """
    if EXCLUSIVE_TIMING:
        return 1
    else:
        return 0

def auto_mapping_memory(wildcards):
    """
    Determine the memory to use for Giraffe mapping, in MB, from subset and tech.
    """
    thread_count = auto_mapping_threads(wildcards)

    base_mb = 25000

    if wildcards["tech"] == "illumina":
        scale_mb = 50000
    elif wildcards["tech"] == "hifi":
        scale_mb = 175000
    else:
        scale_mb = 475000

    # Scale down memory with threads
    return scale_mb / MAPPER_THREADS * thread_count + base_mb



def choose_partition(minutes):
    """
    Get a Slurm partition that can fit a job running for the given number of
    minutes, or raise an error.
    """
    for name, limit in SLURM_PARTITIONS:
        if minutes <= limit:
            return name
    raise ValueError(f"No Slurm partition accepts jobs that run for {minutes} minutes")

def subset_to_number(subset):
    """
    Take a subset like 1m and turn it into a number.
    """

    if subset.endswith("m"):
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
    return os.path.join(REFS_DIR, wildcards["reference"] + "-pansn.repetitive_k15.txt")

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
        "hifi": "map-pb",
        "illumina": "sr" # Only Minimap2 has this one, Winnowmap doesn't.
    }

    return MODE_BY_TECH[wildcards["tech"]]

def minimap2_index(wildcards):
    """
    Find the minimap2 index from reference and tech.
    """
    
    mode_part = minimap_derivative_mode(wildcards)
    return os.path.join(REFS_DIR, wildcards["reference"] + "-pansn." + mode_part + ".mmi")
    
def reference_fasta(wildcards):
    """
    Find the linear reference FASTA from a reference.
    """
    return os.path.join(REFS_DIR, wildcards["reference"] + "-pansn.fa")

def graph_base(wildcards):
    """
    Find the base name for a collection of graph files from reference.
    """
    if wildcards["refgraph"] == "hprc-v1.1-mc":
        return os.path.join(GRAPHS_DIR, "hprc-v1.1-mc-" + wildcards["reference"] + ".d9")
    else:
        return os.path.join(GRAPHS_DIR, wildcards["refgraph"] + "-" + wildcards["reference"])

def gbz(wildcards):
    """
    Find a graph GBZ file from reference.
    """
    return graph_base(wildcards) + ".gbz"

def gfa(wildcards):
    """
    Find a graph GFA file from reference.
    """
    return graph_base(wildcards) + ".gfa"

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
    """
    base = graph_base(wildcards)
    indexes = dist_indexed_graph(wildcards)
    new_indexes = {
        "minfile": base + "." + wildcards["minparams"] + ".withzip.min",
        "zipfile": base + "." + wildcards["minparams"] + ".zipcodes"
    }
    new_indexes.update(indexes)
    return new_indexes

def base_fastq(wildcards):
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
    if len(results) == 0:
        # Can't find it
        if wildcards["trimmedness"] == ".trimmed":
            # Look for an untrimmed version
            untrimmed_pattern = os.path.join(READS_DIR, "{realness}/{tech}/{sample}/*{sample}*.f*q.gz".format(**wildcards))
            results = glob.glob(untrimmed_pattern)
            if len(results) == 1:
                # We can use this untrimmed one to make a trimmed one
                without_gz = os.path.splitext(results[0])[0]
                without_fq, fq_ext = os.path.splitext(without_gz)
                trimmed_base = without_fq + ".trimmed" + fq_ext + ".gz"
                return trimmed_base
        raise FileNotFoundError(f"No files found matching {full_gz_pattern}")
    elif len(results) > 1:
        raise RuntimeError("Multiple files matched " + full_gz_pattern)
    return results[0]

def fastq(wildcards):
    """
    Find a FASTQ from realness, tech, sample, trimmedness, and subset.

    Works even if there is extra stuff in the name besides sample. Accounts for
    being able to make a FASTQ from a GAM.
    """
    import glob
    fastq_by_sample_pattern = os.path.join(READS_DIR, "{realness}/{tech}/{sample}/*{sample}*{trimmedness}[._-]{subset}.f*q".format(**wildcards))
    results = glob.glob(fastq_by_sample_pattern)
    if wildcards["trimmedness"] != ".trimmed":
        # Don't match trimmed files when not trimmed.
        results = [r for r in results if ".trimmed" not in r]
    if len(results) == 0:
        if wildcards["realness"] == "real":
            # Make sure there's a full .fq.gz to extract from (i.e. this doesn't raise)
            full_file = base_fastq(wildcards)
            # And compute the subset name
            without_gz = os.path.splitext(full_file)[0]
            without_fq = os.path.splitext(without_gz)[0]
            return without_fq + ".{subset}.fq".format(**wildcards)
        elif wildcards["realness"] == "sim":
            # Assume we can get this FASTQ.
            # For simulated reads we assume the right subset GAM is there. We
            # don't want to deal with the 1k/1000 difference here.
            return os.path.join(READS_DIR, "{realness}/{tech}/{sample}/{sample}-{realness}-{tech}{trimmedness}-{subset}.fq".format(**wildcards))
        else:
            raise FileNotFoundError(f"No files found matching {fastq_by_sample_pattern}")
    elif len(results) > 1:
        raise AmbiguousRuleException("Multiple files matched " + fastq_by_sample_pattern)
    return results[0]

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

    for value in possible_values:
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
        for with_rest in augmented_with_all(base_dict, rest):
            # Augment with the rest
            for with_first in augmented_with_each(with_rest, first_key, first_values):
                # And augment with this key
                yield with_first


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
    Filter dowen wildcards to just the condition parameters for the experiment in expname.
    
    Raises an error if any variable in the experiment cannot be determined.
    """

    exp_dict = config.get("experiments", {}).get(all_wildcards["expname"], {})
    base_condition = exp_dict.get("control", {})
    to_vary = exp_dict.get("vary", {})
    all_vars = list(base_condition.keys()) + list(to_vary.keys())

    condition = {}

    for var in all_vars:
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
    
    # Paste together all the varied variable values from the condition.
    varied = list(to_vary.keys())
    varied_values = [condition[v] for v in varied]
    return ",".join(varied_values)

def all_experiment(wildcard_values, pattern, filter_function=None, debug=False):
    """
    Produce all values of pattern substituted with the wildcards and the experiment conditions' values, from expname.

    If provided, restricts to conditions passing the filter function.
    
    Needs to be used like:
        lambda w: all_experiment(w, "your pattern")
    """

    for condition in all_experiment_conditions(wildcard_values["expname"], filter_function=filter_function):
        merged = dict(wildcard_values)
        merged.update(condition)
        if debug:
            print(f"Evaluate {pattern} in {merged} from {wildcard_values} and {condition}")
        filename = pattern.format(**merged)
        yield filename

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
        case "fragonly":
            return "--fragment-max-lookback-bases 3000 --max-lookback-bases 0"
        case minfrag_number if minfrag_number[0:7] == "minfrag":
            return "--fragment-score-fraction 0  --fragment-min-score " + minfrag_number[7:]
        case maxminfrag_number if maxminfrag_number[0:10] == "maxminfrag":
            return "--fragment-max-min-score " + maxminfrag_number[10:]
        case maxgap_number if minfrag_number[0:6] == "maxgap":
            return "--max-tail-gap " + maxgap_number[6:]
        case gapscale_number if gapscale_number[0:8] == "gapscale":
            return "--gap-scale " + gapscale_number[8:]
        case "gappy":
            return "--gap-scale 0.05 --max-lookback-bases-per-base 0.3"
        case "gappier":
            return "--gap-scale 0.05 --max-lookback-bases-per-base 0.6 --hard-hit-cap 40000 --max-min 200 --num-bp-per-min 250 --downsample-window-count 1000 --downsample-window-length 10"
        case "noflags":
            return ""
        case unknown:
            #otherwise this is a hash and we get the flags from ParameterSearch
            return PARAM_SEARCH.hash_to_parameter_string(wildcard_flag)
        #raise ValueError(f"Unknown flag set: \"{unknown}\"")

def get_vg_version(wildcard_vgversion):
    if wildcard_vgversion == "default":
        return "vg"
    else:
        return "./vg_"+wildcard_vgversion


def param_search_tsvs(wildcards, statname="time_used.mean", realness="real"):
    """
    Get the combined (i.e. mean) TSVs for the conditions in the parameter search.

    TSVs are in the same order as PARAM_SEARCH.get_hashes().

    Needs to be used like:
        lambda w: param_search_tsv(w, "time_used.mean")
    """

    trimmedness = ".trimmed" if wildcards["tech"] in ("r9", "r10", "q27") and realness == "real" else wildcards["trimmedness"]
    values = dict(wildcards)
    values["trimmedness"] = trimmedness
    values["param_hash"] = PARAM_SEARCH.get_hashes()
    values["realness"] = realness
    values["statname"] = statname
    
    return expand("{root}/stats/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{param_hash}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv", **values)


rule minimizer_index_graph:
    input:
        unpack(dist_indexed_graph)
    output:
        minfile="{graphs_dir}/{refgraph}-{reference}.{d9}k{k}.w{w}{weightedness}.withzip.min",
        zipfile="{graphs_dir}/{refgraph}-{reference}.{d9}k{k}.w{w}{weightedness}.zipcodes"
    wildcard_constraints:
        weightedness="\\.W|",
        k="[0-9]+",
        w="[0-9]+",
        reference="chm13|grch38",
        d9="d9\.|"
    params:
        weighting_option=lambda w: "--weighted" if w["weightedness"] == ".W" else ""
    threads: 16
    resources:
        mem_mb=lambda w: 320000 if w["weightedness"] == ".W" else 80000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "vg minimizer --progress -k {wildcards.k} -w {wildcards.w} {params.weighting_option} -t {threads} -p -d {input.dist} -z {output.zipfile} -o {output.minfile} {input.gbz}"

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
    threads: 4
    resources:
        mem_mb=10000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "seqkit subseq -j {threads} -r 100:-10 {input.fq_gz} -o {output.fq_gz}"
    

rule subset_base_fastq_gz:
    input:
        base_fastq=base_fastq
    output:
        fastq="{reads_dir}/{realness}/{tech}/{sample}/{basename}{trimmedness}.{subset}.fq"
    wildcard_constraints:
        realness="real"
    params:
        lines=lambda w: str(subset_to_number(w["subset"]) * 4)
    threads: 8
    resources:
        mem_mb=10000,
        runtime=120,
        slurm_partition=choose_partition(120)
    shell:
        # We need to account for bgzip getting upset that we close the pipe before it is done writing.
        "(bgzip -d <{input.base_fastq} || true) | head -n {params.lines} >{output.fastq}"

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

rule giraffe_real_reads:
    input:
        unpack(indexed_graph),
        fastq=fastq,
    output:
        # Giraffe can dump out pre-annotated reads at annotation range -1.
        gam="{root}/aligned/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    benchmark: "{root}/aligned/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    wildcard_constraints:
        realness="real"
    threads: auto_mapping_threads
    resources:
        mem_mb=auto_mapping_memory,
        runtime=600,
        slurm_partition=choose_partition(600),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    run:
        vg_binary = get_vg_version(wildcards.vgversion)
        flags=get_vg_flags(wildcards.vgflag)

        shell(vg_binary + " giraffe -t{threads} --parameter-preset {wildcards.preset} --progress --track-provenance -Z {input.gbz} -d {input.dist} -m {input.minfile} -z {input.zipfile} -f {input.fastq} " + flags + " >{output.gam}")

rule giraffe_sim_reads:
    input:
        unpack(indexed_graph),
        gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/annotated-1/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
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

        shell(vg_binary + " giraffe -t{threads} --parameter-preset {wildcards.preset} --progress --track-provenance --set-refpos -Z {input.gbz} -d {input.dist} -m {input.minfile} -z {input.zipfile} -G {input.gam} " + flags + " >{output.gam}")

rule giraffe_sim_reads_with_correctness:
    input:
        unpack(indexed_graph),
        gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/correctness/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
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

        shell(vg_binary + " giraffe -t{threads} --parameter-preset {wildcards.preset} --progress --track-provenance --track-correctness --set-refpos -Z {input.gbz} -d {input.dist} -m {input.minfile} -z {input.zipfile} -G {input.gam} " + flags + " >{output.gam}")

rule winnowmap_sim_reads:
    input:
        reference_fasta=reference_fasta,
        repetitive_kmers=repetitive_kmers,
        fastq=fastq
    params:
        mode=minimap_derivative_mode,
    output:
        sam="{root}/aligned-secsup/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam"
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
        "winnowmap -t {threads} -W {input.repetitive_kmers} -ax {params.mode} {input.reference_fasta} {input.fastq} >{output.sam}"

rule winnowmap_real_reads:
    input:
        reference_fasta=reference_fasta,
        repetitive_kmers=repetitive_kmers,
        fastq=fastq
    params:
        mode=minimap_derivative_mode,
    output:
        sam="{root}/aligned-secsup/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam"
    benchmark: "{root}/aligned-secsup/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    wildcard_constraints:
        realness="real"
    wildcard_constraints:
        # Winnowmap doesn't have a short read preset, so we can't do Illumina reads.
        # So match any string but that. See https://stackoverflow.com/a/14683066
        tech="(?!illumina).+"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "winnowmap -t {threads} -W {input.repetitive_kmers} -ax {params.mode} {input.reference_fasta} {input.fastq} >{output.sam}"

rule minimap2_index_reference:
    input:
        reference_fasta=reference_fasta
    output:
        index=REFS_DIR + "/{reference}-pansn.{preset}.mmi"
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
        sam="{root}/aligned-secsup/{reference}/minimap2-{minimapmode}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam"
    wildcard_constraints:
        realness="sim"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "minimap2 -t {threads} -ax {wildcards.minimapmode} -N 0 {input.minimap2_index} {input.fastq} >{output.sam}"

rule minimap2_real_reads:
    input:
        minimap2_index=minimap2_index,
        fastq=fastq
    output:
        sam="{root}/aligned-secsup/{reference}/minimap2-{minimapmode}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam"
    benchmark: "{root}/aligned-secsup/{reference}/minimap2-{minimapmode}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    wildcard_constraints:
        realness="real"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "minimap2 -t {threads} -ax {wildcards.minimapmode} -N 0 {input.minimap2_index} {input.fastq} >{output.sam}"

# Minimap2 and Winnowmap include secondary alignments in the output by default, and Winnowmap doesn't quite have a way to limit them (minimap2 has -N)
# Also they only speak SAM and we don't want to benchmark the BAM-ification time.
# So drop secondary and supplementary alignments and conver tto BAM.
# TODO: Get the downstream stats tools to be able to do things like measure average softclips when there are supplementary alignments.
rule drop_secondary_and_supplementary:
    input:
        sam="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.sam"
    output:
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    wildcard_constraints:
        mapper="(minimap2-.*|winnowmap)"
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
        gam="{root}/aligned/{reference}/{refgraph}/graphaligner/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        realness="sim"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "GraphAligner -t {threads} -g {input.gfa} -f {input.fastq} -x vg -a {output.gam}"

rule graphaligner_real_reads:
    input:
        gfa=gfa,
        fastq=fastq
    output:
        gam="{root}/aligned/{reference}/{refgraph}/graphaligner/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    benchmark: "{root}/aligned/{reference}/{refgraph}/graphaligner/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    wildcard_constraints:
        realness="real"
    threads: auto_mapping_threads
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600),
        slurm_extra=auto_mapping_slurm_extra,
        full_cluster_nodes=auto_mapping_full_cluster_nodes
    shell:
        "GraphAligner -t {threads} -g {input.gfa} -f {input.fastq} -x vg -a {output.gam}"

rule inject_bam:
    input:
        gbz=gbz,
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        mapper="(minimap2.+|winnowmap)"
    threads: 64
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg inject --threads {threads} -x {input.gbz} {input.bam} >{output.gam}"

rule compare_alignments:
    input:
        gam="{root}/annotated-1/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        truth_gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        tsv="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv",
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    threads: 16
    resources:
        mem_mb=200000,
        runtime=800,
        slurm_partition=choose_partition(800)
    shell:
        "vg gamcompare --threads 16 --range 200 {input.gam} {input.truth_gam} --output-gam {output.gam} -T -a {wildcards.mapper} > {output.tsv} 2>{output.compare}"

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
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.correct.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.compare} | grep -o '[0-9]* reads correct' | cut -f1 -d' ' >{output.tsv}"

rule accuracy_from_comparison:
    input:
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.accuracy.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "cat {input.compare} | grep -o '[0-9%.]* accuracy' | cut -f1 -d' ' >{output.tsv}"

rule wrong_from_comparison:
    input:
        compare="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.wrong.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "echo \"$(cat {input.compare} | grep -o '[0-9]* reads eligible' | cut -f1 -d' ') - $(cat {input.compare} | grep -o '[0-9]* reads correct' | cut -f1 -d' ')\" | bc -l >{output.tsv}"

rule comparison_experiment_stat:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{comparisonstat}.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{comparisonstat}.tsv"
    wildcard_constraints:
        comparisonstat="(wrong|accuracy|correct)"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && cat {input.tsv} >>{output.tsv}"


rule experiment_stat_table:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv", filter_function=has_stat_filter(w["statname"]))
    output:
        table="{root}/experiments/{expname}/results/{statname}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "cat {input} >{output.table}"

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
        "python3 barchart.py {input.tsv} --title '{wildcards.expname} Correctness' --y_label 'Correct Reads' --x_label 'Condition' --x_sideways --no_n --save {output}"

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
        "python3 barchart.py {input.tsv} --title '{wildcards.expname} Wrongness' --y_label 'Wrong Reads' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule compared_named_from_compared:
    input:
        tsv="{root}/compared/{reference}/{refgraph}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv",
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compared.tsv"
    threads: 3
    resources:
        mem_mb=1000,
        runtime=120,
        slurm_partition=choose_partition(120)
    shell:
        "printf 'correct\\tmq\\taligner\\tread\\teligible\\n' >{output.tsv} && cat {input.tsv} | grep -v '^correct' | awk -F '\\t' -v OFS='\\t' '{{ $3 = \"{params.condition_name}\"; print }}' >>{output.tsv}"


rule experiment_compared_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compared.tsv", lambda condition: condition["realness"] == "sim")
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
        "{root}/experiments/{expname}/plots/pr.{ext}"
    threads: 1
    resources:
        mem_mb=10000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "Rscript plot-pr.R {input.tsv} {output}"

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
        "python3 barchart.py {input.tsv} --title '{wildcards.expname} Mapping Rate' --y_label 'Mapped Reads' --x_label 'Condition' --x_sideways --no_n --save {output}"

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
        "python3 barchart.py {input.tsv} --title '{wildcards.expname} Speed' --y_label 'Reads per Second' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule softclips_from_mean_softclips:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips.mean.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips.tsv"
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
        "python3 barchart.py {input.tsv} --title '{wildcards.expname} Softclips' --y_label 'Mean Softclip (bp)' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule softclipped_or_unmapped_from_softclipped_or_unmapped:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclipped_or_unmapped.tsv"
    params:
        condition_name=condition_name
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
        "python3 barchart.py {input.tsv} --title '{wildcards.expname} Softclipped or Unmapped Bases' --y_label 'Total (bp)' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule chain_coverage_from_mean_best_chain_coverage:
    input:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.mean.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.chain_coverage.tsv"
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
        "python3 barchart.py {input.tsv} --title '{wildcards.expname} Chain Coverage' --y_label 'Best Chain Coverage (fraction)' --x_label 'Condition' --x_sideways --no_n --save {output}"


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

rule time_used:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.time_used.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"time_used\" {input.gam} | grep -v \"#\" >{output}"

rule stage_time:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.stage_{stage}_time.tsv"
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

rule length_by_mapping:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "(vg filter --only-mapped {input.gam} -T length | grep -v '^#' | sed 's/^/mapped\t/'; vg filter --only-mapped --complement {input.gam} -T length | grep -v '^#' | sed 's/^/unmapped\t/') >{output}"

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
        "cat {input.tsv} | grep '^unmapped' | cut -f2 >{output}"

rule length_by_correctness:
    input:
        gam="{root}/compared/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_correctness.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"correctness;sequence\" {input.gam} | grep -v \"#\" | awk -v OFS='\t' '{{print $1, length($2)}}' > {output}"

rule softclips_by_name_gam:
    input:
        gam="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips_by_name.tsv"
    wildcard_constraints:
        mapper="(giraffe.*|graphaligner)"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"name;softclip_start;softclip_end\" {input.gam} | grep -v \"#\" > {output}"

rule softclips_by_name_other:
    input:
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips_by_name.tsv"
    wildcard_constraints:
        mapper="(?!giraffe).+"
    threads: 7
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        r"samtools view {input.bam} | cut -f1,2,6 | sed 's/\t\(\([0-9]*\)S\)\?\([0-9]*[IDMH]\|\*\)*\(\([0-9]*\)S\)\?$/\t\2\t\5/g' | sed 's/\t\t/\t0\t/g' | sed 's/\t$/\t0/g' | sed 's/16\t\([0-9]*\)\t\([0-9]*\)/\2\t\1/g' | sed 's/\t[0-9]\+\t\([0-9]*\t[0-9]*\)$/\t\1/g' > {output}"

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

rule softclipped_or_unmapped:
    input:
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.softclips.total.tsv",
        "{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.unmapped_length.total.tsv"
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
        

rule memory_usage_gam:
    input:
        benchmark="{root}/aligned/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_usage.tsv"
    wildcard_constraints:
        realness="real",
        mapper="(giraffe.+|graphaligner)"
    threads: MAPPER_THREADS
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        # max_rss happens to be column 3, but try to check
        "cat {input.benchmark} | cut -f3 | grep -A1 max_rss | tail -n1 >{output.tsv}"

rule memory_usage_sam:
    input:
        benchmark="{root}/aligned-secsup/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.benchmark"
    output:
        tsv="{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.memory_usage.tsv"
    wildcard_constraints:
        realness="real",
        mapper="(minimap2.+|winnowmap)"
    threads: MAPPER_THREADS
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        # max_rss happens to be column 3, but try to check
        "cat {input.benchmark} | cut -f3 | grep -A1 max_rss | tail -n1 >{output.tsv}"

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
        # Input files must be in the same order as STAGES
        expand("{{root}}/stats/{{reference}}/{{refgraph}}/{{mapper}}/{{realness}}/{{tech}}/{{sample}}{{trimmedness}}.{{subset}}.stage_{stage}_time.mean.tsv", stage=STAGES)
    output:
        "{root}/tables/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_stage_time.tsv"
    threads: 1
    resources:
        mem_mb=512,
        runtime=20,
        slurm_partition=choose_partition(20)
    run:
        # Make a TSV of stage name and its average value
        with open(output[0], "w") as out_stream:
            for (stage, filename) in zip(STAGES, input):
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
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {STAGES} --title '{wildcards.tech} {wildcards.realness} Mean Stage Times' --y_label 'Time (s)' --x_label 'Stage' --no_n --save {output}"

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

rule mapping_stats:
    input:
        compared_tsv="{root}/compared/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{param_hash}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv"
    output:
        "{root}/stats/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{param_hash}/sim/{tech}/{sample}{trimmedness}.{subset}.mapping_stats.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=100,
        slurm_partition=choose_partition(100)
    log: "{root}/stats/{reference}/{refgraph}/giraffe-{minparams}-{preset}-{vgversion}-{param_hash}/sim/{tech}/{sample}{trimmedness}.{subset}.mapping_stats.log"
    run:

        correct_count = 0
        mapq60_count = 0
        wrong_mapq60_count = 0

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
        f.close()
        shell("printf \"" + str(correct_count) + "\t" + str(mapq60_count) + "\t" + str(wrong_mapq60_count) + "\" > {output}")
        

rule parameter_search_mapping_stats:
    input:
        times = lambda w: param_search_tsvs(w, "time_used.mean"),
        mapping_stats = expand("{{root}}/stats/{{reference}}/{{refgraph}}/giraffe-{{minparams}}-{{preset}}-{{vgversion}}-{param_hash}/sim/{{tech}}/{{sample}}{{trimmedness}}.{{subset}}.mapping_stats.tsv",param_hash=PARAM_SEARCH.get_hashes())
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
        f.write("#correct\tmapq60\twrong_mapq60\tspeed(r/s/t)\t" + '\t'.join([param.name for param in PARAM_SEARCH.parameters]))
        for param_hash, stats_file, times_file in zip(PARAM_SEARCH.get_hashes(), input.mapping_stats, input.times):

            param_f = open(stats_file)
            l = param_f.readline().split()
            correct_count = l[0]
            mapq60_count = l[1]
            wrong_mapq60_count = l[2]
            param_f.close()

            time_f = open(times_file)
            l = time_f.readline().split()
            speed = str(1/float(l[0]))
            time_f.close()

            parameters = PARAM_SEARCH.hash_to_parameters[param_hash]
            f.write("\n" + correct_count + "\t" + mapq60_count + "\t" + wrong_mapq60_count + "\t" + speed + "\t" + '\t'.join([str(x) for x in parameters])) 
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
        shell("cat <(cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $4}}' | sed 's/^/RPS /g') <(cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $1}}' | sed 's/^/Correct /g') | ./scatter.py --title 'Speed vs Correctness vs " + wildcards.parameter + "' --x_label " + wildcards.parameter + "  --y_per_category --categories 'RPS' 'Correct' --y_label 'Reads per Second' 'Reads Correct' --legend_overlay 'best' --save {output.plot} /dev/stdin")

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

