###########################################################
# Experimental procedure for evaluating Long Read Giraffe #
###########################################################

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
#│   │   └── HiFi
#│   │       ├── HiFi_reads_10k.fq
#│   │       ├── HiFi_reads_1k.fq
#│   │       ├── HiFi_reads_1m.fq
#│   │       └── HiFi_reads.fq.gz
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
# A few more threads will be used for filters
MAPPER_THREADS=32



wildcard_constraints:
    trimmedness="\\.trimmed|",
    sample=".+(?<!\\.trimmed)",
    basename=".+(?<!\\.trimmed)",
    subset="[0-9]+[km]?",
    statname="[a-zA-Z0-9_]+(?<!compared)"

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

def chunk_count(items, chunk_size):
    """
    Return the number of chunks of the given size needed to fit the given number of items.
    """

    # Since integer division rounds negatively, we can work in negative numbers
    # to make it round away from 0. See <https://stackoverflow.com/a/33300093>
    return -(-items // chunk_size)

def each_chunk_of(subset):
    """
    Given a subset string like "10k", produce a collection of all the p[added chunk number strings.
    """
    return [f"{i:06}" for i in range(1, chunk_count(subset_to_number(subset), CHUNK_SIZE) + 1)]

def all_chunk(wildcard_values, pattern, debug=False):
    """
    Produce all values of pattern substituted with the wildcards and the
    0-padded GAM chunk numbers as {chunk}, from subset.
    
    Needs to be used like:
        lambda w: all_chunk(w, "your pattern")
    """

    for chunk in each_chunk_of(wildcard_values["subset"]):
        merged = dict(wildcard_values)
        merged.update(chunk=chunk)
        if debug:
            print(f"Evaluate {pattern} in {merged}")
        filename = pattern.format(**merged)
        yield filename

def repetitive_kmers(wildcards):
    """
    Find the Winnowmap repetitive kmers file from a reference.
    """
    return os.path.join(REFS_DIR, wildcards["reference"] + "-pansn.repetitive_k15.txt")

def minimap_derivative_mode(wildcards):
    """
    Determine the right Minimap2/Winnowmap preset (map-pb, etc.) from tech.
    """

    return {
        "r9": "map-ont",
        "r10": "map-ont",
        "hifi": "map-pb",
        "illumina": "sr" # Only Minimap2 has this one, Winnowmap doesn't.
    }[wildcards["tech"]]

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
    return os.path.join(GRAPHS_DIR, "hprc-v1.1-mc-" + wildcards["reference"] + ".d9")

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
        raise AmbiguousRuleException("Multiple files matched " + full_gz_pattern)
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
    {root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv

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

        if stat_name.startswith("time_used") or stat_name == "mapping_speed":
            # This is a Giraffe time used stat or mean thereof. We need to be a
            # Giraffe condition.
            if not condition["mapper"].startswith("giraffe"):
                return False

        return True

    return filter_function

def get_vg_flags(wildcard_flag):
    if wildcard_flag == "gapExt":
        return "--do-gapless-extension"
    elif wildcard_flag == "mqCap":
        return "--explored-cap"
    elif wildcard_flag[0:10] == "downsample":
        return "--downsample-min " + wildcard_flag[10:]
    elif wildcard_flag == "fragonly":
        return "--fragment-max-lookback-bases 3000 --max-lookback-bases 0"
    else:
        assert(wildcard_flag == "noflags")
        return ""

def get_vg_version(wildcard_vgversion):
    if wildcard_vgversion == "default":
        return "vg"
    else:
        return "./vg_"+wildcard_vgversion

rule minimizer_index_graph:
    input:
        unpack(dist_indexed_graph)
    output:
        minfile="{graphs_dir}/hprc-v1.1-mc-{reference}.d9.k{k}.w{w}{weightedness}.withzip.min",
        zipfile="{graphs_dir}/hprc-v1.1-mc-{reference}.d9.k{k}.w{w}{weightedness}.zipcodes"
    wildcard_constraints:
        weightedness="\\.W|",
        k="[0-9]+",
        w="[0-9]+"
    threads: 16
    resources:
        mem_mb=80000,
        runtime=240,
        slurm_partition=choose_partition(240)
    shell:
        "vg minimizer --progress -k {wildcards.k} -w {wildcards.w} -t {threads} -p -d {input.dist} -z {output.zipfile} -o {output.minfile} {input.gbz}"

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
        runtime=60,
        slurm_partition=choose_partition(60)
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
        gam="{root}/annotated-1/{reference}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        realness="real"
    threads: MAPPER_THREADS
    resources:
        mem_mb=500000,
        runtime=600,
        slurm_partition=choose_partition(600)
    run:
        vg_binary = get_vg_version(wildcards.vgversion)
        flags=get_vg_flags(wildcards.vgflag)

        shell(vg_binary + " giraffe -t{threads} --parameter-preset {wildcards.preset} --progress --track-provenance --set-refpos -Z {input.gbz} -d {input.dist} -m {input.minfile} -z {input.zipfile} -f {input.fastq} " + flags + " >{output.gam}")

rule giraffe_sim_reads:
    input:
        unpack(indexed_graph),
        gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/annotated-1/{reference}/giraffe-{minparams}-{preset}-{vgversion}-{vgflag}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        realness="sim"
    threads: MAPPER_THREADS
    resources:
        mem_mb=500000,
        runtime=600,
        slurm_partition=choose_partition(600)
    run:
        vg_binary = get_vg_version(wildcards.vgversion)
        flags=get_vg_flags(wildcards.vgflag)

        shell(vg_binary + " giraffe -t{threads} --parameter-preset {wildcards.preset} --progress --track-provenance --set-refpos -Z {input.gbz} -d {input.dist} -m {input.minfile} -z {input.zipfile} -G {input.gam} " + flags + " >{output.gam}")

rule winnowmap_reads:
    input:
        reference_fasta=reference_fasta,
        repetitive_kmers=repetitive_kmers,
        fastq=fastq
    params:
        mode=minimap_derivative_mode
    output:
        bam="{root}/aligned/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    wildcard_constraints:
        # Winnowmap doesn't have a short read preset, so we can't do Illumina reads.
        # So match any string but that. See https://stackoverflow.com/a/14683066
        tech="(?!illumina).+"
    threads: MAPPER_THREADS + 4
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "winnowmap -t {MAPPER_THREADS} -W {input.repetitive_kmers} -ax {params.mode} {input.reference_fasta} {input.fastq} | samtools view --threads 4 -h -F 2048 -F 256 --bam - >{output.bam}"

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

rule minimap2_reads:
    input:
        minimap2_index=minimap2_index,
        fastq=fastq
    params:
        mode=minimap_derivative_mode
    output:
        bam="{root}/aligned/{reference}/minimap2/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    threads: MAPPER_THREADS + 4
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "minimap2 -t {MAPPER_THREADS} -ax {params.mode} {input.minimap2_index} {input.fastq} | samtools view --threads 4 -h -F 2048 -F 256 --bam - >{output.bam}"

rule graphaligner_reads:
    input:
        gfa=gfa,
        fastq=fastq
    output:
        gam="{root}/aligned/{reference}/graphaligner/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    threads: MAPPER_THREADS
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "GraphAligner -t {threads} -g {input.gfa} -f {input.fastq} -x vg -a {output.gam}"

rule inject_bam:
    input:
        gbz=gbz,
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        gam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        mapper="(minimap2|winnowmap)"
    threads: 64
    resources:
        mem_mb=300000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg inject --threads {threads} -x {input.gbz} {input.bam} >{output.gam}"

rule compare_alignments:
    input:
        gbz=gbz,
        gam="{root}/annotated-1/{reference}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        truth_gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/compared/{reference}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        tsv="{root}/compared/{reference}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv",
        report="{root}/compared/{reference}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    threads: 16
    resources:
        mem_mb=200000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg gamcompare --threads 16 --range 200 {input.gam} {input.truth_gam} --output-gam {output.gam} -T -a {wildcards.mapper} > {output.tsv} 2>{output.report}"

rule annotate_alignments:
    input:
        gbz=gbz,
        gam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    output:
        gam="{root}/annotated-1/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        # Giraffe pre-annotates output reads
        mapper="(?!giraffe).+"
    threads: 16
    resources:
        mem_mb=100000,
        runtime=600,
        slurm_partition=choose_partition(600)
    shell:
        "vg annotate -t16 -a {input.gam} -x {input.gbz} -m --search-limit=-1 >{output.gam}"

rule correct_from_comparison:
    input:
        report="{root}/compared/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.correct.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && cat {input.report} | grep -o '[0-9]* reads correct' | cut -f1 -d' ' >>{output.tsv}"

rule accuracy_from_comparison:
    input:
        report="{root}/compared/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.accuracy.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && cat {input.report} | grep -o '[0-9%.]* accuracy' | cut -f1 -d' ' >>{output.tsv}"

rule wrong_from_comparison:
    input:
        report="{root}/compared/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.wrong.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && echo \"$(cat {input.report} | grep -o '[0-9]* reads eligible' | cut -f1 -d' ') - $(cat {input.report} | grep -o '[0-9]* reads correct' | cut -f1 -d' ')\" | bc -l >>{output.tsv}"

rule experiment_stat_table:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv", filter_function=has_stat_filter(w["statname"]))
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
        tsv="{root}/compared/{reference}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv",
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compared.tsv"
    threads: 3
    resources:
        mem_mb=1000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "printf 'correct\\tmq\\taligner\\tread\\teligible\\n' >{output.tsv} && cat {input.tsv} | grep -v '^correct' | awk -F '\\t' -v OFS='\\t' '{{ $3 = \"{params.condition_name}\"; print }}' >>{output.tsv}"


rule experiment_compared_tsv:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compared.tsv", lambda condition: condition["realness"] == "sim")
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
        gam="{root}/annotated-1/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        stats="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gamstats.txt"
    threads: 16
    resources:
        mem_mb=10000,
        runtime=90,
        slurm_partition=choose_partition(90)
    shell:
        "vg stats -p {threads} -a {input.gam} >{output.stats}"

rule facts_from_alignments:
    input:
        gam="{root}/annotated-1/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        facts="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.facts.txt"
        facts_dir="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.facts"
    wildcard_constraints:
        mapper="giraffe-.+"
    threads: 2
    resources:
        mem_mb=10000,
        runtime=90,
        slurm_partition=choose_partition(90)
    shell:
        "python3 giraffe-facts.py {input.gam} {output.facts_dir} >{output.facts}"

rule mapping_rate_from_stats:
    input:
        stats="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gamstats.txt"
    params:
        condition_name=condition_name
    output:
        rate="{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_rate.tsv"
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
        tsv="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.time_used.mean.tsv"
    params:
        condition_name=condition_name
    output:
        tsv="{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_speed.tsv"
    wildcard_constraints:
        mapper="giraffe-.+"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
        slurm_partition=choose_partition(5)
    shell:
        "printf '{params.condition_name}\\t' >{output.tsv} && echo \"1 / $(cat {input.tsv})\" | bc -l >>{output.tsv}"

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

for subset in KNOWN_SUBSETS:
    for stage in ["annotated-1", "compared"]:
        # We can chunk reads either before or after comparison.
        # TODO: This is now like 3 copies of the whole GAM.

        # This rule has a variable number of outputs so we need to generate it in a loop.
        rule:
            input:
                gam="{root}/" + stage + "/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}." + str(subset) + ".gam"
            params:
                basename="{root}/" + stage + "/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}." + str(subset) + ".chunk"
            output:
                expand("{{root}}/{stage}/{{reference}}/{{mapper}}/{{realness}}/{{tech}}/{{sample}}{{trimmedness}}.{subset}.chunk{chunk}.gam", stage=stage, subset=subset, chunk=each_chunk_of(subset))
            threads: 1
            resources:
                mem_mb=4000,
                runtime=90,
                slurm_partition=choose_partition(90)
            shell:
                "vg chunk -t {threads} --gam-split-size " + str(CHUNK_SIZE) + " -a {input.gam} -b {params.basename}"

rule chain_coverage:
    input:
        gam="{root}/annotated-1/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"annotation.best_chain_coverage\" {input.gam} | grep -v \"#\" >{output}"

rule time_used:
    input:
        gam="{root}/annotated-1/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.time_used.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"time_used\" {input.gam} | grep -v \"#\" >{output}"

rule stage_time:
    input:
        gam="{root}/annotated-1/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.stage_{stage}_time.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"annotation.stage_{wildcards.stage}_time\" {input.gam} | grep -v \"#\" >{output}"

rule length_by_mapping_chunk:
    input:
        gam="{root}/annotated-1/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.chunk{chunk}.gam",
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.chunk{chunk}.length_by_mapping.tsv"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=30,
        slurm_partition=choose_partition(30)
    shell:
        "vg view -aj {input.gam} | jq -r '[if (.path.mapping // []) == [] then \"unmapped\" else \"mapped\" end, (.sequence | length)] | @tsv' >{output}"

rule length:
    input:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv"
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"sequence\" {input.gam} | grep -v \"#\" | awk '{{print length($1)}}' >{output}"

rule length_by_correctness:
    input:
        gam="{root}/compared/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_correctness.tsv"
    threads: 5
    resources:
        mem_mb=2000,
        runtime=60,
        slurm_partition=choose_partition(60)
    shell:
        "vg filter -t {threads} -T \"correctness;sequence\" {input.gam} | grep -v \"#\" | awk -v OFS='\t' '{{print $1, length($2)}}' > {output}"

rule merge_stat_chunks:
    input:
        lambda w: all_chunk(w, "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.chunk{chunk}.{statname}.tsv")
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=20,
        slurm_partition=choose_partition(20)
    shell:
        "cat {input} >{output}"

rule mean_stat:
    input:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.tsv"
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.{statname}.mean.tsv"
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

rule average_stage_time_table:
    input:
        # Input files must be in the same order as STAGES
        expand("{{root}}/stats/{{reference}}/{{mapper}}/{{realness}}/{{tech}}/{{sample}}{{trimmedness}}.{{subset}}.stage_{stage}_time.mean.tsv", stage=STAGES)
    output:
        "{root}/tables/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_stage_time.tsv"
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

rule chain_coverage_histogram:
    input:
        tsv="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.tsv"
    output:
        "{root}/plots/{reference}/{mapper}/best_chain_coverage-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title '{wildcards.tech} {wildcards.realness} Fraction Covered' --y_label 'Items' --x_label 'Coverage' --no_n --save {output}"

rule time_used_histogram:
    input:
        tsv="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.time_used.tsv",
        mean="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.time_used.mean.tsv"
    output:
        "{root}/plots/{reference}/{mapper}/time_used-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title \"{wildcards.tech} {wildcards.realness} Time Used, Mean=$(cat {input.mean})\" --y_label 'Items' --x_label 'Time (s)' --no_n --save {output}"

rule stage_time_histogram:
    input:
        tsv="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.stage_{stage}_time.tsv",
        mean="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.stage_{stage}_time.mean.tsv"
    output:
        "{root}/plots/{reference}/{mapper}/stage_{stage}_time-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title \"{wildcards.tech} {wildcards.realness} Stage {wildcards.stage} Time, Mean=$(cat {input.mean})\" --y_label 'Items' --x_label 'Time (s)' --no_n --save {output}"

rule average_stage_time_barchart:
    input:
        tsv="{root}/tables/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.average_stage_time.tsv"
    output:
        "{root}/plots/{reference}/{mapper}/average_stage_time-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 barchart.py {input.tsv} --categories {STAGES} --title '{wildcards.tech} {wildcards.realness} Mean Stage Times' --y_label 'Time (s)' --x_label 'Stage' --no_n --save {output}"

rule length_by_mapping_histogram:
    input:
        tsv="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv",
        mean="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length.mean.tsv"
    output:
        "{root}/plots/{reference}/{mapper}/length_by_mapping-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --bins 100 --title \"{wildcards.tech} {wildcards.realness} Read Length, Mean=$(cat {input.mean})\" --y_label 'Items' --x_label 'Length (bp)' --no_n --categories mapped unmapped --category_labels Mapped Unmapped --legend_overlay 'best' --save {output}"


rule length_by_correctness_histogram:
    input:
        tsv="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_correctness.tsv"
    output:
        "{root}/plots/{reference}/{mapper}/length_by_correctness-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
        slurm_partition=choose_partition(10)
    shell:
        "python3 histogram.py {input.tsv} --log_counts --bins 100 --title '{wildcards.tech} {wildcards.realness} Read Length for {wildcards.mapper}' --y_label 'Items' --x_label 'Length (bp)' --no_n --categories correct incorrect off-reference --category_labels Correct Incorrect 'Off Reference' --legend_overlay 'best' --stack --save {output}"


ruleorder: chain_coverage > merge_stat_chunks 
ruleorder: time_used > merge_stat_chunks 
ruleorder: stage_time > merge_stat_chunks
ruleorder: length > merge_stat_chunks 
ruleorder: length_by_correctness > merge_stat_chunks  

