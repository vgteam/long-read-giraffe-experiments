# Long Read Giraffe Experiments

This repository contains scripts for running experiments to evaluate the long-read mode for the `vg giraffe` pangenome sequencing read mapper. It knows how to run Giraffe, Winnowmap, and Minimap2, and how to analyze the results with `vg` and the included scripts.

## Usage

The general flow for using this repository is:

1. Obtain some real long reads for your sequencing technologies.

2. Use `make_pbsim_reads.sh` to simulate reads with plausible base qualities, for a sample in a sample-and-reference graph. Use `simulate_illumina_reads.sh` to make Illumina reads and `simulate_element_reads.sh` to make Element reads. All scripts take their input configuration via environment variables. For example:
```
SAMPLE_FASTQ=/private/groups/patenlab/xhchang/reads/real/element/HG002/HG002.GAT-LI-C044.fq.gz VG="${HOME}/workspace/lr-giraffe/vg_v1.64.1" READ_DIR=/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads GRAPH_DIR=/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs ./simulate_element_reads.sh
```

Note that the Illumina and Element read subsets won't be uniformly sampled from the full set.

3. Lay out your input files in the expected format, including a GBZ graph and distance index, the reads, a linear reference FASTA, and Winnowmap indexes, and any Minimap2 indexes. See `Snakefile` for a description of the required layout for input data. If you need to get a linear reference from a graph, you can run something like:
```
vg paths --extract-fasta --sample CHM13 -x graphs/hprc-v2.0-mc-chm13-eval.gbz > references/chm13-pansn-newY.fa
```

4. Write a YAML configuration file defining one or more experiments in which different variables (sequencing technology, number of reads, read aligner used, etc.) are either held constant or varied. See `lr-config.yaml`, the default config file, for an example and documentation for all the experimental variables and their possible values. You can also use the config file to set `graphs_dir`, `reads_dir`, and `refs_dir` to directories where your graphs, reads, and Winnowmap2/Minimap indexed linear references are stored, if not running with access to the UCSC Genomics Institute `/private/groups` directory hierarchy. 

5. Use Snakemake and the included `Snakefile` to do the work necessary to create the experimental results files you are interested in. For example if you write `lr-config.yaml` to define an experiment named `r10_accuracy_small` to compare Winnowmap and Giraffe on some R10 nanopore reads, you can get a plot of the rate of mapped reads across the two conditions in the experiment with:

    ```
    snakemake --rerun-incomplete --use-singularity --singularity-args "-B /private" output/experiments/r10_accuracy_small/plots/mapping_rate.png
    ```

    Make sure that `output` is a directory with enough free space to store the mapped reads. Also, if any inputs or outputs are not in your home directory, make sure that the Singularity mountpoint (here `/private`) is an absolute path to a parent of all input and output files (including symlink destinations), but is not `/`. (If there is no common parent, you can use multiple `-B` flags in the string.) We use the `--rerun-incomplete` flag to make Snakemake rerun any commands that were interrupted on previous runs.

    You can also run jobs on Slurm, where they will be automatically assigned to the correct partitions for UCSC's Phoenix cluster:

    ```
    snakemake --rerun-incomplete --use-singularity --singularity-args "-B /private" --slurm --latency-wait 120 -j128 output/experiments/r10_accuracy_small/plots/mapping_rate.png
    ```

    This will run on Slurm, waiting up to 120 seconds for files created on the worker nodes to become visible to the head node, and using up to 128 cores total.

    You can also ask for multiple output files in a single run:

    ```
    snakemake --rerun-incomplete --use-singularity --singularity-args "-B /private" --slurm --latency-wait 120 -j128 \
        output/experiments/r10_accuracy_small/plots/mapping_rate.png \
        output/experiments/r10_accuracy_small/plots/correct.png \
        output/plots/chm13/hprc-v1.1-mc-d9/minimap2-lr:hq/length_by_correctness-sim-r10-HG002.trimmed.1k.png
    ```

    And you can use the `--dry-run` flag to see all the rules that will be executed and the files that will be created. (You may need to include `--rerun-triggers mtime` to work around https://github.com/snakemake/snakemake/issues/3675 causing `Params have changed since last execution` reruns in dry run mode.)

    To test the Snakemake:

    ```
    snakemake all_paper_figures --dry-run --debug-dag
    ```

    And to run all the plots for the paper into the configured `all_out_dir` (expected to be group sticky):

    ```
    (umask 0002; snakemake -j128 --rerun-incomplete --use-singularity --singularity-args "-B /private" --latency-wait 120 --slurm all_paper_figures)
    ```

## Dependencies

To run all the rules, you will need to have installations of:

* wget
* curl
* vg
* minimap2
* pbmm2
* winnowmap
* samtools
* Java
* Picard (as `picard.jar` in the current directory)
* bcftools
* Singularity, for running fastqsplitter, [Truvari](https://github.com/ACEnglish/truvari), [mafft](https://mafft.cbrc.jp/alignment/software/installation_without_root.html) which Truvari calls, vcfwave, vcfbub, sniffles, and GraphAligner
* [toil](https://github.com/DataBiosphere/toil) for DeepVariant

## Available Experiment Outputs:

In these file name templates, `{root}` is your base output directory, `{expname}` is the name of the experiment defined in the config file, and `{ext}` is the image format you want the plot in, such as `png`.

* `{root}/experiments/{expname}/plots/correct.{ext}`: A bar chart of the fraction of eligible reads in each condition that are mapped correctly.
* `{root}/experiments/{expname}/plots/wrong.{ext}`: A bar chart of the fraction of eligible reads in each condition that are mapped incorrectly.
* `{root}/experiments/{expname}/plots/softclipped_or_unmapped.{ext}`: A bar chart of the number of bases left either softclipped or unmapped by each condition.
* `{root}/experiments/{expname}/plots/mapping_rate.{ext}`: A bar chart of the fraction of reads in each condition that are mapped.
* `{root}/experiments/{expname}/plots/mapping_speed.{ext}`: A bar chart of the speed of each Giraffe condition, in reads per second per thread.
* `{root}/experiments/{expname}/plots/chain_coverage.{ext}`: A bar chart of the best-chain coverage fraction of each Giraffe condition.
* `{root}/experiments/{expname}/plots/pr.{ext}`: A precision-recall plot for mapping accuracy showing each condition.
* `{root}/experiments/{expname}/plots/qq.{ext}`: A "QQ" plot with error bars, showing the calibration of mapping quality for detecting incorrectly-mapped reads.

## Useful Per-Condition Outputs

In these file name templates, `{root}` is your base output directory, `{ext}` is the image format you want the plot in, such as `png`, and the other placeholders are experimental variables that define the condition being run.

* `{root}/plots/{reference}/{refgraph}/{mapper}/best_chain_coverage-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of the amount of the read covered by the best chain, for Giraffe conditions.
* `{root}/plots/{reference}/{refgraph}/{mapper}/chain_anchor_length-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of the number of bases in seed anchors in the best chain.
* `{root}/plots/{reference}/{refgraph}/{mapper}/chain_anchors-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of the number of seed anchors in the best chain.
* `{root}/plots/{reference}/{refgraph}/{mapper}/time_used-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of the CPU time used to map each read, for Giraffe conditions.
* `{root}/plots/{reference}/{refgraph}/{mapper}/average_stage_time-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A chart of the average CPU time used per mapping stage, for Giraffe conditions.
* `{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_time-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A chart of the average wall clock time used per dynamic programming method, for Giraffe conditions.
* `{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_bases-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A chart of the average bases aligned per dynamic programming method, for Giraffe conditions.
* `{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_invocations-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A chart of the average invocation count per dynamic programming method, for Giraffe conditions.
* `{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_fraction-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A chart of the average read fraction aligned per dynamic programming method, for Giraffe conditions.
* `{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_speed-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A chart of the average bases aligned per second per dynamic programming method, for Giraffe conditions.
* `{root}/plots/{reference}/{refgraph}/{mapper}/average_aligner_probsize-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A chart of the average problem size averaged again over reads, per dynamic programming method, for Giraffe conditions.
* `{root}/plots/{reference}/{refgraph}/{mapper}/length_by_mapping-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of read length, broken out by whether the read mapped or not.
* `{root}/plots/{reference}/{refgraph}/{mapper}/length_by_correctness-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of read length, broken out by whether the read was correct, incorrect, or without a truth position.
* `{root}/stats/{reference}/{refgraph}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.facts.txt`: A Giraffe Facts report about where candidates are filtered, for Giraffe conditions.

# Parameter Search

This repo also includes machinery for doing parameter search experiments, to find the best combinations of parameters for `vg giraffe`, in order to generate Giraffe parameter presets.

To use it, first update `parameter_search_config.tsv` to describe the names of the parameters you want to vary, their Python type (`int` or `float`), their min, max, and current default values, and the sampling strategy to distribute sampled values (`uniform` or `log`). (The default value is used for previously-sampled conditions from before you added the parameter to the file.)

Then, run `parameter_search.py`, with the `--coount` option set to the number of points to sample in multi-dimensional parameter space. This populates `./hash_to_parameters.tsv` with several sample point hash values and their corresponding parameter sets.

Then use Snakemake to request any of the parameter search plots, such as `{root}/parameter_search/plots/chm13/hprc-v1.1-mc-sampled4o/giraffe-k29.w11.W-sr-default/HG002.100000/illumina.correct_speed_vs_hit-cap.png`. This will use the current default `vg` binary, with the `sr` preset, on 100000 Illumina reads, against a haplotype-sampled version of the HPRC v1.1 graph, and make a plot of simulated read correctness and real read mapping speed as a function of the `hit-cap` parameter (assuming it's defined as one of the parameters to vary). There are other Snakemake reles available for files in `parameter_search` that cna make parametric plots of pairs of real or simulated alignment statistics, or plot those statistics against parameter values. 

# Multi-User Operation

To collaborate with multiple people on a single set of intermediate files, first make a directory owned by a group they have in common. Set it to be group-writable and set the group sticky bit:
```
mkdir /private/groups/patenlab/project-lrg
chgrp patenlab /private/groups/patenlab/project-lrg
chmod g+w /private/groups/patenlab/project-lrg
chmod g+s /private/groups/patenlab/project-lrg
```

Then run the Snakemake in a subshell with a group-writable umask set, and ask for a file under that directory:
```
(umask 0002; snakemake -j128 --rerun-incomplete --use-singularity --singularity-args "-B /private" --latency-wait 120 --slurm /private/groups/patenlab/project-lrg/output/plots/chm13/hprc-v1.1-mc-d9/giraffe-k31.w50.W-hifi-default-noflags/time_used-real-hifi-HG002.1k.png)
```

The file and all intermediates will end up owned and writable by the correct group.


