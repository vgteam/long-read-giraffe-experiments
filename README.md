# Long Read Giraffe Experiments

This repository contains scripts for running experiments to evaluate the long-read mode for the `vg giraffe` pangenome sequencing read mapper. It knows how to run Giraffe, Winnowmap, and Minimap2, and how to analyze the results with `vg` and the included scripts.

## Usage

The general flow for using this repository is:

1. Obtain some real long reads for your sequencing technologies.

2. Use `make_pbsim_reads.sh` to simulate reads with similar base qualities for a sample with set genotypes.

3. Lay out your input files in the expected format, including a GBZ graph and distance index, the reads, a linear reference FASTA, and Winnowmap indexes, and any Minimap2 indexes. See `Snakefile` for a description of the required layout for input data.

4. Write a YAML configuration file defining one or more experiments in which different variables (sequencing technology, number of reads, read aligner used, etc.) are either held constant or varied. See `lr-config.yaml`, the default config file, for an example and documentation for all the experimental variables and their possible values. You can also use the config file to set `graphs_dir`, `reads_dir`, and `refs_dir` to directories where your graphs, reads, and Winnowmap2/Minimap indexed linear references are stored, if not running with access to the UCSC Genomics Institute `/private/groups` directory hierarchy. 

5. Use Snakemake and the included `Snakefile` to do the work necessary to create the experimental results files you are interested in. For example if you write `lr-config.yaml` to define an experiment named `r10_accuracy_small` to compare Winnowmap and Giraffe on some R10 nanopore reads, you can get a plot of the rate of mapped reads across the two conditions in the experiment with:

    ```
    snakemake --rerun-incomplete output/experiments/r10_accuracy_small/plots/mapping_rate.png
    ```

    Make sure that `output` is a directory with enough free space to store the mapped reads. We use the `--rerun-incomplete` flag to make Snakemake rerun any commands that were interrupted on previous runs.

    You can also run jobs on Slurm, where they will be automatically assigned to the correct partitions for UCSC's Phoenix cluster:

    ```
    snakemake --rerun-incomplete --slurm --latency-wait 120 -j128 output/experiments/r10_accuracy_small/plots/mapping_rate.png
    ```

    This will run on Slurm, waiting up to 120 seconds for files created on the worker nodes to become visible to the head node, and using up to 128 cores total.

    You can also ask for multiple output files in a single run:

    ```
    snakemake --rerun-incomplete --slurm --latency-wait 120 -j128 \
        output/experiments/r10_accuracy_small/plots/mapping_rate.png \
        output/experiments/r10_accuracy_small/plots/correct.png \
        output/plots/chm13/minimap2/length_by_correctness-sim-r10-HG002.trimmed.1k.png
    ```

    And you can use the `--dry-run` flag to see all the rules that will be executed and the files that will be created.

## Available Experiment Outputs:

In these file name templates, `{root}` is your base output directory, `{expname}` is the name of the experiment defined in the config file, and `{ext}` is the image format you want the plot in, such as `png`.

* `{root}/experiments/{expname}/plots/correct.{ext}`: A bar chart of the fraction of eligible reads in each condition that are mapped correctly.
* `{root}/experiments/{expname}/plots/mapping_rate.{ext}`: A bar chart of the fraction to reads in each condition that are mapped.
* `{root}/experiments/{expname}/plots/pr.{ext}`: A precision-recall plot for mapping accuracy showing each condition.
* `{root}/experiments/{expname}/plots/qq.{ext}`: A "QQ" plot with error bars, showing the calibration of mapping quality for detecting incorrectly-mapped reads.

## Useful Per-Condition Outputs

In these file name templates, `{root}` is your base output directory, `{ext}` is the image format you want the plot in, such as `png`, and the other placeholders are experimental variables that define the condition being run.

* `{root}/plots/{reference}/{mapper}/best_chain_coverage-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of the amount of the read covered by the best chain, for Giraffe conditions.
* `{root}/plots/{reference}/{mapper}/time_used-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of the CPU time used to map each read, for Giraffe conditions.
* `{root}/plots/{reference}/{mapper}/average_stage_time-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A chart of the average CPU time used per mapping stage, for Giraffe conditions.
* `{root}/plots/{reference}/{mapper}/length_by_mapping-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of read length, broken out by whether the read mapped or not.
* `{root}/plots/{reference}/{mapper}/length_by_correctness-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}`: A histogram of read length, broken out by whether the read was correct, incorrect, or without a truth position.

## Missing Features

Right now, Giraffe isn't able to be run with combinations of additional flags under the experiment system.


