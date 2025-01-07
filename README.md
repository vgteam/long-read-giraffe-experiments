# Long Read Giraffe Experiments

This repository contains scripts for running experiments to evaluate the long-read mode for the `vg giraffe` pangenome sequencing read mapper. It knows how to run Giraffe, Winnowmap, and Minimap2, and how to analyze the results with `vg` and the included scripts.

## Usage

The general flow for using this repository is:

1. Obtain some real long reads for your sequencing technologies.

2. Use `make_pbsim_reads.sh` to simulate reads with similar base qualities for a sample with set genotypes. For Illumina reads, you should simulate paired-end reads manually with `vg sim` like this:

```
GRAPH_DIR=/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs
READ_DIR=/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads
mkdir -p "${READ_DIR}/sim/illumina/HG002"
# Prepare a separete GBWT
vg gbwt -o"${GRAPH_DIR}/hprc-chm-hg002-2024-03-25-mc.full.gbwt" -g "${GRAPH_DIR}/hprc-chm-hg002-2024-03-25-mc.full.gg" -Z "${GRAPH_DIR}/hprc-chm-hg002-2024-03-25-mc.full.gbz"
# Prepare an xg graph
vg convert --drop-haplotypes --xg-out "${GRAPH_DIR}//hprc-chm-hg002-2024-03-25-mc.full.gbz" >"${GRAPH_DIR}//hprc-chm-hg002-2024-03-25-mc.full.xg"
# Simulate reads (assuming you already put 3M reads in "${READ_DIR}/real/illumina/HG002/HG002.novaseq.pcr-free.40x.3m.fq")
vg sim -r -n 2500000 -a -s 12345 -p 570 -v 165 -i 0.00029 -x "${GRAPH_DIR}//hprc-chm-hg002-2024-03-25-mc.full.xg" -g "${GRAPH_DIR}//hprc-chm-hg002-2024-03-25-mc.full.gbwt" --sample-name HG002 -F "${READ_DIR}/real/illumina/HG002/HG002.novaseq.pcr-free.40x.3m.fq" --multi-position > "${READ_DIR}/sim/illumina/HG002/HG002-sim-illumina.gam"
# Subset reads
for READ_COUNT in 100 1000 10000 100000 1000000 ; do
    vg filter --interleaved -t1 --max-reads "${READ_COUNT}" "${READ_DIR}/sim/illumina/HG002/HG002-sim-illumina.gam" >"${READ_DIR}/sim/illumina/HG002/HG002-sim-illumina-${READ_COUNT}.gam"
done
```

Note that the Illumina read subsets made like this won't be uniformly sampled from the full set.

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
        output/plots/chm13/hprc-v1.1-mc-d9/minimap2-lr:hq/length_by_correctness-sim-r10-HG002.trimmed.1k.png
    ```

    And you can use the `--dry-run` flag to see all the rules that will be executed and the files that will be created.

    To test the Snakemake:

    ```
    snakemake all_paper_figures --dry-run --debug-dag
    ```

    And to run all the plots for the paper into the configured `all_out_dir` (expected to be group sticky):

    ```
    (umask 0002; snakemake -j128 --rerun-incomplete --latency-wait 120 --slurm all_paper_figures)
    ```

## Dependencies

To run all the rules, you will need to have installations of:

* wget
* curl
* vg
* minimap2
* winnowmap
* GraphAligner
* samtools
* Java
* Picard (as `picard.jar` in the current directory)
* bcftools
* [Truvari](https://github.com/ACEnglish/truvari)
* [mafft](https://mafft.cbrc.jp/alignment/software/installation_without_root.html) for Truvari to call
* [toil](https://github.com/DataBiosphere/toil) for DeepVariant

## Available Experiment Outputs:

In these file name templates, `{root}` is your base output directory, `{expname}` is the name of the experiment defined in the config file, and `{ext}` is the image format you want the plot in, such as `png`.

* `{root}/experiments/{expname}/plots/correct.{ext}`: A bar chart of the fraction of eligible reads in each condition that are mapped correctly.
* `{root}/experiments/{expname}/plots/wrong.{ext}`: A bar chart of the fraction of eligible reads in each condition that are mapped incorrectly.
* `{root}/experiments/{expname}/plots/softclipped_or_unmapped.ext`: A bar chart of the number of bases left either softclipped or unmapped by each condition.
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
(umask 0002; snakemake -j128 --rerun-incomplete --latency-wait 120 --slurm /private/groups/patenlab/project-lrg/output/plots/chm13/hprc-v1.1-mc-d9/giraffe-k31.w50.W-hifi-default-noflags/time_used-real-hifi-HG002.1k.png)
```

The file and all intermediates will end up owned and writable by the correct group.


