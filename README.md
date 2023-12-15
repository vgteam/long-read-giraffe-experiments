# Long Read Giraffe Experiments

This repository contains scripts for running experiments to evaluate the long-read mode for the `vg giraffe` pangenome sequencing read mapper. It knows how to run Giraffe, Winnowmap, and Minimap2, and how to analyze the results with `vg` and the included scripts.

The general flow for using this repository is:

1. Obtain some real long reads for your sequencing technologies.

2. Use `make_pbsim_reads.sh` to simulate reads with similar base qualities for a sample with set genotypes.

3. Write a YAML configuration file defining one or more experiments in which different variables (sequencing technology, number of reads, read aligner used, etc.) are either held constant or varied. See `lr-config.yaml` for an example. You can also use the config file to set `graphs_dir`, `reads_dir`, and `refs_dir` to directories where your graphs, reads, and Winnowmap2/Minimap indexed linear references are stored, if not running with access to the UCSC Genomics Institute `/private/groups` directory hierarchy.

TODO: Document the experiment definition semantics.

TODO: Document the required input file layout.

4. Use Snakemake and the included `Snakefile` to do the work necessary to create the experimental results files you are interested in. For example if you write `lr-config.yaml` to define an experiment named `exp1` to comapre Winnowmap and Giraffe on some HiFi reads, you can get a plot of the rate of mapped reads across the two conditions in the experiment with:

```
snakemake --config lr-config.yaml my/out/dir/experiments/exp1/plots/mapping_rate.png
```

TODO: Document all the available output files.
