#!/usr/bin/env Rscript

# plot-qq.R <stats TSV> <destination image file> [<comma-separated "aligner" names to include> [title]]

list.of.packages <- c("tidyverse", "ggrepel", "svglite", "binom")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require("tidyverse")
require("ggrepel")
require("binom")

# Read in the combined toil-vg stats.tsv, listing:
# correct, mapq, aligner (really graph name), read name, count, eligible
dat <- read.table(commandArgs(TRUE)[1], header=T, colClasses=c("aligner"="factor"))

if (("eligible" %in% names(dat))) {
    # If the eligible column is present, remove ineligible reads
    dat <- dat[dat$eligible == 1, ]
}

if (! ("count" %in% names(dat))) {
    # If the count column is not present, add i
    dat$count <- rep(1, nrow(dat))
}

if (length(commandArgs(TRUE)) > 2) {
    # A set of aligners to plot is specified. Parse it.
    aligner.set <- unlist(strsplit(commandArgs(TRUE)[3], ";"))
    # Subset the data to those aligners
    dat <- dat[dat$aligner %in% aligner.set,]
    # And restrict the aligner factor levels to just the ones in the set
    dat$aligner <- factor(dat$aligner, levels=aligner.set)
}

# Determine title
title <- ''
if (length(commandArgs(TRUE)) > 3) {
    title <- commandArgs(TRUE)[4]
}

# Determine the order of aligners, based on sorting in a dash-separated tag aware manner
aligner.names <- levels(dat$aligner)
#name.lists <- aligner.names %>% (function(name) map(name,  (function(x) as.list(unlist(strsplit(x, "-"))))))
## Transpose name fragments into a list of vectors for each position, with NAs when tag lists end early
#max.parts <- max(sapply(name.lists, length))
#name.cols <- list()
#for (i in 1:max.parts) {
#    name.cols[[i]] <- sapply(name.lists, function(x) if (length(x) >= i) { x[[i]] } else { NA })
#}
#name.order <- do.call(order,name.cols)
#aligner.names <- aligner.names[name.order]
dat$aligner <- factor(dat$aligner, levels=aligner.names)
#name.lists <- name.lists[name.order]

# Determine colors for aligners
colors <- c("#1f78b4","#e31a1c","#33a02c","#6600cc","#ff8000","#5c415d","#458b74","#698b22","#008b8b","#6caed1")

if (length(commandArgs(TRUE)) > 4) {
    # A set of colors to plot is specified. Parse it.
    colors <- unlist(strsplit(commandArgs(TRUE)[5], ","))
}


dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))

x <- as.data.frame(summarize(group_by(dat, bin, aligner), N=n(), mapq=mean(mq), mapprob=mean(1-10^(-mapq/10)), observed=weighted.mean(correct, count), select(binom.confint(sum(correct * count), sum(count), conf.level=0.9, methods="lrt"), c("lower", "upper"))))

print(names(x))
print(x$ci)
print(colors)

# Now plot the points as different sizes, but the error bar line ranges as a consistent size
dat.plot <- ggplot(x, aes(1-mapprob+1e-7, 1-observed+1e-7, color=aligner, size=N, weight=N, label=round(mapq,2))) +
    scale_color_manual(values=colors, guide=guide_legend(title=NULL, ncol=1)) +
    scale_y_log10("measured error", limits=c(1e-7,2), breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0)) +
    scale_x_log10("error estimate", limits=c(1e-7,2), breaks=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0)) +
    scale_size_continuous("number", guide=guide_legend(title=NULL, ncol=4)) +
    geom_point() +
    # Only aesthetics that depend on each point need to be in the aes() mapping
    geom_linerange(aes(x=1-mapprob+1e-7, ymin=1-upper+1e-7, ymax=1-lower+1e-7), linewidth=0.2, position=position_dodge(.05)) +
    geom_smooth() +
    geom_abline(intercept=0, slope=1, linetype=2) +
    theme_bw()
    
if (title != '') {
    # And a title
    dat.plot + ggtitle(title)
}

filename <- commandArgs(TRUE)[2]
ggsave(filename, height=4, width=7)
