#!/usr/bin/env Rscript
#
# plot-roc.R <stats TSV> <destination image file> [<comma-separated "aligner" names to include> [title]]
#
# plots a pseudo-ROC that allows the comparison of different alignment methods and their mapping quality calculations
# the format is clarified in the map-sim script, and should be a table (tab separated) of:
#      correct   mq    score   aligner
# where "correct" is 0 or 1 depending on whether the alignnment is correct or not and "aligner" labels the mapping method
#
# This is not a true ROC because we are not purely plotting the binary classification performance of
# each of the methods' mapping quality calculation over the same set of candidate alignments.
# Rather, we are mixing both the alignment sensitivity of the method with the MQ classification performance.
# As such we do not ever achieve 100% sensitivity, as we have effectively scaled the y axis (TPR) by the total
# sensitivity of each mapper.

list.of.packages <- c("tidyverse", "ggrepel", "svglite")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require("tidyverse")
require("ggrepel")
require("scales") # For squish

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
aligner.names <- levels(factor(dat$aligner))
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
colors <- c("#27A974","#9B8BF4","#33a02c","#C44601","#ff8000","#F2C300","#458b74","#698b22","#008b8b")

if (length(commandArgs(TRUE)) > 4) {
    # A set of colors to plot is specified. Parse it.
    colors <- unlist(strsplit(commandArgs(TRUE)[5], ","))
}

dat$bin <- cut(dat$mq, c(-Inf,seq(0,60,1),Inf))
dat.roc <- dat %>%
    mutate(Positive = (correct == 1) * count, Negative = (correct == 0) * count) %>%
    group_by(aligner, mq) %>%
    summarise(Positive = sum(Positive), Negative = sum(Negative)) %>%
    arrange(-mq) %>%
    mutate(Total=sum(Positive+Negative)) %>%
    mutate(TPR = cumsum(Positive) / Total, FPR = cumsum(Negative) / Total)
    
# We want smart scales that know how tiny a rate of things we can care about
total.reads <- max(dat.roc$Total)
min.log10 <- floor(log10(1/total.reads))
max.log10 <- 0
# Work out a set of bounds to draw the plot on
range.log10 <- min.log10 : max.log10
range.unlogged = 10^range.log10

if (length(commandArgs(TRUE)) > 5) {
    # A set of colors to plot is specified. Parse it.
    zoom_coords <- as.numeric(strsplit(commandArgs(TRUE)[6], ",")[[1]])
    print(zoom_coords)
    print(zoom_coords[1])
    dat.plot <- ggplot(dat.roc, aes( x= FPR, y = TPR, color = aligner, label=mq)) +
        geom_line() +
        #geom_text_repel(data = subset(dat.roc, mq %% 60 == 0), size=3.5, point.padding=unit(0.7, "lines"), segment.alpha=I(1/2.5)) +
        geom_point(aes(size=Positive+Negative)) +
        scale_color_manual(values=colors, guide=guide_legend(title=NULL, ncol=2)) +
        scale_size_continuous("number", guide=guide_legend(title=NULL, ncol=4)) +
        coord_trans(xlim = c(zoom_coords[1], zoom_coords[2]), x = "log",ylim = c(zoom_coords[3], zoom_coords[4])) +
        geom_vline(xintercept=1/total.reads) + # vertical line at one wrong read
        theme_bw() + 
        theme(legend.position=c(1,0), legend.justification=c(1,0)) +
        ggtitle(title)
} else {
 
    dat.plot <- ggplot(dat.roc, aes( x= FPR, y = TPR, color = aligner, label=mq)) +
        geom_line() +
        #geom_text_repel(data = subset(dat.roc, mq %% 60 == 0), size=3.5, point.padding=unit(0.7, "lines"), segment.alpha=I(1/2.5)) +
        geom_point(aes(size=Positive+Negative)) +
        scale_color_manual(values=colors, guide=guide_legend(title=NULL, ncol=2)) +
        scale_size_continuous("number", guide=guide_legend(title=NULL, ncol=4)) +
        scale_x_log10(limits=c(range.unlogged[1],range.unlogged[length(range.unlogged)]), breaks=range.unlogged, oob=squish) +
        geom_vline(xintercept=1/total.reads) + # vertical line at one wrong read
        theme_bw() + 
        theme(legend.position=c(1,0), legend.justification=c(1,0)) +
        ggtitle(title)
}
    
if (title != '') {
    # And a title
    dat.plot + ggtitle(title)
}
    
filename <- commandArgs(TRUE)[2]
pdf(filename, 4.5, 4)
dat.plot + guides(color="none", size="none")
dat.plot
dev.off()
