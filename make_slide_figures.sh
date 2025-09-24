#!/usr/bin/env bash

CATEGORIES=(minimap2-map-hifi,hprc-v2.0-mc-eval-d46 pbmm2,hprc-v2.0-mc-eval-d46 giraffe,primary minigraph,hprc-v2.0-minigraph-eval graphaligner-fast,hprc-v2.0-mc-eval-d46 giraffe,hprc-v2.0-mc-eval-d46 giraffe,hprc-v2.0-mc-eval.ec1M-sampled16o)
COLORS=('#9B8BF4' '#B3C7F7' '#53AEDF' '#6396B0' '#D066A3' '#C44601' '#F57600')
LABELS=("Minimap2" "pbmm2" "$(printf 'Giraffe\n(Linear)')" "Minigraph" "$(printf 'GraphAligner\n(Fast)')" "Giraffe" "$(printf 'Giraffe\n(Personal)')")
grep -E $(echo "${CATEGORIES[@]}" | sed 's/ /|/g') /private/groups/patenlab/project-lrg/experiments/hifi_real_full/results/mapping_stats_real.tsv > hifi-runtime.tsv
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' hifi-runtime.tsv) --title "$(printf 'HiFi Runtime\n(66x coverage, 64 cores)')" --y_label 'Time (hours)' --x_label Mapper --no_n --categories "${CATEGORIES[@]}" --colors "${COLORS[@]}" --category_labels "${LABELS[@]}" --x_sideways --save hifi_runtime_ga.png --width 3 --height 3 --dpi 200 && imgcat hifi_runtime_ga.png --height 20
python3 barchart.py <(tail -n +2 /private/groups/patenlab/project-lrg/experiments/hifi_sim_1m/results/mapping_stats_sim.tsv | cut -f1,5)  --title 'HiFi Wrong Reads' --y_label 'Count' --x_label 'Mapper' --x_sideways --no_n --categories "${CATEGORIES[@]}" --colors "${COLORS[@]}" --category_labels "${LABELS[@]}" --save hifi_incorrect_ga.png --width 3 --height 3 --dpi 200 && imgcat hifi_incorrect_ga.png --height 20

CATEGORIES=(minimap2-map-hifi,hprc-v2.0-mc-eval-d46 pbmm2,hprc-v2.0-mc-eval-d46 giraffe,primary minigraph,hprc-v2.0-minigraph-eval giraffe,hprc-v2.0-mc-eval-d46 giraffe,hprc-v2.0-mc-eval.ec1M-sampled16o)
COLORS=('#9B8BF4' '#B3C7F7' '#53AEDF' '#6396B0' '#C44601' '#F57600')
LABELS=("Minimap2" "pbmm2" "$(printf 'Giraffe\n(Linear)')" "Minigraph" "Giraffe" "$(printf 'Giraffe\n(Personal)')")
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' hifi-runtime.tsv) --title "$(printf 'HiFi Runtime\n(66x coverage, 64 cores)')" --y_label 'Time (hours)' --x_label Mapper --no_n --categories "${CATEGORIES[@]}" --colors "${COLORS[@]}" --category_labels "${LABELS[@]}" --x_sideways --save hifi_runtime_noga.png --width 3 --height 3 --dpi 200 && imgcat hifi_runtime_noga.png --height 20
python3 barchart.py <(tail -n +2 /private/groups/patenlab/project-lrg/experiments/hifi_sim_1m/results/mapping_stats_sim.tsv | cut -f1,5)  --title 'HiFi Wrong Reads' --y_label 'Count' --x_label 'Mapper' --x_sideways --no_n --categories "${CATEGORIES[@]}" --colors "${COLORS[@]}" --category_labels "${LABELS[@]}" --save hifi_incorrect_noga.png --width 3 --height 3 --dpi 200 && imgcat hifi_incorrect_noga.png --height 20


CATEGORIES=(minimap2-lr:hqae,hprc-v2.0-mc-eval-d46 giraffe,primary minigraph,hprc-v2.0-minigraph-eval graphaligner-fast,hprc-v2.0-mc-eval-d46 giraffe,hprc-v2.0-mc-eval-d46 giraffe,hprc-v2.0-mc-eval.ec1M-sampled16o)
COLORS=('#9B8BF4' '#53AEDF' '#6396B0' '#D066A3' '#C44601' '#F57600')
LABELS=("Minimap2" "$(printf 'Giraffe\n(Linear)')" "Minigraph" "$(printf 'GraphAligner\n(Fast)')" "Giraffe" "$(printf 'Giraffe\n(Personal)')")
grep -E $(echo "${CATEGORIES[@]}" | sed 's/ /|/g') /private/groups/patenlab/project-lrg/experiments/r10y2025_real_full/results/mapping_stats_real.tsv > r10-runtime.tsv
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' r10-runtime.tsv) --title "$(printf 'R10 Runtime\n(50x coverage, 64 cores)')" --y_label 'Time (hours)' --x_label Mapper --no_n --categories "${CATEGORIES[@]}" --colors "${COLORS[@]}" --category_labels "${LABELS[@]}" --x_sideways --save r10_runtime_ga.png --width 3 --height 3 --dpi 200 && imgcat r10_runtime_ga.png --height 20
python3 barchart.py <(tail -n +2 /private/groups/patenlab/project-lrg/experiments/r10y2025_sim_1m/results/mapping_stats_sim.tsv | cut -f1,5)  --title 'R10 Wrong Reads' --y_label 'Count' --x_label 'Mapper' --x_sideways --no_n --categories "${CATEGORIES[@]}" --colors "${COLORS[@]}" --category_labels "${LABELS[@]}" --save r10_incorrect_ga.png --width 3 --height 3 --dpi 200 && imgcat r10_incorrect_ga.png --height 20

CATEGORIES=(minimap2-lr:hqae,hprc-v2.0-mc-eval-d46 giraffe,primary minigraph,hprc-v2.0-minigraph-eval giraffe,hprc-v2.0-mc-eval-d46 giraffe,hprc-v2.0-mc-eval.ec1M-sampled16o)
COLORS=('#9B8BF4' '#53AEDF' '#6396B0' '#C44601' '#F57600')
LABELS=("Minimap2" "$(printf 'Giraffe\n(Linear)')" "Minigraph" "Giraffe" "$(printf 'Giraffe\n(Personal)')")
grep -E $(echo "${CATEGORIES[@]}" | sed 's/ /|/g') /private/groups/patenlab/project-lrg/experiments/r10y2025_real_full/results/mapping_stats_real.tsv > r10-runtime.tsv
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' r10-runtime.tsv) --title "$(printf 'R10 Runtime\n(50x coverage, 64 cores)')" --y_label 'Time (hours)' --x_label Mapper --no_n --categories "${CATEGORIES[@]}" --colors "${COLORS[@]}" --category_labels "${LABELS[@]}" --x_sideways --save r10_runtime_noga.png --width 3 --height 3 --dpi 200 && imgcat r10_runtime_noga.png --height 20
python3 barchart.py <(tail -n +2 /private/groups/patenlab/project-lrg/experiments/r10y2025_sim_1m/results/mapping_stats_sim.tsv | cut -f1,5)  --title 'R10 Wrong Reads' --y_label 'Count' --x_label 'Mapper' --x_sideways --no_n --categories "${CATEGORIES[@]}" --colors "${COLORS[@]}" --category_labels "${LABELS[@]}" --save r10_incorrect_noga.png --width 3 --height 3 --dpi 200 && imgcat r10_incorrect_noga.png --height 20

R --no-save <<'EOF'

library(dplyr)
library(ggplot2)
library(tidyr)

#Input is a tsv of condition name, followed by tp, fn, fp, recall, precision, f1
dat <- read.table("/private/groups/patenlab/project-lrg/experiments/dv_calling/results/dv_snp_summary.tsv", header=T)

categories <- c("hifi,minimap2-map-hifi,eval-d46,.nomodel", "hifi,pbmm2,eval-d46,.nomodel", "hifi,giraffe-hifi-mmcs100,eval-d46,.model2025-03-26noinfo", "hifi,giraffe-hifi-mmcs100,eval.ec1M-sampled16o,.model2025-03-26noinfo")
colors <- c("#9B8BF4", "#B3C7F7", "#C44601", "#F57600")
category.names <- c("minimap2", "pbmm2", "Giraffe", "Giraffe\n(Personal)")

# A set of aligners to plot is specified. Parse it.
condition.set <- unlist(categories)
# Subset the data to those aligners
dat <- dat[dat$condition %in% condition.set,]
# And restrict the condition factor levels to just the ones in the set
dat$condition <- factor(dat$condition, levels=condition.set)

# Then rename them to the human-readable names
levels(dat$condition) <- category.names

## save a PNG
png("hifi_dv_snp_plot.png", 3, 3, units="in", res=300)

dat %>% 
  select(condition, FP, FN) %>% 
  pivot_longer(cols=c(FP, FN), names_to='error', values_to='count') %>%
  ggplot(aes(x=condition, fill=error, y=count)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  labs(x="Mapper", y="Errors", title="HiFi SNP Errors", fill="Error Type") + 
  theme_bw() +
  theme(legend.position='bottom') +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))


dev.off()

EOF
imgcat "hifi_dv_snp_plot.png" --height 20

R --no-save <<'EOF'

library(dplyr)
library(ggplot2)
library(tidyr)

#Input is a tsv of condition name, followed by tp, fn, fp, recall, precision, f1
dat <- read.table("/private/groups/patenlab/project-lrg/experiments/dv_calling/results/dv_indel_summary.tsv", header=T)

categories <- c("hifi,minimap2-map-hifi,eval-d46,.nomodel", "hifi,pbmm2,eval-d46,.nomodel", "hifi,giraffe-hifi-mmcs100,eval-d46,.model2025-03-26noinfo", "hifi,giraffe-hifi-mmcs100,eval.ec1M-sampled16o,.model2025-03-26noinfo")
colors <- c("#9B8BF4", "#B3C7F7", "#C44601", "#F57600")
category.names <- c("minimap2", "pbmm2", "Giraffe", "Giraffe\n(Personal)")

# A set of aligners to plot is specified. Parse it.
condition.set <- unlist(categories)
# Subset the data to those aligners
dat <- dat[dat$condition %in% condition.set,]
# And restrict the condition factor levels to just the ones in the set
dat$condition <- factor(dat$condition, levels=condition.set)

# Then rename them to the human-readable names
levels(dat$condition) <- category.names

## save a PNG
png("hifi_dv_indel_plot.png", 3, 3, units="in", res=300)

dat %>% 
  select(condition, FP, FN) %>% 
  pivot_longer(cols=c(FP, FN), names_to='error', values_to='count') %>%
  ggplot(aes(x=condition, fill=error, y=count)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  labs(x="Mapper", y="Errors", title="HiFi Indel Errors", fill="Error Type") + 
  theme_bw() +
  theme(legend.position='bottom') +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))


dev.off()

EOF
imgcat "hifi_dv_indel_plot.png" --height 20

R --no-save <<'EOF'

library(dplyr)
library(ggplot2)
library(tidyr)

#Input is a tsv of condition name, followed by tp, fn, fp, recall, precision, f1
dat <- read.table("/private/groups/patenlab/project-lrg/experiments/dv_calling/results/dv_snp_summary.tsv", header=T)

categories <- c("r10y2025,minimap2-lr:hqae,eval-d46,.nomodel", "r10y2025,giraffe-r10-noflags,eval-d46,.nomodel", "r10y2025,giraffe-r10-noflags,eval.ec1M-sampled16o,.nomodel")
colors <- c("#9B8BF4", "#C44601", "#F57600")
category.names <- c("minimap2", "pbmm2", "Giraffe", "Giraffe\n(Personal)")

# A set of aligners to plot is specified. Parse it.
condition.set <- unlist(categories)
# Subset the data to those aligners
dat <- dat[dat$condition %in% condition.set,]
# And restrict the condition factor levels to just the ones in the set
dat$condition <- factor(dat$condition, levels=condition.set)

# Then rename them to the human-readable names
levels(dat$condition) <- category.names

## save a PNG
png("r10_dv_snp_plot.png", 3, 3, units="in", res=300)

dat %>% 
  select(condition, FP, FN) %>% 
  pivot_longer(cols=c(FP, FN), names_to='error', values_to='count') %>%
  ggplot(aes(x=condition, fill=error, y=count)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  labs(x="Mapper", y="Errors", title="R10 SNP Errors", fill="Error Type") + 
  theme_bw() +
  theme(legend.position='bottom') +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))


dev.off()

EOF
imgcat "r10_dv_snp_plot.png" --height 20

R --no-save <<'EOF'

library(dplyr)
library(ggplot2)
library(tidyr)

#Input is a tsv of condition name, followed by tp, fn, fp, recall, precision, f1
dat <- read.table("/private/groups/patenlab/project-lrg/experiments/dv_calling/results/dv_indel_summary.tsv", header=T)

categories <- c("r10y2025,minimap2-lr:hqae,eval-d46,.nomodel", "r10y2025,giraffe-r10-noflags,eval-d46,.nomodel", "r10y2025,giraffe-r10-noflags,eval.ec1M-sampled16o,.nomodel")
colors <- c("#9B8BF4", "#C44601", "#F57600")
category.names <- c("minimap2", "Giraffe", "Giraffe\n(Personal)")

# A set of aligners to plot is specified. Parse it.
condition.set <- unlist(categories)
# Subset the data to those aligners
dat <- dat[dat$condition %in% condition.set,]
# And restrict the condition factor levels to just the ones in the set
dat$condition <- factor(dat$condition, levels=condition.set)

# Then rename them to the human-readable names
levels(dat$condition) <- category.names

## save a PNG
png("r10_dv_indel_plot.png", 3, 3, units="in", res=300)

dat %>% 
  select(condition, FP, FN) %>% 
  pivot_longer(cols=c(FP, FN), names_to='error', values_to='count') %>%
  ggplot(aes(x=condition, fill=error, y=count)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  labs(x="Mapper", y="Errors", title="R10 Indel Errors", fill="Error Type") + 
  theme_bw() +
  theme(legend.position='bottom') +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))


dev.off()

EOF
imgcat "r10_dv_indel_plot.png" --height 20

