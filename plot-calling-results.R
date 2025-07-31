library(dplyr)
library(ggplot2)
library(tidyr)

# plot-calling-results.R <stats TSV> <destination image file> [<comma-separated "aligner" names to include> [title]]

#Input is a tsv of condition name, followed by tp, fn, fp, recall, precision, f1
dat <- read.table(commandArgs(TRUE)[1], header=T)


## save a PDF
pdf(commandArgs(TRUE)[2], 9, 3)

## comparing F1
dat %>% 
  select(condition, F1, recall, precision) %>% 
  pivot_longer(cols=c(F1, recall, precision), names_to='metric', values_to='value') %>%
  ggplot(aes(x=condition, color=metric, y=value)) +
  geom_hline(yintercept=1, linetype=2) + 
  geom_line(aes(group=condition), linewidth=0.75, color='black', alpha=.8) +
  geom_point( alpha=.8, size=3) +
  theme_bw() +
  scale_color_brewer(palette="Set2") +
  coord_flip() + 
  theme(legend.position='bottom')

## comparing errors (FP and FN)
dat %>% 
  select(condition, FP, FN) %>% 
  pivot_longer(cols=c(FP, FN), names_to='error', values_to='count') %>%
  ggplot(aes(x=condition, fill=error, y=count)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  ylab('variant call') + 
  theme_bw() +
  coord_flip() + 
  theme(legend.position='bottom')


dev.off()
