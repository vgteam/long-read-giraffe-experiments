library(dplyr)
library(ggplot2)
library(tidyr)

# plot-calling-results.R <stats TSV> <destination image file> [<semicolon-separated "condition" names to include> [title]]

#Input is a tsv of condition name, followed by tp, fn, fp, recall, precision, f1
dat <- read.table(commandArgs(TRUE)[1], header=T)

if (length(commandArgs(TRUE)) > 2) {
    # A set of aligners to plot is specified. Parse it.
    condition.set <- unlist(strsplit(commandArgs(TRUE)[3], ";"))
    # Subset the data to those aligners
    dat <- dat[dat$condition %in% condition.set,]
    # And restrict the condition factor levels to just the ones in the set
    dat$condition <- factor(dat$condition, levels=condition.set)
}

## save a PDF
pdf(commandArgs(TRUE)[2], 7, 3)

## comparing F1
## Without name
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
  theme(legend.position='bottom', axis.text.y=element_blank())

## With name
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
## Without name
dat %>% 
  select(condition, FP, FN) %>% 
  pivot_longer(cols=c(FP, FN), names_to='error', values_to='count') %>%
  ggplot(aes(x=condition, fill=error, y=count)) +
  geom_col() +
  scale_fill_brewer(palette="Set2") +
  ylab('variant call') + 
  theme_bw() +
  coord_flip() + 
  theme(legend.position='bottom', axis.text.y=element_blank())

## With name
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
