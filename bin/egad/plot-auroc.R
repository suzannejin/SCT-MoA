# Plot and statistically test differences in median AUC between measures of 
# association.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(lawstat)
library(aggregation)
source("R/theme.R")
source("bin/functions.R")

# set colors
colors = list("rho"         = "#f06c45",
              "pearson"     = "#b294c7",
              "spearman"    = "#f88a89",
              "zi_kendall"  = "#569ea4"
             )

message("plotting main figure")

# read data
dat = read.delim("results/egad/auroc.txt.gz") %>%
  mutate(dataset = gsub("-expr.*$", "", dataset))
med = dat %>%
  group_by(dataset, coefficient, network) %>%
  dplyr::summarise(median_auc = median(auroc)) %>%
  ungroup() %>%
  mutate(dataset = gsub("-expr.*$", "", dataset))

# main plot
labels = levels(with(med, reorder(coefficient, median_auc, median)))
p1 = ggplot(med, aes(x = reorder(coefficient, median_auc, median), 
                      y = median_auc, fill = coefficient)) + 
  facet_wrap(~ network, ncol = 3) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
  scale_y_continuous("AUC", limits = c(0.45, 0.7)) + 
  scale_fill_manual(name = '', values = colors, guide = F) + 
  coord_flip() + 
  clean_theme + 
  theme(legend.position = 'right',
        axis.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major.y = element_line(color = 'grey80', 
                                          linetype = 'dotted'),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color = 'grey50'),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = 'grey90', color = 'white'))
p1
ggsave("results/fig/egad.pdf", p1, width = 18.5, height = 9, units = 'cm')


# write p-values
for (db in c("GO","GOslim","Reactome")){

  message("plot figure for ", db)
  # read data
  file = paste0("results/egad/auroc_", db, ".txt.gz")
  dat = read.delim(file)
  # get median for each dataset
  med = dat %>%
    group_by(dataset, coefficient) %>%
    dplyr::summarise(median_auc = median(auroc)) %>%
    ungroup() %>%
    mutate(dataset = gsub("-expr.*$", "", dataset))
  # clean up coefficeint
  # med$coefficient %<>% clean_metric() %>% as.character()

  # filter to one dataset per publication
  opp = read.delim("data/geo/one-dataset-per-publication.txt")
  med0 = filter(med, dataset %in% opp$Dataset)

  # load other Fig 1 panels
  # load("fig/Rdata/datasets.Rdata", verbose = T)

  # plot
  labels = levels(with(med0, reorder(coefficient, -median_auc, median)))
  p6 = ggplot(med0, aes(x = reorder(coefficient, -median_auc, median), 
                      y = median_auc, fill = coefficient)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
    scale_fill_manual(values = colors, guide = F) +
    scale_y_continuous("AUC", limits = c(0.45, 0.65), expand = c(0, 0)) + 
    theme_sc + 
    theme(axis.title.x = element_blank(),
          axis.line.y = element_line(color = 'grey50'),
          panel.grid.major.x = element_line(color = 'grey80', 
                                            linetype = 'dotted'),
          axis.ticks.x = element_blank(),
          plot.background = element_blank())
  p6
  ggsave(paste0("results/fig/egad_", db, ".pdf"), p6, width = 10, height = 7, units = 'cm')

  # print median AUCs
  med0 %>% 
    group_by(coefficient) %>%
    summarise(median = median(median_auc)) %>%
    arrange(desc(median))

  # run statistical tests 
  message("testing p-value for ", db)
  coefs = labels
  pvals = matrix(NA, nrow = length(coefs), ncol = length(coefs), 
                dimnames = list(coefs, coefs))
  for (i in seq_len(length(coefs))) {
    coef1 = coefs[i]
    message("analyzing coefficient ", coef1, " ...")
    for (j in seq_len(length(coefs))) {
      coef2 = coefs[j]
      if (coef2 == coef1)
        next
      message("  analyzing coefficient ", coef2, " ...")
      # run Brunner--Munzel tests within each dataset
      p = numeric(0)
      for (d in unique(med0$dataset)) {
        x = dat %>% dplyr::filter(dataset == d & coefficient == coef1) %>% 
          pull(auroc)
        y = dat %>% dplyr::filter(dataset == d & coefficient == coef2) %>% 
          pull(auroc)
        test = lawstat::brunner.munzel.test(x, y)
        p[d] = test$p.value
      }
      # Fisher integration 
      fisher = aggregation::fisher(p)
      pvals[coef1, coef2] = fisher
    }
  }
  write.table(pvals, paste0("results/egad/auroc_pvals_", db, ".txt"), quote = F, row.names = T,
              sep = "\t")
}
