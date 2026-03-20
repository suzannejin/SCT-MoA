# Plot overlap between different networks and single-cell coexpression networks. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(patchwork)
source("R/theme.R")
source("bin/functions.R")

# read data
dat = read.delim("results/overlap/overlap.txt.gz") %>%
  mutate(dataset = gsub("-expr.*$", "", dataset))

# set colors
colors = list("rho"         = "#f06c45",
              "pearson"     = "#b294c7",
              "spearman"    = "#f88a89",
              "zi_kendall"  = "#569ea4"
             )

# plot default cutoff (50k edges)
dat1 = filter(dat, cutoff == 5e4)
labeler = as_labeller(function(x) forcats::fct_recode(
  x, Signalling = "OmniPath", "Metabolic pathways" = "Reactome", 
  "Text mining" = "STRING", "PPIs" = "HIPPIE"))
labels = with(dat1, reorder(coefficient, z_score, median)) %>% levels()
p1 = ggplot(dat1, aes(x = reorder(coefficient, z_score, median), 
                      y = z_score, fill = coefficient)) + 
  facet_wrap(~ network, scales = 'free_x', ncol = 4) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous("Z score") +
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
ggsave("results/fig/overlap_1.pdf", p1, width = 18.3, height = 9, units = 'cm')

# plot 20k and 100k edges
dat2 = filter(dat, cutoff == 2e4) 
labels2 = with(dat2, reorder(coefficient, z_score, median)) %>% levels()
p2 = ggplot(dat2, aes(x = reorder(coefficient, z_score, median), 
                      y = z_score, fill = coefficient)) + 
  facet_wrap(~ network, scales = 'free_x', labeller = labeler, ncol = 4) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous("Z score") +
  scale_fill_manual(name = '', values=colors, guide = F) + 
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
p2

dat3 = filter(dat, cutoff == 1e5) 
labels3 = with(dat3, reorder(coefficient, z_score, median)) %>% levels()
p3 = ggplot(dat3, aes(x = reorder(coefficient, z_score, median), 
                      y = z_score, fill = coefficient)) + 
  facet_wrap(~ network, scales = 'free_x', labeller = labeler, ncol = 4) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous("Z score") +
  scale_fill_manual(name = '', values=colors, guide = F) + 
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
p3

# plot together
p4 = plot_grid(p2, p3, labels = letters, ncol = 1, label_size = 10)
ggsave("results/fig/overlap_2.pdf", p4, width = 18.3, height = 18, units = 'cm')


# calculate p-vals
for (db in unique(dat$network)){
  message("compute pvals for ", db)
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
      x = dat %>% dplyr::filter(network == db & coefficient == coef1 & cutoff==5e4) %>% pull(z_score)
      y = dat %>% dplyr::filter(network == db & coefficient == coef2 & cutoff==5e4) %>% pull(z_score)
      test = lawstat::brunner.munzel.test(x, y)
      pvals[coef1, coef2] = test$p.value
    }
  }
  write.table(pvals, paste0("results/overlap/pvals_", db, ".txt"), quote = F, row.names = T,
              sep = "\t")

}