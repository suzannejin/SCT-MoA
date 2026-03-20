# Plot and statistically test differences in median AUC between measures of 
# association.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(lawstat)
library(aggregation)
source("R/theme.R")
source("bin/functions.R")

# set method colors
methods = c("proportionality", "pearson", "spearman", "zi_kendall")
colors =  c("proportionality"       = "#4892C2",
            "pearson"               = "#899D5B",
            "spearman"              = "#52AE43",
            "zi_kendall"            = "#559EA3"
            )
colors1 = list("proportionality" = "#4892C2", 
               "pearson"         = "#899D5B",
               "pearsonCLR"      = "#899D5B",
               "spearman"        = "#52AE43",
               "spearmanCLR"     = "#52AE43",
               "zi_kendall"      = "#559EA3",
               "zi_kendallCLR"   = "#559EA3"
              )
colors2 = c("original"       =  "#0a0a0a",
            "recomputed.clr" =  "#B6B6B6")


# ========= #
# read data #
# ========= #

# load recomputed data
dat1 = read.delim("results/egad/auroc.txt.gz") %>% 
        filter(network == "GOslim" & grepl("CLR", coefficient)) %>% 
        mutate(dataset = gsub("-expr.*$", "", dataset),
               coefficient1=coefficient,
               coefficient=gsub("CLR", "", coefficient1),
               type="recomputed.clr") %>% 
        filter(coefficient %in% methods)

# load original data
dat2 = read.delim("data/function/auroc.txt.gz") %>% 
        mutate(dataset = gsub("-expr.*$", "", dataset), coefficient=gsub("rho_p","proportionality",coefficient)) %>% 
        filter(dataset %in% dat1$dataset & coefficient %in% methods) %>% 
        mutate(coefficient1=coefficient, network="GOslim", type="original")

# merge dataset
dat = rbind(dat1, dat2)

# group data
med0 = dat %>%
  group_by(dataset, coefficient, coefficient1, network, type) %>%
  dplyr::summarise(median_auc = median(auroc)) %>%
  ungroup() 


# plot
labels = levels(with(med0, reorder(coefficient1, -median_auc, median)))
p6 = ggplot(med0, aes(x = reorder(coefficient1, -median_auc, median), 
                    y = median_auc, fill = coefficient, color=type)) + 
        geom_boxplot(outlier.shape = NA) +
        geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
        scale_fill_manual(values = colors) +
        scale_color_manual(name='', values=colors2) +
        scale_y_continuous("AUC", limits = c(0.45, 0.65), expand = c(0, 0)) + 
        theme_sc + 
        theme(axis.title.x = element_blank(),
                axis.line.y = element_line(color = 'grey50'),
                panel.grid.major.x = element_line(color = 'grey80', 
                                                linetype = 'dotted'),
                axis.ticks.x = element_blank(),
                plot.background = element_blank(),
                legend.position="right")
p6
ggsave(paste0("results/fig_nov_2021/egad_GOslim_clr+original.pdf"), p6, width = 14, height = 9, units = 'cm')