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
methods = c("rho", "pearson", "spearman", "zi_kendall", "partialcor")
colors = c("rho"       = "#f06c45",
              "pearson"   = "#b294c7",
              "spearman"  = "#f88a89",
              "zi_kendall"= "#569ea4",
              "partialcor"= "#c7b699"
              )
colors1 = list("rho.recomputed"           = "#f06c45", 
              "rho.recomputed.tpq"        = "#f06c45",
              "rho.original"              = "#f06c45", 
              "pearson.recomputed"        = "#b294c7",
              "pearsonCLR.recomputed"     = "#b294c7",
              "pearson.original"          = "#b294c7",
              "spearman.recomputed"       = "#f88a89",
              "spearmanCLR.recomputed"    = "#f88a89",
              "spearman.original"         = "#f88a89",
              "zi_kendall.recomputed"     = "#569ea4",
              "zi_kendallCLR.recomputed"  = "#569ea4",
              "zi_kendall.original"       = "#569ea4",
              "partialcor.recomputed"     = "#c7b699",
              "partialcorCLR.recomputed"  = "#c7b699"
             )
colors2 = c("original"       =  "#0a0a0a",
            "recomputed"     =  "#d5d6d5",
            "recomputed.tpq" =  "#ee3f30",
            "recomputed.clr" =  "#52af43")

message("plotting main figure")

# read data
dat = read.delim("results/egad/auroc.txt.gz") %>%
  mutate(dataset = gsub("-expr.*$", "", dataset), 
        coefficient1=paste0(coefficient, ".recomputed"), 
        type="recomputed")
dat[grep("CLR",dat$coefficient), "type"] = "recomputed.clr"
dat$coefficient=gsub("CLR","",dat$coefficient)
dat2 = read.delim("data/function/auroc.txt.gz") %>% 
  mutate(dataset = gsub("-expr.*$", "", dataset), coefficient=gsub("rho_p","rho",coefficient)) %>% 
  filter(dataset %in% dat$dataset & coefficient %in% methods) %>% 
  mutate(coefficient1=paste0(coefficient, ".original"), network="GOslim", type="original")
dat3 = read.delim("results/egad/auroc_tpq_GOslim.txt.gz") %>%
  mutate(coefficient=gsub("rhoTPQ","rho",coefficient)) %>%
  filter(coefficient == "rho") %>%
  mutate(coefficient1=paste0(coefficient, ".recomputed.tpq"), type="recomputed.tpq")
dat = rbind(dat, dat2, dat3)
med = dat %>%
  group_by(dataset, coefficient1, coefficient, network, type) %>%
  dplyr::summarise(median_auc = median(auroc)) %>%
  ungroup() %>%
  mutate(dataset = gsub("-expr.*$", "", dataset))


# # main plot
# labels = levels(with(med, reorder(coefficient1, median_auc, median)))
# p1 = ggplot(med, aes(x = reorder(coefficient1, median_auc, median), 
#                       y = median_auc, fill = coefficient, color=type)) + 
#   facet_wrap(~ network, ncol = 3) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
#   scale_y_continuous("AUC", limits = c(0.45, 0.7)) + 
#   scale_fill_manual(name = '', values = colors) + 
#   scale_color_manual(name='', values=colors2) +
#   coord_flip() + 
#   clean_theme + 
#   theme(legend.position = 'right',
#         axis.text.y = element_text(angle = 0, hjust = 1),
#         panel.grid.major.y = element_line(color = 'grey80', 
#                                           linetype = 'dotted'),
#         axis.ticks.y = element_blank(),
#         axis.line.x = element_line(color = 'grey50'),
#         axis.title.y = element_blank(),
#         strip.background = element_rect(fill = 'grey90', color = 'white'))
# p1
# ggsave("results/fig/egad_recomputed+tpq+original.pdf", p1, width = 24.5, height = 12, units = 'cm')



# GOslim
db = "GOslim"
message("plot figure for ", db)

med0 = dat %>%
    filter(network == "GOslim") %>%
    group_by(dataset, coefficient1, coefficient, type) %>%
    dplyr::summarise(median_auc = median(auroc)) %>%
    ungroup() %>%
    mutate(dataset = gsub("-expr.*$", "", dataset))

# plot
p6 = ggplot(med0, aes(x = reorder(coefficient, -median_auc, median), 
                    y = median_auc, fill = coefficient)) + 
        facet_wrap(~type, nrow=1) +
        geom_boxplot(outlier.shape = NA) +
        geom_hline(aes(yintercept = 0.5), linetype = 'dotted') +
        scale_fill_manual(values = colors) +
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
ggsave(paste0("results/fig/egad_", db, "_recomputed+tpq+original_grouped+type.pdf"), p6, width = 18, height = 9, units = 'cm')



p6 = ggplot(med0, aes(x = reorder(type, -median_auc, median), 
                    y = median_auc, fill = coefficient, color=type)) + 
        facet_wrap(~coefficient, nrow=1) +
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
ggsave(paste0("results/fig/egad_", db, "_recomputed+tpq+original_grouped+metric.pdf"), p6, width = 18, height = 9, units = 'cm')
