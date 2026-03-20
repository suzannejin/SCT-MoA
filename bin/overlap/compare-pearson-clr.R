options(stringsAsFactors = F)
library(tidyverse)
library(ggpubr)

files1 = list.files("../SCT-MoA/results/overlap/main/", 
                   pattern="pearson.txt", 
                   full.names=T)
files2 = list.files("../SCT-MoA/results/overlap/main/", 
                   pattern="pearsonCLR.txt", 
                   full.names=T)  
dat1 =  map(files1, read.delim) %>%
  bind_rows()
dat2 = map(files2, read.delim) %>%
  bind_rows() %>% 
  mutate(coefficient=gsub("pearsonCLR", "pearson", coefficient))
dat = merge(dat1, dat2, by=c("dataset","coefficient","network","cutoff"), suffix=c(".pearson",".pearsonCLR"))

info = read.csv("data/one-per-publication/counts.txt")
dat = merge(dat, info, by="dataset")



# ================== #
# color by min count #
# ================== #

g1 = ggplot(dat, aes(x=obs.pearson, y=obs.pearsonCLR, color=min)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Pearson - overlap") +
        ylab("Pearson+CLR - overlap") +
        theme(legend.text=element_text(size=6))
g2 = ggplot(dat, aes(x=rnd_mean.pearson, y=rnd_mean.pearsonCLR, color=min)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Pearson - overlap rnd_mean") +
        ylab("Pearson+CLR - overlap rnd_mean") +
        theme(legend.text=element_text(size=6))
g3 = ggplot(dat, aes(x=rnd_sd.pearson, y=rnd_sd.pearsonCLR, color=min)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Pearson - overlap rnd_sd") +
        ylab("Pearson+CLR - overlap rnd_sd") +
        theme(legend.text=element_text(size=6))
g4 = ggplot(dat, aes(x=z_score.pearson, y=z_score.pearsonCLR, color=min)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Pearson - zscore") +
        ylab("Pearson+CLR - zscore") +
        theme(legend.text=element_text(size=6))
g = ggarrange(g1, g2, g3, g4, nrow=2, ncol=2, common.legend=T, legend="right")
out = paste0("results/overlap/compare_results/pearson+clr_mincount.jpg")
ggsave(out, g, height=6, width=8)

# ================== #
# color by min range #
# ================== #

dat2 = dat
dat2$iscount = "0"
dat2[grep("count", dat2$count), "iscount"] = "1"
dat2$mincount = NA 
dat2[dat2$min <= 1e-10, "mincount"] = "<= 1e-10"
dat2[dat2$min < 1 & dat2$min > 1e-10, "mincount"] = "1e-10 < m < 1"
dat2[dat2$min >= 1, "mincount"] = ">= 1"
dat2$mincount = factor(dat2$mincount, level=c("<= 1e-10", "1e-10 < m < 1", ">= 1"))

g1 = ggplot(dat2, aes(x=obs.pearson, y=obs.pearsonCLR, color=mincount, shape=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Pearson - overlap") +
        ylab("Pearson+CLR - overlap") +
        theme(legend.text=element_text(size=6))
g2 = ggplot(dat2, aes(x=rnd_mean.pearson, y=rnd_mean.pearsonCLR, color=mincount, shape=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Pearson - overlap rnd_mean") +
        ylab("Pearson+CLR - overlap rnd_mean") +
        theme(legend.text=element_text(size=6))
g3 = ggplot(dat2, aes(x=rnd_sd.pearson, y=rnd_sd.pearsonCLR, color=mincount, shape=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Pearson - rnd_sd") +
        ylab("Pearson+CLR - rnd_sd") +
        theme(legend.text=element_text(size=6))
g4 = ggplot(dat2, aes(x=z_score.pearson, y=z_score.pearsonCLR, color=mincount, shape=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Pearson - zscore") +
        ylab("Pearson+CLR - zscore") +
        theme(legend.text=element_text(size=6))
g = ggarrange(g1, g2, g3, g4, nrow=2, ncol=2, common.legend=T, legend="right")
out = paste0("results/overlap/compare_results/pearson+clr_mincountrange.jpg")
ggsave(out, g, height=6, width=8)