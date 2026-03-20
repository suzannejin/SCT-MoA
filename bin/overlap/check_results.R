options(stringsAsFactors = F)
library(tidyverse)
library(ggpubr)

f1 = "results/overlap/overlap.txt.gz"
f2 = "data/networks/overlap.txt.gz"
files = list.files("../SCT-MoA_old/results_new/overlap/main/", 
                   pattern="pearson.txt|spearman.txt|rhoORI.txt|zi_kendall.txt", 
                   full.names=T)
files2 = list.files("../SCT-MoA_old/results_new/overlap/main/", 
                   pattern="pearson.txt|spearman.txt|rhoTPQ.txt|zi_kendall.txt", 
                   full.names=T)
files3 = list.files("../SCT-MoA_old/results_new/overlap/95",
                    pattern="pearson.txt|spearman.txt|rhoORI.txt|zi_kendall.txt", 
                    full.names=T
                   )

dat1 = read.delim(f1)
dat2 = read.delim(f2) %>%
    mutate(dataset=gsub("-expr.*$", "", dataset), coefficient=gsub("^rho_p$","rho",coefficient)) %>% 
    filter(coefficient %in% c("pearson","spearman","rho","zi_kendall") & dataset %in% dat1$dataset)
dat3 = map(files, read.delim) %>%
  bind_rows() %>% 
  mutate(coefficient=gsub("rhoORI","rho",coefficient))
dat4 = map(files2, read.delim) %>%
  bind_rows() %>% 
  mutate(coefficient=gsub("rhoTPQ","rho",coefficient))
dat5 = map(files3, read.delim) %>% 
  bind_rows() %>% 
  mutate(coefficient=gsub("rhoORI","rho",coefficient))

cols = c("dataset","coefficient","network","cutoff")
dats = list(dat1,dat2,dat4,dat4,dat5)
suffix = c(".recomputed",".original",".recomputed.old",".recomputed.tpq",".recomputed.95")
for (i in 1:length(dats)){
  a = names(dats[[i]])
  pos = which(a %in% cols)
  names(dats[[i]])[-pos] = paste0(names(dats[[i]])[-pos], suffix[i])
}
dat = Reduce(function(x,y) merge(x,y,all=T,by=cols), dats)
dat$cutoff = as.character(dat$cutoff)

# add info
info = read.csv("data/one-per-publication/counts.txt")
dat = merge(dat, info, by="dataset")


# ====================== #
# ORIGINAL VS RECOMPUTED #
# ====================== #

# plot by network and coefficient
for (network in unique(dat$network)){
    gdf = dat[dat$network==network,] # & dat$cutoff==1e5,]
    g = ggplot(gdf, aes(x=obs.recomputed, y=obs.original, color=coefficient)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed overlap") +
            ylab("Original overlap") +
            ggtitle(network)  
            # geom_text(aes(label=ifelse(obs.x!=obs.y, dataset, '')),hjust = -0.1,size=1.5)  #, check_overlap = TRUE, hjust = -0.1
            # theme(legend.position = "none") 
    out = paste0("results/overlap/compare_results/compare_recomputed_vs_original_obs/", network, ".jpg")
    ggsave(out, g, height=4, width=6)
}

# obs
g = ggplot(dat, aes(x=obs.recomputed, y=obs.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed overlap") +
            ylab("Original overlap") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_original_obs/obs.jpg")
ggsave(out, g, height=4, width=6)

# a = dat[dat$coefficient=="rho" & dat$obs.x!=dat$obs.y,]
# b = dat[dat$coefficient=="rho" & dat$obs.x==dat$obs.y,]

# rnd_mean
g = ggplot(dat, aes(x=rnd_mean.recomputed, y=rnd_mean.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed mean overlap (vs random ref)") +
            ylab("Original mean overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_original_obs/rnd.jpg")
ggsave(out, g, height=4, width=6)

# rnd_sd
g = ggplot(dat, aes(x=rnd_sd.recomputed, y=rnd_sd.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed sd overlap (vs random ref)") +
            ylab("Original sd overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_original_obs/rnd_sd.jpg")
ggsave(out, g, height=4, width=6)

# z-score
g = ggplot(dat, aes(x=z_score.recomputed, y=z_score.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed z-score") +
            ylab("Original z-score") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_original_obs/zscore.jpg")
ggsave(out, g, height=4, width=6)


# ============================ #
# RECOMPUTED OLD VS RECOMPUTED #
# ============================ #

# plot by network and coefficient
for (network in unique(dat$network)){
    gdf = dat[dat$network==network,] # & dat$cutoff==1e5,]
    g = ggplot(gdf, aes(x=obs.recomputed, y=obs.recomputed.old, color=coefficient)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed overlap") +
            ylab("Recomputed-old overlap") +
            ggtitle(network)  
            # geom_text(aes(label=ifelse(obs.x!=obs.y, dataset, '')),hjust = -0.1,size=1.5)  #, check_overlap = TRUE, hjust = -0.1
            # theme(legend.position = "none") 
    out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedold_obs/", network, ".jpg")
    ggsave(out, g, height=4, width=6)
}

# obs
g = ggplot(dat, aes(x=obs.recomputed, y=obs.recomputed.old, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed overlap") +
            ylab("Recomputed-old overlap") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedold_obs/obs.jpg")
ggsave(out, g, height=4, width=6)

# rnd_mean
g = ggplot(dat, aes(x=rnd_mean.recomputed, y=rnd_mean.recomputed.old, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed mean overlap (vs random ref)") +
            ylab("Recomputed-old mean overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedold_obs/rnd.jpg")
ggsave(out, g, height=4, width=6)

# rnd_sd
g = ggplot(dat, aes(x=rnd_sd.recomputed, y=rnd_sd.recomputed.old, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed sd overlap (vs random ref)") +
            ylab("Recomputed-old sd overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedold_obs/rnd_sd.jpg")
ggsave(out, g, height=4, width=6)

# z-score
g = ggplot(dat, aes(x=z_score.recomputed, y=z_score.recomputed.old, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed z-score") +
            ylab("Recomputed-old z-score") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedold_obs/zscore.jpg")
ggsave(out, g, height=4, width=6)


# ============================ #
# RECOMPUTED TPQ VS RECOMPUTED #
# ============================ #

# plot by network and coefficient
for (network in unique(dat$network)){
    gdf = dat[dat$network==network,] # & dat$cutoff==1e5,]
    g = ggplot(gdf, aes(x=obs.recomputed, y=obs.recomputed.tpq, color=coefficient)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed overlap") +
            ylab("Recomputed-tpq overlap") +
            ggtitle(network)  
            # geom_text(aes(label=ifelse(obs.x!=obs.y, dataset, '')),hjust = -0.1,size=1.5)  #, check_overlap = TRUE, hjust = -0.1
            # theme(legend.position = "none") 
    out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedTPQ_obs/", network, ".jpg")
    ggsave(out, g, height=4, width=6)
}

# obs
g = ggplot(dat, aes(x=obs.recomputed, y=obs.recomputed.tpq, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed overlap") +
            ylab("Recomputed-tpq overlap") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedTPQ_obs/obs.jpg")
ggsave(out, g, height=4, width=6)

# rnd_mean
g = ggplot(dat, aes(x=rnd_mean.recomputed, y=rnd_mean.recomputed.tpq, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed mean overlap (vs random ref)") +
            ylab("Recomputed-tpq mean overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedTPQ_obs/rnd.jpg")
ggsave(out, g, height=4, width=6)

# rnd_sd
g = ggplot(dat, aes(x=rnd_sd.recomputed, y=rnd_sd.recomputed.tpq, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed sd overlap (vs random ref)") +
            ylab("Recomputed-tpq sd overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedTPQ_obs/rnd_sd.jpg")
ggsave(out, g, height=4, width=6)

# z-score
g = ggplot(dat, aes(x=z_score.recomputed, y=z_score.recomputed.tpq, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed z-score") +
            ylab("Recomputed-tpq z-score") 
out = paste0("results/overlap/compare_results/compare_recomputed_vs_recomputedTPQ_obs/zscore.jpg")
ggsave(out, g, height=4, width=6)



# ========================== #
# ORIGINAL VS RECOMPUTED TPQ #
# ========================== #

# plot by network and coefficient
for (network in unique(dat$network)){
    gdf = dat[dat$network==network,] # & dat$cutoff==1e5,]
    g = ggplot(gdf, aes(x=obs.recomputed.tpq, y=obs.original, color=coefficient)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed-tpq overlap") +
            ylab("Original overlap") +
            ggtitle(network)  
            # geom_text(aes(label=ifelse(obs.x!=obs.y, dataset, '')),hjust = -0.1,size=1.5)  #, check_overlap = TRUE, hjust = -0.1
            # theme(legend.position = "none") 
    out = paste0("results/overlap/compare_results/compare_recomputedTPQ_vs_original_obs/", network, ".jpg")
    ggsave(out, g, height=4, width=6)
}

# obs
g = ggplot(dat, aes(x=obs.recomputed.tpq, y=obs.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed-tpq overlap") +
            ylab("Original overlap") 
out = paste0("results/overlap/compare_results/compare_recomputedTPQ_vs_original_obs/obs.jpg")
ggsave(out, g, height=4, width=6)

# rnd_mean
g = ggplot(dat, aes(x=rnd_mean.recomputed.tpq, y=rnd_mean.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed-tpq mean overlap (vs random ref)") +
            ylab("Original mean overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputedTPQ_vs_original_obs/rnd.jpg")
ggsave(out, g, height=4, width=6)

# rnd_sd
g = ggplot(dat, aes(x=rnd_sd.recomputed.tpq, y=rnd_sd.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed-tpq sd overlap (vs random ref)") +
            ylab("Original sd overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputedTPQ_vs_original_obs/rnd_sd.jpg")
ggsave(out, g, height=4, width=6)

# z-score
g = ggplot(dat, aes(x=z_score.recomputed.tpq, y=z_score.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed-tpq z-score") +
            ylab("Original z-score") 
out = paste0("results/overlap/compare_results/compare_recomputedTPQ_vs_original_obs/zscore.jpg")
ggsave(out, g, height=4, width=6)



# ========================= #
# ORIGINAL VS RECOMPUTED 95 #
# ========================= #

# plot by network and coefficient
for (network in unique(dat$network)){
    gdf = dat[dat$network==network,] # & dat$cutoff==1e5,]
    g = ggplot(gdf, aes(x=obs.recomputed.95, y=obs.original, color=coefficient)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed.95 overlap") +
            ylab("Original overlap") +
            ggtitle(network)  
            # geom_text(aes(label=ifelse(obs.x!=obs.y, dataset, '')),hjust = -0.1,size=1.5)  #, check_overlap = TRUE, hjust = -0.1
            # theme(legend.position = "none") 
    out = paste0("results/overlap/compare_results/compare_recomputed95_vs_original_obs/", network, ".jpg")
    ggsave(out, g, height=4, width=6)
}

# obs
g = ggplot(dat, aes(x=obs.recomputed.95, y=obs.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed.95 overlap") +
            ylab("Original overlap") 
out = paste0("results/overlap/compare_results/compare_recomputed95_vs_original_obs/obs.jpg")
ggsave(out, g, height=4, width=6)

# rnd_mean
g = ggplot(dat, aes(x=rnd_mean.recomputed.95, y=rnd_mean.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed.95 mean overlap (vs random ref)") +
            ylab("Original mean overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputed95_vs_original_obs/rnd.jpg")
ggsave(out, g, height=4, width=6)

# rnd_sd
g = ggplot(dat, aes(x=rnd_sd.recomputed.95, y=rnd_sd.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed.95 sd overlap (vs random ref)") +
            ylab("Original sd overlap (vs random ref)") 
out = paste0("results/overlap/compare_results/compare_recomputed95_vs_original_obs/rnd_sd.jpg")
ggsave(out, g, height=4, width=6)

# z-score
g = ggplot(dat, aes(x=z_score.recomputed.95, y=z_score.original, color=network, shape=cutoff)) + 
            facet_wrap(~coefficient, nrow=2) +
            geom_point(alpha=.5, size=1) + 
            geom_abline(linetype = "dashed") + 
            xlab("Recomputed.95 z-score") +
            ylab("Original z-score") 
out = paste0("results/overlap/compare_results/compare_recomputed95_vs_original_obs/zscore.jpg")
ggsave(out, g, height=4, width=6)

# ============================================== #
# COMPARE RESULTS - RHO - ORIGINAL VS RECOMPTUED #
# ============================================== #

dat2 = dat[dat$coefficient=="pearson",]
g1 = ggplot(dat2, aes(x=obs.recomputed, y=obs.original, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed observed overlap") +
        ylab("Original observed overlap") +
        theme(axis.title=element_text(size=8.5))
g2 = ggplot(dat2, aes(x=rnd_mean.recomputed, y=rnd_mean.original, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed random.mean overlap") +
        ylab("Original random.mean overlap") +
        theme(axis.title=element_text(size=8.5))
g3 = ggplot(dat2, aes(x=rnd_sd.recomputed, y=rnd_sd.original, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed random.sd overlap") +
        ylab("Original random.sd overlap") +
        theme(axis.title=element_text(size=8.5))
g4 = ggplot(dat2, aes(x=z_score.recomputed, y=z_score.original, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed z-score") +
        ylab("Original z-score") +
        theme(axis.title=element_text(size=8.5))
g = ggarrange(g1, g2, g3, g4, nrow=2, ncol=2, common.legend=T, legend="right")
out = paste0("results/overlap/compare_results/compare_recomputed_vs_original_obs/pearson.jpg")
ggsave(out, g, height=4, width=6)


# ===================== #
# COMPARE RESULTS - RHO # 
# ===================== #

dat2 = dat[dat$coefficient=="rho",]
# original vs recomputed
g1 = ggplot(dat2, aes(x=obs.recomputed, y=obs.original, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed overlap") +
        ylab("Original overlap") 
# original vs recomputed.tpq
g2 = ggplot(dat2, aes(x=obs.recomputed.tpq, y=obs.original, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed-tpq overlap") +
        ylab("Original overlap") 
# recomputed vs recomputed.tpq
g3 = ggplot(dat2, aes(x=obs.recomputed.tpq, y=obs.recomputed, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed-tpq overlap") +
        ylab("Recomputed overlap") 
# original vs recomputed
g4 = ggplot(dat2, aes(x=z_score.recomputed, y=z_score.original, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed z-score") +
        ylab("Original z-score") 
# original vs recomputed.tpq
g5 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.original, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed-tpq z-score") +
        ylab("Original z-score") 
# recomputed vs recomputed.tpq
g6 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.recomputed, color=network, shape=cutoff)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed-tpq z-score") +
        ylab("Recomputed z-score") 
g = ggarrange(g4, g5, g6, nrow=1, ncol=3, common.legend=T)
out = paste0("results/overlap/compare_results/rho.jpg")
ggsave(out, g, height=3, width=8)


# color by data info
dat2 = dat[dat$coefficient=="rho",]
dat2$iscount = "0"
dat2[grep("count", dat2$Data_processing), "iscount"] = "1"
# original vs recomputed
g4 = ggplot(dat2, aes(x=z_score.recomputed, y=z_score.original, color=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed z-score") +
        ylab("Original z-score") 
# original vs recomputed.tpq
g5 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.original, color=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed-tpq z-score") +
        ylab("Original z-score") 
# recomputed vs recomputed.tpq
g6 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.recomputed, color=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed-tpq z-score") +
        ylab("Recomputed z-score") 
g = ggarrange(g4, g5, g6, nrow=1, ncol=3, common.legend=T)
out = paste0("results/overlap/compare_results/rho_iscount.jpg")
ggsave(out, g, height=3, width=8)


# color by min count
dat2 = dat[dat$coefficient=="rho",]
# original vs recomputed
g4 = ggplot(dat2, aes(x=z_score.recomputed, y=z_score.original, color=min)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Recomputed z-score") +
        ylab("Original z-score") +
        theme(legend.text=element_text(size=6))
# original vs recomputed.tpq
g5 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.original, color=min)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Recomputed-tpq z-score") +
        ylab("Original z-score") +
        theme(legend.text=element_text(size=6))
# recomputed vs recomputed.tpq
g6 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.recomputed, color=min)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Recomputed-tpq z-score") +
        ylab("Recomputed z-score") +
        theme(legend.text=element_text(size=6))
g = ggarrange(g4, g5, g6, nrow=1, ncol=3, common.legend=T, legend="right")
out = paste0("results/overlap/compare_results/rho_mincount.jpg")
ggsave(out, g, height=3, width=10)


# by min count range
dat2 = dat[dat$coefficient=="rho",]
dat2$iscount = "0"
dat2[grep("count", dat2$count), "iscount"] = "1"
dat2$mincount = NA 
dat2[dat2$min <= 1e-10, "mincount"] = "<= 1e-10"
dat2[dat2$min < 1 & dat2$min > 1e-10, "mincount"] = "1e-10 < m < 1"
dat2[dat2$min >= 1, "mincount"] = ">= 1"
dat2$mincount = factor(dat2$mincount, level=c("<= 1e-10", "1e-10 < m < 1", ">= 1"))
g4 = ggplot(dat2, aes(x=z_score.recomputed, y=z_score.original, color=mincount, shape=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed z-score") +
        ylab("Original z-score") +
        theme(legend.text=element_text(size=10))
# original vs recomputed.tpq
g5 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.original, color=mincount, shape=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed-tpq z-score") +
        ylab("Original z-score") +
        theme(legend.text=element_text(size=10))
# recomputed vs recomputed.tpq
g6 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.recomputed, color=mincount, shape=iscount)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        xlab("Recomputed-tpq z-score") +
        ylab("Recomputed z-score") +
        theme(legend.text=element_text(size=10))
g = ggarrange(g4, g5, g6, nrow=1, ncol=3, common.legend=T, legend="right")
out = paste0("results/overlap/compare_results/rho_mincountrange.jpg")
ggsave(out, g, height=3, width=10)


# color by max count
dat2 = dat[dat$coefficient=="rho",]
# original vs recomputed
g4 = ggplot(dat2, aes(x=z_score.recomputed, y=z_score.original, color=max)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Recomputed z-score") +
        ylab("Original z-score") +
        theme(legend.text=element_text(size=6))
# original vs recomputed.tpq
g5 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.original, color=max)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Recomputed-tpq z-score") +
        ylab("Original z-score") +
        theme(legend.text=element_text(size=6))
# recomputed vs recomputed.tpq
g6 = ggplot(dat2, aes(x=z_score.recomputed.tpq, y=z_score.recomputed, color=max)) +
        geom_point(alpha=.5, size=1) + 
        geom_abline(linetype = "dashed") + 
        scale_color_gradientn(colours = rainbow(5), trans = "log10") +
        xlab("Recomputed-tpq z-score") +
        ylab("Recomputed z-score") +
        theme(legend.text=element_text(size=6))
g = ggarrange(g4, g5, g6, nrow=1, ncol=3, common.legend=T, legend="right")
out = paste0("results/overlap/compare_results/rho_maxcount.jpg")
ggsave(out, g, height=3, width=10)

