# Consolidate all network overlap data into a single file and gzip it. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)

# list files
files = list.files("../SCT-MoA_old/results_new/overlap/main", 
                   pattern = "pearson.txt|spearman.txt|rhoTPQ.txt|zi_kendall.txt", 
                   full.names = T)
# read data
dat = map(files, read.delim) %>%
  bind_rows()
# write
output = "results/overlap/overlap_tpq.txt"
write.table(dat, output, quote = F, sep = "\t", row.names = F)
system(paste("gzip --force", output))
