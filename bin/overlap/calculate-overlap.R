# Calculate the enrichment of overlap between single-cell gene co-expression
# networks and PPI, TF/target, or metabolic networks relative to rewired 
# versions of the same. 
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(igraph)
library(this.path)

# basedir
scriptdir = this.dir()
mypath <- function(path){
    scriptdir = this.dir()
    basedir = file.path(scriptdir, "../..")
    file.path(basedir, path)
}

source(mypath("bin/functions.R"))

# usage: calculate-overlap.R <input> <output name> <rewire_dir>
args = commandArgs(trailingOnly = T)
if (length(args) < 3)
  stop("must provide input coexpression file, output filename, and rewired networks directory")
file = args[1]
output = args[2]
rewire_dir = args[3]

# ================== #
# INPUT COEXPRESSION #
# ================== #

# check whether output file already exists
message(".. processing file ", file, " ...")
filename = basename(file)
coef = get_coef(filename)
coef2 = filename2coef(filename, ".Rdata")
idx = unlist(gregexpr(coef, filename))
dataset = substr(filename, 0, idx - 2)

# load coexpression matrix
load(file, verbose = T) ## coexpr
# replace missing values with median (ZI kendall)
coexpr[is.na(coexpr)] = median(coexpr, na.rm = T)
# replace infinite values with the minimum (binomial)
coexpr[is.infinite(coexpr)] = min(coexpr, na.rm = T)

# detect species
species = get_species(colnames(coexpr))

# create ranked coexpression data frame
tri = upper.tri(coexpr)
idxs = which(tri, arr.ind = T)
int = data.frame(id1 = rownames(coexpr)[idxs[,1]], 
                 id2 = rownames(coexpr)[idxs[,2]], 
                 correlation = coexpr[tri])
int %<>% arrange(-correlation)

# ============== #
# REFERENCE DATA #
# ============== #

# read original networks
net_names = c("HIPPIE", "Reactome", "STRING", "OmniPath")
original_files = mypath(paste0("data/networks/", net_names, "/", species, ".txt.gz"))
nets = map(original_files, ~ graph_from_data_frame(drop_na(read.delim(.)))) %>%
  setNames(net_names)

# read rewired networks one at a time
bootstraps = 1e2
results = data.frame()
for (net_name in net_names) {
  message("analyzing network ", net_name, " ...")
  net = nets[[net_name]]
  
  # calculate observed overlap 
  message("  calculating observed overlap with ", net_name, " ...")
  observed = numeric(0)
  cutoffs = c(2e4, 5e4, 1e5)
  for (cutoff in cutoffs) {
    subset = graph_from_data_frame(int[seq_len(cutoff), ])
    obs = length(E(intersection(subset, net, keep.all.vertices = F)))
    observed %<>% c(obs)
  }

  # calculate overlap in rewired networks     
  random = list()
  for (bootstrap in seq_len(bootstraps)) {
    if (bootstrap %% 10 == 0)
      message("  analyzing bootstrap ", bootstrap, " of ", bootstraps, " ...")
    
    # read rewired network
    rewire_file = file.path(rewire_dir, paste0(
      net_name, "-", species, "-", bootstrap, ".txt.gz"))
    rewire = graph_from_data_frame(drop_na(read.delim(rewire_file)))
    
    # calculate rewired network overlap at different cutoffs
    rnd = numeric(0)
    for (cutoff in cutoffs) {
      subset = graph_from_data_frame(int[seq_len(cutoff), ])
      rnd %<>% c(length(
        E(intersection(subset, rewire, keep.all.vertices = F))))
    }
    # record
    random[[bootstrap]] = rnd
  }
  
  # calculate outcomes
  rnd_mean = map_dbl(seq_len(length(cutoffs)), ~ 
                       mean(map_dbl(random, .), na.rm = T))
  rnd_median = map_dbl(seq_len(length(cutoffs)), ~
                         median(map_dbl(random, .), na.rm = T))
  rnd_sd = map_dbl(seq_len(length(cutoffs)), ~
                     sd(map_dbl(random, .), na.rm = T))
  enrs = observed / rnd_median
  z_scores = (observed - rnd_mean) / rnd_sd
  
  # record all
  results %<>% rbind(data.frame(dataset = dataset, coefficient = coef2, 
                                network = net_name, cutoff = cutoffs,
                                obs = observed, rnd_mean = rnd_mean, 
                                rnd_median = rnd_median, rnd_sd = rnd_sd, 
                                enrichment = enrs, z_score = z_scores))
}

# write output 
write.table(results, output, quote = F, sep = "\t", row.names = F)
