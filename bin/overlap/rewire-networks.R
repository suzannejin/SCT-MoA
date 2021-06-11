# Generate 1,000 randomized versions of human and mouse PPI, signalling,
# metabolic, and co-occurrence networks. 
options(stringsAsFactors = F)
library(igraph)
library(tidyverse)
library(magrittr)
library(this.path)

# basedir
scriptdir = this.dir()
mypath <- function(path){
    scriptdir = this.dir()
    basedir = file.path(scriptdir, "../..")
    file.path(basedir, path)
}

# load functions
source(mypath("bin/functions.R"))

# usage: rewire-networks.R <input> <idx> <network> <species> <output_dir>
args = commandArgs(trailingOnly = T)
paste(args, collapse = " ")
if (length(args) < 5)
  stop("must provide the network file, idx file, network name, species name, and output directory")
file = args[1]
idxfile = args[2]
network = args[3]
species = args[4]
output_dir = args[5]
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)

# read idx file, with columns: file, network, species
inputs = read.table(idxfile)
idx = which(inputs$network==network & inputs$species==species)
message("processing database ", network, " - species ", species, " - idx ", idx)
set.seed(idx)

# read network
net = drop_na(read.delim(file))
g = graph_from_data_frame(net, directed = F)
m = length(E(g))

# rewire network
bootstraps = 1e3
for (b in seq_len(bootstraps)) {
  if (b %% 10 == 0) 
    message("bootstrap ", b, " of ", bootstraps, " ...")
  rewired = rewire(g, with = keeping_degseq(
    loops = F, niter = 10 * m)) %>%
    igraph::as_data_frame()
  output_file = paste0(network, "-", species, "-", b, ".txt")
  output_path = file.path(output_dir, output_file)
  write_and_gzip(rewired, output_path)
}
