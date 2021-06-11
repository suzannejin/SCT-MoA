# Analyze the functional connectivity of gene coexpression networks constructed
# using a series of measures of correlation for single-cell RNA-seq datasets.
options(stringsAsFactors = F)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(flavin)
library(ontologyIndex)
library(EGAD)
library(dplyr)
library(data.table)
library(scales)
library(dismay)
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

# read arguments
# usage: calculate-auroc.R <coexpr> <output> <db name>
args = commandArgs(trailingOnly = T)
if (length(args) < 3)
  stop("Must set input coexpression file, output filename, and database name (GO or Reactome)", call. = F)
file = args[1]
output = args[2]
db = args[3]


# ================== #
# INPUT COEXPRESSION #
# ================== #

# check input file 
message(".. processing file ", file, " ...")
filename = basename(file)
coef = get_coef(filename)
coef2 = filename2coef(filename, ".Rdata")
idx = unlist(gregexpr(coef, filename))   # method
dataset = substr(filename, 0, idx - 2)   # dataset

# load coexpression matrix
load(file, verbose = T) ## coexpr

# replace missing values with median (ZI kendall)
if (any(is.na(coexpr)))
  coexpr[is.na(coexpr)] = median(coexpr, na.rm = T)
# replace infinite values with the minimum (binomial)
if (any(is.infinite(coexpr)))
  coexpr[is.infinite(coexpr)] = min(coexpr, na.rm = T)
  
# scale metrics to [-1, 1] range
if (!dismay::is_bounded(coef)) {
  coexpr = scales::rescale(coexpr, to = c(-1, 1))
}

# make EGAD coexpression network
## modified from EGAD::build_coexp_network to accommodate missing values
message("making EGAD coexpression network ...")
n = nrow(coexpr)
genes = rownames(coexpr)
net = matrix(rank(coexpr, na.last = "keep", ties.method = "average"), 
             nrow = n, ncol = n)
rownames(net) = colnames(net) = genes
net = net / max(net, na.rm = T)  
diag(net) = 1

# ============== #
# REFERENCE DATA #
# ============== #

# detect species
species = get_species(colnames(coexpr))

# get annotations for genes
annotations = get_annotations(db, species, genes) 

rm(coexpr)

# ======== #
# RUN EGAD #
# ======== #

# run EGAD
message("running guilt-by-association analysis ...")
gba = run_GBA(net, annotations, min = 0, max = 1e4)

# record AUROCs 
aurocs = gba[[1]][, "auc"]
# get number of proteins to which each term was annotated
all_terms = colSums(annotations)
n_terms = all_terms[names(aurocs)]
# write result
result = data.frame(dataset = dataset, coefficient = coef2, 
                    term = names(aurocs), auroc = aurocs, 
                    n_proteins = n_terms, 
                    pct_proteins = n_terms / length(genes))
output = stringr::str_replace(output, ".gz", "")
message("writing output file ", output)
write.table(result, output, quote = F, sep = "\t", row.names = F)
system(paste("gzip --force", output))
