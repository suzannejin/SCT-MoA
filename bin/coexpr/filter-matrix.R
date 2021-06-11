# Filter coexpression matrices based on the proportion of cells in the dataset
# detectably expressing each gene. 
library(methods)
library(dismay)
library(tidyverse)
library(magrittr)

# usage: filter-matrices.R <expr_matrix> <coexpr_matrix> <output_dir> <output filename>
args = commandArgs(trailingOnly = T)
if (length(args) < 3)
  stop("must provide expression matrix, coexpression matrix, output directory, output filename")
expr = args[1]
coexpr = args[2]
output_dir = args[3]
filename = args[4]

# define thresholds
thresholds = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
# make sure output directories exist
threshold_dirs = file.path(output_dir, thresholds * 100)
for (threshold_dir in threshold_dirs)
  if (!dir.exists(threshold_dir))
    dir.create(threshold_dir, recursive = T)

# expression matrix
message("filtering coexpression networks for file ", basename(expr))
dat = as.matrix(read.delim(expr))
 
# define gene lists
gene_lists = map(thresholds, ~ names(which(colMeans(dat == 0) < .))) %>%
setNames(thresholds)
  
# correlation matrix
message("  filtering ", basename(coexpr), " matrix ...")
load(coexpr, verbose = T) ## coexpr
coexpr_in = coexpr

# filter based on proportion of dropouts
for (threshold in thresholds) {
    genes = gene_lists[[as.character(threshold)]]
    coexpr = coexpr_in[genes, genes]
    
    # write
    output_file = file.path(output_dir, threshold * 100, filename)
    save(coexpr, file = output_file)
}

