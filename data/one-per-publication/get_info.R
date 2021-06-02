#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)
datadir = args[1]   # ./data/one-dataset-per-publication


# list datasets
files = list.files(datadir, pattern='.txt.gz', full.names=T)

for (perc in c(NA, .8)){

    # create empty data frame
    csv = data.frame(matrix(NA, ncol=4, nrow=length(files)))
    colnames(csv) = c('dataset', 'cells', 'genes', 'dropouts')

    for (i in 1:length(files)){

        # read expression data 
        f = files[i]
        dataset = gsub('.txt.gz', '', basename(f))
        expr = read.delim(f)
        print(dataset)

        # filter: remove genes that are zero in > perc of samples
        if (is.na(perc)){
            out = paste0(datadir, "/info1.csv")
        }else{
            cols = which(colMeans(expr == 0) < perc)
            expr = expr[, cols]
            out = paste0(datadir, "/info2.csv")
        }

        # fill data frame with cells, genes, and dropout
        size = dim(expr)
        dropout = round(mean(expr == 0, na.rm = T), 6)
        csv[i,] = c(dataset, size[1], size[2], dropout)
    }

    write.csv(csv, out, row.names=F, quote=F)
}