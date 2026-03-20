# Construct gene correlation matrix for a given dataset and measure
library(methods)
library(dismay)

# usage: write-matrix.R <input> <output> <method>
args = commandArgs(trailingOnly = T)
if (length(args) != 3)
    stop("usage: write-matrix.R <input> <output> <method>", call. = F)
file = args[1]   # .gz
output = args[2] # .Rdata
coef = args[3]

# read input data
filename = gsub("\\.txt.*$", "", basename(file))
message("creating coexpression networks for file ", filename)
dat = as.matrix(read.delim(file))

# create correlation matrix 
message("  calculating ", coef, " matrix ...")
if (coef %in% dismay::metrics()){
    coexpr = dismay::dismay(dat, metric = coef)
}else if( grepl("CLR", coef) ){
    library(propr)
    coef = stringr::str_replace(coef, "CLR", "")
    packageVersion("propr")
    print(coef)
    pro = propr(dat, metric=coef, p=0) 
    coexpr = pro@matrix
}else if( grepl("SJIN", coef) ){
    library(propr, lib.loc=paste0(.libPaths()[1], "/propr_sjin"))
    coef = stringr::str_replace(coef, "SJIN", "")
    packageVersion("propr")
    print(coef)
    pro = propr(dat, metric=coef, p=0) 
    coexpr = pro@matrix
}
save(coexpr, file = output)
