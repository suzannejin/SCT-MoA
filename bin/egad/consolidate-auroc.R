# Filter the complete set of AUROCs to (i) GO terms within the GO slim, and
# (ii) annotated to between 10 and 1,000 proteins in the network. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(ontologyIndex)
library(tidyverse)
library(magrittr)
source("bin/functions.R")

# read GO slim
slim = get_ontology("data/go/goslim_generic.obo.gz")

# process files
dir = "results/egad/main"
message("reading files in directory ", dir, " ...")
files = list.files(dir, pattern = "*.gz", full.names = T)

for (db in c("GO","GOslim","Reactome")){
    all = list()
    if (db %in% c("GO","GOslim")) {
      files2 = files[grepl("GO", files)]
    }else{
      files2 = files[grepl(db, files)]
    }
    pb = progress::progress_bar$new(
      format = "file :what [:bar] :percent eta: :eta",
      clear = F, total = length(files2), width = 100)
    for (i in seq_len(length(files2))) {
      file = files2[i]
      dat = read.delim(file)
      if (db == "GOslim"){
        filtered = dplyr::filter(dat, term %in% slim$id & 
                                  dplyr::between(n_proteins, 10, 1000)) %>%
          dplyr::select(dataset, coefficient, term, auroc)
      }else{
        filtered = dplyr::filter(dat, dplyr::between(n_proteins, 10, 1000)) %>%
          dplyr::select(dataset, coefficient, term, auroc)
      }
      filtered["network"] = db
      all[[i]] = filtered
      # tick progress bar
      pb$tick(tokens = list(what = sprintf(
        paste0("%-", nchar(length(files2)), "s"), i)))
    }
    combined = bind_rows(all)

    # write 
    write_and_gzip(combined, paste0("results/egad/auroc_", db, ".txt"))
}


dats = list()
for (db in c("GO","GOslim","Reactome")){
  file = paste0("results/egad/auroc_", db, ".txt.gz")
  dats[[db]] = read.delim(file)
}
combined = bind_rows(dats)
write_and_gzip(combined, "results/egad/auroc.txt")