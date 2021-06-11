# General-purpose functions

map_genes = function(dat, keys, from = "SYMBOL", to = "ENSEMBL", db = NULL) {
  if (is.null(db) | !"OrgDb" %in% class(db))
    stop("argument `db` must be an object of class OrgDb")
  if (is.null(from) | !from %in% keytypes(db))
    stop("invalid `from` argument: ", from)
  if (is.null(to) | !to %in% keytypes(db))
    stop("invalid `to` argument: ", to)
  map = suppressMessages(AnnotationDbi::select(db, keys = keys, columns = to, 
                                               keytype = from))
  values = map[[to]][match(keys, map[[from]])]
  aggregate(dat, by = list(values), FUN = sum)
}

write_and_gzip = function(data, filepath) {
  write.table(data, filepath, sep = "\t", quote = F, row.names = F)
  system(paste("gzip --force", filepath))
}

clean_metric = function(vec) {
  as.character(fct_recode(vec,
                          "Manhattan distance" = "manhattan",
                          "Euclidean distance" = "euclidean",
                          "Canberra distance" = "canberra",
                          "Zero-inflated Kendall correlation" = "zi_kendall",
                          "Kendall correlation" = "kendall",
                          "Spearman correlation" = "spearman",
                          "Pearson correlation" = "pearson",
                          "Jaccard index" = "jaccard",
                          "Dice coefficient" = "dice",
                          "Hamming distance" = "hamming",
                          "Co-dependency index" = "binomial",
                          "Biweight midcorrelation" = "bicor",
                          "Cosine distance" = "cosine",
                          "Weighted rank correlation" = "weighted_rank",
                          "ϕs" = "phi",
                          "ρp" = "rho"
  ))
}

get_children = function(term) {
  library(MeSH.PCR.db)
  # get term and children
  terms = select(MeSH.PCR.db, keys = term, columns = c("PARENT", "CHILD"), 
                 keytype = "PARENT")
  # get CUI-MeSH map
  map = read.delim(mypath("data/disease/phenopedia/CUI-MeSH-map.txt.gz"))
  # subset
  sub = map %>% filter(mesh %in% terms$CHILD)
  # return children
  return(unique(sub$CUI))
}

get_species = function(vector) {
  coding = read.delim(mypath("data/ensembl/protein_coding_genes.txt.gz"))
  subset = coding %>% dplyr::filter(gene %in% vector)
  names(which.max(table(subset$source)))
}

get_coef = function(filename){
  ## get real coef implemented in dismay. eg pearson, rho, spearman
  coefs = dismay::metrics()
  coef = dplyr::last(coefs[sapply(coefs, function(c) grepl(c, filename))])
  if (is.na(coef) | length(coef) == 0) {
    stop(".. could not match file ", filename, " to coexpression method")
  }
  return(coef)
}

filename2coef = function(filename, extension){
  ## get full coef based on filename. eg pearson, pearsonCLR, rho, spearman, spearmanCLR, etc.
  require(stringr)
  name = str_replace(filename, extension, "")
  split = str_split(name, "-")[[1]]
  coef = split[length(split)]
  return(coef)
}


# ============= #
# GET REFERENCE # 
# ============= #

get_annotations = function(db, species, genes){
  if (db == "GO"){
    go = get_go(species)
    annotations = get_annotations_go(go, genes)
  }else if(db == "Reactome"){
    react = get_reactome(species)
    annotations = get_annotations_reactome(react, genes)
  }
  return(annotations)
}

get_go = function(species){
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(flavin)
  require(ontologyIndex)

  ontology = get_ontology(mypath("data/go/go-basic.obo.gz"))
  if (species == "human"){
    data = mypath("data/go/goa_human.gpa.gz")
    database = org.Hs.eg.db
  }else if(species == "mouse"){
    data = mypath("data/go/goa_mouse.gpa.gz")
    database = org.Mm.eg.db
  }else {
    stop("Erroneous species ", species)
  }
  go = read_gpa(data, accession = "ENSEMBL", database = database, ontology = ontology, propagate = T)
  go = filter_breadth(go, min = 1, max = 1e4)
  return(go)
}

get_annotations_go = function(go, genes){
  require("EGAD")
  go_subset = go %>%
    dplyr::filter(ENSEMBL %in% genes) %>%
    dplyr::select(ENSEMBL, GO.ID) %>%
    as.matrix()
  if (nrow(go_subset) == 0)
    stop("no GO terms annotated to network genes. check right GO file")
  go_terms = unique(go_subset[, "GO.ID"])
  annotations = make_annotations(go_subset, genes, go_terms) 
  return(annotations) 
}

get_reactome = function(species){
  react = readr::read_tsv(mypath("data/networks/Reactome/Ensembl2Reactome.txt.gz"),
                 col_names = c("ensembl", "pathway", "url", "name", "ec", "sp"))
  if (species == "mouse"){
    react = react[which(react$sp == "Mus musculus"),]
  }else if(species == "human"){
    react = react[which(react$sp == "Homo sapiens"),]
  }else{
    stop("species ", species, " not available")
  }
  return(react)
}

get_annotations_reactome = function(react, genes){
  require("EGAD")
  react_subset = react %>%
    dplyr::filter(ensembl %in% genes) %>%
    dplyr::select(ensembl, pathway) %>%
    as.matrix()
  if (nrow(react_subset) == 0)
    stop("no Reactome pathways annotated to network genes. check right Reactome file")
  react_terms = unique(react_subset[, "pathway"])
  annotations = make_annotations(react_subset, genes, react_terms) 
  return(annotations) 
}
