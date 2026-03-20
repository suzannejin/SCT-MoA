setwd("~/git/SCT-MoA")

# list files
files = list.files("data/one-per-publication/", ".txt.gz", full.names=T)

# info
info = read.csv("data/one-per-publication/one-dataset-per-publication.txt")[,c("Dataset","Data_processing")]
counts = info$Data_processing
names(counts) = info$Dataset

# create matrix
cols = c("dataset","count","min","max","mean","median","sd")
df = data.frame(matrix(NA, nrow=length(files), ncol=length(cols)))
colnames(df) = cols

for(i in 1:length(files)){
    file = files[i]
    filename = basename(file)
    dataset = gsub(".txt.gz","",filename)
    count = counts[dataset]
    dat = read.delim(file)
    values = unlist(c(dat))
    df[i,] = c(dataset,
               as.character(count),
               min(values[values!=0]),
               max(values),
               mean(values),
               median(values),
               sd(values)
               )
}

out = "data/one-per-publication/counts.txt"
write.csv(df, file=out, quote=F, row.names=F)