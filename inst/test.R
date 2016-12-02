library(Biostrings)
file = file.path(path.package("prozor"),"extdata/shortfasta.fasta" )
fasta = read.fasta(file = file, as.string = TRUE, seqtype="AA")
data(pepdata)
head(pepdata)


biofasta <- lapply(fasta, function(x){Biostrings::AAString(as.character(x))})
peptides <- lapply(pepdata[,1],function(x){Biostrings::AAString(as.character(x))})

biofasta <- AAStringSet(biofasta)

length(biofasta)
length(peptides)

Biostrings::matchPattern(Biostrings::AAString(pepdata[1,1]),biofasta)

library(Biostrings)
biofasta <- lapply(fasta, function(x){Biostrings::AAString(as.character(x))})
Biostrings::matchPattern(Biostrings::AAString(pepdata[1,1]),biofasta)
tt <- AAStringSet(biofasta)
Biostrings::matchPattern(Biostrings::AAString(pepdata[1,1]),tt)
Biostrings::vmatchPattern(Biostrings::AAString(pepdata[1,1]),tt)
Biostrings::vmatchPattern(Biostrings::AAString(pepdata[2,1]),tt)
pepdata[2,1]
pepdata[1,1]
tmp <- Biostrings::vmatchPattern(Biostrings::AAString(pepdata[1,1]),tt)
tmp[[1]]
length(tmp)
