library(Biostrings)
library(prozor)


library(seqinr)
file = file.path(path.package("prozor"),"extdata/shortfasta.fasta" )
fasta = read.fasta(file = file, as.string = TRUE, seqtype="AA")

class(fasta)
names(fasta)
as.character(fasta[[1]])
attributes(fasta[[1]])
x<-"blub"
class(x)

attributes(x)$newattrib <- "special"
attributes(x)$class <- "numeric"
class(x)

data(pepdata)
head(pepdata)


biofasta <- lapply(fasta, function(x){Biostrings::AAString(as.character(x))})
peptides <- lapply(pepdata[,1],function(x){Biostrings::AAString(as.character(x))})

biofasta <- AAStringSet(biofasta)

length(biofasta)
length(peptides)

Biostrings::vmatchPattern(Biostrings::AAString(pepdata[1,1]),biofasta)


res <-  vector(lenght())

length(peptides)
tmp <-Biostrings::matchPDict(AAStringSet(peptides),biofasta[[1]])

xx <- (lapply(tmp, function(x){x@start}))

unlist(xx[1:10])
