[![Build Status](https://travis-ci.org/protViz/prozor.svg?branch=master)](https://travis-ci.org/protViz/prozor)
[![Project Stats](https://www.ohloh.net/p/prozor/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/prozor)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/prozor)](https://cran.r-project.org/package=prozor)
[![](http://cranlogs.r-pkg.org/badges/prozor)](https://cran.r-project.org/package=prozor)
[![](http://cranlogs.r-pkg.org/badges/grand-total/prozor)](https://cran.r-project.org/package=prozor)


# prozor
Determine minimal Protein set given list of peptide protein mappings. Various weights can be assigned to peptides (i.e. inverse peptide frequencies).
Generate reverse decoy sequences

 
## How to install:
For CRAN version (not the newest). Please use the github version.

```r
install.packages("prozor")
```

This is how you install the github version.

```r
install.packages("devtools")
library(devtools)
install_github("protviz/prozor")
```

### for Developers

downlod git repo. Use roxygenize2 to document new functions. Than run these 2 commands to update namespace and Rd files:

```r
library(devtools)
document()
```


Example fo creating a fasta file with the `fgcz_create_fasta.R` script

```bash
ls ./fasta_db/fgcz_3071_Chlorella
more ./fasta_db/fgcz_3071_Chlorella/annotation.txt
more ./fasta_db/fgcz_3071_Chlorella/uniprot-taxonomy_3071.fasta
clear
/home/wolski/R/x86_64-pc-linux-gnu-library/3.5/prozor/script/fgcz_create_fasta.R -h
/home/wolski/R/x86_64-pc-linux-gnu-library/3.5/prozor/script/fgcz_create_fasta.R ./fasta_db/fgcz_3071_Chlorella -o /srv/www/htdocs/FASTA/

cat fgcz_3071_Chlorella_d.txt
cat fgcz_3071_Chlorella_d.txt | bfabric_save_fasta.py 3071  /srv/www/htdocs/FASTA/fgcz_3071_Chlorella_d_20200604.fasta
```
