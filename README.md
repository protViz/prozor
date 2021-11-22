[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/prozor)](https://cran.r-project.org/package=prozor)
[![](http://cranlogs.r-pkg.org/badges/prozor)](https://cran.r-project.org/package=prozor)
[![](http://cranlogs.r-pkg.org/badges/grand-total/prozor)](https://cran.r-project.org/package=prozor)


# prozor

- Determine minimal protein set explaining peptide spectrum matches. 
- Utility functions for creating fasta amino acid databases with decoys and contaminants.
- Peptide false discovery rate estimation for target decoy search results on psm, precursor, peptide and protein level. 
- Computing dynamic swath window sizes based on MS1 or MS2 signal distributions.
 
An HTML version of the package documentation can be found here:

https://protviz.github.io/prozor

## How to install:
For CRAN version (not the newest). Please use the github version.

```r
install.packages("prozor")
```

This is how you install the github version.

```r
install.packages("remotes")
remotes::install_github("protviz/prozor")
```

### for Developers

downlod git repo. Use roxygenize2 to document new functions. Than run these 2 commands to update namespace and Rd files:

```r
library(devtools)
document()
```


Example for creating a fasta file with the `fgcz_create_fasta.R` script

Go to fgcz-r-035.uzh.ch


```bash
ls ./fasta_db/p3071_Chlorella
more ./fasta_db/p3071_Chlorella/annotation.txt
more ./fasta_db/p3071_Chlorella/uniprot-taxonomy_3071.fasta
clear
/usr/local/lib/R/site-library/prozor/script/fgcz_create_fasta.R -h

/usr/local/lib/R/site-library/prozor/script/fgcz_create_fasta.R nodecoy ./fasta_db/p3071_Chlorella -o /srv/www/htdocs/FASTA/
/usr/local/lib/R/site-library/prozor/script/fgcz_create_fasta.R ./fasta_db/p3071_Chlorella -o /srv/www/htdocs/FASTA/

cat p3071_Chlorella_d.txt
cat p3071_Chlorella_d.txt | bfabric_save_fasta.py 3071  /srv/www/htdocs/FASTA/fgcz_3071_Chlorella_d_20200604.fasta
# 3071 is the bfabric project name
```
