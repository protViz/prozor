[![Build Status](https://travis-ci.org/wolski/prozor.svg?branch=master)](https://travis-ci.org/wolski/prozor)
[![Project Stats](https://www.ohloh.net/p/prozor/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/prozor)

# prozor
 Determine minimal Protein set given list of peptide protein mappings. Various weights can be assigned to peptides (i.e. inverse peptide frequencies).

## How to install:
For CRAN version (not the newest)

```r
install.packages("prozor")
```

For development version from github

```r
install.packages("devtools")
library(devtools)
install_github("protviz/prozor")
```

### for developers

downlod git repo. Use roxygenize2 to document new functions. Than run these 2 commands to update namespace and Rd files:

```r
library(prozor)
document()
```
