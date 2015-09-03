[![Build Status](https://travis-ci.org/wolski/prozor.svg?branch=master)](https://travis-ci.org/wolski/prozor)
[![Project Stats](https://www.ohloh.net/p/prozor/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/prozor)

# prozor
 Determine minimal Protein set given list of peptide protein mappings. Weights can be assigned to peptides or proteins 
 (inverse peptide frequencies) and various weighting functions can be assigned.

## How to install:
for CRAN version

```r
install.packages("prozor")
```

for development version from github

```r
install.packages("devtools")
library(devtools)
install_github("wolski/prozor")
```

### for developers

downlod git repo. Use roxygenize2 to document new functions. Than run these 2 commands to update namespace and Rd files:

```r
library(prozor)
document()
```
