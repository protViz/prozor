[![Build Status](https://travis-ci.org/wolski/prozor.svg?branch=master)](https://travis-ci.org/wolski/prozor)
[![Project Stats](https://www.ohloh.net/p/prozor/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/prozor)

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
library(prozor)
document()
```
