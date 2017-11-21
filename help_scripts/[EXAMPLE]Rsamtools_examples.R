############## Testing Rsamtools : General information ##################
# This script is to test the library

########### [INFO] Previous requirement (libraries) ###############
# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsamtools")

########### Load library #############
library(Rsamtools)

########## Indextabix ###############
from <- system.file("extdata", "ex1.sam", package="Rsamtools",
                    mustWork=TRUE)
to <- tempfile()
zipped <- bgzip(from, to)
idx <- indexTabix(zipped, "sam")

tab <- TabixFile(zipped, idx) 

