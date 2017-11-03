# When running script 3
#  Error: Required package curl not found. Please run: install.packages('curl') 
install.packages('curl')


# Script 1
library(VariantAnnotation)
# Error: package or namespace load failed for 'VariantAnnotation' in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]):
#  there is no package called 'GenomicFeatures'
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
