############## How to install a package in R ################
# It depends where the package is hosted

############## From CRAN (Standard) #####################
# https://cran.r-project.org/

# Example with the package "jsonlite"
install.packages("jsonlite")

############## From Bioconductor #####################
# https://www.bioconductor.org/



source("https://bioconductor.org/biocLite.R")
## try http:// if https:// URLs are not supported
biocLite("cellbaseR")

############## From GITHUB #####################
# You have to install devtools to install libraries from Github
install.packages("devtools") 
library(devtools)
# Example with the package "cellbaseR" develloped by user "melsiddieg"
devtools::install_github("melsiddieg/cellbaseR")
