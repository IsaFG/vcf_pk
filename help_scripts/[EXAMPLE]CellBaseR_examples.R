############## Testing cellbaseR : General information ##################
# This script is to test the library

############# [FACULTATIVE] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [FACULTATIVE] Install libraries ###################
source("https://bioconductor.org/biocLite.R")
## try http:// if https:// URLs are not supported
biocLite("cellbaseR")

############# [CODE] Load libraries ###################
library(cellbaseR)
cb <- CellBaseR()

############ getCellBase ##############
# getCellBase 8/26
# Description
# The generic method for querying CellBase web services.
cb <- CellBaseR()
res <- getCellBase(object=cb, category="feature", subcategory="gene",
                   ids="TET1", resource="info")

############ getGene ##############
cb <- CellBaseR()
genes <- c("TP73","TET1")
res <- getGene(object = cb, ids = genes, resource = "info")
str(res,2)

############ getVariant ##############
cb <- CellBaseR()
res2 <- getVariant(object=cb, ids="1:169549811:A:G", resource="annotation")
# to get the data 
# res2 <- cbData(res2) #this line throw an error
str(res2, 1)

########### getClinical ##############
cb <- CellBaseR()
cbParam <- CellBaseParam(feature=c("TP73","TET1"), limit=100)
res3 <- getClinical(object=cb,param=cbParam)
