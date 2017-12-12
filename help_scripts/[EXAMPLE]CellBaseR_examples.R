############## Testing cellbaseR : General information ##################
# This script is to test the library

############# [FACULTATIVE] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [FACULTATIVE] Install libraries ###################
# source("https://bioconductor.org/biocLite.R")
# ## try http:// if https:// URLs are not supported
# biocLite("cellbaseR")

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

################ bach ###############
cb <- CellBaseR()
res <- getVariant(object=cb, ids="19:45411941:T:C", resource="annotation")


###### getCellBaseResourceHelp #######
cb <- CellBaseR()
# Get help about what resources are available to the getGene method
getCellBaseResourceHelp(cb, subcategory="gene")
# Get help about what resources are available to the getRegion method
getCellBaseResourceHelp(cb, subcategory="region")
# Get help about what resources are available to the getXref method
getCellBaseResourceHelp(cb, subcategory="id")

getCellBaseResourceHelp(cb, subcategory="variant")
getCellBaseResourceHelp(cb, subcategory="clinical")

res <- getMeta(object=cb, resource="species")
res[1,]
