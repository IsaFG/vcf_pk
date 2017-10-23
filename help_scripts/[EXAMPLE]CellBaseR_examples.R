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

########### Annotate VCF ############
library(VariantAnnotation)

fl <- system.file("extdata", "hapmap_exome_chr22_500.vcf.gz",package = "cellbaseR" )
res <- AnnotateVcf(object=cb, file=fl, BPPARAM = bpparam(workers=2))
vcf <- readVcf(fl, "hg19")
samples(header(vcf))

########## Test Annotate VCF #############
my_vcf <- "C:/Users/FollonIn/Documents/GitHub/vcf_pk/R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
res3 <- AnnotateVcf(object=cb, file=my_vcf, BPPARAM = bpparam(workers=2))
vcf3 <- readVcf(fl3, "hg19")
samples(header(vcf3))

bgzip(from, to)

fl
