########### [INFO] Script 0 General information ############
# Script to generate indexed VCF file from a non indexed VCF file

########### [INFO] Input and Output ############################
# Inputs : a non indexed VCF file
# Output : Instance file for a VCF file and his tabix file 
# + indexed bgzip file
# + tabix file

########### [INFO] Previous requirement (libraries) ###############
# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsamtools")

########### [INFO] Problems ############################
# In TIBCO, problem :
# TIBCO Enterprise Runtime for R returned warnings
# In function(e) : package 'Rsamtools' not found: using .GlobalEnv instead

########### [Code] Main method ########
indexVCFfile <- function(path_to_vcf, path_to_indexed_files) {
  library(Rsamtools)

  from <- path_to_vcf
  to <- tempfile(pattern = "indexedVCF.", tmpdir = path_to_indexed_files, fileext = "")
  zipped <- bgzip(from, to)
  idx <- indexTabix(file = zipped, format = "vcf", comment = "#")
  tab <- TabixFile(zipped, idx)
  return (tab)
}

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  ########### [TIBCO] Load RinR library ###################
  library(RinR)
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
  
  ########### [TIBCO] Get the input variables ############
  # path_to_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\samples\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
  path_to_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\samples\\NA12878.chr.22.vcf"
  path_to_indexed_files <-"C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\test_files\\indexedVCFfiles"
  
  ########### [TIBCO] Create the REvaluate object to execute main method ########
  locIndexedVCF <- REvaluate({
    loc_indexed_VCF <- indexVCFfile(path_to_vcf, path_to_indexed_files)
  }
  , data = list(indexVCFfile = indexVCFfile, path_to_vcf = path_to_vcf, path_to_indexed_files = path_to_indexed_files)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
  ########### [TIBCO] Convert the instance file to a Blob Object ########
  locIndexedVCFBlob <- SObjectToBlob(locIndexedVCF)
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Get the input variables ############
  # path_to_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
  path_to_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\samples\\NA12878.chr.22.vcf"
  path_to_indexed_files <-"C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\test_files\\indexedVCFfiles"
  
  ########### [RStudio] Execute main method ###########
  locIndexedVCF <- indexVCFfile(path_to_vcf, path_to_indexed_files)

  ########### [RStudio] Print the ouput in a txt file ###########
  sink("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\test_files\\indexedVCFfiles\\locIndexedVCF")
  str(locIndexedVCF)
  sink()
}