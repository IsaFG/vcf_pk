########### [INFO] Script 0 General information ############
# Script to generate indexed VCF file from a non indexed VCF file

########### [INFO] Input and Output ############################
# Inputs : a non indexed VCF file. Please ensure to edit the path if you run this script in R.
# Output : Instance file for a VCF file and his tabix file 
# + indexed bgzip file
# + tabix file

########### [INFO] Previous requirement (libraries) ###############
# Library "Rsamtools"

# Working directory :
# Please uncoment the next line changing the working directory by the correct one:
# working_dir <- "C:\\..."

# R version : please uncoment the next line indicating the location of your current R version:
# r_version <- "C:/Program Files/R/R-3.4.1"

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
  Rversion <- makeREvaluator("R", RHome = r_version)
  
  ########### [TIBCO] Get the input variables ############
  # PLEASE EDIT THE PATH 
  path_to_vcf <- "samples\\NA12878.chr.22.vcf"
  path_to_indexed_files <-"test_files\\indexedVCFfiles"
  
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
  setwd(working_dir)
  
  ########### [RStudio] Get the input variables ############
  path_to_vcf <- "samples\\NA12878.chr.22.vcf"
  path_to_indexed_files <-"test_files\\indexedVCFfiles"
  
  ########### [RStudio] Execute main method ###########
  locIndexedVCF <- indexVCFfile(path_to_vcf, path_to_indexed_files)

  ########### [RStudio] Print the ouput in a txt file ###########
  sink("test_files\\indexedVCFfiles\\locIndexedVCF")
  str(locIndexedVCF)
  sink()
}
