############## [INFO] SCRIPT 1 General information ############
# PROJECT : Create an R script capable of reading in a vcf formatted file
# and extracting a list of variants.
# IMPORTANT NOTE : this script has been written with functions, is not linear

############## [INFO] Input and Output ############################
# INPUT : a VCF file. 
# OUTPUT : a list of variants. The list will be stored as a dataframe.
## NOTE 1 : For the first version of the script, the list will have only
## the cromosome and the range of each variant
## NOTE 2 : For this TESTING version, only the first 5 variants are retrieved

############# [FACULTATIVE] Install libraries ###################
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
biocLite("GenomicAlignments")
biocLite("snpStats")
# install.packages("tibble") # could come with VariantAnnotation package ?

############# [FACULTATIVE] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [CODE] Load libraries ###################
library(GenomicAlignments)
library(VariantAnnotation)

############# [CODE] [PENDING] Define clases ? ###################
# I dont know if this would be better with a defined class ?
# Define the class - NOT DONE YET
# reading_vcf <- 

############# [CODE] Extract the VCF ###################
# Extract the VCF to an object type "CollapsedVCF" with VariantAnnotation
extractVCF <- function(vcfname) {
  my_vcf <- vcfname
  vcf2 <- readVcf(my_vcf, "hg19")
  return (vcf2)}

############# [CODE] [TEST] Print info from the VCF extracted ###################
# ## NO FUNCIONA ?????
# get_header <- function(vcf_extracted) {
#   vcfHeader <- header(vcf_extracted)
#   print (vcfHeader)}

############# [CODE] Extract a GRanges object from the VCF ###################
# NOTE : in the example file, only the first 5 variants will be extracted
get_GRanges <- function (vcf_extracted){
  # Variable to get the desired size of the Granges
  # In this example, I will set the number as 5
  GRsize <- 5
  rowRangeVCF <- head(rowRanges(vcf_extracted), GRsize)
  return(rowRangeVCF)}

############# [CODE] Build the dataframe with cromosomes and ranges ###################
get_dataframe <- function(vcfGRanges) {
  # Extract chromosomes and ranges from the GRanges object and store them in both vectors
  # tip got from https://www.biostars.org/p/137765/
  GRranges <- as.character(ranges(vcfGRanges))
  GRchroms_chr <- as.character(seqnames(vcfGRanges))
  # Next to that, remove the "chr" characters in the chromosome vector, in order to only maintain the number
  # store it in a new vector
  string_to_remove <- "chr"
  GRchroms <- gsub(x = GRchroms_chr, pattern = paste(string_to_remove, collapse = "|"), replacement = "")
  # Build the dataframe
  variant.data <- data.frame(GRchroms, GRranges)
  return(variant.data)}

############# [CODE] Send the dataframe to a txt file ###################
get_datafile <- function(vcfDataframe) {
  write.table(vcfDataframe,"test_files\\Crom_Loc_Aleles_table.txt",sep="\t",row.names=FALSE)
}

############# [CODE] Calls to all the previews method ###################
# extract the vcf
vcf_extracted <- extractVCF("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf")
# print info header
# vcfHeader <- get_header # NO FUNCIONA
# get the Granges
vcfGRanges <- get_GRanges(vcf_extracted)
# get de dataframe
vcfDataframe <- get_dataframe(vcfGRanges)
# generate a text file with the dataframe
vcfDatafile <- get_datafile(vcfDataframe)

############ [TEST] Testing block ###########################
# The following block is fot testing
# # print info header
# header(vcf_extracted)
# # print the GRange object and the obtained dataframe
# vcfGRanges
# vcfDataframe