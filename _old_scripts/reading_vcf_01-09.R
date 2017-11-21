# Create an R script capable of reading in a vcf formatted file and extracting a list of variants.
# The list has to be extracted to be used in this format :
# http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/{category}/{subcategory}/{id}/{resource}?{filters}
# For example, we have these links :
# "bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
# "bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/1:18149476-18149476/gene"
# these REST call will get all the INFO of the gene BRCA2,BRCA1,KRAS of human in the latest version.

##########################################
## OPTIONAL : TO INSTALL THE LIBRARIES TO BE USE
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
biocLite("GenomicAlignments")
biocLite("snpStats")
install.packages("tibble")
############################################

# load libraries
library(GenomicAlignments)
library(VariantAnnotation)

############################################
# METHOD : EXTRACT THE VCF
# Extract the VCF to an object type "CollapsedVCF" with VariantAnnotation
extractVCF <- function(vcfname) {
  my_vcf <- vcfname
  vcf2 <- readVcf(my_vcf, "hg19")
  return (vcf2)}

############################################
# METHOD : Print info from the VCF extracted
## NO FUNCIONA ?????
get_header <- function(vcf_extracted) {
  vcfHeader <- header(vcf_extracted)
  print (vcfHeader)}

############################################
# METHOD : Extract a GRanges object from the VCF, with 5 first variants
get_GRanges <- function (vcf_extracted){
  # Variable to get the desired size of the Granges
  # In this example, I will set the number as 5
  GRsize <- 5
  rowRangeVCF <- head(rowRanges(vcf_extracted), GRsize)
  return(rowRangeVCF)}

############################################
# METHOD : Build the dataframe with cromosomes and ranges
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

##########################################
# CALLS
# extract the vcf
vcf_extracted <- extractVCF("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf")
# print info header
# vcfHeader <- get_header # NO FUNCIONA
# get the Granges
vcfGRanges <- get_GRanges(vcf_extracted)
# get de dataframe
vcfDataframe <- get_dataframe(vcfGRanges)

# print info header
header(vcf_extracted)

vcfGRanges
vcfDataframe

