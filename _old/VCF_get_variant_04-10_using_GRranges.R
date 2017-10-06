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
my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
vcf_extracted <- readVcf(my_vcf, "hg19")

############# [CODE] [TEST] Print info from the VCF extracted ###################
# ## NO FUNCIONA ?????
# get_header <- function(vcf_extracted) {
#   vcfHeader <- header(vcf_extracted)
#   print (vcfHeader)}

############# [CODE] Extract a GRanges object from the VCF ###################
# NOTE : in the example file, only the first 5 variants will be extracted
# Variable to get the desired size of the Granges
# In this example, I will set the number as 5
GRsize <- 5
vcfGRanges <- head(rowRanges(vcf_extracted), GRsize)
GRangesdf <- data.frame(vcfGRanges)

# This block is for testing
sink("test_files\\vcfGRanges.txt")
vcfGRanges
sink()

############# [CODE] Build the dataframe with cromosomes and ranges ###################
# Extract chromosomes, ranges and ALT and REF alleles from the GRanges object
# and store them in vectors
# tip got from https://www.biostars.org/p/137765/
GRranges <- as.character(ranges(vcfGRanges))
GRchroms_chr <- as.character(seqnames(vcfGRanges))
GRREFal <- as.character(REF(vcfGRanges))
GRALTal <- as.character(ALT(vcfGRanges))
# Next to that, remove the "chr" characters in the chromosome vector, in order to only maintain the number
# store it in a new vector
string_to_remove <- "chr"
GRchroms <- gsub(x = GRchroms_chr, pattern = paste(string_to_remove, collapse = "|"), replacement = "")
# Build the dataframe
vcfDataframe <- data.frame(GRchroms, GRranges, GRREFal, GRALTal)

############# [CODE] Send the dataframe to a txt file ###################
write.table(vcfDataframe,"test_files\\Crom_Loc_Aleles_table.txt",sep="\t",row.names=FALSE)


############ [TEST] Testing block ###########################
# The following block is fot testing
# # print info header
# header(vcf_extracted)
# # print the GRange object and the obtained dataframe
# vcfGRanges
# vcfDataframe