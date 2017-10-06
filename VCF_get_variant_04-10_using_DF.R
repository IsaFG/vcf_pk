############## [INFO] SCRIPT 1 General information ############
# PROJECT : Create an R script capable of reading in a vcf formatted file
# and extracting a list of variants.
# NOTE : This script is incomplete because of a problem:
# [PROBLEM] Extract the ALTERNATE alleles of the variants.
## See if in cellbase there is not a tip for that 

############## [INFO] Input and Output ############################
# INPUT : a VCF file. 
# OUTPUT : a list of variants. The list will be stored as a dataframe. For each variant,
#         the dataframe  will contain the chromosome, the range, the reference allele and the alterate allele
# NOTE : For this TESTING version, only the first 5 variants are retrieved

############# [FACULTATIVE] Install libraries ###################
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
biocLite("GenomicAlignments")
biocLite("snpStats")
# install.packages("tibble") # could come with VariantAnnotation package ?

############# [FACULTATIVE] Set Working directory ############
# Change it by the one you want
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [CODE] Load libraries ###################
library(GenomicAlignments)
library(VariantAnnotation)

############# [CODE] Extract the VCF ###################
# Extract the VCF as an object type "CollapsedVCF" (VariantAnnotation package)
my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
vcf_extracted <- readVcf(my_vcf, "hg19")

############# [CODE] Extract a GRanges object from the VCF ###################
# REMEMBER : in this testing version of the script, only the first 5 variants will be extracted

GRsize <- 5 # number of variant to be extracted

# Extract the variants to a Formal class GRanges
vcfGRanges <- head(rowRanges(vcf_extracted), GRsize)

# This block is for testing
# Pass the Formal class GRanges to a txt file
sink("test_files\\vcfGRanges.txt")
vcfGRanges
sink()

# The following block is for visualisation
# Pass the Formal class GRanges to a txt file
sink("test_files\\vcfGRanges_str.txt")
str(vcfGRanges)
sink()

# Pass the Formal class GRanges to a data frame
GRangesdf <- data.frame(vcfGRanges)

# There is a PROBLEM here, 
# the ALT alleles are stored as result from S4 object.
# See below for more information about this issue.

# The following block is for visualisation
# Pass the data.frame to a txt file
sink("test_files\\GRangesdf.txt")
write.table(GRangesdf,"test_files\\GRangesdf.txt",sep="\t",row.names=FALSE)
sink()

# The following block is for visualisation
sink("test_files\\ALTdf.txt")
str(ALTdf)
sink()

############# [CODE] Extract the chromosomes and the ranges of the variants ###################
# Extract chromosomes, ranges from the Formal class GRanges
# and store them in vectors
# tip got from https://www.biostars.org/p/137765/
GRranges <- as.character(ranges(vcfGRanges))
GRchroms_chr <- as.character(seqnames(vcfGRanges))

# Next to that, remove the "chr" characters in the chromosome vector, in order to only maintain the number
# store it in a new vector
string_to_remove <- "chr"
GRchroms <- gsub(x = GRchroms_chr, pattern = paste(string_to_remove, collapse = "|"), replacement = "")

############# [CODE] Extract the REFERENCE alleles of the variants ###################
# The REF allele comes from the dataframe
GRREFal <- as.character(GRangesdf$REF)

############# [PROBLEM] Extract the ALTERNATE alleles of the variants ###################
# Since the ALT allele is the result of a S4 method, I didnt figure out yet how to retrive it
GRangesdf$ALT
class(GRangesdf$ALT)

ALTdf <- data.frame(GRangesdf$ALT)
print(ALTdf[2,1])
class(ALTdf)
str(ALTdf)
ALTdf

# Just for testing :
GRALTal <- c("A", "T", "A", "G", "A")
  
############# [CODE] Build the dataframe with the extracted information ###################
# Build the dataframe
vcfDataframe <- data.frame(GRchroms, GRranges, GRREFal)
vcfDataframe <- data.frame(GRchroms, GRranges, GRREFal, GRALTal)

############# [CODE] Send the dataframe to a txt file ###################
write.table(vcfDataframe,"test_files\\variants_table.txt",sep="\t",row.names=FALSE)

############ [TEST] Testing block ###########################
# The following block is for testing
# # print info header
# header(vcf_extracted)
# # print the GRange object and the obtained dataframe
# vcfGRanges
# vcfDataframe