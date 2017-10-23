############## [INFO] SCRIPT 1 General information ############
# PROJECT : Create an R script capable of reading in a vcf formatted file
# and extracting a list of variants.
# NOTE : For this TESTING version, only the first 5 variants are retrieved

############## [INFO] Input and Output ############################
# INPUT : a VCF file. 
# OUTPUT : a list of variants. The list will be stored as a dataframe. For each variant,
#         the dataframe  will contain the chromosome, the range, the reference allele and the alterate allele

############# [FACULTATIVE] Install libraries ###################
# source("https://bioconductor.org/biocLite.R")
# biocLite("VariantAnnotation")
# biocLite("GenomicAlignments")
# biocLite("snpStats")
# # install.packages("tibble") # could come with VariantAnnotation package ?

########### Load RinR library #############
library(RinR)

########### Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### Determinate path of the VCF ##############
# Ask user to introduce the path to the VCF
# PENDING CODE OR TIBCO FUNCTION

my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"

########### REvaluate use #############
vcfDataframe <- REvaluate({
  ############# [CODE] Load libraries ###################
  library(GenomicAlignments)
  library(VariantAnnotation)
  
  ############# [CODE] Extract the VCF ###################
  # Extract the VCF as an object type "CollapsedVCF" (VariantAnnotation package)
  vcf_extracted <- readVcf(my_vcf, "hg19")
  
  ############# Determinate a sample for testing #########
  # Get the samples number
  samples_num <- (header(vcf_extracted))@samples
  
  # Get a little sample for testing
  my_sample <- runif(40, min=0, max=samples_num)
  
  ############# [CODE] Extract a GRanges object from the VCF ###################
  # REMEMBER : in this testing version of the script, only the first 5 variants will be extracted
  
  GRsize <- 5 # number of variant to be extracted
  
  # Extract the variants to a Formal class GRanges
  vcfGRanges <- head(rowRanges(vcf_extracted), GRsize)
  
  # Pass the Formal class GRanges to a data frame
  GRangesdf <- data.frame(vcfGRanges)
  
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
  
  ############# Extract the REFERENCE alleles of the variants ###################
  # The REF allele comes from the dataframe
  GRREFal <- as.character(GRangesdf$REF)
  
  ############# Extract the ALTERNATE alleles of the variants ###################
  # ALT allele is the result of a S4 method
  # the data will be extracted directly from cvf_extracted
  # as a DNAStringSetList which will be passed to a dataframe
  GRALTal_df <- data.frame(alt(vcf_extracted)[1:GRsize])
  GRALTal <- GRALTal_df[3][[1]]
  
  ############# [CODE] Build the dataframe with the extracted information ###################
  # Build the dataframe
  vcfDataframe <- data.frame(GRchroms, GRranges, GRREFal, GRALTal)
  vcfDataframe
  }
  ,
  data = list(my_vcf = my_vcf)
)

############# [Rstudio only] Send the dataframe to a txt file ###################
write.table(vcfDataframe,"test_files\\variants_table.txt",sep="\t",row.names=FALSE)