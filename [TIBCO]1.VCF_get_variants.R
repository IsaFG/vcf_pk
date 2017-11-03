############## [INFO] SCRIPT 1 General information ############
# PROJECT : Create an R script capable of reading in a vcf formatted file
# and extracting a list of variants.

############## [INFO] Input and Output ############################
# INPUT : a VCF file. 
# OUTPUT : a list of variants. The list will be stored as a dataframe. For each variant,
#         the dataframe  will contain the chromosome, the range, the reference allele and the alterate allele

############# [INFO] Previous requirement (libraries) ###############
# source("https://bioconductor.org/biocLite.R")
# biocLite("VariantAnnotation")
# biocLite("GenomicAlignments")
# biocLite("snpStats")
# # install.packages("tibble") # could come with VariantAnnotation package ?

########### [INFO] PROBLEMS ########
# Warning message:
#   In .bcfHeaderAsSimpleList(header) :
#   duplicate keys in header will be forced to unique rownames

########### [Code] Main method ########
getvcfDataframe <- function(my_vcf) {
  ############# [CODE] Load libraries ###################
  library(GenomicAlignments)
  library(VariantAnnotation)
  
  ############# [CODE] Extract the VCF ###################
  # Extract the VCF as an object type "CollapsedVCF" (VariantAnnotation package)
  vcf_extracted <- readVcf(my_vcf, "hg19")
  
  ############# [CODE] Extract a GRanges object from the CollapsedVCF ###################
  vcfGRanges <- rowRanges(vcf_extracted)
  
  # ############# [TEST] Determinate a sample for testing #########
  # # Get a random sample from the GRanges object for testing
  # my_sample_num <- sort(sample(length(vcfGRanges), 40))
  # vcf_sample <- vcfGRanges[c(my_sample_num)]
  # 
  # # assign the vcf_sample to the main variable (only for testing)
  # vcfGRanges <- vcf_sample
  
  ############# [CODE] Pass the Formal class GRanges to a data frame ###################
  GRangesdf <- data.frame(vcfGRanges)
  
  ranges_start <- start(ranges(vcfGRanges))
  ranges_end <- end(ranges(vcfGRanges))
  chroms_vec_chr <- as.character(seqnames(vcfGRanges))
  
  # Next to that, remove the "chr" characters in the chromosome vector,
  # in order to only maintain the number, then store it in a new vector
  string_to_remove <- "chr"
  chroms_vec <- gsub(x = chroms_vec_chr, pattern = paste(string_to_remove, collapse = "|"), replacement = "")
  
  ############# [CODE] Extract the REFERENCE alleles of the variants ###################
  # The REF allele comes from the dataframe
  REFal <- as.character(vcfGRanges$REF)
  
  ############# [CODE] Extract the ALTERNATE alleles of the variants ###################
  # ALT allele is the result of a S4 method
  # the data will be extracted directly from vcf_extracted
  # as a DNAStringSetList which will be passed to a dataframe
  # Then the duplicated value (which display in different rows) will be unificate in an unique row
  
  ALTal_df <- data.frame(vcfGRanges$ALT)
  i <- 1
  for (boolresult in duplicated(ALTal_df[,1])) {
    if (boolresult == TRUE) {
      ALTal_df[,3][[i-1]] <- paste(ALTal_df[,3][[i-1]], ALTal_df[,3][[i]], sep=",")
      ALTal_df <- ALTal_df[-c(i),]
    }
    i <- i + 1
  }
  ALTal <- ALTal_df[3][[1]]
  
  ############# [CODE] Build the dataframe with the extracted information ###################
  # Build the dataframe
  vcfDataframe <- data.frame(chroms_vec, ranges_start, ranges_end, REFal, ALTal)
  colnames(vcfDataframe) <- c("Chromosome", "Range start", "Range end", "Ref allele", "Alt allele")
  
  return(vcfDataframe)
}

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  ########### [TIBCO] Load RinR library #############
  library(RinR)
  
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
  
  ########### [TIBCO] Determinate path of the VCF ##############
  # Ask user to introduce the path to the VCF
  # PENDING CODE OR TIBCO FUNCTION
  my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
  
  ########### [TIBCO] Create the REvaluate object ########
  vcfDataframe <- REvaluate({
    vcfDataframe <- getvcfDataframe(my_vcf)
    vcfDataframe
    }
    ,
    data = list(getvcfDataframe = getvcfDataframe, my_vcf = my_vcf)
  )
  
  } else {
    
    # ############# [RStudio] Set Working directory ############
    # Please uncoment the next line changin the working directory by the correct one:
    setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
    
    ########### [RStudio] Determinate path of the VCF ##############
    # Please edit the following code line to put the correct path to your VCF
    my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
    
    ########### [RStudio] Get the dataframe ########################
    vcfDataframe <- getvcfDataframe(my_vcf)
    
    ########### [Rstudio] Send the dataframe to a txt file ###################
    write.table(vcfDataframe,"test_files\\variants_table.txt",sep="\t",row.names=FALSE)
  }
