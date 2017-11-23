########### [INFO] Script [num] General information ############
# Split a VCF file into small VCF.

########### [INFO] Input and Output ############################
# Inputs : VCF file
# Output : Various VCF files, as parts of the input

########### [INFO] Previous requirement (libraries) ###############

########### [INFO] Problems ############################

########### [Code] Main method ########
sampleVCF <- function(locIndexedVCF) {
  library(Rsamtools)
  library(VariantAnnotation)
  
  file.gz <- locIndexedVCF$path
  stopifnot(file.exists(file.gz))
  file.gz.tbi <- paste(file.gz, ".tbi", sep="")
  if(!(file.exists(file.gz.tbi)))
    indexTabix(file.gz, format="vcf")

  # original
  start.loc <- 55000000
  end.loc <- 56000000
  chr7.gr <- GRanges("7", IRanges(start.loc, end.loc))
  params <- ScanVcfParam(which=chr7.gr)
  # original
  
  # testeando
  vcf1 <- readVcf(file.gz, "hg19")
  GRangesVcf1 <- rowRanges(vcf1)
  my_sample_num <- sort(sample(length(GRangesVcf1), 40))
  vcf_sample <- GRangesVcf1[c(my_sample_num)]
  params <- ScanVcfParam(which=vcf_sample)
  # testeando
  
  vcf <- readVcf(TabixFile(file.gz), "hg19", params)
  
  writeVcf(vcf, "chr7-sub.vcf")
  bgzip("chr7-sub.vcf", overwrite=TRUE)
  indexTabix("chr7-sub.vcf.gz", format="vcf")
  
  
  # ############# [TEST] Determinate a sample for testing #########
  # # Get a random sample from the GRanges object for testing
  # my_sample_num <- sort(sample(length(vcfGRanges), 40))
  # vcf_sample <- vcfGRanges[c(my_sample_num)]
  # 
  # # assign the vcf_sample to the main variable (only for testing)
  # vcfGRanges <- vcf_sample
  
  
}

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  ########### [TIBCO] Load RinR library ###################
  library(RinR)
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
  
  ########### [RStudio] Get the input variables ############
  global_input <- "whatever"
  
  ########### [TIBCO] Create the REvaluate object to execute main method ########
  mainVariable <- REvaluate({
    mainVariable <- mainMethod(local_input)
  }
  , data = list(mainMethod = mainMethod, local_input = global_input)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Get the input variables ############
  global_input <- "whatever"
  
  ########### [RStudio] Execute main method ###########
  mainVariable <- mainMethod(global_input)
  
  ########### [RStudio] Print the ouput in a txt file ###########
}

########### [TEST block] ###############
