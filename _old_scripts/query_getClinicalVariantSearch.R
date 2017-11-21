########### [INFO] Script 4 General information ############
# Clinical info extraction for a list of variant

########### [INFO] Input and Output ############################
# Inputs : variant list (rsXXXX format) available in variable basicTable (produced by script 3)
# As a test, you can use the variant rs2296035
# Output : clinical annotation for the variant

########### [INFO] Previous requirement (libraries) ###############

########### [INFO] Problems ############################

########### [Code] Main method ########
getClinicalVarInfo <- function(basicTable) {
  
  ############# [Code] Load libraries ############# 
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  library(dplyr)
  library(tidyr)
  
  ############# [Code] Build the GET URL and query CellBase (CellBaseR) ############
  variants_vct <- basicTable$rsID
  cb <- CellBaseR()
  
  # Initialize the annotation table
  clinicalVariants_table <- data.frame()
  # Initialize a table to record any problem about dimensions in the annotated table
  problems <- data.frame(variant=integer(0), column=integer(0), dimension=integer(0), class=character(0)) # testing line

  
  cbParam <- CellBaseParam(rsid="rs2296035")
  cbParam <- CellBaseParam(rsid=variants_vct)
  
  print (paste("Processing variant number:", i)) # this line is for testing
  
  # Get clinical info with cellbaseR package
  # the call will return as a data.frame
  clinicalVariant <- getClinical(object=cb,param=cbParam)

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
