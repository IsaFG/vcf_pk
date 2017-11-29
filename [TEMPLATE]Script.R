########### [INFO] Script [num] General information ############


########### [INFO] Input and Output ############################
# Inputs : 
# Output : 

########### [INFO] Previous requirement (libraries) ###############

########### [INFO] Problems ############################

########### [Code] Main method ########
mainMethod <- function(input) {
  
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
