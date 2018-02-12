########### [INFO] SCRIPT x General information ############
# 

########### [INFO] Input and Output ############################
# Inputs : 
# Output : 

########### [INFO] Previous requirement (libraries) ###############

########### [INFO] Problems ############################

########### [Code] Main methods #############
mainMethod <- function(input) {
  ############# [CODE] Load libraries ###################

  test <- merge(clinvarDisease, cytobandTable[, c("rsID", "end")], by="rsID", all.x=TRUE) 
  
  test <- merge(clinvarDisease, cytobandTable[, c("rsID", "end")], by="rsID", all.x=TRUE) 
  
  names(test)[names(test) == "rsID"] <- "rs"
  
  ############# [CODE] Method ###############
  
  return(output)
}

########### [Code] Set main variables (if needed) ###############


########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  
  ########### [TIBCO] Load RinR library ###################
  library(RinR)
  
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
  
  ########### [TIBCO] Create the REvaluate object to execute main method ########
  globalOutput <- REvaluate({
    
    localOutput <- mainMethod(localInput)

  }
  , data = list(localMethod = mainMethod, localInput = globalInput)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Execute main method ###########
  globalOutput <- mainMethod(globalInput)
  
  ########### [RStudio] Print the ouput in a txt file ###########

}

########### [TEST block] ###############
