########### [INFO] Script 4 General information ############
# Clinical info extraction for a list of variant

########### [INFO] Input and Output ############################
# Inputs : variant list (rsXXXX format) available in variable basicTable (produced by script 3)
# Output : clinical annotation for the variant

########### [INFO] Previous requirement (libraries) ###############

########### [INFO] Problems ############################
# I need to know if there is a method available in cellbaseR

########### [Code] Main method ########
getClinicalVarInfo <- function(basicTable) {
  variants_vct <- basicTable$rsID
  
  ############# BUILDING QUERY : common part of the GET URL ############
  # Part of the get url which is common
  # NOTE : maybe the common URL could suffer some changes
  common_URL <- "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/"
  category_URL <- "clinical/variant/search"
  subcategory_URL <- N/A
  parameters_URL <- "limit=-1&skip=-1&skipCount=false&count=false&Output%20format=json"
  URL_vector <- character()
  
  for (variantID in variants_vct) {
    
    # build a query like this one :
    # http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/clinical/variant/search?limit=-1&skip=-1&skipCount=false&count=false&Output%20format=json&id=rs6025
    
    # Get the Swagger link
    URL_base <- paste(common_URL,category_URL,"?",parameters_URL,"&id=",variantID,sep = "")
    URL_vector <- c(URL_vector, URL_base)
    
  }
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
