########### [INFO] SCRIPT 2.1 General information ############
# This script will extract structural information from CellBase SWAGGER CELL BASE,
# in order to obtain categories and link's structure (main tags and paths)
# of the API.

########### [INFO] Input and Output ############################
# Inputs : cellbase URL
# Output : paths list in dataframe and character format
# + tags list in character format

########### [INFO] Previous requirement (libraries) ###############
# install.packages("RCurl")
# install.packages("jsonlite")
# install.packages('purrr')

########### [INFO] Problems ############################
# Some problems can appear depending on the libraries that are loaded.
# But so far no problems

########### [Code] Main methods #############
getCBinfo <- function(cellbase_URL) {
  ############# [CODE] Load libraries ###################
  library(RCurl)
  library(jsonlite)
  
  # For map function
  library(purrr)
  # library(tidyverse) # not necessary in the last versions and can create bugs
  # library(rjson) # not necessary in the last versions and can create bugs
  
  ############# [CODE] Get the main Swagger.json from CELLBASE ###############
  swagger_URL <- getURL(cellbase_URL)
  swagger_JSON <- fromJSON(swagger_URL)

  ############# [CODE] Declare the list info list #####
  JSON_info <- list()
  
  ############# [CODE] Get the main categories (tags) ###############
  # The tags will be adden to a vector of tags
  extract_tags <- swagger_JSON$tags
  # extract_tags # this line is for testing
  tags <- flatten(extract_tags)
  tags_list <- character()
  for (i in tags){
    tag <- tags[i]
    # print(paste("The element is: ", i)) # this line is for testing
    tags_list <- c(tags_list, i)
  }
  JSON_info$tagList <- tags_list
  
  ############# [CODE] Get the GET paths ###############
  level1_keys <- map(swagger_JSON, names)
  
  # Create a vector with all the paths
  # This creates problems when you pass it to a txt file in order to use it in TIBCO
  paths_list <- level1_keys$paths
  
  JSON_info$pathsList <- paths_list
  
  ############# [CODE] Get the GET paths ###############
  level1_keys <- map(swagger_JSON, names)
  
  # Create a dataframe from the vector
  # This creates problems when you pass it to a txt file in order to use it in TIBCO
  paths_df <- data.frame(level1_keys$paths)
  
  JSON_info$pathsDF <- paths_df
  
  return(JSON_info)
}

########### [Code] Set the cellbase URL ###############
cellbase_URL <- "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/swagger.json"

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  
  ########### [TIBCO] Load RinR library ###################
  library(RinR)
  
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
  
  ########### [TIBCO] Create the REvaluate object to execute main method ########
  CBinfo <- REvaluate({
    JSONinfo <- getCBinfo(cellbase_URL)
    JSONinfo
  }
  , data = list(getCBinfo = getCBinfo, cellbase_URL = cellbase_URL)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
  CBinfoBlob <- SObjectToBlob(CBinfo)
  tagList <- CBinfo$tagList
  pathsList <- CBinfo$pathsList
  pathsDF <- CBinfo$pathsDF
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Execute main method ###########
  CBinfo <- getCBinfo(cellbase_URL)
  
  tagsList <- CBinfo$tagList
  pathsList <- CBinfo$pathsList
  pathsDF <- CBinfo$pathsDF
  
  ########### [RStudio] Print the ouput in a txt file ###########
  sink("test_files\\tags_list.txt")
  tagsList
  sink()
  
  sink("test_files\\paths_df.txt")
  pathsDF
  sink()
  
  sink("test_files\\paths_list.txt")
  pathsList
  sink()
}

########### [TEST block] ###############
