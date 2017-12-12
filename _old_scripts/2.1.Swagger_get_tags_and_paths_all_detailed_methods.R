########### [INFO] Script [num] General information ############


########### [INFO] Input and Output ############################
# Inputs : 
# Output : 

########### [INFO] Previous requirement (libraries) ###############

########### [INFO] Problems ############################

########### [Code] Main methods #############
getCBjson <- function(cellbase_URL) {
  ############# [CODE] Load libraries ###################
  library(RCurl)
  library(jsonlite)
  
  ############# [CODE] Get the main Swagger.json from CELLBASE ###############
  swagger_URL <- getURL(cellbase_URL)
  swagger_JSON <-fromJSON(swagger_URL)
  return(swagger_JSON)
}

getJSONinfo <- function(swagger_JSON) {
  ############# [CODE] Load libraries ###################
  library(RCurl)
  library(jsonlite)
  
  # For map function
  library(purrr)
  library(tidyverse)
  library(rjson)
  
  ############# [CODE] Declare the list info list #####
  JSONinfo <- list()
  
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

getTagsList <- function(swagger_JSON) {
  ############# [CODE] Load libraries ###################
  library(RCurl)
  library(jsonlite)
  
  # For map function
  library(purrr)
  library(tidyverse)
  library(rjson)
  
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
  return(tags_list)
}
  
getPathsList <- function(swagger_JSON) {
  ############# [CODE] Load libraries ###################
  library(RCurl)
  library(jsonlite)
  
  # For map function
  library(purrr)
  library(tidyverse)
  library(rjson)
  
  ############# [CODE] Get the GET paths ###############
  level1_keys <- map(swagger_JSON, names)
  
  # Create a vector with all the paths
  # This creates problems when you pass it to a txt file in order to use it in TIBCO
  paths_list <- level1_keys$paths
  
  return(paths_list)
  }

getPathsDF <- function(swagger_JSON){
  ############# [CODE] Load libraries ###################
  library(RCurl)
  library(jsonlite)
  
  # For map function
  library(purrr)
  library(tidyverse)
  library(rjson)
  
  ############# [CODE] Get the GET paths ###############
  level1_keys <- map(swagger_JSON, names)
  
  # Create a dataframe from the vector
  # This creates problems when you pass it to a txt file in order to use it in TIBCO
  paths_df <- data.frame(level1_keys$paths)
  
  return(paths_df)
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
    swagger_JSON <- getCBjson(cellbase_URL)
    JSONinfo <- getJSONinfo(swagger_JSON)
    JSONinfo
  }
  , data = list(getCBjson = getCBjson, getJSONinfo = getJSONinfo, cellbase_URL = cellbase_URL)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
  tagList <- CBinfo$tagList
  pathsList <- CBinfo$pathsList
  pathsDF <- CBinfo$pathsDF
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Execute main method ###########
  swagger_JSON <- getCBjson(cellbase_URL)
  CBinfo <- getJSONinfo(swagger_JSON)
  
  tagList <- CBinfo$tagList
  pathsList <- CBinfo$pathsList
  pathsDF <- CBinfo$pathsDF
  
  ########### [RStudio] Print the ouput in a txt file ###########
}

########### [TEST block] ###############
