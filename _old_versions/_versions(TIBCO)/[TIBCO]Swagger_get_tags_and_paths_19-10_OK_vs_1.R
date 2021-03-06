############## [INFO] SCRIPT 2.1 General information ############
# This script works on TIBCO
# but is not as efficient that [TIBCO]Swagger_get_tags_and_paths_19-10_NOT_OK

# TIBCO script
# This script will extract structural information from SWAGGER CELL BASE,
# in order to obtain categories and link's structure (main tags and paths)
# of the API.

############## [INFO] Input and Output ############################
# Input : We will work with the resource list, which is the main JSON of the CellBase
# "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/swagger.json

# Output :
# 1) a TIBCO table with all the main tags
# 2) a TIBCO table with all the availables paths

# ############# [FACULTATIVE] Install libraries ###################
# install.packages("RCurl")
# install.packages("jsonlite")
# 
# # The follwing packages are for 
# install.packages('purrr')
# install.packages('tidyverse')
# install.packages('rjson')

# ############# [FACULTATIVE] Set Working directory ############
# setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

########### Load RinR library #############
library()

########### Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### REvaluate use #############
swagger_JSON <- REvaluate(
  {
    ############# [CODE] Load libraries ###################
    library(RCurl)
    library(jsonlite)
    
    # For map function
    library(purrr)
    library(tidyverse)
    library(rjson)
    
    ############# [CODE] Get the main Swagger.json from CELLBASE ###############
    swagger_URL <- getURL("http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/swagger.json")
    swagger_JSON <-fromJSON(swagger_URL)
    swagger_JSON
  }
  # ,
  # data = list(swagger_JSON = swagger_JSON)
  # ,
  # REvaluator = Rversion
  # ,
  # verbose	= TRUE
)

########### REvaluate use #############
tagsList <- REvaluate(
  {
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
    tags_list
  }
  ,
  data = list(swagger_JSON = swagger_JSON)
  # ,
  # REvaluator = Rversion
  # ,
  # verbose	= TRUE
)

pathsList <- REvaluate(
  {
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
    
    paths_list
  }
  ,
  data = list(swagger_JSON = swagger_JSON)
  # ,
  # REvaluator = Rversion
  # ,
  # verbose	= TRUE
)

pathsDF <- REvaluate(
  {
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
    
    paths_df
  }
  ,
  data = list(swagger_JSON = swagger_JSON)
  # ,
  # REvaluator = Rversion
  # ,
  # verbose	= TRUE
)

############# [INFO] [NOT STARTED] Ask the user to choose a particular GET (path) ###############
# This part will probably be done in Spotify (TIBCO)