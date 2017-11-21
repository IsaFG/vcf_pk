############## [INFO] SCRIPT 2.1 General information ############
# This script will extract structural information from SWAGGER CELL BASE,
# in order to obtain categories and link's structure (main tags and paths) of the API.
# With this information, the goal is to create another script for TIBCO
# to ask the user which options he wants to use from CELLBASE.

############## [INFO] Input and Output ############################
# Input : We will work with the resource list, which is the main JSON of the CellBase
# "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/swagger.json

# Output :
# 1) a txt file with all the main tags
# 2) a txt file with all the availables paths
# (PROBLEM when passing the data to a txt file that will be used in TIBCO)

# ############# [FACULTATIVE] Install libraries ###################
# install.packages("RCurl")
# install.packages("jsonlite")
# 
# # The follwing packages are for 
# install.packages('purrr')
# install.packages('tidyverse')
# install.packages('rjson')

############# [FACULTATIVE] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

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

# store the tags list in a txt file
sink("test_files\\tags_list.txt")
tags_list
sink()

# tags_list # this line is for testing

############# [CODE] Get the GET paths ###############
level1_keys <- map(swagger_JSON, names)

# Create a vector with all the paths
# This creates problems when you pass it to a txt file in order to use it in TIBCO
paths_list <- level1_keys$paths

# Create a dataframe from the vector
# This creates problems when you pass it to a txt file in order to use it in TIBCO
paths_df <- data.frame(level1_keys$paths)

sink("test_files\\paths_df.txt")
paths_df
sink()

sink("test_files\\paths_list.txt")
paths_list
sink()

############# [INFO] [NOT STARTED] Ask the user to choose a particular GET (path) ###############
# This part will probably be done in Spotify (TIBCO)