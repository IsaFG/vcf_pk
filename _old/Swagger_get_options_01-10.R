############## [INFO] SCRIPT 2 General information ############
# This script will extract structural information from SWAGGER CELL BASE,
# in order to obtain categories, link's sructure and parameters fields of the API.
# With this information, the goal is to create another script
# to ask the user which options he wants to retrieve from CELLBASE.

############## [INFO] Input and Output ############################
# Input : We will work with the resource list, which is the main JSON of the CellBase
# "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/swagger.json

# Output :
# 1) a txt file with all the main tags
# 2) a txt file with all the availables paths
# 3) a txt file with a matrix with the paramaters of a particular chosen path

############## [INFO] Path example used in this script ############################
# Retreive info from VARIANT ANNOTATION GET LINK
# This is the link :
# GET /{version}/{species}/genomic/variant/{variants}/annotation

# Implementation Notes :
# Include and exclude lists take values from the following set:
# {variation, clinical, conservation, functionalScore, consequenceType,
# expression, geneDisease, drugInteraction, populationFrequencies, repeats}

# Response Class (Status 200) :
# For informative purpose, all the responses sections of the example link
# are in the following json :
# variant_annotation_Response_Class.json

############# [FACULTATIVE] Install libraries ###################
install.packages("RCurl")
install.packages("jsonlite")

# The follwing packages are for 
install.packages('purrr')
install.packages('tidyverse')
install.packages('rjson')

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

sink("test_files\\tags_list.txt")
tags_list
sink()

# tags_list # this line is for testing

############# [CODE] Get the GET paths ###############
level1_keys <- map(swagger_JSON, names)
paths_list <- level1_keys$paths
# paths_list # this line is for testing

sink("test_files\\paths_list.txt")
paths_list
sink()

############# [INFO] [NOT STARTED] Ask the user to choose a particular GET (path) ###############
# This part will probably be done in Spotify (TIBCO)

############# [CODE] [IN PROGRESS] Get the filter parameters from a particular path ###############
# The idea would be to reach the "parameters" level
# to retrieve all the paramater for a particular given link
# In this DRAFT, we will use a specific path as an example.
# But the idea would be that the USER itself would choose the path he wants
# and then the script will return the availables parameters.

# Idealy, the path chosen by User would be stored in a variable :
# In this DRAFT, we will use the path "variant annotation" as an example.
chosen_path <- '/{version}/{species}/genomic/variant/{variants}/annotation'

# Create a list with the paths and theis content:
paths_JSON <- swagger_JSON$paths
# paths_JSON # this line is for testing

# Create a list with all the content of the choosen path (try 1)
# The following code does not work but why ?
chosen_path <-  paths_JSON$chosen_path
# chosen_path # this line is for testing

# Create a list with all the content of the choosen path (try 2)
# This code works but it is limited because the path is pre-chosen
var_annot_path <- paths_JSON$`/{version}/{species}/genomic/variant/{variants}/annotation`
# map(var_annot_path, names) # this line is for testing
parameters_list <- var_annot_path$get$parameters
# parameters_list # this line is for testing
param_list_length <- length(parameters_list)

param_name_vector <- character()
param_descr_vector <- character()

for (i in 1:param_list_length){
  one_param <- parameters_list[[i]]
  param_name <- one_param$name
  param_description <- one_param$description
  param_name_vector <- c(param_name_vector, param_name)
  param_descr_vector <- c(param_descr_vector, param_description)
}

param_mtrx <- cbind(param_name_vector, param_descr_vector)
write.table(param_mtrx,"test_files\\param_mtrx.txt",sep="\t",row.names=FALSE)
