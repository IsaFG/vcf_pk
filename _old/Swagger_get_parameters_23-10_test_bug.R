# ########### Load RinR library #############
# library(RinR)
# 
# chosenPath<-SelectedPath

# ########### Determinate R interpreter location ########
# Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
# 
# param_mtrx <- REvaluate({

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


parameters_list[[i]]

class(parameters_list)

for (i in 1:param_list_length){
  one_param <- parameters_list[[i]]
  param_name <- one_param$name
  param_description <- one_param$description
  param_name_vector <- c(param_name_vector, param_name)
  param_descr_vector <- c(param_descr_vector, param_description)
}

param_name_vector
param_descr_vector

param_mtrx <- cbind(param_name_vector, param_descr_vector)
param_mtrx

# }
# # ,
# # data = list(variable_to_use_in_REvaluate = variable_in_main_script)
# # ,
# # REvaluator = Rversion
# # ,
# # verbose	= TRUE
# )

paramTable <- param_mtrx