############## [INFO] SCRIPT 2.2 General information ############
# This script will extract structural information from SWAGGER CELL BASE,
# in order to obtain parameters fields (filters) of a particular path (query).

# With this information, the goal is to create another script for TIBCO
# to ask the user which information he wants to retrieve from CELLBASE
# for a particular query.

############## [INFO] Input and Output ############################
# Input :
# 1) The path chosen by the user (the path would be provided by TIBCO)
# 2) We will work with the resource list, which is the main JSON of the CellBase
# "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/swagger.json

# Output : a table with the paramaters of a the chosen path

# ############# [INFO] Previous requirement (libraries) ###################
# install.packages("RCurl")
# install.packages("jsonlite")
# 
# # The follwing packages are for 
# install.packages('purrr')
# install.packages('tidyverse')
# install.packages('rjson')

# ############# [RStudio] Set Working directory ############
# Please uncoment the next line changin the working directory by the correct one:
# setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

########### [TIBCO] Load RinR library #############
library(RinR)

########### [TIBCO] Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### [TIBCO] Provide the chosen path ########
chosenPath <- SelectedPath

########### [TIBCO] Create the REvaluate object ########
param_mtrx <- REvaluate({
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
  
  ############# [Rstudio] Provide the chosen path ###############
  # If you are using RStudio, please provide the path
  # Below, we provide the path "variant annotation" as default.
  # Please uncoment the next line:
  # chosen_path <- '/{version}/{species}/genomic/variant/{variants}/annotation'
  
  ############# [CODE] Get the filter parameters from provided path ###############
  # The idea would be to reach the "parameters" level
  # to retrieve all the paramater for a particular given link
  
  # Create a list with the paths and their content:
  paths_JSON <- swagger_JSON$paths
  # paths_JSON # this line is for testing
  
  # Create a list with all the content of the chosen path
  var_annot_path <- paths_JSON[[paste(chosen_path)]]
  
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
  param_mtrx
}
,
# provide the path chosen by the user
data = list(chosen_path = chosenPath)
# ,
# REvaluator = Rversion
# ,
# verbose	= TRUE
)

paramTable <- param_mtrx

############# [Rstudio] Send the parameters table to a txt file ###############
write.table(param_mtrx,"test_files\\param_mtrx.txt",sep="\t",row.names=FALSE)