############## [INFO] Building query ############
# NOTE: UNFINISHED SCRIPT
# NOTE : this script has been written only for the GET Annotated variant
# This script depends on several other scripts
# which will provide the path and the parameters list

############## [INFO] Input and Output ############################
# Input 1 : path chosen by user
## For testing, we will user the path Annotated variant

# Input 2 : Parameters chosen by user
## For testing, we will use the table param_mtrx.txt
## coming from script Swagger_get_parameters

# Input 3 : The variable(s) to be queried (variant, gene, etc...)
## For testing, we will user a variant table (variants_table.txt)
## coming from a VCF file (script VCF_get_variant)

# Output : builded QUERY

############# [FACULTATIVE] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [FACULTATIVE] Install libraries ###################
install.packages("RCurl")
install.packages("jsonlite")

# ############# [FACULTATIVE] Install library cellbaseR ##########
# # From Bioconductor :
# source("https://bioconductor.org/biocLite.R")
# ## try http:// if https:// URLs are not supported
# biocLite("cellbaseR")
# 
# # From the develloper github (Mohammed Elsiddieg <melsiddieg@gmail.com>) :
# install.packages("devtools") # you may want to install devtools to install libraries from Github
# library(devtools)
# devtools::install_github("melsiddieg/cellbaseR")

############# [CODE] Load libraries ###################
library(RCurl)
library(jsonlite)
# library(cellbaseR)

############# [CODE] Get the main Swagger.json from CELLBASE ###############
swagger_URL <- getURL("http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/swagger.json")
swagger_JSON <-fromJSON(swagger_URL)

############# [TIBCO] Let the user choose the path he wants ##########
# The USER itself would choose the path he wants
# and then the script will return the availables parameters.
# The path chosen by User would be stored in a variable.

chosen_path <- SelectedPath

# In this DRAFT, we will use a specific path as an example.
chosen_path <- '/{version}/{species}/genomic/variant/{variants}/annotation'

# Another way to build the URL would be divide the path in different parts ?
category_URL <- "genomic/variant/"
subcategory_URL <- "/annotation"

############# [PENDING] BUILDING QUERY with Parameters chosen by user (TIBCO) ###############
## NOT IN USE YET
## The user will choose the parameters he want with TIBCO
# chosen_parameters # table with parameters chose by user
## Process the parameters to build the parameters_URL
parameters_URL <- ""

############# [FACULTATIVE] Parameters by default ###############
parameters_URL <- "limit=-1&skip=-1&skipCount=false&count=false&Output%20format=json&normalize=false&phased=false&useCache=false&imprecise=true&svExtraPadding=0&cnvExtraPadding=0"

############# BUILDING QUERY : common part of the GET URL ############
# Part of the get url which is common
# NOTE : maybe the common URL could suffer some changes
common_URL <- "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/"
URL_vector <- character()

############# BUILDING QUERY for getVariant [Ony for GET VARIANT QUERIES] ############
# If the user has chosen the path Annotated variant
variant_table <- read.table("test_files\\variants_table.txt", header=TRUE)
# results_list <- list()
var_number <- nrow(variant_table)

for (i in 1:var_number) { 
  var_chrom <- variant_table[i,1]
  var_range <- variant_table[i,2]
  var_refAl <- variant_table[i,3]
  var_altAl <- variant_table[i,4]
  
  # Get the Swagger link
  variant_URL <- paste(var_chrom,"%3A",var_range,"%3A",var_refAl,"%3A",var_altAl,sep = "")
  URL_base <- paste(common_URL,category_URL,variant_URL,subcategory_URL,"?",parameters_URL,sep = "")
  URL_vector <- c(URL_vector, URL_base)

  # # Testing block
  # print (paste("Processing variant number:", i)) # this line is for testing
  # print (URL_vector)
}

# You get a vector with the getVariant URL for the 5 variants
URL_vector
