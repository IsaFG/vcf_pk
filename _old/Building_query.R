############## [INFO] Building query ############
# NOTE: UNFINISHED SCRIPT

############## [INFO] Input and Output ############################
# Inputs :
# 1) A variant table (variants_table.txt) coming from a VCF file (script VCF_get_variant)
# 2) A filter parameters table (param_mtrx.txt) coming from script Swagger_get_options

# The user will have to choice filter parameters so that the QUERY can be build
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

############# Set the different parts of the GET URL ############
# Ideally category and subcategory would be chosen by the user
# NOTE : maybe the common URL could suffer some changes
common_URL <- "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/"
category_URL <- "genomic/variant/"
subcategory_URL <- "/annotation"
# param_start_URL <- "?" # probably unnecesary
# param_separator_URL <- "&" # probably unnecesary

############# Get the parameters table ############
parameters_table <- read.table("test_files\\param_mtrx.txt", header=TRUE)
param_number <- nrow(parameters_table)

############# [FACULTATIVE] Parameters by default ###############
parameters_URL <- "limit=-1&skip=-1&skipCount=false&count=false&Output%20format=json&normalize=false&phased=false&useCache=false&imprecise=true&svExtraPadding=0&cnvExtraPadding=0"

############# [UNFINISHED] BUILDING QUERY: Parameters chosen by user (TIBCO?) ###############
## NOT IN USE YET
## I will probably use TIBCO for this
# param_coded_table <- parameters_table

############# [Ony for GET VARIANT QUERIES] Get the variants table ############
variant_table <- read.table("test_files\\variants_table.txt", header=TRUE)

############# Build the GET URL ############
# query_vector <-character()
URL_vector <- character()
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
