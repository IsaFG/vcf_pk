############## [INFO] SCRIPT 3 (TIBCO) General information ############
# This script is a draft for TIBCO
# IMPORTANT : this script will use the library "cellbaseR" 

# Bioconductor link : https://bioconductor.org/packages/release/bioc/html/cellbaseR.html
# To avoid this package, please refer to a previous version of this script (possibly unfinished)

# Este script va a hacer uso :
# 1) De la tabla de variantes CromAndLocation_table.txt (generado por script VCF_get_variant)
# 2) De las tablas de parametros de filtrado param_mtrx.txt (generado por script Swagger_get_options)
# 3) Del GET de Varian Annotation
#   (idealmente debería cambiarse a un GET cualquiera elegido por el usuario)
#   GET /{version}/{species}/genomic/variant/{variants}/annotation

############## [INFO] Input and Output ############################
# Input :
# 1) A variant table coming from a VCF file previously processed
# 2) tablas de parametros de filtrado param_mtrx.txt
# 2) pide al usuario elegir paramtetros de filtrado
# Codigo : con el resultado construye el enlace GET Variant annotation completo
# Output 1 : el enlace de GET 
# Codigo : manda el enlace a CellBase (otro script?)
# Output 2 : recibe la tabla de anotaciones (otro script?)

############# [FACULTATIVE] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [FACULTATIVE] Install libraries ###################
install.packages("RCurl")
install.packages("jsonlite")

############# Install library cellbaseR ##########
# From Bioconductor :
source("https://bioconductor.org/biocLite.R")
## try http:// if https:// URLs are not supported
biocLite("cellbaseR")

# From the develloper github (Mohammed Elsiddieg <melsiddieg@gmail.com>) :
install.packages("devtools") # you may want to install devtools to install libraries from Github
library(devtools)
devtools::install_github("melsiddieg/cellbaseR")

############# [CODE] Load libraries ###################
library(RCurl)
library(jsonlite)
library(cellbaseR)

############# BUILDING QUERY: Set the different parts of the GET URL ############
# Ideally category and subcategory would be chosen by the user
# NOTE : maybe the common URL could suffer some changes
common_URL <- "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/"
category_URL <- "genomic/variant/"
subcategory_URL <- "/annotation"
# param_start_URL <- "?" # probably unnecesary
# param_separator_URL <- "&" # probably unnecesary

############# BUILDING QUERY: Get the parameters table ############
parameters_table <- read.table("test_files\\param_mtrx.txt", header=TRUE)
param_number <- nrow(parameters_table)

############# BUILDING QUERY: Parameters by default ###############
parameters_URL <- "limit=-1&skip=-1&skipCount=false&count=false&Output%20format=json&normalize=false&phased=false&useCache=false&imprecise=true&svExtraPadding=0&cnvExtraPadding=0"

############# BUILDING QUERY: Parameters chosen by user (TIBCO?) ###############
# I will probably use TIBCO for this
param_coded_table <- parameters_table

############# [Ony for GET VARIANT] Get the variants table ############
variant_table <- read.table("test_files\\CromAndLocation_table.txt", header=TRUE)
variant_table <- read.table("test_files\\variants_table.txt", header=TRUE)
var_number <- nrow(variant_table)

############# Build the GET URL (cellbaseR use) ############
query_vector <-character()
URL_vector <- character()
results_list <- list()
cb <- CellBaseR()

for (i in 1:var_number) { 
  var_chrom <- variant_table[i,1]
  var_range <- variant_table[i,2]
  var_refAl <- variant_table[i,3]
  var_altAl <- variant_table[i,4]
  variant <- paste(var_chrom, ":", var_range, ":", var_refAl, ":", var_altAl, sep = "")
  # ANIADIR AQUI EL CODIGO PARA EXTRAER LOS ALELOS DE REF Y ALT
  # colon_separator <- "%3A" # not necessary
  # variant_URL <- paste(var_chrom,"%3A",var_range,"%3A",var_refAl,"%3A",var_altAl,sep = "")
  # URL_base <- paste(common_URL,category_URL,variant_URL,subcategory_URL,"?",parameters_URL,sep = "")
  # URL_vector <- c(URL_vector, URL_base)
  # query_v <- getURL(URL_base)
  # query_vector <- c(query_vector, query_v)
  try(res2 <- getVariant(object=cb, ids=variant, resource="annotation"))
  # to get the data 
  try(res2 <- cbData(res2))
  
  if (i==1) {
    testvar1 <- res2
  } else if ( i ==2) {
    testvar2 <- res2
  } else if ( i ==3) {
    testvar3 <- res2
  # } else {
  #   getVariant_table <- data.frame(rbind(as.matrix(getVariant_table), as.matrix(res2)))
  }
  results_list <- append(results_list, res2)
  try(write.table(res2,"test_files\\CB_variants_table.txt", append = TRUE, sep="\t",row.names=FALSE))
  print (paste("Processing variant number:", i)) # this line is for testing
  print(" ") # this line is for testing
}

# Testing block:
table1 <- testvar1[c("chromosome", "start", "reference", "alternate", "id")]
table2 <- testvar2[c("chromosome", "start", "reference", "alternate", "id")]
table3 <- testvar3[c("chromosome", "start", "reference", "alternate", "id")]

# getVariant_table <- data.frame(rbind(as.matrix(table1), as.matrix(table2)))
getVariant_table <- (rbind(as.matrix(table1), as.matrix(table2)))
getVariant_table <- (rbind(as.matrix(getVariant_table), as.matrix(table3)))

# add the URL to the dataframe
variant_table$URL <- URL_vector

query_vector[1]
