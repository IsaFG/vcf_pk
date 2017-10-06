############## [INFO] SCRIPT 3 (TIBCO) General information ############
# This script is a draft for TIBCO
# Este script va a hacer uso :
# 1) De la tabla de variantes CromAndLocation_table.txt (generado por script VCF_get_variant)
# 2) De las tablas de parametros de filtrado param_mtrx.txt (generado por script Swagger_get_options)
# 3) Del GET de Varian Annotation
#   (idealmente debería cambiarse a un GET cualquiera elegido por el usuario)
#   GET /{version}/{species}/genomic/variant/{variants}/annotation

############## [INFO] Input and Output ############################
# Input :
# 1) tabla de variantes CromAndLocation_table.txt
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

############# [CODE] Load libraries ###################
library(RCurl)
library(jsonlite)

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

############# Parameters by default ###############
parameters_URL <- "limit=-1&skip=-1&skipCount=false&count=false&Output%20format=json&normalize=false&phased=false&useCache=false&imprecise=true&svExtraPadding=0&cnvExtraPadding=0"

############# Parameters chosen by user (TIBCO?) ###############
# I will probably use TIBCO for this
param_coded_table <- parameters_table

############# [Ony for GET VARIANT] Get the variants table ############
variant_table <- read.table("test_files\\CromAndLocation_table.txt", header=TRUE)
variant_table <- read.table("test_files\\variants_table.txt", header=TRUE)
var_number <- nrow(variant_table)

############# Build the GET URL ############

query_vector <-character()
URL_vector <- character()

for (i in 1:var_number) { 
  var_chrom <- variant_table[i,1]
  var_range <- variant_table[i,2]
  # ANIADIR AQUI EL CODIGO PARA EXTRAER LOS ALELOS DE REF Y ALT
  # colon_separator <- "%3A" # not necessary
  variant_URL <- paste(var_chrom,"%3A",var_range,"%3A",ref_al,"%3A",alt_al,sep = "")
  URL_base <- paste(common_URL,category_URL,variant_URL,subcategory_URL,"?",parameters_URL,sep = "")
  URL_vector <- c(URL_vector, URL_base)
  query_v <- getURL(URL_base)
  query_vector <- c(query_vector, query_v)
  print (paste("This is the i value: ", i)) # this line os for testing
}

# add the URL to the dataframe
variant_table$URL <- URL_vector

