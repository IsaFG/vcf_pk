############## [INFO] General information of the SCRIPT ############
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

############# Set the category part of the GET URL ############
# Ideally this part would be chosen by the user
category_URL <- "bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/variant/"
subcategory_URL <- "/annotation"

############# Insert the variant in the GET URL ############
variant_table <- read.table("test_files\\CromAndLocation_table.txt", header=TRUE)
chroms <- variant_table[1] # could be unnecesary 
ranges <- variant_table[2] # could be unnecesary 
var_number <- nrow(variant_table)

for (i in nrow) { # OJO ESTO ES INVALIDO > corregir
  var_chrom <- variant_table[i,1]
  var_range <- variant_table[i,2]
  colon_separator <- ":"
  URL_base <- paste(category_URL,var_chrom,colon_separator,var_range,subcategory_URL, sep = "")
  query_v <- getURL(URL_base)
  print(query_v)
  # TEST NOTE : in the next line the URL is replaced by some URL that work better
  # query_v <-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/1:18149476-18149476/gene")
  # query_v <-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
  return (query_v)
}

# add a new column...
variant_table$geturl <- c()
