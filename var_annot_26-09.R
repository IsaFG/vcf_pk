############## General information of the SCRIPT ############
# Get annotation of a variant found in a VCF file
# Annotation information is extracted from CellBase

############## Example used in this script ############################
#### First variant of the VCF file
# R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf
#### -> chr1 - 18149476 - Ref:G and Alt:A

# The link that will be created will be this one:
# http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/1%3A18149476%3AG%3AA/annotation?limit=-1&skip=-1&count=false&Output%20format=json&normalize=false
# short version:
# http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/1%3A18149476%3AG%3AA/annotation

############# Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# Load libraries ###################
library(RCurl)
library(jsonlite)

############# Get the json from cellbase #######
query<-getURL("http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/1%3A18149476%3AG%3AA/annotation")
# esta direccion se puede pegar en el navegador para ver la pinta que tiene
# en esta query tenemos lo que es el JSON. Tipo : Large character

# Almacenamos el contenido json de la url en un objeto tipo lista (este caso de cinco)
# este objeto contiene pues toda la informacion de los genes solicitados
querydf<-fromJSON(query)
# querydf # CUIDADO esto cuelga R con datos masivos

########### Extraer la respuesta del Json ############################
## NOTA : warning, error y response estan en la primera linea de la URL
## NOTA : el contenido response es todo el resto

#Get any warnings that may be worth noting
query_warnings<-querydf$warning
#Get any errors produced by the query
query_errors<-querydf$error

# conseguir la respuesta a la query
# nos da un query_result de 1 observaciones de 8 variables
query_result<-querydf$response
# query_result  # CUIDADO esto cuelga R con datos masivos

############# Extraer el numero de respuestas del Json ##############
# Create a table with the number of responses from each of the queries:
# crear una tabla con el numero de respuestas para cada una de las queries
query_stats<-data.frame(query_result[c("id","numTotalResults","numResults")])
names(query_stats)[1]<-"query"
# query_stats

############# Construir una tabla con todos los resultados #################
# Filter from the result list the output of those queries that produced no result
# filtrar de la lista de resultados los outputs de aquellas queries que no devuelven resultado
queries_with_results<- rowSums(query_stats[,c(2,3)])>0
# queries_with_results
retrieved_results<-query_result$result[queries_with_results]
# retrieved_results  # CUIDADO esto cuelga R con datos masivos

############# Construir una tabla con las anotaciones ConsecuenceTypes #################
consequenceTypes_list <- data.frame(retrieved_results[[1]]$consequenceTypes)
consequenceTypes_list
sink("test_files\\var_annot_consequenceTypes.txt")
consequenceTypes_list
sink()

clinicalSignificance_list <- data.frame(retrieved_results[[1]]$clinicalSignificance)
clinicalSignificance_list
sink("test_files\\var_annot_clinicalSignificance.txt")
clinicalSignificance
sink()

############# Extraer el nombre de gen #################
geneName_o <- data.frame(consequenceTypes_list[1])
gene_found <- geneName_o[1,]
# gene_found

