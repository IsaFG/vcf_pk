################# SCRIPT 2 DESCRIPTION #############################
# PROJECT : Create an R script that queries Cellbase
## in order to retrieve the available annotation information.
# INPUTS :
## - a cromosome and a range
## - the information type that the user wants to retrieve
# OUTPUT :
## the information which corresponds to the cromosome-range submited and the chosen info type

################# SOME NOTES #############################
# Example of a query to cellbase :
## query<-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
## This REST call will get all the INFO of the genes
## BRCA2,BRCA1,KRAS of human in the latest version.

# EXAMPLES USED TO TEST THE SCRIPT :
## Input example 1 : does not produce any information
#### first variant found in the VCF file used as example
#### -> chr1 - 18149476 - Ref:G and Alt:A

## Input example 2 : breaks at the last step :
#### "bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/1:18149476-18149476/gene")
#### Category : genomics - subcategory region
#### The information of this URL type is here :
#### http://docs.bioinfo.cipf.es/projects/cellbase/wiki/Genomic_rest_ws_api#Region

## Input example 3 : works until the last step
#### ("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")

## Input example 4 :- DOES NOT WORK SO FAR
#### URL THAT RETREIVE VARIANT INFO 
#### para variantes aun estoy buscando un ejemplo valido
#### http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/variant/1:18149476:G/snp_phenotype
#### Guardar el enlace en el objeto query. Para eso se usa la funcion geturl de Rcurl, que descarga URLs
#### query_v <-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/variant/10:52575931:G/mutation_phenotype")

################# INSTALL LIBRARIES #############################
## OPTIONAL : INSTALL THE LIBRARIES NEEDED FOR THE SCRIPT
# install.packages("RCurl")
# install.packages("jsonlite")

################# LOAD LIBRARIES #############################
# Load libraries
library(RCurl)
library(jsonlite)

################# METHODS AND FUNCTIONS ############################# 

# METHOD get all the categories from CellBase
categories_URL <- "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/#/"
categories <- getURL(categories_URL)
fileConn<-file("CB_categories.txt")
writeLines(categories, fileConn)
close(fileConn)


# METHOD 1 : insert crom and range in the URL
# NOTE : the crom and the range will come from another file,
# but in this version, they are preconfigurate
# So the variables of the function have to be defined later by other funtions
get_RESTURL <- function() {
  variant_crom <- "1"
  variant_range <- "18149476"
  colon_separator <- ":"
  category_URL <- "bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/"
  subcategory_URL <- "/gene"
  URL_base <- paste(category_URL,variant_crom,colon_separator,variant_range,subcategory_URL, sep = "")
  query_v <- getURL(URL_base)
  # TEST NOTE : in the next line the URL is replaced by some URL that w¡work better
  query_v <-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/1:18149476-18149476/gene")
  query_v <-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
  query_v <-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2/info")
  return (query_v)}

# METHOD 2 : convert URL to a JSON object
get_json <- function(restfull_URL) {
  # pass the URL content to a JSON
  querydf_v<-fromJSON(restfull_URL)
  return (querydf_v)}

# METHOD 3 : get the results from the query
get_results <- function(cellbase_json){
  # Get any warnings that may be worth noting
  query_warnings_v<-cellbase_json$warning # NOT USED SO FAR
  # Get any errors produced by the query
  query_errors_v<-cellbase_json$error # NOT USED SO FAR
  # get the results
  query_result_v<-cellbase_json$response
  return(query_result_v)
}

# METHOD : Create a table with the number of responses from each of the queries.
get_num_responses <- function(query_results){
  query_stats_v<-data.frame(query_results[c("id","numTotalResults","numResults")])
  names(query_stats_v)[1]<-"query"
  return(query_stats_v)
}

# METHOD : produce the data table and send it to a txt file
# get_results_table <- function(num_responses, query_results) {
#   # Filter from the result list the output of those queries that produced no result
#   queries_with_results_v<- rowSums(num_responses[,c(2,3)])>0
#   write.table(queries_with_results_v,"4.queries_with_results_GENE_1.txt",sep="\t",row.names=FALSE)
#   # get the results
#   retrieved_results_v<-query_results$result[queries_with_results_v]
#   # write.table(retrieved_results_v,"3.retrieved_results.txt",sep="\t",row.names=FALSE)
#   # get the table
#   results_table <-data.frame(retrieved_results_v[[1]]$id,retrieved_results_v[[1]]$name,retrieved_results_v[[1]]$description,retrieved_results_v[[1]]$drugInteractions)
#   # TEST PRINT IN TXT
#   write.table(results_table,"5.RESTFUL_call_GENE_1.txt",sep="\t",row.names=FALSE)
#   return (results_table)
# }

################# CALLS TO FUNCTIONS #############################
# create and obtain the URL
restfull_URL <- get_RESTURL()
# Generate the JSON from the URL
cellbase_json <- get_json(restfull_URL)
# get resuts from JSON
query_results <- get_results(cellbase_json)
# get the table of response
num_responses <- get_num_responses(query_results)

# get the data table
# results_dataframe <- get_results_table(num_responses, query_results)

################# CODE TO TEST ####################################
# Filter from the result list the output of those queries that produced no result
queries_with_results_v<- rowSums(num_responses[,c(2,3)])>0
write.table(queries_with_results_v,"4.queries_with_results_GENE_1.txt",sep="\t",row.names=FALSE)
# get the results
retrieved_results_v<-query_results$result[queries_with_results_v]
# write.table(retrieved_results_v,"3.retrieved_results.txt",sep="\t",row.names=FALSE)
# get the table
results_table <-data.frame(retrieved_results_v[[1]]$id,retrieved_results_v[[1]]$name,retrieved_results_v[[1]]$description,retrieved_results_v[[1]]$drugInteractions)
# TEST PRINT IN TXT
write.table(results_table,"5.RESTFUL_call_GENE_1.txt",sep="\t",row.names=FALSE)
# return (results_table)

################# SOME TEST PRINTS #############################
# TEST PRINT IN TXT
write.table(num_responses,"3.num_of_responses_GENE_1.txt",sep="\t",row.names=FALSE)
# TEST PRINT IN TXT
write.table(query_results,"2.query_results_GENE_1.txt",sep="\t",row.names=FALSE)
# ESTE DA PROBLEMAS, PONERLO AL FINAL
write.table(cellbase_json,"1.query_from_json_GENE_1.txt",sep="\t",row.names=FALSE)

