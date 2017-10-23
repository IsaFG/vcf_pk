# Example of ho to use REvaluate from RinR library
# This script will query a RESTfull database and retrive information
# This script works so far

########### Load RinR library #############
library(RinR)

########### Testing block #############
# The following block is for testing
variant_table <- variants_table.txt
TestingTable01 <- variant_table[c(1,4)]

########### Load URL #############
# JSON URL
jsonURL <- URL_string

########### Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### REvaluate use #############
preJSONtable <- REvaluate({
  library(RCurl)
  library(jsonlite)
  
  # Query the database
  query<-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
  querydf<-fromJSON(query)
  
  # #Get any warnings that may be worth noting
  # query_warnings<-querydf$warning
  # #Get any errors produced by the query
  # query_errors<-querydf$error
  
  # conseguir la respuesta a la query
  # nos da un query_result de 3 observaciones (los genes) de 9 variables
  query_result<-querydf$response
  
  # Crear una tabla con el numero de respuestas para cada una de las queries
  query_stats<-data.frame(query_result[c("id","numTotalResults","numResults")])
  names(query_stats)[1]<-"query"
  
  # Filtrar de la lista de resultados los outputs de aquellas queries que no devuelven resultado
  queries_with_results<- rowSums(query_stats[,c(2,3)])>0
  retrieved_results<-query_result$result[queries_with_results]
  david_test_table <-data.frame(retrieved_results[[1]]$id,retrieved_results[[1]]$name,retrieved_results[[1]]$description,retrieved_results[[1]]$drugInteractions)
  david_test_table
}
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)

########### Convert the REvaluate Result #############
JSONtable <- preJSONtable

########### Testing block #############
# This line is for testing is this table changes correctly
TestingTable02 <- variant_table[c(3,3,3,3)]