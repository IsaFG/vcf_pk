library("RinR")
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
# variant <- testRinR # only for TIBCO

REvaluate({
  library(RCurl)
  library(jsonlite)
  "bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info"
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

}
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)

# to get the data
