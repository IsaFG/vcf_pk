# EJEMPLO 1 DEL API
# ESTE EJEMPLO FUNCIONA !!!!

library(RCurl)
library(jsonlite)

# guardar el enlace en el objeto query. Para eso se usa la funcion geturl de Rcurl, que descarga URLs
query<-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
# this REST call will get all the INFO of the gene BRCA2,BRCA1,KRAS of human in the latest version.
# esta direccion se puede pegar en el navegador para ver la pinta que tiene
query

# almacenamos el contenido json de la url en un objeto
# este objeto contiene pues toda la informacion de los genes solicitados
querydf<-fromJSON(query)

## NOTA : warning, error y response estan en la primera linea de la URL
## NOTA : el contenido response es todo el resto
#Get any warnings that may be worth noting
query_warnings<-querydf$warning
#Get any errors produced by the query
query_errors<-querydf$error
# conseguir la respuesta a la query
query_result<-querydf$response

# nos da un query_result de 3 observaciones (los genes) de 9 variables

# query_result #esto cuelga R

# Create a table with the number of responses from each of the queries:
# crear una tabla con el numero de respuestas para cada una de las queries
query_stats<-data.frame(query_result[c("id","numTotalResults","numResults")])
names(query_stats)[1]<-"query"
query_stats


# Filter from the result list the output of those queries that produced no result
# filtrar de la lista de resultados los outputs de aquellas queries que no devuelven resultado
queries_with_results<- rowSums(query_stats[,c(2,3)])>0
queries_with_results
retrieved_results<-query_result$result[queries_with_results]

test<-data.frame(retrieved_results[[1]]$id,retrieved_results[[1]]$name,retrieved_results[[1]]$description,retrieved_results[[1]]$drugInteractions)
test
