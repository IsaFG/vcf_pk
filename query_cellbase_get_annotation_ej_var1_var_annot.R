# EJEMPLO 1 DEL API
# Ejemplo usando la primera variante
# que encontramos en el VCF file usado como ejemplo
#### first variant found in the VCF file used as example
#### -> chr1 - 18149476 - Ref:G and Alt:A

# COn este link obtenemos anotación para la variante seleccionada
# "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/1%3A18149476%3AG%3AA/annotation?limit=-1&skip=-1&count=false&Output%20format=json&normalize=false
# version corta :
# http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/1%3A18149476%3AG%3AA/annotation

library(RCurl)
library(jsonlite)
# guardar el enlace en el objeto query. Para eso se usa la funcion geturl de Rcurl, que descarga URLs
query<-getURL("http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/1%3A18149476%3AG%3AA/annotation")
# esta direccion se puede pegar en el navegador para ver la pinta que tiene

# TEST CON OTRA URL
# query <-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/1:18149476-18149476/gene")
# en esta query tenemos lo que es el JSON. Tipo : Large character
query # podemos mostrar el JSON en pantalla
sink("query.txt")
query
sink()



# almacenamos el contenido json de la url en un objeto tipo lista
# este objeto contiene pues toda la informacion de los genes solicitados
querydf<-fromJSON(query)
querydf # CUIDADO esto cuelga R con datos masivos
sink("querydf.txt")
querydf
sink()

# Convertir un JSON a dataframe
asFrame <- do.call("rbind.fill", lapply(querydf, as.data.frame))

asDataFrame <- fromJSON(query) %>% as.data.frame

## NOTA : warning, error y response estan en la primera linea de la URL
## NOTA : el contenido response es todo el resto
#Get any warnings that may be worth noting
query_warnings<-querydf$warning
#Get any errors produced by the query
query_errors<-querydf$error

# conseguir la respuesta a la query
# nos da un query_result de 3 observaciones (los genes) de 9 variables
query_result<-querydf$response
query_result  # CUIDADO esto cuelga R con datos masivos
sink("query_result.txt")
query_result
sink()

# Create a table with the number of responses from each of the queries:
# crear una tabla con el numero de respuestas para cada una de las queries
query_stats<-data.frame(query_result[c("id","numTotalResults","numResults")])
names(query_stats)[1]<-"query"
query_stats
sink("query_stats.txt")
query_stats
sink()


# Filter from the result list the output of those queries that produced no result
# filtrar de la lista de resultados los outputs de aquellas queries que no devuelven resultado
queries_with_results<- rowSums(query_stats[,c(2,3)])>0
queries_with_results
retrieved_results<-query_result$result[queries_with_results]
retrieved_results  # CUIDADO esto cuelga R con datos masivos
sink("retrieved_results.txt")
retrieved_results
sink()

class(retrieved_results)
results_DF <- do.call(rbind.data.frame, retrieved_results)
sink("results_DF")
results_DF
sink()
class(results_DF)

# test<-data.frame(retrieved_results[[1]]$id,retrieved_results[[1]]$name,retrieved_results[[1]]$description,retrieved_results[[1]]$drugInteractions)
test<-data.frame(retrieved_results[[1]]$id,retrieved_results[[1]]$name,retrieved_results[[1]]$description)
test
sink("tessst.txt")
test
sink()
