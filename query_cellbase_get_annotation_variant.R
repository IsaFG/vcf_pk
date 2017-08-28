library(RCurl)
library(jsonlite)

# guardar el enlace en el objeto query. Para eso se usa la funcion geturl de Rcurl, que descarga URLs
# query<-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
# this REST call will get all the INFO of the gene BRCA2,BRCA1,KRAS of human in the latest version.

# ESTRUCTURA A RESPETAR :
# http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/13:32972105-32973105/gene
# EJEMPLO QUE USAREMOS :
# la primera variante del fichero que voy a usar como prueba -> chr1	18149476 Ref:G y Alt:A

# Empezemos por anotar el posible gen localizado en la region de interés
# DA RESULTADO
query_v <-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/1:18149476-18149476/gene")

# ejemplo con variantes
# para variantes aun estoy buscando un ejemplovalido
# http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/variant/1:18149476:G/snp_phenotype
# query_v <-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/variant/10:52575931:G/mutation_phenotype")

# almacenamos el contenido json de la url en un objeto
querydf_v<-fromJSON(query_v)

#Get any warnings that may be worth noting
query_warnings_v<-querydf_v$warning
#Get any errors produced by the query
query_errors_v<-querydf_v$error

query_result_v<-querydf_v$response

# Create a table with the number of responses from each of the queries:
query_stats_v<-data.frame(query_result_v[c("id","numTotalResults","numResults")])
names(query_stats_v)[1]<-"query"

query_stats_v

# Filter from the result list the output of those queries that produced no result
queries_with_results_v<- rowSums(query_stats_v[,c(2,3)])>0
queries_with_results_v
retrieved_results_v<-query_result_v$result[queries_with_results_v]

test_v<-data.frame(retrieved_results_v[[1]]$id,retrieved_results_v[[1]]$name,retrieved_results_v[[1]]$description,retrieved_results_v[[1]]$drugInteractions)
test_v
