library(RCurl)
library(jsonlite)

query<-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
querydf<-fromJSON(query)

#Get any warnings that may be worth noting
query_warnings<-querydf$warning
#Get any errors produced by the query
query_errors<-querydf$error

query_result<-querydf$response

# Create a table with the number of responses from each of the queries:
query_stats<-data.frame(query_result[c("id","numTotalResults","numResults")])
names(query_stats)[1]<-"query"

# Filter from the result list the output of those queries that produced no result
queries_with_results<- rowSums(query_stats[,c(2,3)])>0
retrieved_results<-query_result$result[queries_with_results]

test<-data.frame(retrieved_results[[1]]$id,retrieved_results[[1]]$name,retrieved_results[[1]]$description,retrieved_results[[1]]$drugInteractions)

