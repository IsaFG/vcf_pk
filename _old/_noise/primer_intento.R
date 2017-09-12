library(RCurl)
library(jsonlite)
library(curl)

# query<-getURL("bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/3:1000-200000,X:35-459000,4:2334-5555/pathway")

query<-getURL("ws.bioinfo.cipf.es/cellbase/rest/latest/hsa/feature/gene/BRCA2/transcript")

# text <- readLines(curl(query))

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
