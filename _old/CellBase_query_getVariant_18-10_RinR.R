############## [INFO] SCRIPT 3 (TIBCO) General information ############
# This script is a draft for TIBCO

############# [FACULTATIVE] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [FACULTATIVE] Install libraries ###################
# install.packages("RCurl")
# install.packages("jsonlite")

############# [FACULTATIVE] Install library cellbaseR ##########
# # From Bioconductor :
# source("https://bioconductor.org/biocLite.R")
# ## try http:// if https:// URLs are not supported
# biocLite("cellbaseR")
# 
# # From the develloper github (Mohammed Elsiddieg <melsiddieg@gmail.com>) :
# install.packages("devtools") # you may want to install devtools to install libraries from Github
# library(devtools)
# devtools::install_github("melsiddieg/cellbaseR")

############# HERE STARTS TIBCO SCRIPT ###########
############# PREVIOUS REQUIREMENT ###############
# LIBRARIES from CRAN : RCurl, jsonlite
# LIBRARY from Bioconductor : "cellbaseR":
# https://bioconductor.org/packages/release/bioc/html/cellbaseR.html

############## [INFO] Input and Output ############################
# Inputs :
# 1) A variant table (variants_table.txt) coming from a VCF file (script VCF_get_variant)
# 2) A filter parameters table (param_mtrx.txt) coming from script Swagger_get_options

# The builded query will be sent to CellBase
# Output :
# -> annotation table

############# Load RinR library ###################
library(RinR)

########### Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

############# [CODE] Load libraries ###################
library(RCurl)
library(jsonlite)
library(cellbaseR)

############# [CODE] Get the variants table ############
variant_table <- read.table("test_files\\variants_table.txt", header=TRUE) # for RStudio

############# [CellBaseR] Build the GET URL and query CellBase ############
preVariable <- REvaluate({
  var_number <- nrow(variant_table)
  query_vector <- character()
  URL_vector <- character()
  results_list <- list()
  cb <- CellBaseR()
  
  for (i in 1:var_number) { 
    var_chrom <- variant_table[i,1]
    var_range <- variant_table[i,2]
    var_refAl <- variant_table[i,3]
    var_altAl <- variant_table[i,4]
    
    # Get variant cellbase info with cellbaseR package
    variant <- paste(var_chrom, ":", var_range, ":", var_refAl, ":", var_altAl, sep = "")
    res2 <- getVariant(object=cb, ids=variant, resource="annotation")
    res2table <- res2[c("chromosome", "start", "reference", "alternate", "id", "displayConsequenceType")]
    if (i==1) {
      getVariant_table <- res2table
    } else {
      getVariant_table <- (rbind(as.matrix(getVariant_table), as.matrix(res2table)))
    }
    
    # THE FOLLOWING BLOCK DEPENDS ON ANOTHER SCRIPT
    # Get the Swagger link
    # variant_URL <- paste(var_chrom,"%3A",var_range,"%3A",var_refAl,"%3A",var_altAl,sep = "")
    # URL_base <- paste(common_URL,category_URL,variant_URL,subcategory_URL,"?",parameters_URL,sep = "")
    # URL_vector <- c(URL_vector, URL_base)
    # query_v <- getURL(URL_base)
    # query_vector <- c(query_vector, query_v)
    
    print (paste("Processing variant number:", i)) # this line is for testing
  }
  
},
data = list(variant_table = variants_table.txt) # not sure if it is correctly written
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)

# To print the table in a txt file
try(write.table(getVariant_table,"test_files\\CB_variants_table.txt", append = TRUE, sep="\t",row.names=FALSE))
