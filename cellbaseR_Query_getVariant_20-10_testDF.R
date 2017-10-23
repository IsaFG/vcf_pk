############## [INFO] SCRIPT 3 General information ############
# NOTE : UNFINISHED WORK

# We provide a table with some variants, then the script will query the RESTfull DB CellBase
# to retreive information about these variants

############# PREVIOUS REQUIREMENT ###############
# LIBRARIES from CRAN : RCurl, jsonlite
# LIBRARY to use this script in TIBCO : RinR
# LIBRARY from Bioconductor : "cellbaseR":
# Bioconductor link for library "cellbaseR":
# https://bioconductor.org/packages/release/bioc/html/cellbaseR.html

############## [INFO] Input and Output ############################
# Inputs : A variant table (variants_table.txt) coming from a VCF file (script VCF_get_variant)
# Output : annotated table

# ############# [FACULTATIVE] Install libraries ###################
# install.packages("RCurl")
# install.packages("jsonlite")
# 
# ############# [FACULTATIVE] Install library cellbaseR ##########
# # From Bioconductor :
# source("https://bioconductor.org/biocLite.R")
# ## try http:// if https:// URLs are not supported
# biocLite("cellbaseR")
# 
# # From the develloper github (Mohammed Elsiddieg <melsiddieg@gmail.com>) :
# install.packages("devtools") # you may want to install devtools to install libraries from Github
# library(devtools)
# devtools::install_github("melsiddieg/cellbaseR")

############# [FACULTATIVE] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [CODE] Load libraries ###################
library(RCurl)
library(jsonlite)
library(cellbaseR)

############# [Ony for GET VARIANT QUERIES] Get the variants table ############
variant_table <- read.table("test_files\\variants_table.txt", header=TRUE)

############# [cellbaseR] Build the GET URL and query CellBase ############
query_vector <-character()
URL_vector <- character()
results_list <- list()
cb <- CellBaseR()
var_number <- nrow(variant_table)

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
    
    # testing lines
    completeTable1 <- res2

  } else if (i==2) {
    completeTable2 <- res2
    getVariant_table <- (rbind(as.matrix(getVariant_table), as.matrix(res2table)))
    
  } else if (i==3) {
    completeTable3 <- res2
    getVariant_table <- (rbind(as.matrix(getVariant_table), as.matrix(res2table)))
    
  } else if (i==4) {
    completeTable4 <- res2
    getVariant_table <- (rbind(as.matrix(getVariant_table), as.matrix(res2table)))
    
  } else {
    getVariant_table <- (rbind(as.matrix(getVariant_table), as.matrix(res2table)))
    
    # testing line
    completeTable5 <- res2
  }
  
  print (paste("Processing variant number:", i)) # this line is for testing
}

########### Print the results by console ############

########### Testing block ############
# Podemos comprobar que obtenemos tablas con numeros diferentes
c1 <- colnames(completeTable1)
c5 <- colnames(completeTable5)

completeTable1$geneDrugInteraction # class : list
completeTable1$hgvs # class : list
completeTable1$consequenceTypes # class: list
completeTable1$consequenceTypes[[1]] # class: data.frame

c1
c5
# La linea siguiente muestra que campos son diferentes entre ambas tablas
setdiff(c5, c1)

# La linea siguiente muestra los campos comunes entre ambas tablas
intersect(c1,c5)

class(res2)
class(res2table)
getVariant_table

########### Print the table in a txt file ###########
try(write.table(getVariant_table,"test_files\\CB_variants_table.txt", append = TRUE, sep="\t",row.names=FALSE))
