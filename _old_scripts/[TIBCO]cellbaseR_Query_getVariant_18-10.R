############## [INFO] SCRIPT 3 (TIBCO) General information ############
# This script only works in TIBCO Spotfire
# We provide a table with some variants, then the script
# will query the RESTfull DB CellBase to retreive information about these variants

############# PREVIOUS REQUIREMENT ###############
# LIBRARIES from CRAN : RCurl, jsonlite
# LIBRARY to use this script in TIBCO : RinR
# LIBRARY from Bioconductor : "cellbaseR":
# https://bioconductor.org/packages/release/bioc/html/cellbaseR.html

############## [INFO] Input and Output ############################
# Inputs : A variant table coming from a VCF file (script VCF_get_variant)
# Output : annotated table

############# [FACULTATIVE RStudio] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

########### [RStudio] Get the variants table ############
variant_table <- read.table("test_files\\variants_table.txt", header=TRUE) # for RStudio

############# [TIBCO] Load RinR library ###################
library(RinR)

########### [TIBCO] Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

############# [CellBaseR] Build the GET URL and query CellBase ############
preVariable <- REvaluate({
  # Load libraries
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  library(dplyr)
  library(tidyr)
  
  var_number <- nrow(variant_table)
  query_vector <- character()
  URL_vector <- character()
  results_list <- list()
  cb <- CellBaseR()

  
  # test <- variant_table[5,]$Alt.allele
  # test2 <- variant_table[1,]$Alt.allele
  # test2
  # substring(test2, 1, 1)
  
  annotVariants_table <- data.frame()
  TESTannotVariants_table <- data.frame()
  for (i in 1:var_number) { 
    print (paste("Processing variant number:", i)) # this line is for testing
    # extract the chromosome
    var_chrom <- variant_table[i,1]
    
    # extract the range
    var_range <- variant_table[i,2]
    
    # extract the ref and alt alleles
    # WARNING:you could have more than one allele in each field
    # so that this formula extract only the fiest one
    var_refAl <- substring((variant_table[i,4]), 1, 1)
    var_altAl <- substring((variant_table[i,5]), 1, 1)
    
    # Get variant cellbase info with cellbaseR package
    variant <- paste(var_chrom, ":", var_range, ":", var_refAl, ":", var_altAl, sep = "")
    annotVariant <- getVariant(object=cb, ids=variant, resource="annotation")

    if (nrow(annotVariant)==0) {
      print (paste("WARNING : result call for variant number:", i, "is EMPTY"))
      annotVariant <- data.frame(var_chrom, as.integer(var_range), var_refAl, var_altAl)
      TESTcolumnas <- annotVariant
      colnames(annotVariant) <- c("chromosome", "start", "reference", "alternate")
      # annotVariants_table <- bind_rows(annotVariants_table, annotVariant)
      reducedTable <- annotVariant
      TESTannotVariants_table <- bind_rows(TESTannotVariants_table, reducedTable)
      
    } else {
      if (i==2) {
        # annotVariants_table <- bind_rows(annotVariants_table, annotVariant)
        reducedTable <- annotVariant[c("chromosome", "start", "reference", "alternate", "id", "displayConsequenceType", "displayConsequenceType", "geneDrugInteraction", "functionalScore")]
        TESTannotVariants_table <- bind_rows(TESTannotVariants_table, reducedTable)
        completeTable2 <- annotVariant
        reducedTable2 <- reducedTable
        variantTraitAssociation_df2 <- annotVariant$variantTraitAssociation
        # unnest_test2 <- unnest(reducedTable2) # does not work since the nested list are different
        
        
      } else if (i==13) {
        # annotVariants_table <- bind_rows(annotVariants_table, annotVariant)
        reducedTable13 <- ""
        completeTable13 <- annotVariant
        reducedTable <- annotVariant[c("chromosome", "start", "reference", "alternate", "displayConsequenceType", "displayConsequenceType", "geneDrugInteraction", "functionalScore")]
        reducedTable13 <- reducedTable
        TESTannotVariants_table <- bind_rows(TESTannotVariants_table, reducedTable)
        completeTable13
        #### end of testing block ####
        
      } else {
        # annotVariants_table <- bind_rows(annotVariants_table, annotVariant)
        completeTable_last <- annotVariant
        reducedTable <- annotVariant[c(1:9)]
        TESTannotVariants_table <- bind_rows(TESTannotVariants_table, reducedTable)
        
        # annotVariants_table <- (rbind(as.matrix(annotVariants_table), as.matrix(res2table)))
        # annotated_table <- (rbind(as.matrix(annotated_table), as.matrix(res2)))
      }
    } 
  }
  annotVariants_table
},
data = list(variant_table = variants_table.txt)
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)
getVariantTable <- preVariable

########### Testing block ############

completeTable2$variantTraitAssociation

str(completeTable13)
# Testing loop to analyze what there is inside the table 
j <- 0
for (column in completeTable13) {
  print (paste("Column", j, "has a class: ", class(column)))
  print (paste("Length of the column:", length(column)))
  if (class(column) == "data.frame") {
    print ("WARNING, THE FOLLWING COLUMN IS A DATAFRAME")
    print (paste("Column", j))
    print ("---------------------------")
    print (column)
  }
  print ("###################")
  j <- j + 1
}
column

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
annotVariants_table

########### [RStudio] Print the table in a txt file ###########
try(write.table(annotVariants_table,"test_files\\CB_variants_table.txt", append = TRUE, sep="\t",row.names=FALSE))
