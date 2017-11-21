############## [INFO] SCRIPT 3 General information ############
# We provide a table with some variants, then the script
# will query the RESTfull DB CellBase to retreive CADD annotations about these variants

############## [INFO] Input and Output ############################
# Inputs : A variant table coming from a VCF file (script VCF_get_variant)
# Output : CADD annotated table

############# [INFO] Previous requirement (libraries) ###############
# LIBRARIES from CRAN : RCurl, jsonlite
# LIBRARY to use this script in TIBCO : RinR
# LIBRARY from Bioconductor : "cellbaseR":
# https://bioconductor.org/packages/release/bioc/html/cellbaseR.html

############## [INFO] Problems ############################
# Problem 1:  Some variants have more than 1 REF or ALT allele
# so that this script do the annotation call with the first allele of each field.
# I will change that just deleting these variants

# Problem 2: warning message : 
# In bind_rows_(x, .id) :
#   binding character and factor vector, coercing into character vector

########### [Code] Main method ########
getCaddScoresTable <- function(variants_table) {
  ############# [Code] Load libraries ############# 
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  library(dplyr)
  library(tidyr)
  
  ############# [Code] Build the GET URL and query CellBase (CellBaseR) ############
  var_number <- nrow(variants_table)
  cb <- CellBaseR()
  
  # Initialize the annotation table
  caddScores_table <- data.frame()
  
  i <- 1

  for (i in 1:var_number) { 
    print (paste("Processing variant number:", i)) # this line is for testing
    
    # Extract the first columns from the variant table, in order to build the CADD annotation table
    table_first_columns <- variants_table[i,-3]

    # extract the chromosome
    var_chrom <- variants_table[i,1]
    
    # extract the range
    var_range <- variants_table[i,2]
    
    # extract the ref and alt alleles
    # WARNING: you could have more than one allele in each field
    # so that this formula extract only the first one to do the annotation call
    # IS THAT A GOOD APPROACH ?
    var_refAl <- substring((variants_table[i,4]), 1, 1)
    var_altAl <- substring((variants_table[i,5]), 1, 1)
    
    # Get variant cellbase info with cellbaseR package
    # the call will return as a data.frame
    variant <- paste(var_chrom, ":", var_range, ":", var_refAl, ":", var_altAl, sep = "")
    caddScores_row <- getCaddScores(object=cb, id=variant, param = NULL)
    
    # Edit the table to make it readable (transpose, colname and remove unusefull column)
    caddScores_row <- as.data.frame(t(caddScores_row[,2:1]))
    caddScores_row <- caddScores_row[2,]

    annotVariant <- cbind(table_first_columns, caddScores_row)
    caddScores_table <- bind_rows(caddScores_table, annotVariant)
    # warnings()
  }
  # Change the name of the "id" column to avoid repetition
  colnames(caddScores_table) <- (c("chromosome", "start", "reference", "alternate", "cadd_raw_score", "cadd_scaled_score"))
  caddScores_table[,5:6] <- sapply(caddScores_table[,5:6], as.numeric)
  return(caddScores_table)
}

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  ########### [TIBCO] Load RinR library ###################
  library(RinR)
  
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
  
  ########### [TIBCO] Create the REvaluate object to execute main method ########
  caddScoresTable <- REvaluate({
    caddScores_table <- getCaddScoresTable(variants_table)
    caddScores_table
  }
  , data = list(getCaddScoresTable = getCaddScoresTable, variants_table = variants_table)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Get the variants table ############
  variants_table <- read.table("test_files\\variants_table.txt", header=TRUE) # for RStudio
  
  ########### [RStudio] Execute main method ###########
  caddScoresTable <- getCaddScoresTable(variants_table)
  
  ########### [RStudio] Print the basic table in a txt file ###########
  # Works only with basic table
  try(write.table(caddScoresTable,"test_files\\CB_caddScores_table.txt", append = FALSE, sep="\t",row.names=FALSE))
}