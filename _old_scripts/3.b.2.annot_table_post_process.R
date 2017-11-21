############## [INFO] SCRIPT 3 General information ############
# Post processing of the pre-processed annotated table

############## [INFO] Input and Output ############################
# Inputs : a preprocessed annotated table (script 3.b.1.annotateIndexedVCF)
# Output : a post-processed annotated table 

# IMPORTANT NOTE : the OUTPUT table cannot be send to a txt file neither to a TIBCO table
# The table should be keeped in as an R object in R environment or in a TIBCO Blob object

############# [INFO] Previous requirement (libraries) ###############
# LIBRARIES from CRAN : RCurl, jsonlite
# LIBRARY to use this script in TIBCO : RinR
# LIBRARY from Bioconductor : "cellbaseR":
# https://bioconductor.org/packages/release/bioc/html/cellbaseR.html

############## [INFO] Problems ############################
# Problem 1:  Some variants have more than 1 REF or ALT allele
# so that this script do the annotation call with the first allele of each field.
# I will change that just deleting these variants

# Problem 2: The annotated table has sometimes columns with dataframe inside
# this has been fixed with additional code, but it could slow down the script...

# Problem 3 : Some warnings messages when running the entire scripts :
# "There were 50 or more warnings (use warnings() to see the first 50)
# 1: In bind_rows_(x, .id) : Unequal factor levels: coercing to character
# 2: In bind_rows_(x, .id) :
#   binding character and factor vector, coercing into character vector

########### [Code] Main method ########
getAnnotVariantsTable <- function(preAnnotVCFTable) {
  ############# [Code] Load libraries ############# 
  library(RCurl)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  
  ############# [Code] ############
  var_number <- nrow(preAnnotVCFTable)
  
  class_vec <- lapply(preAnnotVCFTable, class)
  dim_vec <- lapply(preAnnotVCFTable, dim)

  for (i in 1:length(preAnnotVCFTable)) { 

    col_name <- names(class_vec)[i]
    col_class <- class_vec[[col_name]]
    col_dim <- dim_vec[[col_name]]
    
    print (paste("Processing column", i, ":", col_name, ", class: ", col_class, ", dimension", col_dim)) # this line is for testing

    # if (class_chr=="data.frame")  { 
    
    if (is.null(col_dim) == FALSE) {  # anormal situation, the column has a dimension > 0
      print (paste("WARNING for column", col_name, ": MULTIDIMENSIONAL column(s)"))
      
      # col_name <- "variantTraitAssociation" # TESTING LINE
      # i <- 16 # TESTING LINE
      
      # assign the data.frame to a new variable
      nested_df <- data.frame(preAnnotVCFTable[[col_name]])
      
      # Scan the data.frame in order to find the problematic column(s)
      print (paste("Analyzing columns of column", col_name)) # testing line
    
      # Build the table with all the annotation
      # Subset the original table
      preAnnotVCFTable_sub <- bind_cols(preAnnotVCFTable[1:(i-1)], preAnnotVCFTable[(i+1):length(preAnnotVCFTable)])
      annotVCFTable <- bind_cols(preAnnotVCFTable_sub, nested_df)
    }
    # warnings()
    i <- i + 1
  }
  return(annotVCFTable)
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
  annotVCFTable <- REvaluate({
    annotVariants_table <- getAnnotVariantsTable(preAnnotVCFTable)
    annotVariants_table
  }
  , data = list(getAnnotVariantsTable = getAnnotVariantsTable, preAnnotVCFTable = preAnnotVCFTable)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
  ########### [TIBCO] Comvert the table to a Blob Object ########
  annotVCFBlob <- SObjectToBlob(annotVariants_table)
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Execute main method ###########
  annotVCFTable <- getAnnotVariantsTable(preAnnotVCFTable)
  
  ########### [RStudio] Print the basic table and available annotations in a txt file ###########
  # Works only with basic table
  try(write.table(annotVariants_table[,1:5],"test_files\\VCF_basic_table.txt", append = FALSE, sep="\t",row.names=FALSE))
  try(write.table(colnames(annotVariants_table[6:length(annotVariants_table)]),"test_files\\available_annotations.txt", append = FALSE, sep="\t",row.names=FALSE))
}
