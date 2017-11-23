############## [INFO] SCRIPT 3 General information ############
# We provide an indexed VCF (file with his tabix file),
# then the script will query the RESTfull DB CellBase to retreive information about the variants

############## [INFO] Input and Output ############################
# Inputs : an instance file for a VCF file and his tabix file
# -> instance file comes from the script 1.b.index_VCF_file

# Output : a table with all the variant with their annotations
# + a table with only the rsid annotation
# + a table with all the available annotations type

# IMPORTANT NOTE : the OUTPUT table with the annotations cannot be send to a txt file neither to a TIBCO table
# The table should be kept as an R object in R environment or in a TIBCO Blob object

############# [INFO] Previous requirement (libraries) ###############
# LIBRARIES from CRAN : dplyr
# LIBRARY to use this script in TIBCO : RinR
# LIBRARY from Bioconductor : "cellbaseR":
# https://bioconductor.org/packages/release/bioc/html/cellbaseR.html

############## [INFO] Problems ############################
# No problem detected so far (finger crossed)

########### [Code] Method 1 : query CellBase with the VCF file ########
annotIndexedVCF<- function(path_to_indexed_files) {
  ############# [Code] Load libraries ############# 
  library(cellbaseR)
  
  ############# [Code] Annotate the VCF file with cellbase ############
  path_to_indexed_files <- locIndexedVCF$path # TEST LINE
  
  cb <- CellBaseR()
  fl <- path_to_indexed_files
  annotVCF_table <- AnnotateVcf(object=cb, file=fl, BPPARAM = bpparam(workers=2))

  # Change the name of the "id" column to avoid repetition
  colnames(annotVCF_table)[5] <- "rsID"
  return(annotVCF_table)
}

########### [Code] Method 2 : post-process the annotation table ########
getAnnotVariantsTable <- function(preAnnotVCFTable) {
  ############# [Code] Load libraries ############# 
  library(dplyr)

  ############# [Code] Post-process the table ############
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
  
  path_to_indexed_files <- locIndexedVCF$path
  
  ########### [TIBCO] Create the REvaluate object to execute method 1 ########
  preAnnotVCFTable <- REvaluate({
    annotVCF_table <- annotIndexedVCF(path_to_indexed_files)
    annotVCF_table
  }
  , data = list(annotIndexedVCF = annotIndexedVCF, path_to_indexed_files = path_to_indexed_files)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
  ########### [TIBCO] Create the REvaluate object to execute method 2 ########
  annotVCFTable <- REvaluate({
    annotVariants_table <- getAnnotVariantsTable(preAnnotVCFTable)
    annotVariants_table
  }
  , data = list(getAnnotVariantsTable = getAnnotVariantsTable, preAnnotVCFTable = preAnnotVCFTable)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
  ########### [TIBCO] Convert the tables to Blob Objects ########
  preAnnotVCFTableBlob <- SObjectToBlob(preAnnotVCFTable)
  annotVCFTableBlob <- SObjectToBlob(annotVCFTable)
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Execute method 1 ###########
  # Note: the variable locIndexedVCF is the instance file of the indexed files
  # and comes from a previous script
  preAnnotVCFTable <- annotIndexedVCF(locIndexedVCF$path)
  
  ########### [RStudio] Execute method 2 ###########
  annotVCFTable <- getAnnotVariantsTable(preAnnotVCFTable)
  
  ########### [RStudio] Print the basic table and available annotations in a txt file ###########
  # Works only with basic table
  try(write.table(preAnnotVCFTable[,1:5],"test_files\\annotated_indVCF.txt", append = FALSE, sep="\t",row.names=FALSE))
  try(write.table(colnames(preAnnotVCFTable[6:length(preAnnotVCFTable)]),"test_files\\available_annotations_indVCF.txt", append = FALSE, sep="\t",row.names=FALSE))
}

########### [Code] Build a basic table that could be loaded whithout problems ########
# This table has no nested data.frame neither nested list 
basicTable <- annotVCFTable[,1:5]

########### [Code] Get the available annotations ############
availableAnnots <- colnames(annotVCFTable[5:length(annotVCFTable)])
