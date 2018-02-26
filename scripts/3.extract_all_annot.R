############## [INFO] Script 3 General information ############
# This method will process complex columns coming from a variants annotated table.
# The script will take the dataframe containing all the annotation obtained after querying cellbase
# and extract all annotation types stored to build one new table per annotation

############## [INFO] Input and Output ############################
# Input 1: an annotated variant table, variable "annotVCFTable".
# -> In TIBCO, this table is not readable and is stored as a Document propertie : annotVCFTableBlob
# -> In R console, this table should be previously loaded in global environment as a data.frame
# Input 2: List of available annotation : "availableAnnots"
# -> TIBCO: stored as a column
# -> R console: stored as a variable

# NOTE: Both input tables come from the script : 2.annotateIndexedVCF

# Output : a dataframe per annotation + a list storing all these dataframes

############# [INFO] Previous requirement ###############
# This script will only work with specifics variables stored as R objects
# and with some libraries (see input section)

# Working directory :
# Please uncoment the next line changing the working directory by the correct one:
# working_dir <- "C:\\..."

# R version : please uncoment the next line indicating the location of your current R version:
# r_version <- "C:/Program Files/R/R-3.4.1"

############## [INFO] Problems ############################
# Running this script throw the following error :
# Warning messages:
#   1: In data.frame(..., check.names = FALSE) :
#   row names were found from a short variable and have been discarded

########### [Code] Method to create a new table from a specific annotation ########
getSpecAnnotTable <- function(specific_annot, annotVariants_table) {
  ############# [Code] Load libraries ############# 
  # library(RCurl)
  library(jsonlite) #
  library(dplyr)
  # library(tidyr)
  
  ############# [Code] Sub-Method to process cells containing a dataframe #############
  # This function will convert any cell containing a dataframe in readable information
  # cell_to_convert <- annot_cell # For TESTING purpose
  simplifyDFcell <- function (cell_to_convert, loaded_annotations) {
    # print ("Processing a dataframe") # testing line
    
    # Create a new dataframe from the cell
    cell_df <- data.frame(cell_to_convert)
    
    # To know what classes are inside the dataframe
    class_info <- lapply(cell_df, class)
    classes_vec <- as.character(class_info)
    
    # Check if there still is a dataframe inside the new dataframe
    if ('data.frame' %in% classes_vec) {
      cell_df <- flatten(cell_df)
      class_info <- lapply(cell_df, class)
      classes_vec <- as.character(class_info)
    } 
    
    # Check if there is a list inside the new dataframe
    i <- 1
    for (class in classes_vec) {
      if (class == "list") {
        cell_df[,i] <- as.character(cell_df[,i])
        }
      i <- i + 1
    }
    
    # The result will be a new row created from the cell
    return(cell_df)
  }

  # Get the index of the chosen annotation
  index_annot <- grep(specific_annot, colnames(annotVariants_table))
  
  # Extract, as a single variable, the chosen annotation from the table
  annot_column <- annotVariants_table[,index_annot]
  
  # Declare the basic table
  basicTable <- annotVariants_table[,1:5]
  
  # Declare the dataframe that will contain all the info of the chosen annotation
  specific_annot_table <- data.frame()
  
  # Determinate the class of all the cells inside the column
  column_classes <- lapply(annot_column, class)
  
  # For each cell of the annotation column, loop to extract the content of the cell
  i = 1
  for (annot_cell in annot_column) {
    # # Print info of the cell being analyzed
    # print (paste("Variant", i, "with class", class(annot_cell), "and lenght:",length(annot_cell)))
    # print (paste("dimension",(dim(annot_cell))))
    
    # Assign variables to start the analyze
    data_class <- class(annot_cell)
    data_length <- length(annot_cell)
    basic_row <- basicTable[i,]
    
    # print (str(annot_column[i])) # testing line
    
    if (data_length == 0) { # In case the cell would have no result for this annotation
      # print ("Processing empty cell")
      annotated_row <- basic_row
      # annotated_row[, specific_annot] <- "No result"
      
    } else if (data_class == "data.frame") { # In case the cell would contain a dataframe
      # Apply the function  to produce a row/df from the cell
      annotation_row <- simplifyDFcell(annot_cell)
      
      # Build a row with the chosen annotation of the current variant
      # and bind it to basic table
      annotated_row <- cbind(basic_row, annotation_row)
      
    } else if (data_class == "list") { # In case the cell would contain a list
      annotation_row <- as.character(annot_cell)
      
      # Build a row with the chosen annotation of the current variant
      # and bind it to basic table
      annotated_row <- cbind(basic_row, annotation_row)
      
    } else {# In case the cell is not problematic
      annotation_row <- annot_cell
      
      # Build a row with the chosen annotation of the current variant
      # and bind it to basic table
      annotated_row <- cbind(basic_row, annotation_row)
    }
    
    # print ("Building the annotated table") # testing line
    specific_annot_table <- bind_rows(specific_annot_table, annotated_row)
    
    i = i + 1
    # print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")# testing line
  }
  print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")# testing line
  
  # Eliminate possible duplicate rows in table
  specific_annot_table <- unique(specific_annot_table)
  
  return(specific_annot_table)
}

########### [Code] Method to create a list with all the tables ########
getTablesList <- function(availableAnnots, annotVariants_table){
  tables_list <- list()
  i <- 1
  for (annotation in availableAnnots){
    print (paste("Building table", annotation))
    annotation_table <- getSpecAnnotTable(annotation, annotVariants_table)
    tables_list[[i]] <- annotation_table
    i <- i + 1
  }
  return(tables_list)
}

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {

    ########### [TIBCO] Load RinR library ###################
    library(RinR)
    ########### [TIBCO] Determinate R interpreter location ########
    Rversion <- makeREvaluator("R", RHome = r_version)
    
    ############ [TIBCO] Load variables ###############
    # Load the annotated table, stored in memory as a Blob Object created by script "cellbaseR_Query_getVariant"
    annotVCFObject <- BlobToSObject(annotVCFTableBlob)
    annotVCFTable <- annotVCFObject
    # For testing prupose: you can pre-set the annotation "geneDrugInteraction" 
    # specific_annot <- availableAnnots[12]
    
    ########### [TIBCO] Create the REvaluate object ########
    tablesList <- REvaluate({
      tablesList <- getTablesList(availableAnnots, annotVariants_table)
      tablesList
    },
    data = list(getSpecAnnotTable = getSpecAnnotTable, getTablesList = getTablesList, annotVariants_table = annotVCFTable, availableAnnots = availableAnnots)
    # ,
    # REvaluator = Rversion,
    # verbose	= TRUE
    )
    
    for (i in (1:length(tablesList))) {
      # Load each table in a specific variable in global environment
      table_name <- paste(availableAnnots[i],"Table", sep="")
      table_extracted <- tablesList[[i]]
      assign(table_name, table_extracted)
      }
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd(working_dir)
  
  # Check if the directory for output files exists, and create it if needed
  if (file.exists(file.path(working_dir, "output_files/annotatedTables"))){
    # do nothing
  } else {
    dir.create(file.path(working_dir, "output_files/annotatedTables"))
  }

  ########### [RStudio] Execute main method ###########
  tablesList <- getTablesList(availableAnnots, annotVCFTable)
  
  ########### [RStudio] Save the information ###########
  for (i in (1:length(tablesList))) {
    # Load each table in a specific variable in global environment
    table_name <- paste(availableAnnots[i],"Table", sep="")
    table_extracted <- tablesList[[i]]
    assign(table_name, table_extracted)
    # Print the table in a txt file
    # Works only with basic tables
    file_path <- paste("output_files\\annotatedTables\\annotated_table_",availableAnnots[i], ".txt",sep = "")
    try(write.table(table_extracted,file_path, append = FALSE, sep="\t",row.names=FALSE))
  }
}
########### [TESTING] Control classes of the chosen annotation table #######
# Just to be sure this table is suitable for TIBCO
TEST_classes <- lapply(specific_annot_table, class)
View(annot_cell)