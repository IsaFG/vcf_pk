############## [INFO] Script 4 General information ############
# This method will process complex columns coming from a variants annotated table.
# The user must choose the annotation he wishes to be displayed
# and the script will extract this annotation and build a new table
# with the chosen annotation.

############## [INFO] Input and Output ############################
# Input 1: an annotated variant table, variable "annotVCFTable".
# -> In TIBCO, this table is not readable and is stored as a "Document propertie"
# -> In R console, this table should be present in global environment as a data.frame
# Input 2:
# -> TIBCI: the annotation type chosen by the user
# -> R console:List of avaiblable annotations, variable "availableAnnots"

# NOTE: Both input tables come from the script : 3.annotateIndexedVCF

# Output : expanded annotated variant tables, according to the annotation chosen by user.

############# [INFO] Previous requirement ###############
# This script will only work with specifics variables stored as R objects
# and with some libraries (see input section)

############## [INFO] Problems ############################
# Running this script throw the following error :
# Warning messages:
#   1: In data.frame(..., check.names = FALSE) :
#   row names were found from a short variable and have been discarded

########### [Code] Main method ########
getChosenAnnotTable <- function(chosen_annot, annotVariants_table) {
  ############# [Code] Load libraries ############# 
  # library(RCurl)
  library(jsonlite) #
  library(dplyr)
  # library(tidyr)
  
  ########### [Code] [Method] Verify if the annotation has not been already loaded ######
  if (is.element(chosen_annot, loaded_annotations) == FALSE) {
    
    loaded_annotations <- c(loaded_annotations, chosen_annot)
  }
  
  ############# [Code] Method to process cells containing a dataframe #############
  # This function will convert any cell containing a dataframe in readable information
  # cell_to_convert <- annot_cell # For TESTING purpose
  simplifyDFcell <- function (cell_to_convert) {
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
  
  ############# [Code] Main method ############
  # Get the index of the chosen annotation
  index_annot <- grep(chosen_annot, colnames(annotVariants_table))
  
  # Extract, as a single variable, the chosen annotation from the table
  annot_column <- annotVariants_table[,index_annot]
  
  # Declare the basic table
  basicTable <- annotVariants_table[,1:5]
  
  # Declare the dataframe that will contain all the info of the chosen annotation
  chosen_annot_table <- data.frame()
  
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
      print ("Processing empty cell")
      annotated_row <- basic_row
      # annotated_row[, chosen_annot] <- "No result"
      
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
    
    print ("Building the annotated table") # testing line
    chosen_annot_table <- bind_rows(chosen_annot_table, annotated_row)
    
    i = i + 1
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")# testing line
  }
  return(chosen_annot_table)
}

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  ########### [TIBCO] Load RinR library ###################
  library(RinR)
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
  
  ############ [TIBCO] Load variables ###############
  # Load the annotated table, stored in memory as a Blob Object created by script "cellbaseR_Query_getVariant"
  annotVCFObject <- BlobToSObject(annotVCFTableBlob)
  annotVariants_table <- annotVCFObject
  # For testing prupose: you can pre-set the annotation "geneDrugInteraction" 
  # chosen_annot <- availableAnnots[12]
  
  ########### [TIBCO] Create the REvaluate object ########
  AnnotatedTable <- REvaluate({
    chosen_annot_table <- getChosenAnnotTable(chosen_annot, annotVariants_table)
    chosen_annot_table
  },
  data = list(getChosenAnnotTable = getChosenAnnotTable, annotVariants_table = annotVCFTable, chosen_annot = chosen_annot)
  # ,
  # REvaluator = Rversion,
  # verbose	= TRUE
  )
  
  ########### [TIBCO] [PENDING] Save the annotated table as a Blob Object #######
  # annotationsBlob <- SObjectToBlob(AnnotatedTable)
  # assign(paste(chosen_annot,"Table",sep=""), annotationsBlob)
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Get the input variables ############
  chosen_annot = availableAnnots[14]
  # availableAnnots, annotVariants_table
  # They should already be loaded in global environment and come from script 3.cellbaseR_Query_getVariant
  
  ########### [RStudio] Execute main method ###########
  chosen_annot_table <- getChosenAnnotTable(chosen_annot, annotVCFTable)
  
  ########### [RStudio] Print the table in a txt file ###########
  # Works only with basic table
  file_path <- paste("test_files\\annotated_table_",chosen_annot, ".txt",sep = "")
  try(write.table(chosen_annot_table,file_path, append = FALSE, sep="\t",row.names=FALSE))
}

########### [TESTING] Control classes of the chosen annotation table #######
# Just to be sure this table is suitable for TIBCO
TEST_classes <- lapply(chosen_annot_table, class)
View(annot_cell)