############## [INFO] Script 3.2 General information ############
# This methods will process complex columns coming from a variants annotated table.
# The user must choose the annotation he wishes to be displayed
# and the script will extract this annotation and build a new table with the chosen annotation.

############## [INFO] Input and Output ############################
# Input 1: an annotated variant table "annotVariants_table".
# -> In TIBCO, this table is stored as a "Document propertie"
# Input 2: Avaiblable annotations, variable "availableAnnots"

# NOTE: Both input tables come from the script : 3.cellbaseR_Query_getVariant

# Output : expanded annotated variant tables, according to the annotation chosen by user.

############# [INFO] Previous requirement ###############
# This script will only work with specifics variables stored as R objects and with some libraries.
# If you use RStudio, you need to run script 3.cellbaseR_Query_getVariant firstly
# and keep in the R environment the variables described in input section above.

############## [INFO] Problems ############################
# I still have to test this scripts with other annotations.
# Running this script throw the following error :
# Warning messages:
#   1: In data.frame(..., check.names = FALSE) :
#   row names were found from a short variable and have been discarded

# ############# [RStudio] Set Working directory ############
# Please uncoment the next line changing the working directory by the correct one:
# setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [TIBCO] Load RinR library ###################
library(RinR)

########### [TIBCO] Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

############ [TIBCO] Load variables ###############
# Load the annotated table, stored in memory as a Blob Object created by script "cellbaseR_Query_getVariant"
annotVariantsSObject <- BlobToSObject(annotVariantsBlob)
annotVariants_table <- annotVariantsSObject
# Load annotation chosen by user
# in this case the annotation "geneDrugInteraction" is pre-chosen for testing
chosen_annot <- availableAnnots[12]

########### [TIBCO] Create the REvaluate object ########
AnnotatedTable <- REvaluate({
  ############# [Code] Load libraries ############# 
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  library(dplyr)
  library(tidyr)
  
  # Testing block
  cell_to_convert <- annotVariants_table[1,12]
  row_to_verify <- cell_df[1,]
  TEST_verifyCellClasses <- verifyCellClasses(row_to_verify)
  TEST_verifyCellClasses
  
  TEST_classes <- lapply(chosen_annot_table, class)
  # Testing block
  
  ############# [Code] Method to verify the classes in a row #############
  # This method returns a vector with the classes of all the cells included in a row
  # It returns also vectors with the index where there are dataframes and lists
  verifyCellClasses <- function (row_to_verify) {
    classes_vector <- character()
    list_index <- numeric()
    df_index <- numeric()
    i <- 1
    for (any_cell in row_to_verify) {
      cell_class <- class(any_cell)
      if (cell_class == "list") {
        list_index <- c(list_index,i)
      } else if (cell_class == "data.frame") {
        df_index <- c(df_index,i)
      }
      classes_vector <- c(classes_vector, cell_class)
      i <- i + 1
    }
    classes_info <- list(cells_with_list = list_index, cells_with_dataframe = df_index, classes_vector = classes_vector)
    return(classes_info)
  }
  
  

  # ############# [Code] Method to process cells containing a list #############
  # # This function will convert any cell containing a list in readable information
  # simplifyLISTcell <- function (cell_to_convert) {
  #   print ("Processing a list") # testing line
  #   
  #   # Unlist the list. This will create a vector
  #   # This vector is supposed to be included in dataframes as a factor
  #   cell_unlisted <- unlist(cell_to_convert)
  #   
  #   # Convert the cell in a usefull format
  #   # complete_row <- data.frame(cell_unlisted) # not in use
  #   complete_row <- cell_unlisted
  #   
  #   # The result will be a new row created from the cell
  #   return(complete_row)
  # }
  
  
  ############# [Code] Method to process cells containing a dataframe #############
  # This function will convert any cell containing a dataframe in readable information
  simplifyDFcell <- function (cell_to_convert) {
    print ("Processing a dataframe") # testing line
    
    # Create a new dataframe from the cell
    cell_df <- data.frame(cell_to_convert)
    cell_dfBIS <- cell_df

    # To know what classes are inside the dataframe
    class_info <- lapply(cell_df, class)
    classes_vec <- as.character(class_info)

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
  
  annotVariants_tableBIS <- annotVariants_table # TESTING !!
  # [RStudio] Declare the chosen annotation ############
  # In the next line, we have chosen the annotation "geneDrugInteraction" for testing
  chosen_annot <- availableAnnots[12]
  
  # Get the index of the chosen annotation
  index_annot <- grep(chosen_annot, colnames(annotVariants_table))
  
  # Get the classes inside the annotation table
  annotVariants_table_cl <- lapply(annotVariants_table, class) # May be unnecessary

  # Extract, as a single variable, the chosen annotation from the table
  annot_column <- annotVariants_table[[paste(chosen_annot)]]
  # or
  annot_column <- annotVariants_table[,index_annot]
  
  # Declare the dataframe that will contain the chosen annotation
  chosen_annot_table <- data.frame()
  
  # Determinate the class of all the cells inside the column
  column_classes <- lapply(annot_column, class)

  # For each cell of the annotation column, loop to extact the content of the cell
  i = 1
  for (annot_cell in annot_column) {
    # Print info of the cell being analyzed
    print (paste("Variant", i, "with class", class(annot_cell), "and lenght:",length(annot_cell)))
    print (paste("dimension",(dim(annot_cell))))
    
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
  # chosen_annot_table # for TIBCO only
},
data = list(annotVariants_table = annotVariants_table, chosen_annot = chosen_annot, basicTable = basicTable)
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)

########### [RStudio] Print the table in a txt file ###########
# Works only with basic table
file_path <- paste("test_files\\annotated_table_",chosen_annot, ".txt",sep = "")
try(write.table(chosen_annot_table,file_path, append = FALSE, sep="\t",row.names=FALSE))

########### [RStudio] Testing block ###########
class(chosen_annot_table$associationTypes)
class(chosen_annot_table$associationTypes[1])
unlist(chosen_annot_table$associationTypes[1])
