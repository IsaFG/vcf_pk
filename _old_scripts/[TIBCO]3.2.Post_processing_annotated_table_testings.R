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

########### [RStudio] Declare the chosen annotation ############
# In the next line, we have chosen the annotation "geneDrugInteraction" for testing
chosen_annot <- availableAnnots[12]

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
chosen_annot <- availableAnnots[13]

########### [TIBCO] Create the REvaluate object ########
AnnotatedTable <- REvaluate({
  ############# [Code] Load libraries ############# 
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  library(dplyr)
  library(tidyr)
  
  ############# [Code] Main method ############
  # Extract, as a single column, the chosen annotation from the table
  annot_column <- annotVariants_table[[paste(chosen_annot)]]
  
  # Get the index of the chosen annotation
  index_annot <- grep(chosen_annot, colnames(annotVariants_table))
  
  # Declare the new dataframe
  chosen_annot_table <- data.frame()
  
  # For each cell of the annotation column, loop to extact the content of the cell
  i = 1
  for (annot_cell in annot_column) {
    # Print info of the cell
    print (paste("Variant", i, "with class", class(annot_cell), "and lenght:",length(annot_cell)))
    print (paste("dimension",(dim(annot_cell))))
    
    # Assign variables
    data_class <- class(annot_cell)
    data_length <- length(annot_cell)
    basic_row <- basicTable[i,]
    
    # print (str(annot_column[i])) # testing line
    
    # Loop level 1 : analyze first level of the cell
    if (data_length == 0) {
      print ("Processing empty cell")
      annoted_row <- (basicTable[i,])
      annoted_row[, chosen_annot] <- "No result"
      
    } else if (data_class == "data.frame") {
      print ("Processing dataframe") # testing line
      
      # Create a data.frame from the cell
      cell_df <- data.frame(annot_column[i])
      
      # Loop level 2 : analyze cells contained inside the 1st level cell
      j = 1
      for (lev2_annot_cell in cell_df) {
        # Assign variables
        lev2_data_class <- class(lev2_annot_cell)
        lev2_data_length <- length(lev2_annot_cell)
        
        # analyze the level 2 cells 
        if (lev2_data_class == "data.frame") {
          # Print info of the cell
          print (paste("Cell", j, "with class", class(lev2_annot_cell), "and lenght:",length(lev2_annot_cell)))
          print (paste("dimension",(dim(lev2_annot_cell))))
          print ("Processing dataframe") # testing line
          
          # Create a data.frame from the cell
          lev2_cell_df <- data.frame(lev2_annot_cell)
          
          # Build a row with the chosen annotation of the current variant and the basic table
          lev2_annoted_row <- cbind(basicTable[i,], cell_df)
          
        } else if (lev2_data_class == "list") {
          # Print info of the cell
          print (paste("Cell", j, "with class", class(lev2_annot_cell), "and lenght:",length(lev2_annot_cell)))
          print (paste("dimension",(dim(lev2_annot_cell))))
          print ("Processing list") # testing line
          cell_df <- data.frame(annot_column[i])
          lev2_annoted_row <- cbind(basicTable[i,], cell_df)
        }
        j = j + 1
      }
      # Build a row with the chosen annotation of the current variant and the basic table
      annoted_row <- cbind(basicTable[i,], cell_df)
      
      
    } else if (data_class == "list") {
      print ("Processing list") # testing line
      cell_df <- data.frame(annot_column[i])
      annoted_row <- cbind(basicTable[i,], cell_df)
    }

    print ("Building the annotated table") # testing line
    chosen_annot_table <- bind_rows(chosen_annot_table, annoted_row)
    
    i = i + 1
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")# testing line
  }
  chosen_annot_table
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
