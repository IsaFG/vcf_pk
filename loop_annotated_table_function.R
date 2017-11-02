############## [INFO] General information ############
# Function to extract annotation from an annotated table
# created by getVariant method from cellbaseR.

############## [INFO] Input and Output ############################
# Input 1: an annotated variant table "annotVariants_table".
# -> In TIBCO, this table is stored as a "Document propertie"
# Input 2: chosen annotation

# Output : expanded annotated variants table

############# ANALYZE A ROW ##################
test_row <- annotated_row[1,]

test_class_vector <- character()
for (cell in test_row) {
  cell_class <- class(cell)
  test_class_vector <- c(test_class_vector, cell_class)
}

test_class_vector

test_row

############# ANALYZE A CELL ##################
# Declare the cell to process and the row where it belongs
cell_to_convert <- annotVariants_table[1,12]
basic_row <- basicTable[1,]

# Print the cell
cell_to_convert

# Identify the class of the cell
data_class <- class(cell_to_convert)
data_class

# Identify the structure of the cell
str(cell_to_convert)

############# [Code] Load libraries ############# 
library(RCurl)
library(jsonlite)
library(cellbaseR)
library(dplyr)
library(tidyr)

############# [Code] Method to process cells containing a list #############
# This function will convert any cell containing a list in readable information
simplifyLISTcell <- function (cell_to_convert) {
  print ("Processing a list") # testing line
  
  # Unlist the list. This will create a vector
  # This vector is supposed to be included in dataframes as a factor
  cell_unlisted <- unlist(cell_to_convert)
  
  # Convert the cell in a usefull format
  # complete_row <- data.frame(cell_unlisted) # not in use
  complete_row <- cell_unlisted
  
  # The result will be a new row created from the cell
  return(complete_row)
}

############# [Code] Method to process cells containing a dataframe #############
# This function will convert any cell containing a dataframe in readable information
# This function depends on the function "simplifyLISTcell
simplifyDFcell <- function (cell_to_convert) {
  print ("Processing a dataframe") # testing line
  
  # Create a new dataframe from the cell
  cell_df <- data.frame(cell_to_convert)
  
  # # Now, verify if none of the new cell in the new df has a list inside
  # j <- 1
  # rebuild_df <- cell_df
  # for (lev2_cell in cell_df) {
  #   print ("rebuild row:")
  #   print (rebuild_df)
  #   print (paste("Processing cell:", j, "inside annotation cell of variant", i)) # testing line
  #   lev2_data_class <- class(lev2_cell)
  #   if (data_class == "list") {
  #     rebuild_df <- cell_df[1:(j-1)]
  #     lev2_row <- simplifyLISTcell(lev2_cell)
  #     rebuild_df <- cbind(rebuild_df, lev2_row)
  #   } else {
  #     print ("Processing a list inside a dataframe") # testing line
  #     rebuild_df <- rebuild_df[1:j]
  #   }
  #   j <- j + 1
  # }
  
  complete_row <- rebuild_df
  
  # The result will be a new row created from the cell
  return(complete_row)
}

############# [Code] Method for list and for dataframes ############# 
convertCell <- function (cell_to_convert) {
  # In case the cell would contain a dataframe
  if (data_class == "data.frame") { 
    print ("Processing a dataframe") # testing line
    
    # Create a data.frame from the cell
    cell_df <- data.frame(cell_to_convert)
    complete_row <- cell_df
    
    # In case the cell would contain a list
  } else if (data_class == "list") { 
    print ("Processing a list") # testing line
    
    # Unlist the list and convert the new cell in a dataframe
    cell_unlisted <- unlist(cell_to_convert)
    
    # Convert the cell in a usefull format
    # complete_row <- data.frame(cell_unlisted) # not in use
    
    complete_row <- cell_unlisted
  }
  # The result will be a new row created from the cell
  return(complete_row)
}
