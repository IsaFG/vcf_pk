############## [INFO] General information ############
# This method will help usto analyze the composition of the cells in a dataframe

############# ANALYZE ##################
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

############# CONVERTION ##################
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
  # complete_row <- data.frame(cell_unlisted)
  
  complete_row <- cell_unlisted
}
# The result will be a new row created from the cell
return(complete_row)

annotation_row <- complete_row
annotated_row <- cbind(basic_row, annotation_row)


class(annotated_row$annotation_row)

############# View results #################
View(complete_row)
View(cell_unlisted)
str(cell_unlisted)
class(cell_unlisted)
