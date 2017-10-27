############## [INFO] Script 3.2 General information ############
# This scripts WORKS !!

# This methods will process complex columns coming from a variants annotated table.
# Table coming from script : 3.cellbaseR_Query_getVariant.

############## [INFO] Input and Output ############################
# Input 1: an annotated variant table "annotVariants_table" (R object created by script : 3.cellbaseR_Query_getVariant)
# -> NOTE: the table will be an R object, not a txt file, neither a TIBCO table.
# Input 2: Avaiblable annotations, vaiable "availableAnnots"

# Output : expanded annotated variant tables, according to the annotation chosen by user.

############# [INFO] Previous requirement ###############
# This script will only work with R objects, not with any file
# You need to run script 3.cellbaseR_Query_getVariant firstly
# and keep the variables described in input section above 

############## [INFO] Problems ############################


# ############# [RStudio] Set Working directory ############
# Please uncoment the next line changing the working directory by the correct one:
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")


# ############# [TIBCO] Load RinR library ###################
# library(RinR)
# 
# ########### [TIBCO] Determinate R interpreter location ########
# Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")


########### [TIBCO] Annotation chosen by the user ############
# For TIBCO only
chosen_annot <- availableAnnots[13] # example
# Get the index of the chosen annotation
index_annot <- grep(chosen_annot, colnames(annotVariants_table))

# ########### [TIBCO] Create the REvaluate object ########
# annotatedTable <- REvaluate({

############# [Code] Load libraries ############# 

############# [Code] Main method ############
# TEST NOTE : the TEST_table will be used to test the code on a few number of variant
TEST_table <- annotVariants_table[1:3,]
TEST_table <- annotVariants_table
# Extract the chosen annotation from the table as a single column
annot_column <- TEST_table[[paste(chosen_annot)]]
# Preset the new dataframe
chosen_annot_table <- data.frame()

# Loop to extract the information from the MAIN dataframe
i = 1
for (annot_cell in annot_column) {
  # Print info of the cell
  print (paste("Variant", i, "with class", class(annot_cell), "and lenght:",length(annot_cell)))
  print (paste("dimension",(dim(annot_cell))))
  
  # Assign variables
  data_class <- class(annot_cell)
  data_length <- length(annot_cell)
  basic_row <- basicTable[i,]
  
  ## Multiply rows of the basic table by the length of the dataframe / list in the data
  ## This will be usefull in order to create a readable table according to the chosen table
  # prechosen_table <- basic_row[rep(seq_len(nrow(basic_row)), each=data_length),]
  
  print (str(annot_column[i]))
  
  if (data_length == 0) {
    print ("Processing empty cell")
    new_row <- (basicTable[i,])
    new_row[, chosen_annot] <- "No result"

    } else if (data_class == "data.frame") {
      print ("Processing dataframe")
      
      # Create a data.frame from the chosen annotation of the current variant
      column_df <- data.frame(annot_column[i])
      
      # Build a row with the chosen annotation of the current variant and the basic table
      new_row <- cbind(basicTable[i,], column_df)
      # add_annot_columns <- cbind(prechosen_table, annot_cell)
  
    } else if (data_class == "list") {
      print ("Processing list")
      column_df <- data.frame(annot_column[i])
      new_row <- cbind(basicTable[i,], column_df)
  }
  
  print ("Building the annotated table")
  chosen_annot_table <- bind_rows(chosen_annot_table, new_row)
  
  ## Bind multiplied columns
  # chosen_annot_table <- bind_rows(chosen_annot_table, prechosen_table)
  
  # Testing block
  if (i == 1) {
    TEST_OK_column1 <- column_df
    TEST_OK_row1 <- new_row
    TEST_OK_table1 <- chosen_annot_table
  } else if (i == 2) {
    TEST_NOK_column2 <- column_df
    TEST_NOK_row2 <- new_row
    TEST_NOK_table2 <- chosen_annot_table
  }
  # end of Testing block
  
  i = i + 1
  print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
}

TEST_df <- data.frame(annot_column)
View(TEST_df[1,])
View(TEST_df[2,])


str(annot_cell)

for (field in annot_cell) {
  print 
}