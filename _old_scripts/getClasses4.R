# This little script will enable us to determinate all the classes inside the dataframe 
# created by cellbase.

# The point is: the dataframe has dataframe and lists nested inside his cells.
# So we have to understand the structure.

# Input: the annotated variants table: annotVCFTable

##### Libraries #####
library(purrr)

##### METHOD FOR LIST ######
# extract_annot is a list of list
# In case we have a column with one list per cell
# We submit the column as extract_annot, which is a list of elements
getClassesInList <- function(extract_annot) { 
  annotClass <- character()
  
  # To get the classes inside the list,
  # we have to loop until we get one element with information 
  for (sub_extract_annot in extract_annot) { # for each element in list
    if ((length(annotClass)) == 0) {
      if (length(sub_extract_annot) != 0) { # element thaht contains information
        # Get the level class and add it to the vector
        annotClass <- class(sub_extract_annot)
        # annotClass_vec <- c(annotClass_vec, annotClass)
        extract_list <- list(sub_extract_annot, annotClass)
      }
    }
  }
  print(annotClass)
  return(extract_list)
}

#### METHOD 1 FOR DATAFRAME ###########
# Input: a dataframe
getClassesInDataframe <- function(extracted_df) { 
  # Get the header of the dataframe
  header_df <- colnames(extracted_df)
  
  print("Dataframe headers and classes:")
  # Get the list of classes
  columnClass <- lapply(extracted_df, class)
  print(columnClass)

  # extract_list contains header column and classes
  extract_list <- list(header_df, columnClass)
  return (extract_list)
}

##### METHOD 2 FOR DATAFRAME ######
# extract_annot is a list of dataframe
getClassesInDataframeFromList <- function(extract_annot) { 
  # All the cells should have the same structure,
  # So we will work only in the first row
  extracted_df <- extract_annot[[1]]
    
  # Get the header of the dataframe
  header_df <- colnames(extracted_df)
  
  print("Dataframe headers and classes:")
  # Get the list of classes
  annotClass <- lapply(extracted_df, class)
  print(annotClass)
  
  extract_list <- list(extracted_df, annotClass)
  return (extract_list)
}
  

############# Main method ########
annotClass_1 <- getClassesInDataframe(annotVCFTable)

for (annotation in annotClass_1[[1]]) {
  print(annotation)
  print(annotClass_1[[2]][annotation])
}

test <- annotClass_1[[2]]
test_df <- data.frame(test)

# Look through each column to get the classes inside
for (annotation in annotClass_1[[1]]) {
  print(annotation)
  print(annotClass_1[[2]][annotation])
  
  if (annotClass == "data.frame") { # As far as I saw, there are no dataframe in Level 1
    print("NOPE")
    
  } else if (annotClass == "list") { # Each element in list correspond to a variant + annotations
    # Extract the whole column to a new variable
    annotation <- 
    extract_annot <- annotVCFTable[[paste(annotation)]]
    
    # Get the second level class and add it to the vector
    # Since the extracted column was a list,
    # each element of the list corresponds to a variant whith his annotations
    # So we will extract only first variant to see what is the class inside
    annotClass2 <- class(extract_annot[[1]])
    annotClass2_vec <- c(annotClass2_vec, annotClass2)

    print(paste("Second level class:",class(extract_annot[[1]])))
    
    # Set a vector for third level classes
    annotClass3_vec <- character()
    
    if (annotClass2 == "data.frame") { # For each variant, there is a dataframe
      print("Third level:")
      extract_list3 <- getClassesInDataframe(extract_annot)
      
      for (element in extract_list3[2]) {
        if (element == "list") {
          print(element)
          annotClass <- getClassesInList(extract_annot)
        } else if (element == "data.frame") {
          print(element)
          }
      }

      
    } else if (annotClass2 == "list")  { # for each variant, there is a list
      
      annotClass3 <- character()
      for (sub_extract_annot in extract_annot) {
        if ((length(annotClass3)) == 0) {
          if (length(sub_extract_annot) != 0) { # no annotation for this variant
            # Get the third level class and add it to the vector
            annotClass3 <- class(sub_extract_annot)
            # annotClass3_vec <- c(annotClass3_vec, annotClass3)
          }
        }
      }
      print(paste("Third level class:", annotClass3))
    }
  }
  print ("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  }

  
  
  
