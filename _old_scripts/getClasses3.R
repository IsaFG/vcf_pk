# This little script will enable us to determinate all the classes inside the dataframe 
# created by cellbase.

# The point is: the dataframe has dataframe and lists nested inside his cells.
# So we have to understand the structure.

# Input: the annotated variants table: annotVCFTable

##### METHOD FOR LIST ######
getClassesInList <- function(extract_annot) { # extract_annot is a list
  annotClass <- character()
  
  for (sub_extract_annot in extract_annot) { # for each element in list
    if ((length(annotClass)) == 0) {
      
      if (length(sub_extract_annot) != 0) { # the element should contain information
        # Get the level class and add it to the vector
        annotClass <- class(sub_extract_annot)
        # annotClass_vec <- c(annotClass_vec, annotClass)
      }
    }
  }
  print(annotClass)
  return(annotClass)
}

##### METHOD FOR DATAFRAME ######
getClassesInDataframe <- function(extract_annot) { # extract_annot is a dataframe
  # Get the header of the dataframe
  header_df <- colnames(extract_annot[[1]])
  
  print("Dataframe headers and classes:")
  annotClass <- lapply(extract_annot[[1]], class)
  
  print(annotClass)
  return (annotClass)
}

# Get the header of the dataframe
annotations_in_table <- colnames(annotVCFTable)

# Set a vector for first level classes
annotClass1_vec <- character()

# annotation <- "consequenceTypes"

for (annotation in annotations_in_table){
  # Loop inside the Level 1
  # As far as I know, there are only lists as complex variable
  
  # Print annotation name and his class
  print(paste("Annotation",annotation))
  
  # Get the first level class
  annotClass <- class(annotVCFTable[[paste(annotation)]])
  # add the first level class to a vector of classes
  annotClass1_vec <- c(annotClass1_vec, annotClass)
  
  # print class in level 1
  print (paste("First level class:", annotClass))
  
  # Set a vector for second level classes
  annotClass2_vec <- character()
  
  # Look through each column to get the classes inside 
  if (annotClass == "data.frame") { # As far as I saw, there are no dataframe in Level 1
    print("NOPE")
    
  } else if (annotClass == "list") { # Each element in list correspond to a variant + annotations
    # Extract the whole column to a new variable
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
      annotClass3 <- getClassesInDataframe(extract_annot)
      
      for (element in annotClass3) {
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

  
  
  
