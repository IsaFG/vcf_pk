# This little script will enable us to determinate all the classes inside the dataframe 
# created by cellbase.

# The point is: the dataframe has dataframe and lists nested inside his cells.
# So we have to understand the structure.

# Input: the annotated variants table: annotVCFTable

##### Libraries #####
library(purrr)

##### METHOD FOR LIST ######
# extracted_list is a list of list
# In case we have a column with one list per cell
# We submit the column as extracted_list, which is a list of elements
getClassesInList <- function(extracted_list) { 
  elementClass <- character()
  
  # To get the classes inside the list,
  # we have to loop until we get one element with information 
  for (element in extracted_list) { # for each element in list
    if ((length(elementClass)) == 0) {
      if (length(element) != 0) { # element thaht contains information
        # Get the level class and add it to the vector
        elementClass <- class(element)
        # annotClass_vec <- c(annotClass_vec, annotClass)
        extract_list <- list(elementClass, element)
      }
    }
  }
  return(extract_list)
}

#### METHOD FOR DATAFRAME ###########
# Input: a dataframe
getClassesInDataframe <- function(extracted_df) { 
  
  # Get the list of classes
  columnClass <- lapply(extracted_df, class)

  # extract_list contains header column and classes
  extractClasses <- data.frame(columnClass)
  return (extractClasses)
}

############# Main method ########
# LEVEL 1: Get classes and headers from the dataframe
print("Looping LEVEL 1")
annotClassDF_1 <- getClassesInDataframe(annotVCFTable)

# Look through each column to get the classes inside
print(" + + Looping LEVEL 2")
for (index in 1:length(annotClassDF_1)) {
  annotation <- colnames(annotClassDF_1)[index]
  class_1 <- annotClassDF_1[index]

  print (paste(" + Looping", annotation, "with following class:"))
  print(class_1)
  
  if (class_1 == "data.frame") { # As far as I saw, there are no dataframe in Level 1
    print ("ERROR: Second level class: dataframe")
    
  } else if (class_1 == "list") { # Each element in list correspond to a variant + annotations
    # Extract the whole column to a new variable which will be a list of elements
    print (" + + Looping SECOND level class LIST")
    extract_2 <- annotVCFTable[[paste(annotation)]]
    if (annotation == "repeat.") {
      extract_2 <- annotVCFTable["repeat"]
    }
    
    # Get the second level class and add it to the vector
    # Since the extracted column was a list,
    # each element of the list corresponds to a variant whith his annotations
    # So we will extract only first variant to see what is the class inside
    class_2 <- getClassesInList(extract_2)
    
    if (class_2[[1]] == "data.frame") { # For each variant, there is a dataframe
      print(" + + + Looping third level class DATAFRAME")
      the_dataframe <- class_2[[2]]
      class_3 <- getClassesInDataframe(the_dataframe)
      
      for (index in 1:length(class_3)) {
        if (class_3 == "list"){
          print("---> + + + + LIST: Need more loops in level 4")
        } else if (class_3 == "data.frame") {
          print("---> + + + + DATAFRAME : Need more loops in level 4")
        }
       }

    } else if (class_2[[1]] == "list")  { # for each variant, there is a list
      print(" + + + Looping Third level class LIST")
      class_3 <- getClassesInList(class_2[[2]])
      
      if (class_3[[1]] == "list"){
        print("---> + + + + LIST: Need more loops in level 4")
      } else if (class_3[[1]] == "data.frame") {
        print("---> + + + + DATAFRAME : Need more loops in level 4")
      } else {
        print("---> + + + + OK, finish for this one")
      }
      
    } else {
      print ("Third level is OK")
    }
      
  } else {
    print ("Second level is OK")
  }
  
  print ("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print ("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  }


  
  
  
