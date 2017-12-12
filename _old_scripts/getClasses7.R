# This little script will enable us to determinate all the classes inside the dataframe 
# created by cellbase.

# The point is: the dataframe has dataframe and lists nested inside his cells.
# So we have to understand the structure.

# Input: the annotated variants table: annotVCFTable

# ########## test variables #########
annotation <- "traitAssociation"
subAnnotation <- "genomicFeatures"

annotation <- "consequenceTypes"
subAnnotation <- "sequenceOntologyTerms"

# LEVEL 1
# Get the index of the chosen annotation
index_annot <- grep(annotation, colnames(annotVCFTable))
# Extract, as a single variable, the chosen annotation from the table
element1 <- annotVCFTable[,index_annot]
# NOTE: The class of the element will ALWAYS be "list"

# LEVEL 2
element2_class <- character()
basic_df <- annotVCFTable[,1:5]
subAnnotation_df <- data.frame()
j <- 1
for (element2 in element1) { # for each element (variant) in list
  element2_class <- class(element1[[1]])
  
  print("LEVEL 2")
  print(element2_class)
  
  # LEVEL 3
  # Note : level 2 class will almost always be a dataframe (except for repeat)
  if (element2_class == "data.frame") {
    i <- 1
    for (element3 in element2){
      element3_class <- class(element3)
      
      # print("LEVEL 3")
      # print(element3_class)
      
      if (colnames(element2)[i] == subAnnotation) {
        print ("Subannotation FOUND")
        element4_df <- data.frame()
        for (element4 in element3) {
          element4_df <- cbind(basic_df[j,], element4)
          lapply(element4_df, )
          subAnnotation_df <- bind_rows(subAnnotation_df, element4_df[1,])
        }
      }
      i <- i + 1
    }
    
    # NOT DONE FOR LIST SINCE ONLY ONE ANNOTATION TABLE (REPEAT) HAS THIS SITUATION
  } else if (element2_class == "list")  {
    for (element3 in element2){
      element3_class <- class(element3)
      print (element3_class)
    }
  }
  j <- j + 1      
}

