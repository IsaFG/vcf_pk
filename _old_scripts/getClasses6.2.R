############## [INFO] Script 6 General information ############
# Some annotation tables show much information.
# Even if all these tables can be loaded in TIBCO,
# they can have lists inside them, which make the info hard to read.
# This script generate additional sub-tables for some interest information.

############## [INFO] Input and Output ############################
# Input 1: an table with a specific anntotion.
# -> In TIBCO, this table is not readable and is stored as a Document propertie
# -> In R console, this table should be present in global environment as a data.frame
# Input 2: List of annotations of interest  ????????????????
# -> TIBCO: stored as a column ????????????????
# -> R console: stored as a variable ?????????????????

# NOTE: Both input tables come from the script 4.c.

# Output : a dataframe with an annotation of interest

############# [INFO] Previous requirement ###############
# This script will only work with specifics variables stored as R objects / TIBCO tables
# and with some libraries (see input section)

############## [INFO] Problems ############################
# Running this script throw the following error :
# Warning messages:
#   1: In data.frame(..., check.names = FALSE) :
#   row names were found from a short variable and have been discarded

##### Method to loop a column containing a LIST ######
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

####  Method to loop a column containing a DATAFRAME  ###########
# Input: a dataframe
getClassesInDataframe <- function(extracted_df) { 
  
  # Get the list of classes
  columnClass <- lapply(extracted_df, class)
  
  # extract_list contains header column and classes
  extractClasses <- data.frame(columnClass)
  return (extractClasses)
}

########## test variables #########
annotation <- "consequenceTypes"
subAnnotation <- "sequenceOntologyTerms"
specificAnnot_column <- annotVCFTable$consequenceTypes

########### [Code] Method to create a table from a specific annotation ########
getSpecAnnotTable <- function(subAnnotation, specificAnnot_column) {
  ############# [Code] Load libraries ############# 
  # library(RCurl)
  library(jsonlite) #
  library(dplyr)
  # library(tidyr)
  
  class_1 <- class(specificAnnot_column)
  
  if (class_1 == "data.frame") { # As far as I saw, there are no dataframe in Level 1
    print ("ERROR: Second level class: dataframe")
    
  } else if (class_1 == "list") { # Each element in list correspond to a variant + annotations
    print (" + + Looping SECOND level class LIST")
    
    if (annotation == "repeat.") {
      specificAnnot_column <- annotVCFTable["repeat"]
    }
    
    # Since the extracted column was a list,
    # each element of the list corresponds to a variant whith his annotations
    # So we will extract only first variant to see what is the class inside
    class_2 <- getClassesInList(specificAnnot_column)
    
    if (class_2[[1]] == "data.frame") { # For each variant, there is a dataframe
      print(" + + + Looping THIRD level class DATAFRAME")
      dataframe_2 <- class_2[[2]]
      class_3 <- getClassesInDataframe(dataframe_2)
      
      for (index in 1:length(class_3)) {
        if (class_3[index] == "list"){
          print(paste("---> + + + + LIST: Need more loops in level 4 for", colnames(class_3[index])))
          class_4 <- getClassesInList(dataframe_2[index])
          
          
        } else if (class_3[index] == "data.frame") {
          print(paste("---> + + + + DATAFRAME: Need more loops in level 4 for", colnames(class_3[index])))
          class_4 <- getClassesInList(dataframe_2[index])
          
        }
      }
      
    } else if (class_2[[2]] == "list")  { # for each variant, there is a list
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



################### TEST ??? ######
tablesToProcess <- list()
annotation <- "sequenceOntologyTerms"

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {

    ########### [TIBCO] Load RinR library ###################
    library(RinR)
    ########### [TIBCO] Determinate R interpreter location ########
    Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
    
    ############ [TIBCO] Load variables ###############
    # Load the annotated table that will be post-processed
    annotVCFTable <- consequenceTypesTable
    # add the table to post-process to a list 
    # tablesToProcess <-
    
    ########### [TIBCO] Create the REvaluate object ########
    tablesList <- REvaluate({
      # tablesList <- getTablesList(tablesToProcess, specificAnnot_column)
      # tablesList
      annotation_table <- getSpecAnnotTable(annotation, specificAnnot_column)
    },
    data = list(getSpecAnnotTable = getSpecAnnotTable, getTablesList = getTablesList, specificAnnot_column = annotVCFTable, tablesToProcess = tablesToProcess)
    # ,
    # REvaluator = Rversion,
    # verbose	= TRUE
    )
    
    # for (i in (1:length(tablesList))) {
    #   # Load each table in a specific variable in global environment
    #   table_name <- paste(tablesToProcess[i],"Table", sep="")
    #   table_extracted <- tablesList[[i]]
    #   assign(table_name, table_extracted)
    #   }
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Set the annoation which needs a post-processing ############
  # Load the annotated table that will be post-processed
  specificAnnot_column <- annotVCFTable$consequenceTypes
  subAnnotation <- "sequenceOntologyTerms"

  ########### [RStudio] Execute main method ###########
  # tablesList <- getTablesList(tablesToProcess, annotVCFTable)
  annotation_table <- getSpecAnnotTable(subAnnotation, specificAnnot_column)
  
  # ########### [RStudio] Save the information ###########
  # for (i in (1:length(tablesList))) {
  #   # Load each table in a specific variable in global environment
  #   table_name <- paste(tablesToProcess[i],"Table", sep="")
  #   table_extracted <- tablesList[[i]]
  #   assign(table_name, table_extracted)
  #   # Print the table in a txt file
  #   # Works only with basic tables
  #   file_path <- paste("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\test_files\\annotatedTables\\annotated_table_",tablesToProcess[i], ".txt",sep = "")
  #   try(write.table(table_extracted,file_path, append = FALSE, sep="\t",row.names=FALSE))
  # }
  
}
########### [TESTING] Control classes of the chosen annotation table #######
# Just to be sure this table is suitable for TIBCO
TEST_classes <- lapply(subAnnotation_table, class)
View(annot_cell)