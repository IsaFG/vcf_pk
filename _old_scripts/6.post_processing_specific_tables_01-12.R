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

############ TEST ############
# subAnnotation <- "sequenceOntologyTerms"
# specificAnnot_column <- consequenceTypesTable

subAnnotation <- "sequenceOntologyTerms"
specificAnnot_column <- annotVCFTable$consequenceTypes

########### [Code] Method to create a table from a specific annotation ########
getSpecAnnotTable <- function(subAnnotation, specificAnnot_column) {
  ############# [Code] Load libraries ############# 
  # library(RCurl)
  library(jsonlite) #
  library(dplyr)
  # library(tidyr)
  
  # The specific column can be a list or a dataframe
  # So the action should be different depending on the case :
  
  if (class(specificAnnot_column) == "data.frame") { # In case the cell would contain a dataframe
    print ("dataframe")
    # Extract the header of the dataframe
    allSubAnnnotations <- colnames(specificAnnot_column)
    # Get the class of the subAnnotation you wish to extract
    subAnnotation_class <- class(specificAnnot_column[[paste(subAnnotation)]])
    
  } else if (class(specificAnnot_column) == "list") { # In case the cell would contain a list
    print("list")
    
    # Extract first element of the list to determinate of it is a dataframe or a list
    annot_column <- specificAnnot_column[[1]]
    if (class(annot_column) == "data.frame") {
      # Extract the header of the dataframe
      allSubAnnnotations <- colnames(annot_column)
      # Get the class of the subAnnotation you wish to extract
      subAnnotation_class <- class(annot_column[[paste(subAnnotation)]])
    }
    
    # Loop through the specificAnnot_column as a list
    for (annotation in specificAnnot_column){
      # Get the wished list
      subAnnotation_column <- (annotation[[paste(subAnnotation)]])
    }
    
    
  } else {# In case the cell is not problematic
    # PENDING
    print("It is not a list nor a dataframe")
  }
  print(subAnnotation_class)
  
  # Declare the basic table
  basicTable <- specificAnnot_column[,1:5]
  
  # Declare the dataframe that will contain all the info of the chosen annotation
  subAnnotation_table <- data.frame()
  
  # For each cell of the annotation column, loop to extract the content of the cell
  i = 1
  for (annot_cell in annot_column) {
    # # Print info of the cell being analyzed
    # print (paste("Variant", i, "with class", class(annot_cell), "and lenght:",length(annot_cell)))
    # print (paste("dimension",(dim(annot_cell))))
    
    # Assign variables to start the analyze
    data_class <- class(annot_cell)
    data_length <- length(annot_cell)
    basic_row <- basicTable[i,]
    
    # print (str(annot_column[i])) # testing line
    
    if (data_length == 0) { # In case the cell would have no result for this annotation
      # print ("Processing empty cell")
      annotated_row <- basic_row
      # annotated_row[, specific_annot] <- "No result"
      
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
    
    # print ("Building the annotated table") # testing line
    subAnnotation_table <- bind_rows(subAnnotation_table, annotated_row)
    
    i = i + 1
    # print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")# testing line
  }
  print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")# testing line
  return(subAnnotation_table)
}

########### [Code] Method to create a list with all the tables ########
getTablesList <- function(tablesToProcess, specificAnnot_column){
  tables_list <- list()
  i <- 1
  for (annotation in tablesToProcess){
    print (paste("Building table", annotation))
    annotation_table <- getSpecAnnotTable(annotation, specificAnnot_column)
    tables_list[[i]] <- annotation_table
    i <- i + 1
  }
  return(tables_list)
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