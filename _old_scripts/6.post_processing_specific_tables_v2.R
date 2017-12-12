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


########## test variables #########
annotation <- "consequenceTypes"
subAnnotation <- "sequenceOntologyTerms"
subAnnotation_class <- class(annotVCFTable$consequenceTypes)
element1 <- annotVCFTable$consequenceTypes

# Level 1 
# -> column of a dataframe, where each cells (variant) contains one list
# So class is simply "list"
print(class(element1))

# Level 2
# -> each variant list info will contain or a dataframe either a list
for (element2 in element1){
  print(class(element2))
}

# Level 3
# For dataframe inside a single element (variant)
for (element3 in element2) {
  print (class(element3))
  # print (head(element3))
}

for (element4 in element3) {
  print (class(element4))
  # print (head(element4))
}


########### [Code] Method to create a table from a specific annotation ########
getSpecAnnotTable <- function(annotation, subAnnotation, annotVCFTable) {
  ############# [Code] Load libraries ############# 
  # library(RCurl)
  library(jsonlite) #
  library(dplyr)
  # library(tidyr)
  
  # LEVEL 1
  element1 <- annotVCFTable$consequenceTypes
  # The class of the element will ALWAYS be "list"
  
  # LEVEL 2
  element2_class <- character()
  for (element2 in element1) { # for each element in list
    if ((length(element2_class)) == 0) {
      if (length(element2) != 0) { # element that contains information
        element2_class <- class(element1[[1]])
      }
    }
  }
  
  
  # LEVEL 3
  if (element2_class == "data.frame") {
    
    
    
  } else if (element2_class == "list")  {
    
    
  }
  
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