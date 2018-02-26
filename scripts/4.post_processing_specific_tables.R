############## [INFO] Script 4 General information ############
# Some annotation tables show much information.They can be loaded in TIBCO,
# but use to have nested lists, which make the info hard to read.
# This script generate additional sub-tables for some interesting information.

############## [INFO] Input and Output ############################
# Input 1: a variant annotated table that needs post-processing
# -> So far, this method has only been used for "consequenceTypes" and "traitAssociation"
# -> In TIBCO, this table is not readable and is stored as a Document propertie
# -> In R console, this table should be present in global environment as a data.frame
# Input 2: Annotation and his sub-annotation of interest

# NOTE: The input tables come from the script 3.extract_all_annot.R

# Output : a dataframe with the sub-annotation of interest

############# [INFO] Previous requirement ###############
# This script will only work with specifics variables stored as R objects / TIBCO tables
# and with some libraries (see input section)

# Working directory :
# Please uncoment the next line changing the working directory by the correct one:
# working_dir <- "C:\\..."

# R version : please uncoment the next line indicating the location of your current R version:
# r_version <- "C:/Program Files/R/R-3.4.1"

############## [INFO] Problems ############################
# Running this script throw the following error :
# Warning messages:
#   1: In data.frame(..., check.names = FALSE) :
#   row names were found from a short variable and have been discarded

########### [Code] Method to create a table from a specific annotation ########
getSpecAnnotTable <- function(annotation, subAnnotation, annotVCFTable) {
  ############# [Code] Load libraries ############# 
  # library(RCurl)
  library(jsonlite) #
  library(dplyr)
  # library(tidyr)
  
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
    # LEVEL 3
    # Note : level 2 class will almost always be a dataframe (except for repeat)
    if (element2_class == "data.frame") {
      i <- 1
      for (element3 in element2){
        element3_class <- class(element3)
        if (colnames(element2)[i] == subAnnotation) {
          element4_df <- data.frame()
          for (element4 in element3) {
            if (is.null(element4)) {
              element4_df <- basic_df[j,]
            } else {element4_df <- cbind(basic_df[j,], element4)
            }
            element4_df_class <- lapply(element4_df, class)
            k <- 1
            for (element5 in element4_df) {
              if (class(element5) == "data.frame"){
                element4_df[7] <- as.character(element4_df[7])
              }
              k <- k + 1
            }
            subAnnotation_df <- bind_rows(subAnnotation_df, element4_df)
          }
        }
        i <- i + 1
      }
      
      # NOT DONE FOR LIST SINCE ONLY ONE ANNOTATION TABLE HAS THIS SITUATION
    } else if (element2_class == "list")  {
      for (element3 in element2){
        element3_class <- class(element3)
        print (element3_class)
      }
    }
    j <- j + 1      
  }
  return(subAnnotation_df)
}

########### [Code] Set the annotation which needs a post-processing ############
# Pre-set the sub-annotations that you want to post-process
annotation <- "consequenceTypes"
subAnnotation <- "sequenceOntologyTerms"

## Other table that can be post-processed with this method
# annotation <- "traitAssociation"
# subAnnotation <- "additionalProperties"
# subAnnotation <- "heritableTraits"

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {

    ########### [TIBCO] Load RinR library ###################
    library(RinR)
    ########### [TIBCO] Determinate R interpreter location ########
    Rversion <- makeREvaluator("R", RHome = r_version)
    
    ############ [TIBCO] Load variables ###############
    # Load the annotated table that will be post-processed
    # Load the annotated table, stored in memory as a Blob Object created by script "cellbaseR_Query_getVariant"
    annotVCFObject <- BlobToSObject(annotVCFTableBlob)
    annotVCFTable <- annotVCFObject
    
    ########### [TIBCO] Create the REvaluate object ########
    subTable <- REvaluate({
      # tablesList <- getTablesList(tablesToProcess, specificAnnot_column)
      # tablesList
      annotation_table <- getSpecAnnotTable(annotation, subAnnotation, annotVCFTable)
    },
    data = list(getSpecAnnotTable = getSpecAnnotTable, annotation = annotation, subAnnotation = subAnnotation, annotVCFTable = annotVCFTable)
    # ,
    # REvaluator = Rversion,
    # verbose	= TRUE
    )
    
    # Load each table in a specific variable in global environment
    table_name <- paste(subAnnotation,"Table", sep="")
    assign(table_name, subTable)

} else {
  ########### [RStudio] Set Working directory ############
  setwd(working_dir)

  ########### [RStudio] Execute main method ###########
  # tablesList <- getTablesList(tablesToProcess, annotVCFTable)
  subAnnotation_df <- getSpecAnnotTable(annotation, subAnnotation, annotVCFTable)
  
 ########### [RStudio] Save the information ###########

  # Load each table in a specific variable in global environment
  table_name <- paste(annotation,"_", subAnnotation,"Table", sep="")
  assign(table_name, subAnnotation_df)
  # Print the table in a txt file
  # Works only with basic tables
  file_path <- paste("output_files\\annotatedTables\\annotated_table_",annotation,subAnnotation, ".txt",sep = "")
  try(write.table(subAnnotation_df,file_path, append = FALSE, sep="\t",row.names=FALSE))

}
########### [TESTING] Control classes of the chosen annotation table #######
# Just to be sure this table is suitable for TIBCO
TEST_classes <- lapply(subAnnotation_df, class)