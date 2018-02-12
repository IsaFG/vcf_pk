########### [INFO] Script 6 General information ############
# Script to unify variable formats inside an annotation column.
# This script applies to clinvarTable, where the values are strings
# with different cases, unnecesary comas, unnecessary underscores and quotation marks.
# Passing this script aims to unify formating

########### [INFO] Input and Output ############################
# Inputs : clinvarTable
# Output : clinvarSignificance  + clinvarDisease
## where format has been unified and duplicated values have been removed

########### [INFO] Previous requirement (libraries) ###############
# clinvarTable loaded in R environment as a datatable

########### [INFO] Problems ###############
# Unable to find a solution to remove quotation marks
# when running the script in TERR (Spotfire)

########### [Code] Method to remove quotation marks ########
# NOTE: removing quotation marks does not work in TERR 
cleanQuotationMarks <- function(tableToClean, columnToClean) {
  ########### [Code] Determinate if running in TERR or standard R version #############
  isTERR<-R.Version()
  Rversion<-NULL
  if (!is.null(isTERR[["TERR.version"]])) {
    
    # Pending to find a solution
    tableToClean
    
  } else {
    ########### [RStudio] Execute method ###########
    # Uncomment when running in standard R version (no TERR)
    # tableToClean[[paste(columnToClean)]] <- gsub( "\"", "", tableToClean[[paste(columnToClean)]])
  }
  return(tableToClean)
}

######## Method whith all the steps to unify formating and remove duplicates ###########
cleanUselessChr <- function(tableToClean, columnToClean){
  tableToClean[[paste(columnToClean)]] <- gsub( ", ", ",", tableToClean[[paste(columnToClean)]])
  tableToClean[[paste(columnToClean)]] <-  gsub( "_", " ", tableToClean[[paste(columnToClean)]])
  tableToClean[[paste(columnToClean)]] <- gsub( ",", ";", tableToClean[[paste(columnToClean)]])
  tableToClean <- cleanQuotationMarks(tableToClean, columnToClean)
  tableToClean[[paste(columnToClean)]] <- gsub("^c\\(|\\)$", "", tableToClean[[paste(columnToClean)]])
  tableToClean[[paste(columnToClean)]] <- tolower(tableToClean[[paste(columnToClean)]])
  tableToClean <- unique(tableToClean)
  return(tableToClean)
}

######## Create 2 new tables with unified format ###################
# These 2 tables will be used for different purpose
# Table 1 : clinvarSignificance
clinvarSignificance <- clinvarTable
clinvarSignificance$reviewStatus <- gsub( ",", "", clinvarSignificance$reviewStatus)
clinvarSignificance <- cleanUselessChr(clinvarSignificance, "reviewStatus")

# Table 2 : clinvarDisease
# Note: input clinvarTable_dis just exists in Spotfire.
clinvarDisease <- clinvarTable_dis # Comment if running script in standard R version (no TERR)
# clinvarDisease <- clinvarTable # Uncomment when running in standard R version (no TERR)
clinvarDisease$reviewStatus <- gsub( ",", "", clinvarDisease$reviewStatus)
clinvarDisease <- cleanUselessChr(clinvarDisease, "reviewStatus")
clinvarDisease <- cleanUselessChr(clinvarDisease, "geneNames")
clinvarDisease <- cleanUselessChr(clinvarDisease, "traits")