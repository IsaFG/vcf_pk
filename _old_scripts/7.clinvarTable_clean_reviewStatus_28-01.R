########### [INFO] Script 7 General information ############
# Script to unify the values inside an annotation column.
# This script applies to clinvarTable, specificaly to the column "reviewStatus"
# where the values are strings with different case, unnecesary comas and unnecessary underscores.

########### [INFO] Input and Output ############################
# Inputs : clinvarTable
# Output : clinvarSignificance  + clinvarDisease
## where case have been unified to lower, comas and underscore have been removed
# and duplicated have been removed

########### [INFO] Previous requirement (libraries) ###############
# None

########### [Code] Main method ########

cleanUselessChr <- function(tableToClean, columnToClean){
  tableToClean[[paste(columnToClean)]] <- gsub( ", ", ",", tableToClean[[paste(columnToClean)]])
  tableToClean[[paste(columnToClean)]] <-  gsub( "_", " ", tableToClean[[paste(columnToClean)]])
  tableToClean[[paste(columnToClean)]] <- gsub( ",", ";", tableToClean[[paste(columnToClean)]])
  # tableToClean[[paste(columnToClean)]] <- gsub( "\"", "", tableToClean[[paste(columnToClean)]]) # esto no funciona en spotfire
  tableToClean[[paste(columnToClean)]] <- gsub("^c\\(|\\)$", "", tableToClean[[paste(columnToClean)]])
  tableToClean[[paste(columnToClean)]] <- tolower(tableToClean[[paste(columnToClean)]])
  tableToClean <- unique(tableToClean)
  return(tableToClean)
}

clinvarSignificance <- clinvarTable
clinvarSignificance$reviewStatus <- gsub( ",", "", clinvarSignificance$reviewStatus)
clinvarSignificance <- cleanUselessChr(clinvarSignificance, "reviewStatus")

clinvarDisease <- clinvarTable_dis
clinvarDisease$reviewStatus <- gsub( ",", "", clinvarDisease$reviewStatus)
clinvarDisease <- cleanUselessChr(clinvarDisease, "reviewStatus")
clinvarDisease <- cleanUselessChr(clinvarDisease, "geneNames")
clinvarDisease <- cleanUselessChr(clinvarDisease, "traits")