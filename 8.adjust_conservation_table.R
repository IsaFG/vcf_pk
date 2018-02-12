########### [INFO] Script 9 General information ############
# This script will adjust the conservation scores table
# to unify data in order to have one row per variant

########### [INFO] Input and Output ############################
# Inputs : conservationTable, with one row per conservation score type,
# # which means three rows per variant
# Output : adjusted conservationTable, with only one row per variant

########### [INFO] Previous requirement (libraries) ###############
# Library : reshape2

# Working directory :
# Please uncoment the next line changing the working directory by the correct one:
# working_dir <- "C:\\..."

# R version : please uncoment the next line indicating the location of your current R version:
# r_version <- "C:/Program Files/R/R-3.4.1"

########### [INFO] Problems ############################
# No problem found so far

########### [Code] Main method ########
adjustConservationTable <- function(conservationTable) {
  
  library(reshape2)
  
  conScores_table <- dcast(conservationTable, chromosome + start + reference + alternate + rsID ~ source, value.var = "score", fun.aggregate = mean)

  return(conScores_table)
}

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  ########### [TIBCO] Load RinR library ###################
  library(RinR)
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = r_version)
  
  ########### [TIBCO] Create the REvaluate object to execute main method ########
  predictionScoresTable <- REvaluate({
    conScores_table <- adjustConservationTable(conservationTable)
  }
  , data = list(adjustConservationTable = adjustConservationTable, conservationTable = conservationTable)
  )
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd(working_dir)
  
  ########### [RStudio] Execute main method ###########
  predictionScoresTable <- adjustConservationTable(conservationTable)
  
  ########### [RStudio] Print the ouput in a txt file ###########
  
  print ("Done!")
}
