########### [INFO] Script 5.c General information ############
# This script will adjust the CADD score table to unify data for one variant in only one row

########### [INFO] Input and Output ############################
# Inputs : Variant table with CADD annotacion scores : functionalScoreTable
# # with one row per CADD score type, which means three rows per variant
# Output : caddScoresTable

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
adjustcaddScoresTable <- function(functionalScoreTable) {
  library(reshape2)
  
  caddScores_table <- dcast(functionalScoreTable, chromosome + start + reference + alternate + rsID ~ source, value.var = "score", fun.aggregate = mean)

  return(caddScores_table)
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
  caddScoresTable <- REvaluate({
    caddScores_table <- adjustcaddScoresTable(functionalScoreTable)
  }
  , data = list(adjustcaddScoresTable = adjustcaddScoresTable, functionalScoreTable = functionalScoreTable)
  )
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd(working_dir)
  
  ########### [RStudio] Execute main method ###########
  caddScoresTable <- adjustcaddScoresTable(functionalScoreTable)
  
  ########### [RStudio] Print the ouput in a txt file ###########
  try(write.table(caddScoresTable,"test_files\\annotatedTables\\annotated_table_caddScores.txt", append = FALSE, sep="\t",row.names=FALSE))
  print ("Done!")
}
