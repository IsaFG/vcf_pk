########### [INFO] Script 8 General information ############
# Subset the populationFrecuenciesTable according to study and population chosen by user
# NOTE : in R standard version, the study and population are pre-set (ALL and EXAC)
# In Spotfire, the user can choose them through the platform

########### [INFO] Input and Output ############################
# Inputs : populationFrecuenciesTable + chosen pop + chosen study
# Output : a subset of the populationFrecuenciesTable

########### [INFO] Previous requirement (libraries) ###############
# populationFrecuenciesTable loaded in R environment as a datatable

########### [INFO] Problems ############################
# None

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {

  ########### [TIBCO] Main method ########
  filteredPopFrec <- populationFrequenciesTable[populationFrequenciesTable$population == chosenPop,]
  filteredPopFrec <- filteredPopFrec[filteredPopFrec$study == chosenPopStudy,]
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Get the input variables ############
  chosenPop <- "ALL"
  chosenPopStudy <- "EXAC"

  ########### [RStudio] Execute main method ###########
  filteredPopFrec <- populationFrequenciesTable[populationFrequenciesTable$population == chosenPop,]
  filteredPopFrec <- filteredPopFrec[filteredPopFrec$study == chosenPopStudy,]
  
}
