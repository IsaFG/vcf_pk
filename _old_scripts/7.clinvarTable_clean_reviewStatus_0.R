########### [INFO] Script 7 General information ############
# Script to unify the values inside an annotation column.
# This script applies to clinvarTable, specificaly to the column "reviewStatus"
# where the values are strings with different case, unnecesary comas and unnecessary underscores.

########### [INFO] Input and Output ############################
# Inputs : clinvarTable
# Output : clinvarTableRS,
## where case have been unified to lower and comas and underscore have been removed

########### [INFO] Previous requirement (libraries) ###############
# None

########### [Code] Main method ########

clinvarTableRS <- clinvarTable

clinvarTableRS$reviewStatus <- gsub( "[_]", " ", clinvarTableRS$reviewStatus)
clinvarTableRS$reviewStatus <- gsub( "[,]", "", clinvarTableRS$reviewStatus)
clinvarTableRS$reviewStatus <- tolower(clinvarTableRS$reviewStatus)

clinvarSignificance <- clinvarTableRS
clinvarSignificance$accession <- NULL
clinvarSignificance <- unique(clinvarSignificance)
