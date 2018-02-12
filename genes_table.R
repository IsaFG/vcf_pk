########### [INFO] SCRIPT x General information ############
# 

########### [INFO] Input and Output ############################
# Inputs : 
# Output : 

########### [INFO] Previous requirement (libraries) ###############

########### [INFO] Problems ############################

########### [Code] Main methods #############
genesCrude <- merge(consequenceTypesTable, variants_table, by="variantID", all.x=TRUE)
genesCrude <- merge(genesCrude, mostSevereConsequence, by="variantID", all.x=TRUE)
genesMutationsCount <- unique(genesCrude)
genesMutationsCount <- genesMutationsCount[!(is.na(genesMutationsCount$geneName) | genesMutationsCount$geneName==""), ]
genesMutationsCount <- unique(genesMutationsCount)