########### [INFO] Script [num] General information ############


########### [INFO] Input and Output ############################
# Inputs : Annotated Table + chosen annotation
# Output : Annotated table renamed

########### [INFO] Previous requirement (libraries) ###############

########### [INFO] Problems ############################

########### [Code] Main method ########
saveAnnotatedTable <- function(AnnotatedTable, loadedAnnotList, chosen_annot) {
  
  ########### [Code] Verify if the annotation has not been already loaded ######
  # Load the vector with all the annotations already loaded
  loaded_annotations <- loadedAnnotList
  
  if (is.element(chosen_annot, loaded_annotations) == FALSE) {
    
    loaded_annotations <- c(loaded_annotations, chosen_annot)
    
    # Rename the annotated table according the chosen annotation
    assign(chosen_annot, AnnotatedTable)
    
    if (chosen_annot == "rsID"){
                            
    } else if (chosen_annot == "hgvs") {
      
    } else if (chosen_annot == "displayConsequenceType") {
      
    } else if (chosen_annot == "consequenceTypes") {
      
    } else if (chosen_annot == "populationFrequencies") {

    } else if (chosen_annot == "conservation") {
      
    } else if (chosen_annot == "geneExpression") {
    
    } else if (chosen_annot == "geneTraitAssociation") {
      
    } else if (chosen_annot == "geneDrugInteraction") {
    
    } else if (chosen_annot == "traitAssociation") {
      
    } else if (chosen_annot == "functionalScore") {
      
    } else if (chosen_annot == "cytoband") {
    
    } else if (chosen_annot == "repeat"  ) {
      
    } else if (chosen_annot == "clinvar") {
      
    } else if (chosen_annot == "cosmic") {
      
    }
  }
}

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  ########### [TIBCO] Load RinR library ###################
  library(RinR)
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
  
  ########### [TIBCO] Get the input variables ############
  global_input <- "whatever"
  
  ########### [TIBCO] Create the REvaluate object to execute main method ########
  mainVariable <- REvaluate({
    mainVariable <- mainMethod(local_input)
  }
  , data = list(mainMethod = mainMethod, local_input = global_input)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  ########### [Code] Method to create a list with table of interest ########
  getTablesList <- function(availableAnnots, annotVariants_table){
    
    rsIDTable <- getSpecAnnotTable("rsID", annotVariants_table)
    
    hgvsTable <- getSpecAnnotTable("hgvs", annotVariants_table)
    
    displayConsequenceTypeTable <- getSpecAnnotTable("displayConsequenceType", annotVariants_table)
    
    consequenceTypesTable <- getSpecAnnotTable("consequenceTypes", annotVariants_table)
    
    populationFrequenciesTable <- getSpecAnnotTable("populationFrequencies", annotVariants_table)
    
    conservationTable <- getSpecAnnotTable("conservation", annotVariants_table)
    
    geneExpressionTable <- getSpecAnnotTable("geneExpression", annotVariants_table)
    
    geneTraitAssociationTable <- getSpecAnnotTable("geneTraitAssociation", annotVariants_table)
    
    geneDrugInteractionTable <- getSpecAnnotTable("geneDrugInteraction", annotVariants_table)
    
    traitAssociationTable <- getSpecAnnotTable("traitAssociation", annotVariants_table)
    
    functionalScoreTable <- getSpecAnnotTable("functionalScore", annotVariants_table)
    
    cytobandTable <- getSpecAnnotTable("cytoband", annotVariants_table)
    
    repeatTable <- getSpecAnnotTable("repeat", annotVariants_table)
    
    clinvarTable <- getSpecAnnotTable("clinvar", annotVariants_table)
    
    cosmicTable <- getSpecAnnotTable("cosmic", annotVariants_table)
    
    tables_list <- list()
  }
  
  
  ########### [Code] Verify if the annotation has not been already loaded ######
  # Load the vector with all the annotations already loaded
  loaded_annotations <- loadedAnnotList
  
  if (is.element(chosen_annot, loaded_annotations) == FALSE) {
    
    loaded_annotations <- c(loaded_annotations, chosen_annot)
    
    # Rename the annotated table according the chosen annotation
    assign(chosen_annot, AnnotatedTable)
    
    if (chosen_annot == "rsID"){
      rsIDTable <- AnnotatedTable
      
    } else if (chosen_annot == "hgvs") {
      hgvsTable <- AnnotatedTable
      
    } else if (chosen_annot == "displayConsequenceType") {
      displayConsequenceTypeTable <- AnnotatedTable
      
    } else if (chosen_annot == "consequenceTypes") {
      consequenceTypesTable <- AnnotatedTable
      
    } else if (chosen_annot == "populationFrequencies") {
      popFrequenciesTable <- AnnotatedTable
      
    } else if (chosen_annot == "conservation") {
      conservationTable <- AnnotatedTable
      
    } else if (chosen_annot == "geneExpression") {
      geneExpressionTable <- AnnotatedTable
      
    } else if (chosen_annot == "geneTraitAssociation") {
      geneTraitAssociationTable <- AnnotatedTable
      
    } else if (chosen_annot == "geneDrugInteraction") {
      geneDrugInteractionTable <- AnnotatedTable
      
    } else if (chosen_annot == "traitAssociation") {
      traitAssociationTable <- AnnotatedTable
      
    } else if (chosen_annot == "functionalScore") {
      functionalScoreTable <- AnnotatedTable
      
    } else if (chosen_annot == "cytoband") {
      cytobandTable <- AnnotatedTable
      
    } else if (chosen_annot == "repeat"  ) {
      repeatTable <- AnnotatedTable
      
    } else if (chosen_annot == "clinvar") {
      clinvarTable <- AnnotatedTable
      
    } else if (chosen_annot == "cosmic") {
      cosmicTable <- AnnotatedTable
    }
  }
  
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Get the input variables ############
  global_input <- "whatever"
  
  ########### [RStudio] Execute main method ###########
  mainVariable <- mainMethod(global_input)
  
  ########### [RStudio] Print the ouput in a txt file ###########
}

########### [TEST block] ###############
