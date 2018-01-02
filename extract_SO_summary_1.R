############## [INFO] Script 4 General information ############
# This method will create a dataframe with all the variants
# that are splice site, missense, nonsense, frameshift. 
# Data come from Sequence ontology annotation.

############## [INFO] Input and Output ############################
# Input : an annotated variant table, variable "sequenceOntology" table
# Output : a dataframe that includes only variants with specific site types

########### [Code] Main method ########
SO_getSummary <- function(SO_table) {
  
  site_type_df <-  data.frame()
  site_type_list <- list()
  
  for(i in 1:nrow(SO_table)) {
    row <- SO_table[i,]
    print (i) # TEST
    # do stuff with row
    if (row$site_type != "") {
      site_type_list$variantID <- as.character(row$variantID)
      
      if (row$is_splice_site == 1) {
        site_type_list$is_splice_site_2 <- 1
        
        site_type_list$is_missense_2 <- NA
        site_type_list$is_nonsense_2 <- NA
        site_type_list$is_frameshift_2 <- NA
        
      } else if (row$is_missense == 1) {
        site_type_list$is_missense_2 <- 1
        
        site_type_list$is_splice_site_2 <- NA
        site_type_list$is_nonsense_2 <- NA
        site_type_list$is_frameshift_2 <- NA
        
      } else if (row$is_nonsense == 1) {
        site_type_list$is_nonsense_2 <- 1
        
        site_type_list$is_splice_site_2 <- NA
        site_type_list$is_missense_2 <- NA
        site_type_list$is_frameshift_2 <- NA
        
      } else if (row$is_frameshift == 1) {
        site_type_list$is_frameshift_2 <- 1
        
        site_type_list$is_splice_site_2 <- NA
        site_type_list$is_missense_2 <- NA
        site_type_list$is_nonsense_2 <- NA
      }
      site_type_df <- rbind(site_type_df, as.data.frame(site_type_list))
    }
  }
  # Keep only unique couples
  uni_site_type_df <- unique(site_type_df)
  return(uni_site_type_df)
}

########### [Code] Determinate if running in TERR or standard R version #############
isTERR<-R.Version()
Rversion<-NULL
if (!is.null(isTERR[["TERR.version"]])) {
  ########### [TIBCO] Load RinR library ###################
  library(RinR)
  ########### [TIBCO] Determinate R interpreter location ########
  Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
  
  ########### [TIBCO] Create the REvaluate object to execute main method ########
  SOSummary <- SO_getSummary(SO_table)
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Get the input variables ############
  loc_file <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\test_files\\annotatedTables\\sequenceOntologyTerms_TIBCOmodified.txt"
  SO_table <- read.table(loc_file, sep="\t", header = T)
  
  ########### [RStudio] Execute main method ###########
  SO_summary <- SO_getSummary(SO_table)
  
  ########### [RStudio] Print the ouput in a txt file ###########
  # Nah
}

