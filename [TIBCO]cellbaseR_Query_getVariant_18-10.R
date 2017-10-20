############## [INFO] SCRIPT 3 (TIBCO) General information ############
# This script only works in TIBCO Spotfire
# We provide a table with some variants, then the script will query the RESTfull DB CellBase
# to retreive information about these variants

############# PREVIOUS REQUIREMENT ###############
# LIBRARIES from CRAN : RCurl, jsonlite
# LIBRARY to use this script in TIBCO : RinR
# LIBRARY from Bioconductor : "cellbaseR":
# https://bioconductor.org/packages/release/bioc/html/cellbaseR.html

############## [INFO] Input and Output ############################
# Inputs : A variant table coming from a VCF file (script VCF_get_variant)
# Output : annotated table

############# Load RinR library ###################
library(RinR)

########### Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

############# [CellBaseR] Build the GET URL and query CellBase ############
preVariable <- REvaluate({
  # Load libraries
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  
  var_number <- nrow(variant_table)
  query_vector <- character()
  URL_vector <- character()
  results_list <- list()
  cb <- CellBaseR()
  
  for (i in 1:var_number) { 
    var_chrom <- variant_table[i,1]
    var_range <- variant_table[i,2]
    var_refAl <- variant_table[i,3]
    var_altAl <- variant_table[i,4]
    
    # Get variant cellbase info with cellbaseR package
    variant <- paste(var_chrom, ":", var_range, ":", var_refAl, ":", var_altAl, sep = "")
    res2 <- getVariant(object=cb, ids=variant, resource="annotation")
    res2table <- res2[c("chromosome", "start", "reference", "alternate", "id", "displayConsequenceType")]
    if (i==1) {
      getVariant_table <- res2table
    } else {
      getVariant_table <- (rbind(as.matrix(getVariant_table), as.matrix(res2table)))
    }
    
  }
  getVariant_table
},
data = list(variant_table = variants_table.txt)
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)
getVariantTable <- preVariable