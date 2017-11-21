############## [INFO] SCRIPT 3 General information ############
# We provide a table with some variants, then the script
# will query the RESTfull DB CellBase to retreive information about these variants

############## [INFO] Input and Output ############################
# Inputs : A variant table coming from a VCF file (script VCF_get_variant)
# Output : a preprocessed annotated table + a table with rs annotation
# + table with available annotations

# IMPORTANT NOTE : the OUTPUT table cannot be send to a txt file neither to a TIBCO table
# The table needs to be post-processed.

############# [INFO] Previous requirement (libraries) ###############
# LIBRARIES from CRAN : RCurl, jsonlite
# LIBRARY to use this script in TIBCO : RinR
# LIBRARY from Bioconductor : "cellbaseR":
# https://bioconductor.org/packages/release/bioc/html/cellbaseR.html

############## [INFO] Problems ############################
# Problem 1:  Some variants have more than 1 REF or ALT allele
# so that this script do the annotation call with the first allele of each field.
# I will change that just deleting these variants

# Problem 2: The annotated table has sometimes columns with dataframe inside
# this has been fixed with additional code, but it could slow down the script...

# Problem 3 : Some warnings messages when running the entire scripts :
# "There were 50 or more warnings (use warnings() to see the first 50)
# 1: In bind_rows_(x, .id) : Unequal factor levels: coercing to character
# 2: In bind_rows_(x, .id) :
#   binding character and factor vector, coercing into character vector

########### [Code] Main method ########
getAnnotVariantsTable <- function(variants_table) {
  ############# [Code] Load libraries ############# 
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  library(dplyr)
  library(tidyr)
  
  ############# [Code] Build the GET URL and query CellBase (CellBaseR) ############
  var_number <- nrow(variants_table)
  cb <- CellBaseR()
  
  # Initialize the annotation table
  annotVariants_table <- data.frame()
  # Initialize a table to record any problem about dimensions in the annotated table
  problems <- data.frame(variant=integer(0), column=integer(0), dimension=integer(0), class=character(0)) # testing line
  
  for (i in 1:var_number) { 
    print (paste("Processing variant number:", i)) # this line is for testing
    # extract the chromosome
    var_chrom <- variants_table[i,1]
    
    # extract the range
    var_range <- variants_table[i,2]
    
    # extract the ref and alt alleles
    # WARNING: you could have more than one allele in each field
    # so that this formula extract only the first one to do the annotation call
    # IS THAT A GOOD APPROACH ?
    var_refAl <- substring((variants_table[i,4]), 1, 1)
    var_altAl <- substring((variants_table[i,5]), 1, 1)
    
    # Get variant cellbase info with cellbaseR package
    # the call will return as a data.frame
    variant <- paste(var_chrom, ":", var_range, ":", var_refAl, ":", var_altAl, sep = "")
    annotVariant <- getVariant(object=cb, ids=variant, resource="annotation")
    
    # set the vector of column dimension to 0
    # This is important to control calls that return data.frame with dimension problems
    columns_dim <- integer()
    
    if (nrow(annotVariant)==0) { # if the call returns an empty data.frame (no annotations)
      print (paste("WARNING for variant", i, "EMPTY RESULT!"))
      
      # build the data.frame with the available information
      annotVariant <- data.frame(var_chrom, as.integer(var_range), var_refAl, var_altAl)
      colnames(annotVariant) <- c("chromosome", "start", "reference", "alternate")
      annotVariants_table <- bind_rows(annotVariants_table, annotVariant)
      
    } else { # the call returns annotations
      
      # control the dimension of each column in the dataframe
      for (column in annotVariant) { 
        columns_dim <- c(columns_dim, dim(column))
      }
      
      # normal situation, no column with dimension > 1
      if (length(columns_dim) == 0) { 
        annotVariants_table <- bind_rows(annotVariants_table, annotVariant)
        
      } else {
        # anormal situation, there is one or more columns with dimension > 0
        print (paste("WARNING for variant", i, "MULTIDIMENSIONAL column(s)"))
        
        # assign the data.frame to a new variable
        annotVariant_sub <- annotVariant
        
        # Scan the data.frame in order to find the problematic column(s)
        print (paste("Analyzing columns of variant", i)) # testing line
        k = 1
        for (column in annotVariant) {
          # extract class and dimension of the problematic column
          the_class <- class(column)
          the_dim <- dim(column)
          if (length(the_dim) != 0){ # only for columns with problem
            print (paste("Column", k, "has a dimensional problem.")) # testing line
            
            # Remove the problematic column, split it and add it again to the dataframe
            annotVariant_sub <- cbind(annotVariant[,-k], column[1], column[2])
            
            # To record variants with dimensional problem
            variant = i
            column = k
            class = the_class
            dimension = the_dim
            one_problem <- data.frame(variant, column, dimension, class)
            problems <- bind_rows(problems, one_problem)
            
          }
          k = k + 1
        } 
        # Build the table with all the annotation
        annotVariants_table <- bind_rows(annotVariants_table, annotVariant_sub)
      }
    }
    # warnings()
  }
  # Change the name of the "id" column to avoid repetition
  colnames(annotVariants_table)[5] <- "rsID"
  return(annotVariants_table)
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
  annotVariants_table <- REvaluate({
    annotVariants_table <- getAnnotVariantsTable(variants_table)
    annotVariants_table
  }
  , data = list(getAnnotVariantsTable = getAnnotVariantsTable, variants_table = variants_table)
  # , REvaluator = Rversion
  # , verbose	= TRUE
  )
  
  ########### [TIBCO] Comvert the table to a Blob Object ########
  annotVariantsBlob <- SObjectToBlob(annotVariants_table)
  
} else {
  ########### [RStudio] Set Working directory ############
  setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")
  
  ########### [RStudio] Get the variants table ############
  variants_table <- read.table("test_files\\variants_table.txt", header=TRUE) # for RStudio
  
  ########### [RStudio] Execute main method ###########
  annotVariants_table <- getAnnotVariantsTable(variants_table)
  
  ########### [RStudio] Print the basic table and available annotations in a txt file ###########
  # Works only with basic table
  try(write.table(annotVariants_table[,1:5],"test_files\\CB_variants_table.txt", append = FALSE, sep="\t",row.names=FALSE))
  try(write.table(colnames(annotVariants_table[6:length(annotVariants_table)]),"test_files\\available_annotations.txt", append = FALSE, sep="\t",row.names=FALSE))
}

########### [Code] Build a basic table that could be loaded whithout problems ########
# This table has no nested data.frame neither nested list 
basicTable <- annotVariants_table[,1:5]

########### [Code] Get the available annotations ############
availableAnnots <- colnames(annotVariants_table[6:length(annotVariants_table)])
