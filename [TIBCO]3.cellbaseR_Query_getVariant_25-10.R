############## [INFO] SCRIPT 3 General information ############
# TIBCO script
# We provide a table with some variants, then the script
# will query the RESTfull DB CellBase to retreive information about these variants

############## [INFO] Input and Output ############################
# Inputs : A variant table coming from a VCF file (script VCF_get_variant)
# Output : annotated table

############# [INFO] Previous requirement (libraries) ###############
# LIBRARIES from CRAN : RCurl, jsonlite
# LIBRARY to use this script in TIBCO : RinR
# LIBRARY from Bioconductor : "cellbaseR":
# https://bioconductor.org/packages/release/bioc/html/cellbaseR.html

############## [INFO] Problems ############################
# Some variants have more than 1 REF or ALT allele
# so that this script do the annotation call with the first allele of each field.
# IS THAT A GOOD APPROACH?

# The annotated table has sometimes columns with dataframe inside
# this has been fixed with additional code, but it could slow down the script...

# ############# [RStudio] Set Working directory ############
# Please uncoment the next line changing the working directory by the correct one:
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

########### [RStudio] Get the variants table ############
variant_table <- read.table("test_files\\variants_table.txt", header=TRUE) # for RStudio

############# [TIBCO] Load RinR library ###################
library(RinR)

########### [TIBCO] Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### [TIBCO] Create the REvaluate object ########
annotatedTable <- REvaluate({
  ############# [Code] Load libraries ############# 
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  library(dplyr)
  library(tidyr)
  
  ############# [Code] Build the GET URL and query CellBase (CellBaseR) ############
  var_number <- nrow(variant_table)
  cb <- CellBaseR()
  
  # Initialize the annotation table
  annotVariants_table <- data.frame()
  # Initialize a table to record any problem about dimensions in the annotated table
  problems <- data.frame(variant=integer(0), column=integer(0), dimension=integer(0), class=character(0)) # testing line
  
  for (i in 1:var_number) { 
    print (paste("Processing variant number:", i)) # this line is for testing
    # extract the chromosome
    var_chrom <- variant_table[i,1]
    
    # extract the range
    var_range <- variant_table[i,2]
    
    # extract the ref and alt alleles
    # WARNING: you could have more than one allele in each field
    # so that this formula extract only the first one to do the annotation call
    # IS THAT A GOOD APPROACH ?
    var_refAl <- substring((variant_table[i,4]), 1, 1)
    var_altAl <- substring((variant_table[i,5]), 1, 1)
    
    # Get variant cellbase info with cellbaseR package
    # the call will return as a data.frame
    variant <- paste(var_chrom, ":", var_range, ":", var_refAl, ":", var_altAl, sep = "")
    annotVariant <- getVariant(object=cb, ids=variant, resource="annotation")
    
    # set the vector of column dimension to 0
    # This is important to control calls that return data.frame with dimension problems
    columns_dim <- integer()
    
    if (nrow(annotVariant)==0) { # if the call returns an empty data.frame
      print (paste("WARNING for variant", i, "EMPTY RESULT!"))
      
      # build the data.frame with the available information
      annotVariant <- data.frame(var_chrom, as.integer(var_range), var_refAl, var_altAl)
      colnames(annotVariant) <- c("chromosome", "start", "reference", "alternate")
      annotVariants_table <- bind_rows(annotVariants_table, annotVariant)
      
    } else {
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
          # extract class and dimesion of the problematic column
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
        # end of testing block 2
        annotVariants_table <- bind_rows(annotVariants_table, annotVariant_sub)
      }
    }
    # warnings()
  }
},
data = list(variant_table = variants_table.txt)
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)

########### [RStudio] Print the table in a txt file ###########
# Gives an error :
# Error in write.table(annotVariants_table, "test_files\\CB_variants_table.txt",  : 
# unimplemented type 'list' in 'EncodeElement'
try(write.table(annotVariants_table,"test_files\\CB_variants_table.txt", append = TRUE, sep="\t",row.names=FALSE))




