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

# There is another mystery in a tesing loop thath I do not understand.
# Please jump to testing block 2

# ############# [RStudio] Set Working directory ############
# Please uncoment the next line changin the working directory by the correct one:
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

########### [RStudio] Get the variants table ############
variant_table <- read.table("test_files\\variants_table.txt", header=TRUE) # for RStudio

############# [TIBCO] Load RinR library ###################
library(RinR)

########### [TIBCO] Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

############# [Code] Build the GET URL and query CellBase (CellBaseR) ############
########### [TIBCO] Create the REvaluate object ########
preVariable <- REvaluate({
  # Load libraries
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  library(dplyr)
  library(tidyr)
  
  var_number <- nrow(variant_table)
  cb <- CellBaseR()
  
  
  annotVariants_table <- data.frame()
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
    variant <- paste(var_chrom, ":", var_range, ":", var_refAl, ":", var_altAl, sep = "")
    annotVariant <- getVariant(object=cb, ids=variant, resource="annotation")
    
    # set the vector of column dimension to 0
    columns_dim <- integer()
    
    if (nrow(annotVariant)==0) { # the result is an empty data.frame
      print (paste("WARNING for variant", i, "EMPTY RESULT!"))
      annotVariant <- data.frame(var_chrom, as.integer(var_range), var_refAl, var_altAl)
      colnames(annotVariant) <- c("chromosome", "start", "reference", "alternate")
      annotVariants_table <- bind_rows(annotVariants_table, annotVariant)
      
    } else {
      
      for (column in annotVariant) {
        columns_dim <- c(columns_dim, dim(column))
      }
      
      if (length(columns_dim) == 0) { # normal situation, no column with dimension > 0
        annotVariants_table <- bind_rows(annotVariants_table, annotVariant)
        
      } else { # anormal situation, there is one or more columns with dimension > 0
        print (paste("WARNING for variant", i, "MULTIDIMENSIONAL element(s)"))
        # annotVariants_table <- bind_rows(annotVariants_table, annotVariant)

        reducedTable <- annotVariant[c(1:9)] # testing line
        TEST_completeTableDIM <- annotVariant # testing line
        # variantTraitAssociation_df2 <- annotVariant$variantTraitAssociation
        annotVariants_table <- bind_rows(annotVariants_table, reducedTable) # testing line
        
        ############ testing block 2 #############
        print (paste("Analyzing columns of variant", i))
        k = 1
        columns_class <- character()
        whatthefuck <- character()
        columns_dim_problem <- character()
        for (column in annotVariant) {
          
          whatthefuck <- c(columns_class, dim(column))
          columns_dim_problem <- c(columns_dim_problem, dim(column))
          columns_class <- c(columns_class, class(column))
          
          if (length(dim(column)) != 0){
            print (paste("Dimension for column", k, "is the following:", dim(column)))
            print (paste("Class for column", k, "is the following:", class(column)))
            print (paste("Column", k, "has a problem !!!!!!!!!!!!"))
            variant = i
            column = k
            class = class(column)
            dimension = dim(column)
            # one_problem <- c(variant, column, dimension, class)
            one_problem <- data.frame(variant, column, dimension, class)
            # one_problem <- data.frame(one_problem)
            # one_problem <- t(one_problem)
            # colnames(one_problem) <- c("variant", "column", "dimension", "class")
            print(one_problem)
            problems <- bind_rows(problems, one_problem)
          }
          
          k = k + 1
          
        } 
        print (paste("All the dimensions here: ", columns_dim_problem))
        print (paste("All the class here: ", columns_class))
        # end of testing block 2
      }
    }
    # warnings()
  }
  annotVariants_table
},
data = list(variant_table = variants_table.txt)
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)
getVariantTable <- preVariable

########### Testing block 2 ############
# Testing loop to analyze what there is inside the table 

columns_dim_t <- character()
j <- 0
for (column in completeTable2) {
  print (paste("Column", j, "has a class: ", class(column)))
  print (paste("Dimension of the column:", dim(column)))
  columns_dim_t <- c(columns_dim_t, dim(column))
  if (class(column) == "data.frame") {
    print ("WARNING, THE FOLLWING COLUMN IS A DATAFRAME")
    print (paste("Column", j))
    print ("---------------------------")
    print (column)
  }
  print ("###################")
  j <- j + 1
}

columns_dim_t

if (length(columns_dim_t) == 0) {
  print ("empty")
} else {
  print ("not empty")
}
  
column

########### [RStudio] Print the table in a txt file ###########
try(write.table(annotVariants_table,"test_files\\CB_variants_table.txt", append = TRUE, sep="\t",row.names=FALSE))
