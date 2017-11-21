############## [INFO] Script 3.2 General information ############
# We will try to expand the variant table comming from getVariant method
# coming from script : 3.cellbaseR_Query_getVariant

# THIS VERSION CONTAINS SOME EXPERIMENTS FOR UNNESTING THE TABLE, SO I KEEP IT

############## [INFO] Input and Output ############################
# Inputs : an annotated variant table (created by script : 3.cellbaseR_Query_getVariant)
# Output : expanded annotated variant table

############# [INFO] Previous requirement ###############
# You need to run first 3.cellbaseR_Query_getVariant.
# This script will only work with R object, not with any file

############## [INFO] Problems ############################


# ############# [RStudio] Set Working directory ############
# Please uncoment the next line changing the working directory by the correct one:
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

########### [RStudio] Get the input file ############


# ############# [TIBCO] Load RinR library ###################
# library(RinR)
# 
# ########### [TIBCO] Determinate R interpreter location ########
# Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
# 
# ########### [TIBCO] Create the REvaluate object ########
# annotatedTable <- REvaluate({
  ############# [Code] Load libraries ############# 

  
############# [Code] Main method ############
########### [TEST block] ###############

# The TEST_table comes from a table with several columns containing nested list
# TEST_table <- annotVariant # for faster testing
TEST_table <- annotVariants_table[1:3,1:9]

# Add an index column to the TEST_table and move it to the start 
TEST_table$ID <- seq.int(nrow(TEST_table))
TEST_table <- TEST_table %>%
  select(ID, everything())

j = 1 # column iterator
m = 1 # variant iterator
  
# Create a new empty dataframe
definitive_table <- data.frame()
# colnames(definitive_table) <- ("ID")

for (var_row in TEST_table[m,]) {
  j = 1
  # Create a new dataframe from the TEST_table with only ID column
  unnested_Table <- data.frame(TEST_table[m,1])
  colnames(unnested_Table) <- ("ID")
  current_row <- TEST_table[m,]
  print (paste("PROCESSING VARIANT", m))
  for (field in current_row) {
    print (paste("Analising column", j, "with class", class(current_row[,j])))
    one_class <- class(current_row[,j])
    if (j < 10) {
      print(j)
      if (one_class == "list" ) {
        print (paste("Unnesting List in column", j)) # testing line
        
        # Copy only first (key ID) column and last (nested) column
        print ("Extracting rows")
        processingTable <- current_row[,c(1,j)]
        processingTable1 <- subset(current_row, select=c(1, j))
        print ("unnesting rows)")
        # Unnest the table so the table will expand according to the last column
        unnesting <- unnest(processingTable)
        class(processingTable1[2])
        print ("Merge")
        # Merge the unnested table with the main table, according to the key ID
        unnested_Table <- merge(unnested_Table, unnesting, by="ID")
        
        print (paste("State of the table for column", j, "in variant", m))
        print (unnested_Table) # testing line
      } else {
        print (paste("Pasting column", j)) # testing line
        the_column <- current_row[,c(1,j)]
        unnested_Table <- merge(unnested_Table, the_column, by="ID")
        print (paste("State of the table for column", j, "in variant", m))
        print (unnested_Table)
      }
    }
    j = j + 1
    print ("+++++++++++ neXt colmun +++++++++++++++")
  }
  definitive_table <- bind_rows(definitive_table, unnested_Table)
  m = m + 1
  print ("###################### NEXT VARIANT ##################################")
  
}


flatTable <- flatten(TEST_table) # does not work so far

View(TEST_table[1:6])

annotVariant[, length(annotVariant)]

reducedTable <- TEST_table[1:6]

TEST_unn <- unnest(reducedTable)
  
  # Error: All nested columns must have the same number of elements.

  
# }
# # ,
# # data = list(REvaluate_local_var = global_var)
# # ,
# # REvaluator = Rversion,
# # verbose	= TRUE
# )

########### [RStudio] Print the table in a txt file ###########


########### [TEST block] ###############
