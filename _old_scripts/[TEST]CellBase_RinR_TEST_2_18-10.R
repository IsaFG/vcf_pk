############# Load RinR library ###################
library(RinR)

########### Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### REvaluate use #############
variant <- "1:169549811:A:G"

preVariable <- REvaluate({
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  
  cb <- CellBaseR()
  res2 <- getVariant(object=cb, ids=REvariant, resource="annotation")
  
  res2
},
data = list(REvariant = variant)
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)

########### Convert the REvaluate Result #############
cbVariantsTable <- preVariable[c(1,2,3,4,6)]

########### Testing block #############
# The following block is for testing
variant_table <- variants_table.txt
TestingTableCellRinR <- variant_table[c(4,4,1,1,1,1,1,1)]