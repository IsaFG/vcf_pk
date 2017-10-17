# How to use REvaluate from RinR package

########### Load RinR library #############
library(RinR)

########### Testing block #############
# The following block is for testing
variant_table <- variants_table.txt
TestingTable01 <- variant_table[c(1,4)]

########### Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### REvaluate use #############
preVariable <- REvaluate({
  # Here comes libraries
  
  # code
  
  # result
  
}
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)

########### Convert the REvaluate Result #############
JSONtable <- preVariable

########### Testing block #############
# This line is for testing is this table changes correctly
TestingTable02 <- variant_table[c(3,3,3,3)]