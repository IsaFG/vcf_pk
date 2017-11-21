# How to use REvaluate from RinR package

########### Load RinR library #############
library(RinR)

########### Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### REvaluate use #############
preVariable <- REvaluate({
  # Here comes libraries
  
  # code
  
  # result
  
}
# ,
# data = list(variable_to_use_in_REvaluate = variable_in_main_script)
# ,
# REvaluator = Rversion
# ,
# verbose	= TRUE
)

########### Convert the REvaluate Result (if needed) #############
processedVariable <- preVariable

########### Testing block #############
# This line is for testing is this table changes correctly
TestingTable02 <- variant_table[c(3,3,3,3)]