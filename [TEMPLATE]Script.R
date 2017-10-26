############## [INFO] Script [num] General information ############


############## [INFO] Input and Output ############################
# Inputs : 
# Output : 

############# [INFO] Previous requirement (libraries) ###############


############## [INFO] Problems ############################


# ############# [RStudio] Set Working directory ############
# Please uncoment the next line changing the working directory by the correct one:
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

########### [RStudio] Get the input file ############


############# [TIBCO] Load RinR library ###################
library(RinR)

########### [TIBCO] Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### [TIBCO] Create the REvaluate object ########
annotatedTable <- REvaluate({
  ############# [Code] Load libraries ############# 

  
  ############# [Code] Main method ############

  
},
# data = list(REvaluate_local_var = global_var)
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)

########### [RStudio] Print the table in a txt file ###########


########### [TEST block] ###############
