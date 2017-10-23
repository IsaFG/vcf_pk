############## [INFO] An example of CellBase with RinR for TIBCO ############
# This script is a draft for TERR TIBCO
# This script works
# IMPORTANT : this script will use the library "cellbaseR" 

############## [INFO] Input and Output ############################
# Input : a variant with which a CellBaseR method will query CellBase
# Output : the table retrieving the information from CellBase

############# [FACULTATIVE] Set Working directory ############
setwd("C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk")

############# [FACULTATIVE] Install libraries ###################
install.packages("RCurl")
install.packages("jsonlite")

############# [FACULTATIVE] Install library cellbaseR ##########
# From Bioconductor :
source("https://bioconductor.org/biocLite.R")
## try http:// if https:// URLs are not supported
biocLite("cellbaseR")

# From the develloper github (Mohammed Elsiddieg <melsiddieg@gmail.com>) :
install.packages("devtools") # you may want to install devtools to install libraries from Github
library(devtools)
devtools::install_github("melsiddieg/cellbaseR")

############# [HERE STARTS TIBCO SCRIPT] ###########
############# Load RinR library ###################
library(RinR)

########### Determinate R interpreter location ########
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

########### REvaluate use #############
preVariable <- REvaluate({
  library(RCurl)
  library(jsonlite)
  library(cellbaseR)
  
  cb <- CellBaseR()
  res2 <- getVariant(object=cb, ids="1:169549811:A:G", resource="annotation")
  
  res2
}
# ,
# REvaluator = Rversion,
# verbose	= TRUE
)

########### Convert the REvaluate Result #############
cbVariantsTable <- preVariable[c(1,2,3,4,6)]
