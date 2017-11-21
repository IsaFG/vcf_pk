library("RinR")
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")
# variant <- testRinR # only for TIBCO

REvaluate({
  library(cellbaseR)
  try(cb <- CellBaseR())
  variant <- "1:169549811:A:G"
  try(res2 <- getVariant(object=cb, ids=variant, resource="annotation"))
  try(res2 <- cbData(res2))
},
REvaluator = Rversion,
verbose	= TRUE
)

# to get the data 
str(res2, 1)
testRinRresult <- res2