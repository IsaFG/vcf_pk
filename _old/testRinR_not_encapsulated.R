library(cellbaseR)
try(cb <- CellBaseR())
# variant <- testRinR # only for TIBCO
variant <- "1:169549811:A:G"
try(res2 <- getVariant(object=cb, ids=variant, resource="annotation"))
# to get the data 
res2 <- cbData(res2)
str(res2, 1)
testRinRresult <- res2
