library("RinR")
Rversion <- makeREvaluator("R", RHome = "C:/Program Files/R/R-3.4.1")

getVariant <- function(sql) {
  REvaluate(
    
    library(cellbaseR)
    cb <- CellBaseR()
    res2 <- getVariant(object=cb, ids="1:169549811:A:G", resource="annotation")
    # to get the data 
    res2 <- cbData(res2)
    str(res2, 1)
    
    
  )
}