library(cellbaseR)

# getCellBase 8/26
# Description
# The generic method for querying CellBase web services.
cb <- CellBaseR()
res <- getCellBase(object=cb, category="feature", subcategory="gene",
                   ids="TET1", resource="info")

cb <- CellBaseR()
genes <- c("TP73","TET1")
res <- getGene(object = cb, ids = genes, resource = "info")
str(res,2)


cb <- CellBaseR()
res2 <- getVariant(object=cb, ids="1:169549811:A:G", resource="annotation")
# to get the data 
res2 <- cbData(res2)
str(res2, 1)
