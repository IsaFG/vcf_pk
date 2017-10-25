# como construir un objeto Grange
# linea 1 : seqnames : definir los nombres (Rle) (en este caso "chrx")
## y determinar donde aparecen y cuantas veces (segundo parentesis)
# linea 2 : ranges : define los rangos (Iranges) y la etiqueta que reciben esos rangos (letra)
# linea 3 : strand : el sentido de la cadena, tambien se determina por donde y cuantas veces aparece
# linea 4 : score : es obvio y GC viene determinado por una distribucion

gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
score = 1:10, GC = seq(1, 0, length=10))

gr

genomicRange <- granges(gr)
genomicRange


