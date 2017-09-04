# Create an R script capable of reading in a vcf formatted file and extracting a list of variants.
# The list has to be extracted to be used in this format :
# http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/{category}/{subcategory}/{id}/{resource}?{filters}
# For example, we have this link :
# "bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
# "bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/genomic/region/1:18149476-18149476/gene"
# this REST call will get all the INFO of the gene BRCA2,BRCA1,KRAS of human in the latest version.

# Vamos a probar un par de paquetes para leer vcf: vcfR y VariantAnnotation de bioconductor
# en este script sera VariantAnnotation

##########################################
## OPCIONAL, PARA INSTALAR LA LIBRERIA
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
biocLite("GenomicAlignments")
biocLite("snpStats")
install.packages("tibble")
##########################################
# cargar las librerias
library(GenomicAlignments)
library(VariantAnnotation)

##########################################
# Extraemos mi vcf y lo pasamos a un objeto "CollapsedVCF" mediante VariantAnnotation
my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
vcf2 <- readVcf(my_vcf, "hg19")
vcf2

# Para obtener info del cabezal del VCF
header(vcf2)

# Podemos extraer un objecto GRange con la primera variante del fichero VCF
RowRangeVCF2_v1 <- head(rowRanges(vcf2), 1)
RowRangeVCF2_v1

# De esa primera variante podemos extraer el cromosoma y el rango.
# truco sacado de https://www.biostars.org/p/137765/
# esta información nos podria servir como base para hacer una lista de variantes
as.character(seqnames(RowRangeVCF2_v1))
as.character(ranges(RowRangeVCF2_v1))
as.numeric(ranges(RowRangeVCF2_v1))
###########################################################################

##### TEST ################
RowRangeVCF2_v1a5 <- head(rowRanges(vcf2), 5)
RowRangeVCF2_v1a5

# extraer los cromosomas y rangos en listas de caracteres
# creamos pues un vector con los cromosomas y otro con los rangos
cromosomas_v1a5_chr <- as.character(seqnames(RowRangeVCF2_v1a5))
rangos_v1a5 <- as.character(ranges(RowRangeVCF2_v1a5))

# Para hacerlo tratable por el recurso RESTFul, habrá que retirarle el "chr"
string_to_remove <- "chr"
cromosomas_v1a5 <- gsub(x = cromosomas_v1a5_chr, pattern = paste(string_to_remove, collapse = "|"), replacement = "")


# Ahora construimos un dataframe con cromosomas y rangos de las variantes
variant.data <- data.frame(cromosomas_v1a5, rangos_v1a5)
variant.data

