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

##################################
# EJEMPLO SACADO DEL MANUAL DE VARIANT ANNOTATION
# Ver en el manual y en el archivo VariantAnnotation_1.R
##################################

# intentando mi vcf
my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
vcf2 <- readVcf(my_vcf, "hg19")
vcf2
# info del cabezal del VCF
header(vcf2)

# objecto GRange con las tres primeras variantes
headRowRangeVCF2 <- head(rowRanges(vcf2), 3)
headRowRangeVCF2

# objecto GRange con la primera variante
head(rowRanges(vcf2), 1)

#intento de extraccion del cromosoma + la posicion:
# cromosoma :
seqnames_ex1 <- seqnames(head(rowRanges(vcf2), 1))
seqnames_ex1
# determinamos que clase es
sapply(seqnames_ex1, class)
# es un rle object

# posicion:
ranges_ex1 <- ranges(head(rowRanges(vcf2), 1))
ranges_ex1
start(ranges_ex1)

############# ESTO PUEDE SER UTIL #############
# extraer el rango genomico (cromosoma + posicion)
Grangesvcf2 <- granges(headRowRangeVCF2)
Grangesvcf2
as.character(seqnames(Grangesvcf2))
as.character(ranges(Grangesvcf2))
###############################################


# lista de todos los alelos de referencia
ref(vcf2)[1:5]

# lista de todos los alelos alterados
alt(vcf2)[1:5]

# vector numerico con las cualidades QUAL
qual(vcf2)[1:5]

# numero y lista de etiquetas que salen en el campo FORMAT
geno(vcf2)
# podemos averiguar la clase de esas etiquetas
sapply(geno(vcf2), class)

# crear una matriz de genotipos
vcfmat2 <- genotypeToSnpMatrix(vcf2)
vcfmat2

sapply(granges(headRowRangeVCF2), class)
seqnames(granges(headRowRangeVCF2))
