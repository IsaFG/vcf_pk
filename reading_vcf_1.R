# Create an R script capable of reading in a vcf formatted file and extracting a list of variants.
# The list has to be extracted to be used in this format :
# http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/{category}/{subcategory}/{id}/{resource}?{filters}
# For example, we have this link :
# "bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/BRCA2,BRCA1,KRAS/info")
# this REST call will get all the INFO of the gene BRCA2,BRCA1,KRAS of human in the latest version.

# Vamos a probar un par de paquetes para leer vcf: vcfR y VariantAnnotation de bioconductor
# en este script sera VariantAnnotation

## OPCIONAL, PARA INSTALAR LA LIBRERIA
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
biocLite("GenomicAlignments")

# cargar la libreria
library(GenomicAlignments)
library(VariantAnnotation)

# EJEMPLO SACADO DEL MANUAL
# la funcion system.file no pertenece a este paquete.
# tan solo se dedica a buscar la verdadera ubicacion del archivo que se usa como ejemplo
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
# Data are parsed into a VCF object with readVcf.
vcf1 <- readVcf(fl, "hg19")
vcf1
fl
header(vcf1)
head(rowRanges(vcf1), 3)

# intentando mi vcf
my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
vcf2 <- readVcf(my_vcf, "hg19")
vcf2
header(vcf2)
head(rowRanges(vcf2), 3)
