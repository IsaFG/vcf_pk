########### Annotate VCF ############
library(VariantAnnotation)

fl <- system.file("extdata", "hapmap_exome_chr22_500.vcf.gz",package = "cellbaseR" )
res <- AnnotateVcf(object=cb, file=fl, BPPARAM = bpparam(workers=2))
vcf <- readVcf(fl, "hg19")
samples(header(vcf))

########## Test Annotate VCF #############
my_vcf <- "C:/Users/FollonIn/Documents/GitHub/vcf_pk/R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
res3 <- AnnotateVcf(object=cb, file=my_vcf, BPPARAM = bpparam(workers=2))
vcf3 <- readVcf(fl3, "hg19")
samples(header(vcf3))

bgzip(from, to)

fl