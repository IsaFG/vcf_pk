install.packages("WhopGenome")

library(WhopGenome)

my_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\R_2015_01_27_14_48_49_user_XEN-66-Haloplex_316_Nefro_pool24_124_2305.vcf"
compressed_vcf <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\compressedVCF.vcf"
bgzf_compress(my_vcf, compressed_vcf)
