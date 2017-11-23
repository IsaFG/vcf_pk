############ Information #############
# How to read a vcf ?

############ References used #############
# "An Introduction to the GenomicRanges Package" by Carlson et al.
# VariantAnnotation.pdf : "1. Introduction to VariantAnnotation" by Valerie Obenchain
# https://www.biostars.org/p/137765/

############# [FACULTATIVE] Install libraries ###################
# source("https://bioconductor.org/biocLite.R") 
# biocLite("VariantAnnotation")
# biocLite("GenomicAlignments")
# biocLite("snpStats")
# install.packages("tibble")

############# Load libraries ###################
library(GenomicAlignments)
library(VariantAnnotation)

######### Import the VCF file into R #############
# EXAMPLE 1:
# La funcion system.file no pertenece a este paquete.
# tan solo se dedica a buscar la verdadera ubicacion del archivo que se usa como ejemplo
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
# fl # esto para ver cual es el link por pantalla

# EXAMPLE 2 (with the tabix file) 
file.gz     <- system.file("extdata", "chr7-sub.vcf.gz", 
                           package="VariantAnnotation")

file.gz.tbi <- system.file("extdata", "chr7-sub.vcf.gz.tbi", 
                           package="VariantAnnotation")

######### Extract the VCF as an object type "CollapsedVCF" (VariantAnnotation package) #######
# Data are parsed into a VCF object with readVcf.
# EXAMPLE 1:
vcf1 <- readVcf(fl, "hg19")
# EXAMPLE 2 :
vcf_chr7 <- readVcf(file.gz, "hg19")

######### Examinate the VCF objects ##################
samples(header(vcf_chr7))

# Have a look into the VCF object
vcf1

######### VCF exploration of header information ######
# Header information : gives samples names or number, info fields, geno fields...  
header(vcf1) # class VCFHeader

# Extract samples names or number
samples(header(vcf1)) # class character

# Extract "geno" information,
# es decir lo que significan las etiquetas en la columna FORMAT
geno(header(vcf1))

# numero y lista de etiquetas que salen en el campo FORMAT
geno(vcf1) # list

sapply(geno(vcf1), class)

########### Convert extracted vcf to a GRanges object ############
# rowRanges contains information from the CHROM, POS, and ID fields of the VCF file,
# represented as a GRanges.
rowRanges(vcf1) # class GenomicRanges

# Pass the vc object to a GRanges object
GRangesVcf1 <- rowRanges(vcf1)
GRangesVcf1

########### Extract only the first samples form the GRanges object ############
headRowRangeVCF1 <- head(rowRanges(vcf1), 3)
headRowRangeVCF1

# The paramRangeID column is meaningful when reading subsets of
# data and is discussed further below (see pdf)

########### Extracting random records from the vcf ############
my_sample_num <- runif(40, min=0, max=length(GRangesVcf1))
vcf_sample <- GRangesVcf1[c(my_sample_num)]

########### Extract Genomic ranges ############
# The genomic ranges can be extracted without corresponding metadata with granges
GenomRanges <- granges(GRangesVcf1)
GenomRanges

########### Extract metadata ############
# Annotations for these coordinates can be extracted
# as a DataFrame object using the mcols accessor.
metadata_df <- mcols(GRangesVcf1)
metadata_df

############# Extract the chromosomes and the ranges of the variants ###################
# Extract chromosomes, ranges from the Formal class GRanges
# and store them in vectors

# Extract a range vector
# ranges_col <- as.character(ranges(GRangesVcf1)) # This options does not work in some situation

# Extract ranges
ranges_start <- start(ranges(GRangesVcf1))
ranges_end <- end(ranges(GRangesVcf1))
ranges_df <- data.frame(ranges_start, ranges_end)

# Extract a chrom vector
# tip got from https://www.biostars.org/p/137765/
chroms_col <- as.character(seqnames(GRangesVcf1)) # works

########### Exploration of alleles information ############
# Extract the first reference alleles
ref(vcf1)[1:5]


