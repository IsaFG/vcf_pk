######################################################
###### HOW TO USE THE DXP FILE #######################
######################################################
The dxp file is the aplication, and can only be loaded in TIBCO Spotfire.
A free trial of the plaftorm can be downloaded from https://spotfire.tibco.com/es/trial
Once installed in the yout machine, you can execute it and open the dxp.

######################################################
###### HOW TO USE THE R SCRIPTS ######################
######################################################
The zip file must be uncompressed in a folder which should be set at working directory before running the scripts.

It is important to note that almost all these scripts haven been develloped to integrate them in the TIBCO Spotfire platform. Some of them can be run in standard R version, but some others do not make sense outside Spotfire.
Rversion used in TIBCO Spotfire is named "TERR".

######################################################
###### Important note about the scripts structure ####
######################################################
TIBCO Spotfire integrates a specific version of R, called TERR and depending on the R version, the script execution can be slightly different.
In order to use several scripts either in TERR and R standard version, the structure of several scripts includes some code lines that determinate the R version being runned before running any other code line.

################################################
###### Input file ##############################
################################################
We suggest to use the input file named "NA12878.chr.22.vcf"

################################################
###### Scripts description and workflow ########
################################################
0.index_VCF_file.R : generates an indexed VCF file from a non indexed VCF file
# Inputs : a non indexed VCF file
# Output : Instance file for a VCF file and his tabix file 

1.fromVCFgetTable.R : Read in a vcf formatted file and extract a list of variants
# INPUT : a VCF file (not indexed)
# OUTPUT : a list of variants. The list will be stored as a dataframe. 

2.annotateIndexedVCF : Provide an indexed VCF (file with his tabix file), in order to query the RESTfull DB CellBase to retreive information about the variants
# INPUT : an instance file for a VCF file and his tabix file (from the script 0.index_VCF_file.R)
# Output : a table "annotVCFTable" with all the variants with their annotations + List of available annotation : "availableAnnots"

3.extract_all_annot.R : extract all annotations stored in the table "annotVCFTable" to build one new table per annotation
# Input 1: an annotated variant table, variable "annotVCFTable". In R standard version, this table should be previously loaded in global environment as a data.frame
# Input 2: List of available annotation : "availableAnnots"
# Output : a dataframe per annotation, named after the annotation + a list storing all these dataframes "tableList"

#### IMPORTANT NOTE ####
In order to speed up the process, you can load all these tables from the RData file provided in the folder : "NA12878.chr.22_allTables.RData"

4.post_processing_specific_tables.R : generates additional sub-tables for some interesting information. So far, this method has only been used for "consequenceTypes" and "traitAssociation"
# Input 1: a variant annotated table that needs post-processing. In R standard version, this table should be present in global environment as a data.frame. 
# Input 2: Annotation and his sub-annotation of interest
# Output : a dataframe with the sub-annotation of interest

5.adjust_cadd_table.R : adjust "functionalScoreTable" to unify data for one variant in only one row
# Inputs : Variant table with CADD annotacion scores : functionalScoreTable
# Output : caddScoresTable

7.clinvarTable_clean.R : unifies variable formats inside "clinvarTable"
# Inputs : clinvarTable
# Output : clinvarSignificance  + clinvarDisease

7.subsetPop.R : Subset the "populationFrecuenciesTable" according to study and population chosen by user
# Inputs : populationFrecuenciesTable + chosen pop + chosen study 
# Output : a subset of the populationFrecuenciesTable

8.adjust_conservation_table.R : adjust the conservation scores table to unify data in order to have one row per variant
# Inputs : conservationTable, with one row per conservation score type
# Output : adjusted conservationTable, with only one row per variant
