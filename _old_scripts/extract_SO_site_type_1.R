loc_file <- "C:\\Users\\FollonIn\\Documents\\GitHub\\vcf_pk\\test_files\\annotatedTables\\sequenceOntologyTerms_TIBCOmodified.txt"

my_table <- read.table(loc_file, sep="\t", header = T)

site_type_df <-  data.frame()
site_type_list <- list()

for(i in 1:nrow(my_table)) {
  row <- my_table[i,]
  # do stuff with row
  if (row$site_type != "") {
    print (i)
    site_type_list$variantID <- as.character(row$variantID)
    site_type_list$site_type <- as.character(row$site_type)
    site_type_df <- rbind(site_type_df, as.data.frame(site_type_list))
    print (site_type_list)
  }
}

# Keep only unique couples
uni_site_type_df <- unique(site_type_df)

for(i in 1:nrow(uni_site_type_df)) {
  row <- uni_site_type_df[i,]
  # do stuff with row
  if (uni_site_type_df$site_type != "splice") {
    print (i)
    site_type_list$variantID <- as.character(row$variantID)
    site_type_list$site_type <- as.character(row$site_type)
    site_type_df <- rbind(site_type_df, as.data.frame(site_type_list))
    print (site_type_list)
    
  } else if (uni_site_type_df$site_type != "missense") {
    
  } else if (uni_site_type_df$site_type != "nonsense") {
    
  } else if (uni_site_type_df$site_type != "missense") {
    
  } else if (uni_site_type_df$site_type != "frameshift") {

    
  }
}
