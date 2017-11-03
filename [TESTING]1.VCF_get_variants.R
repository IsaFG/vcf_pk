scanVcf(vcf_extracted)


View(data.frame(REFal))

ALTal_df_TEST <- ALTal_df

i <- 1
for (boolresult in duplicated(ALTal_df[,1])) {
  if (boolresult == TRUE) {
    print (i)
    print (ALTal_df[,1][[i]])
    ALTal_df_TEST[,3][[i-1]] <- paste(ALTal_df_TEST[,3][[i-1]], ALTal_df_TEST[,3][[i]], sep=",")
    ALTal_df_TEST <- ALTal_df_TEST[-c(i),]
  }
  i <- i + 1
}
  

ALTal_df_TEST[-c(4,)]

