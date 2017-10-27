############### Iterate a dataframe ####################

TEST_table <- annotVariants_table[1:3,]

# Iterate by column
j = 1
for (any_column in TEST_table) {
  print (paste("Column", j))
  print (any_column)
  j = j + 1
  print ("+++++++++++++++++++++++++")
}


df <- data.frame(a=1:2, b=letters[1:2]) 

df

df[rep(seq_len(nrow(df)), each=2),]


n.times <- c(2,4)
df[rep(seq_len(nrow(df)), n.times),]
