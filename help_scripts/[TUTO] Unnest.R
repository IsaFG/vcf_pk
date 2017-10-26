####### [INFO] How to unnest a data.frame ##########

# Build a dataframe 1 
resp1 <- c("A", "B; A", "B", NA, "B")
resp2 <- c("C; D; F", NA, "C; F", "D", "E")
resp3 <- c(NA, NA, "G; H; I", "H; I", "I")
data <- data.frame(resp1, resp2, resp3, stringsAsFactors = F)
data

# Build a dataframe 2
df1 <- data_frame(
  x = 1,
  y = list(a = 1, b = 2)
)
df1
unnest(df1)

# Build a dataframe 3


# Transform to a dataframe
data1 <- data %>% transform(resp1 = strsplit(resp1, "; "),
                            resp2 = strsplit(resp2, "; "),
                            resp3 = strsplit(resp3, "; "))

unnest(data1)


flatten(data1)

TEST_unn <- unnest(data)

TEST_unn <- flatten(data)

TEST_unn


try1 <- strsplit(resp1, "; ")
try1

unnest(try1)

# This trown an error :
data %>%
  transform(resp1 = strsplit(resp1, "; "),
            resp2 = strsplit(resp2, "; "),
            resp3 = strsplit(resp3, "; ")) %>%
  unnest(resp1) %>% unnest(resp2) %>% unnest(resp3)


