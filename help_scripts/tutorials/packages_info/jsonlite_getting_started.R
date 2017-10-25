library(jsonlite)
all.equal(mtcars, fromJSON(toJSON(mtcars)))

# A JSON array of primitives
json <- '["Mario", "Peach", null, "Bowser"]'

json

# Simplifies into an atomic vector
fromJSON(json)
# No simplification:
fromJSON(json, simplifyVector = FALSE)

# create a json with a dataframe
json <-
  '[
{"Name" : "Mario", "Age" : 32, "Occupation" : "Plumber"}, 
{"Name" : "Peach", "Age" : 21, "Occupation" : "Princess"},
{},
{"Name" : "Bowser", "Occupation" : "Koopa"}
]'
# convert to a dataframe
mydf <- fromJSON(json)
mydf

# adding a new data at each record in the dataframe: Ranking
mydf$Ranking <- c(3, 1, 2, 4)
# convert it to a json again
toJSON(mydf, pretty=TRUE)
