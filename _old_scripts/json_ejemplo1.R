library(RCurl)
library(jsonlite)

# asignar el archivo json a un objeto
# NOTA IMPORTANTE : si nos da Error: '\U' used without hex digits in character string starting ""C:\U",
# tenemos que aniadir dos baras en vez de una
exampleJson1 <- fromJSON("C:\\Users\\FollonIn\\Documents\\R\\json_array1.json")

exampleJson2 <- fromJSON("C:\\Users\\FollonIn\\Documents\\R\\json_object1.json")

exampleJson3 <- fromJSON("C:\\Users\\FollonIn\\Documents\\R\\json_array2.json")

exampleJson3
