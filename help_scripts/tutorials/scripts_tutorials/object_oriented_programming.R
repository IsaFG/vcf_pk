# LEARNING CLASSES
# example from the book R Object-Oriented Programming.pdf p99
oneDie <- list(trials=character(0))
class(oneDie) <- "Die"
oneCoin <- list(trials=character(0))
class(oneCoin) <- "Coin"

reset <- function(theObject)
{
  UseMethod("reset",theObject)
  print("Reset the Trials")
}
reset.default <- function(theObject)
{
  print("Uh oh, not sure what to do here!\n")
  return(theObject)
}
reset.Die <- function(theObject)
{
  theObject$trials <- character(0)
  print("Reset the die\n")
  return(theObject)
}
reset.Coin <- function(theObject)
{
  theObject$trials <- character(0)
  print("Reset the coin\n")
  return(theObject)
}