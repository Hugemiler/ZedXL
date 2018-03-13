restriction.distance.factor <- function(restrictionVector) {

  restrictionVector <- gsub("CB", "", restrictionVector)
  primaryDistance <- sapply(restrictionVector, function(x) { eval(parse(text=x))})
  primaryDistance < abs(primaryDistance)

  rdf <- (primaryDistance/(sum(primaryDistance) / length(primaryDistance)) - 1)

  return(rdf)

}
