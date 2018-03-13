create.avgplot <- function(x) {
  tempAvg <- NULL
  for (i in (1:length(x))) {
    tempAvg[i] <- mean(x[1:i])
  }
  return(tempAvg)
}
