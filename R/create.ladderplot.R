create.ladderplot <- function(x) {
  for (i in (2:length(x))) {
    if (x[i] >= x[i-1]) { x[i] <- x[i-1] }
  }
  return(x)
}
