imperialist.score <- function(models = NULL,
                             similarityTable = NULL) {

  imperialistMatrix <- matrix(nrow = length(models), ncol = length(models))

  for (i in 1:length(models)) {
    for (j in 1:length(models)) {
      imperialistMatrix[i, j] <- similarityTable[which(colnames(similarityTable) == models[i]), which(rownames(similarityTable) == models[j])]
      print(i, j, imperialistMatrix[i, j])
    }
  }

  imperialistScore <- vector()

  for (l in 1:length(models)) {
    imperialistScore[l] <- sum(imperialistMatrix[l,]) / (length(models) - 1)
  }

  return(imperialistScore)

}
