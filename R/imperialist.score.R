imperialist.score <- function(models = NULL,
                              similarityTable = NULL,
                              mode = "cutoff") {

  if (mode == "consensus") {

    tempMatrix <- similarityTable + t(similarityTable)
    imperialistScore <- apply(tempMatrix, 1, function(x) {sum(x)/(length(models) - 1)})
    return(imperialistScore)

  }

  else if (mode == "cutoff") {

    imperialistMatrix <- matrix(nrow = length(models), ncol = length(models))

    for (i in 1:length(models)) {
      for (j in 1:length(models)) {
        imperialistMatrix[i, j] <- similarityTable[which(colnames(similarityTable) == models[i]), which(rownames(similarityTable) == models[j])]
        #print(i, j, imperialistMatrix[i, j])
      }
    }

    imperialistScore <- vector()

    for (l in 1:length(models)) {
      imperialistScore[l] <- sum(imperialistMatrix[l,]) / (length(models) - 1)
    }

    return(imperialistScore)

  }

}
