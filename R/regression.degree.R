regression.degree <- function(scoreTable, similarityMatrix) {

  rsquareMatrix <- matrix(nrow = nrow(scoreTable), ncol = nrow(scoreTable))

  progressBar <- txtProgressBar(min = 0, max = nrow(scoreTable), style = 3)
  for (i in 1:nrow(scoreTable)) {
    for (j in 1:nrow(scoreTable)) {
      rsquareMatrix[i,j] <- summary(lm(similarityMatrix[,rownames(scoreTable)[i]] ~ similarityMatrix[,rownames(scoreTable)[j]]))$r.squared
    }

    setTxtProgressBar(progressBar, i)
  }
  regressionScore <- apply(rsquareMatrix, 2, mean)
  return(regressionScore)
}
