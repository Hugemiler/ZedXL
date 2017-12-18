create.dissimilarity.matrix <- function(mode = "dissimilarity",
                                        diagonal = 0,
                                        alignlistPath,
                                        compactlogPath,
                                        nmodels){

  #######################################
  # This function reads a compactLog and computes its contents as a
  # similarity matrix or a dissimilarity matrix.
  #
  # Dissimilarity is calculated as 1 - TMScore
  #
  # Similarity Matrix is used as an input to BetaK histograms
  #
  # Dissimilarity matrix is used as an input to ForceScheme projections.
  #
  #######################################

  logList <- readLines(alignlistPath)

  for (i in 1:length(logList)){
    logList[i] <- strsplit(logList[i], "/")[[1]][length(strsplit(logList[i], "/")[[1]])]
  }

  # Read COMPACTLOG

  similarityMatrix <- read.table(compactlogPath,
                                 skip = nmodels + 5,
                                 fill = TRUE)

  # Adjust COMPACTLOG matrix (dimension nmodels-1 to dimensionnmodels and upper triangular matrix)

  similarityMatrix[length(logList),length(logList)] <- NA

  similarityMatrix[]<-t(apply(similarityMatrix,1,function(x){
    c(x[is.na(x)], x[!is.na(x)])}))

  similarityMatrix[is.na(similarityMatrix)] <- 0

  # Naming by model list

  dimnames(similarityMatrix)[[1]] <- logList
  dimnames(similarityMatrix)[[2]] <- logList

  # Hard symmetrizing matrix

  similarityMatrix <- t(similarityMatrix) + similarityMatrix

  # Setting diagonal value

  diag(similarityMatrix) <- diagonal

  ## Compute dissimilarity from similarity

  dissimilarityMatrix <- 1 - similarityMatrix
  diag(dissimilarityMatrix) <- diagonal

  # Return desired matrix mased on "MODE"

  if (mode == "similarity") {
    return(as.matrix(similarityMatrix))
  } else if (mode == "dissimilarity") {
    return(as.matrix(dissimilarityMatrix))
  }

}
