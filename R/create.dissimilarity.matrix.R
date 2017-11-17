create.dissimilarity.matrix <- function(alignlistPath,
                                        compactlogPath,
                                        nmodels){

  logList <- readLines(alignlistPath)

  for (i in 1:length(logList)){
    logList[i] <- strsplit(logList[i], "/")[[1]][length(strsplit(logList[i], "/")[[1]])]
  }

  similarityMatrix <- read.table(compactlogPath,
                                 skip = nmodels + 5,
                                 fill = TRUE)

  dissimilarityMatrix <- similarityMatrix
  dissimilarityMatrix[length(logList),length(logList)] <- NA

  dissimilarity.matrix[]<-t(apply(dissimilarityMatrix,1,function(x){
    c(x[is.na(x)], x[!is.na(x)])}))

  diag(dissimilarityMatrix) <- 1
  dissimilarityMatrix <- 1 - dissimilarityMatrix
  dissimilarityMatrix[is.na(dissimilarity.matrix)] <- 0
  dissimilarityMatrix <- t(dissimilarityMatrix) + dissimilarityMatrix

  dimnames(dissimilarityMatrix)[[1]] <- logList
  dimnames(dissimilarityMatrix)[[2]] <- logList

  return(dissimilarityMatrix)

}
