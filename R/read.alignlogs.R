read.alignlogs <- function(alignlist.location = "/home/gbottino/Documents/datasets/SALB3.BISTRIVIAL.002/alignlist.txt",
                           compactlog.location = "/home/gbottino/Documents/datasets/SALB3.BISTRIVIAL.002/gscore/compactlog-TMscore.dat",
                           nmodels = 5000){

  model.list <- readLines(alignlist.location)

  for (i in 1:length(model.list)){
    model.list[i] <- strsplit(model.list[i], "/")[[1]][length(strsplit(model.list[i], "/")[[1]])]
  }

  similarity.matrix <- read.table(compactlog.location,
                                  skip = nmodels + 5,
                                  fill = TRUE)

  dissimilarity.matrix <- similarity.matrix
  dissimilarity.matrix[length(model.list),length(model.list)] <- NA

  dissimilarity.matrix[]<-t(apply(dissimilarity.matrix,1,function(x){
    c(x[is.na(x)], x[!is.na(x)])}))

  diag(dissimilarity.matrix) <- 1
  dissimilarity.matrix <- 1 - dissimilarity.matrix
  dissimilarity.matrix[is.na(dissimilarity.matrix)] <- 0
  dissimilarity.matrix <- t(dissimilarity.matrix) + dissimilarity.matrix

  dimnames(dissimilarity.matrix)[[1]] <- model.list
  dimnames(dissimilarity.matrix)[[2]] <- model.list

  return(dissimilarity.matrix)

}
