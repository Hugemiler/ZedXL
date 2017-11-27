attribute.cur.and.rec <- function(restrictionScores){

  restriction.list <- rownames(restrictionScores)

  # CURATED list
  curated.list <- readLines("~/datasets/SALB3.22CUR.001/utility/curlist.txt")
  curated.factor <- restriction.list %in% curated.list
  restrictionScores$cur <- as.logical(curated.factor)

  # RECOVERED list
  recovered.list <- readLines("~/datasets/SALB3.22CUR.001/utility/reclist.txt")
  recovered.factor <- restriction.list %in% recovered.list
  restrictionScores$rec <- as.logical(recovered.factor)

  #CRYSTALLOGRAPHIC list
  crystallographic.list <- readLines("~/datasets/SALB3.22CUR.001/utility/theorlist.txt")
  crystallographyc.factor <- restriction.list %in% crystallographic.list
  restrictionScores$crys <- as.logical(crystallographyc.factor)

  #OOPTIMAL list
  optimum.list <- readLines("~/datasets/SALB3.22CUR.001/utility/optlist.txt")
  optimum.factor <- restriction.list %in% optimum.list
  restrictionScores$opt <- as.logical(optimum.factor)

  return(restrictionScores)
}
