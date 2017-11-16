attribute.cur.and.rec <- function(){

  restriction.list <- rownames(restriction.parameters)

  # CURATED list
  curated.list <- readLines("~/Documents/datasets/SALB3.22CUR.001/utility/curlist.txt")
  curated.factor <- restriction.list %in% curated.list
  restriction.parameters$cur <- as.logical(curated.factor)

  # RECOVERED list
  recovered.list <- readLines("~/Documents/datasets/SALB3.22CUR.001/utility/reclist.txt")
  recovered.factor <- restriction.list %in% recovered.list
  restriction.parameters$rec <- as.logical(recovered.factor)

  #CRYSTALLOGRAPHIC list
  crystallographic.list <- readLines("~/Documents/datasets/SALB3.22CUR.001/utility/theorlist.txt")
  crystallographyc.factor <- restriction.list %in% crystallographic.list
  restriction.parameters$crys <- as.logical(crystallographyc.factor)

  #OOPTIMAL list
  optimum.list <- readLines("~/Documents/datasets/SALB3.22CUR.001/utility/optlist.txt")
  optimum.factor <- restriction.list %in% optimum.list
  restriction.parameters$opt <- as.logical(optimum.factor)

  return(restriction.parameters)
}
