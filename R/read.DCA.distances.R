read.DCA.distances <- function(pdblist.location = "~/Documents/datasets/ABCt.DCA1.001/pdblist.txt",
                               distance.folder = "~/Documents/datasets/ABCt.DCA1.001/distances/"){

#  pdblist.location = "~/Documents/datasets/ABCt.DCA1.001/pdblist.txt"
#  distance.folder = "~/Documents/datasets/ABCt.DCA1.001/distances/"

  pdb.list <- readLines(pdblist.location)
  pdb.list <- gsub(".pdb", "", pdb.list)
  DCA.distances.list <- list()

  for (i in 1:length(pdb.list)){

    temporary.DCA.table <- read.table(paste0(distance.folder,
                                             pdb.list[i],
                                             "_truedistances"))
    DCA.distances.list[[pdb.list[i]]] <- temporary.DCA.table$V3
    names(DCA.distances.list[[pdb.list[i]]]) <- paste0(temporary.DCA.table$V1,
                                                   "-",
                                                   temporary.DCA.table$V2)
  }

  names(DCA.distances.list) <- pdb.list
  return(DCA.distances.list)
}
