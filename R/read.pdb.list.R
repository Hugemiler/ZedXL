read.pdb.list <- function(pdblist.location = "~/Documents/datasets/ABCt.DCA1.001/pdblist.txt",
                          pdb.folder = "~/Documents/datasets/ABCt.DCA1.001/pdb/",
                          first.line = 785,
                          last.line = 940){

  model.list <- readLines(pdblist.location)
  log.table.list <- list()
  for (i in 1:length(model.list)){
    log.table.list[[model.list[i]]] <- read.table(paste0(output.location, model.list[i]),
                                                  skip = (first.line - 1),
                                                  nrows = (last.line - first.line + 1))
  }
  names(log.table.list) = model.list
  return(log.table.list)
}
