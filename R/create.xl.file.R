create.xl.file <- function(restriction.parameters = "restriction.parameters",
                           type,
                           FASTA = "~/Documents/datasets/SALB3.22CUR.001/loglist.txt") {
  #fasta <- readLines(FASTA)
  fasta <- as.vector(strsplit("MQDEQKRKEIVAEYFRKVNEGDVDAIVEMFTENATIEDPVGKDVREGRAAQREYFNSNVTAEVTIEPGHLSAGQDGKSVAVALAAEMTNILDPNRTRVKINAVDVFTLTPEGKIDSMRVFWGMTDIGVWNSSSV", "", fixed = TRUE))

  restriction.column <- subset(restriction.parameters, select = type)
  selected.restriction.names <- rownames(subset(restriction.column, restriction.column == TRUE))

  selected.restriction.names <- gsub("CB", "", selected.restriction.names)
  strsplit(selected.restriction.names, "-")
}
