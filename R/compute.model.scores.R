compute.model.scores <- function(alignlogPath = "~/datasets/SALB3.OPTIMUM.001/gscore/gscore-TMscore-050.dat",
                                 gscorePath = "~/datasets/SALB3.OPTIMUM.001/lovoalign.log"){

  ######
  #
  # This script relies on logs obtained by TopoLink Software. TopoLink is open-source software.
  # Whenever you use TopoLink software, please cite:
  #
  # Ferrari, Allan & Martínez, Leandro & Cesar Gozzo, Fabio. (2017).
  # TopoLink: A software to validate structural models using chemical crosslinking constraints.
  # Nature Protocol Exchange. 10.1038/protex.2017.035.
  #
  # Obtain TopoLink software at http://leandro.iqm.unicamp.br/topolink/home.shtml
  #
  ######

  tmScoreTable <- read.table(grep(pattern = "#",
                             x = readLines(alignlogPath),
                             invert = TRUE))

  gscoreTable <- read.table(grep(pattern = "#",
                             x = readLines(gscorePath),
                             invert = TRUE))

  modelScores <- data.frame(TMScore = tmScoreTable[,"TM-Score"],
                            gscore = gscoreTable[,"gscore"])

  return(logTableList)
}
