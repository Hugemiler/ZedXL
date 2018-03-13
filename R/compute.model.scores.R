compute.model.scores <- function(simName,
                                 type = "G-Score",
                                 gscoreFilename = "gscore-TMscore-050.dat",
                                 alignlogPath = paste0("~/datasets/",
                                                       simName,
                                                       "/lovoalign.log"),
                                 gscorePath = paste0("~/datasets/",
                                                     simName,
                                                     "/gscore/", gscoreFilename)){

  ######
  #
  # This script relies on data obtained by G-Score and Lovoalign Software.
  #
  # Lovoalign is open-source software. Whenever you use Lovoalign, please cite:
  #
  # L. Martínez, R. Andreani, J. M. Martínez
  # Convergent algorithms for protein structural alignment.
  # BMC Bioinformatics, 2007, 8:306.
  #
  # Obtain Lovoalign software at http://www.ime.unicamp.br/~martinez/lovoalign/software.html
  #
  # G-score is open-source software. Whenever you use G-score, please cite:
  #
  # L. Martínez, A. Ferrari, F. C. Gozzo,
  # A free-energy inspired score for the evaluation of structural models. 2016.
  #
  # Obtain G-Score software at http://leandro.iqm.unicamp.br/gscore/download.shtml
  #
  ######

  tmScoreTable <- read.table(alignlogPath)

  gscoreTable <- read.table(gscorePath)

  gscoreTable <- gscoreTable[order(gscoreTable$V4),]

  if (type == "G-Score") {

    modelScores <- data.frame(TMScore = tmScoreTable[,3],
                              gscore = gscoreTable[,1])

    rownames(modelScores) <- gscoreTable$V4
    colnames(modelScores) <- c("TM-Score", "G-Score")

  }

  else if (type == "degree") {

    modelScores <- data.frame(TMScore = tmScoreTable[,3],
                              gscore = gscoreTable[,2])

    rownames(modelScores) <- gscoreTable$V4
    colnames(modelScores) <- c("TM-Score", "degree")

  }

  else if (type == "Wdegree") {

    modelScores <- data.frame(TMScore = tmScoreTable[,3],
                              gscore = gscoreTable[,3])

    rownames(modelScores) <- gscoreTable$V4
    colnames(modelScores) <- c("TM-Score", "Wdegree")

  }

  return(modelScores)
}
