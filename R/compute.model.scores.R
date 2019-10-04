compute.model.scores <- function(lovoalignLogPath,
                                 computegscore = F,
                                 gscoretype = "gscore",
                                 gscoreLogPath = "../gscore/gscore-TMscore-060.dat",
                                 computeproq3 = F,
                                 proq3listLocation = "../proq3list.txt"){

  ######
  #
  # This function parses output files from lovoalign, G-Score and Proq3D Software.
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
  # ProQ3D is open-source software. Whenever you use ProQ3, please cite:
  #
  # Karolis Uziela, David Menéndez Hurtado, Nanjiang Shu, Björn Wallner and Arne Elofsson,
  # ProQ3D: Improved model quality assessments using Deep Learning.
  # Bioinformatics. 2017 May 15;33(10):1578-1580
  #
  # Obtain ProQ3D software at http://proq3.bioinfo.se/pred/download/
  #
  ######

  modelScores <- data.frame(tmscore = read.table(lovoalignLogPath)[,3])

  ##
  # Computiong G-Score if

  if (computegscore == T) {

    gscoreTable <- read.table(gscoreLogPath)

    gscoreTable <- gscoreTable[order(gscoreTable$V4),]

    if (type == "gscore") {

      modelScores <- cbind(modelScores,
                           gscore = gscoreTable[,1])

      rownames(modelScores) <- gscoreTable$V4

    }

    else if (type == "degree") {

      modelScores <- cbind(modelScores,
                           gscore = gscoreTable[,2])

      rownames(modelScores) <- gscoreTable$V4

    }

    else if (type == "Wdegree") {

      modelScores <- cbind(modelScores,
                           gscore = gscoreTable[,3])

      rownames(modelScores) <- gscoreTable$V4

    }

  }

  ## Computing ProQ3 Scores

  if (computeproq3 == T) {

    proq3list <- readLines(proq3listLocation)
    modelNames <- do.call(rbind, strsplit(proq3list, split = "/"))
    modelNames <- as.character(modelNames[, ncol(modelNames)])
    proq3list <- gsub(".pdb", ".pdb.proq3d.tmscore.global", proq3list,fixed = T)
    proq3Table <- do.call(rbind, lapply(proq3list, function(x) {
      read.table(x, header = T)}))
    rownames(proq3Table) <- modelNames

    ## Compute tmscore

    proq3list <- gsub("sscore", "tmscore", proq3list,fixed = T)
    tmscoreTable <- do.call(rbind, lapply(proq3list, function(x) {
      read.table(x, header = T)}))
    rownames(tmscoreTable) <- modelNames

    ## Appending all modelScores computed

    modelScores <- cbind(modelScores, proq3Table)

  }

  return(modelScores)
}
