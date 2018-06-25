select.constraints.CLI <- function(inputfile, outputfile) {

  source(inputfile)

  ###############
  # Script for the command-line selection methodology
  # May 2018
  #
  # Guilherme Fahur Bottino
  ###############

  # Initial Package Load and Environment Set

  ## descriptive analysis
  require(ZedXL)

  ## Set random seed
  set.seed(102528)

  # Parsing Data

  ## Parsing ModelScores

  modelScores <- compute.model.scores(type = "gscore",
                                      gscoreLogPath,
                                      lovoalignLogPath,
                                      proq3listLocation,
                                      computeproq3)

  ## Computing XlinkMirtTable from topolink logs

  optimumXlinkMirttable <- create.xlink.mirttable(
    prepare.topolink.logs(
      read.topolink.output(mode,
                           loglistLocation,
                           topolinkLogsDirectory)
    )
  )

  freqlist <- colnames(optimumXlinkMirttable)[order(-colSums(optimumXlinkMirttable))][1:nconst]

  ## Computing optimumSimilarityTable from nxn alignments

  optimumSimilarityTable <- create.dissimilarity.matrix(mode = "similarity",
                                                        diagonal = 1,
                                                        alignlistPath,
                                                        compactlogPath,
                                                        nmodels = nrow(optimumXlinkMirttable))

  ## Appendind DavisConsensus to modelScores table

  modelScores$davisconsensus <- imperialist.score(mode = 'consensus',
                                                  models = rownames(optimumSimilarityTable),
                                                  similarityTable = optimumSimilarityTable)

  ## Calculating RegressionScore

  regressionTable <- modelScores[order(-modelScores[, 'davisconsensus']), ][1:100, ]

  regressionTable$regressionScore <- regression.degree(regressionTable,
                                                       optimumSimilarityTable)

  ## Digesting the optimumXlinkMirttable

  optimumXlinkMirttable <- digest.xlink.mirttable(optimumXlinkMirttable)

  ### Classical Analysis

  optimumDescript <- ltm::descript(optimumXlinkMirttable)

  ## Computation of restrictionScores

  globalCronbachAlpha <- optimumDescript$alpha[1]

  ### Creating the RestrictionScores table

  restrictionScores <- data.frame(bis = optimumDescript$bisCorr,
                                  right = optimumDescript$perc[,2],
                                  logit = optimumDescript$perc[,3],
                                  alpha = optimumDescript$alpha[-1],
                                  deltaAlpha = optimumDescript$alpha[2:length(optimumDescript$alpha)] - globalCronbachAlpha)

  restrictionScores$rscore <- restriction.differential.scores(type = "G-Score",
                                                              optimumXlinkMirttable,
                                                              modelScores)$rscore

  restrictionScores$biscore <- -apply(optimumXlinkMirttable, 2, function(x) {
    ltm::biserial.cor(optimumSimilarityTable[which(modelScores$davisconsensus == max(modelScores$davisconsensus)), ], x)})

  restrictionScores$biscore_best <- -apply(optimumXlinkMirttable, 2, function(x) {
    ltm::biserial.cor(optimumSimilarityTable[which(modelScores$`TM-Score` == max(modelScores$`TM-Score`)), ], x)})

  restrictionScores$biscore_native <- -apply(optimumXlinkMirttable, 2, function(x) {
    ltm::biserial.cor(modelScores$`TM-Score`, x)})

  restrictionScores$biscore_regression <- -apply(optimumXlinkMirttable, 2, function(x) {
    ltm::biserial.cor(optimumSimilarityTable[
      which(colnames(optimumSimilarityTable) == rownames(regressionTable)[
        which(regressionTable$regressionScore == max(regressionTable$regressionScore))]), ], x)})

  restrictionScores <- attribute.cur.and.rec(restrictionScores)

  # Correlations and Charts

  ## Index Assignment

  maxIndex <- which(modelScores$`TM-Score` == max(modelScores$`TM-Score`))
  daviesIndex <- which(modelScores$`TM-Score` ==
                         modelScores$`TM-Score`[which(modelScores$davisconsensus == max(modelScores$davisconsensus))])
  regressionIndex <- which(rownames(modelScores) == rownames(regressionTable)[which(regressionTable$regressionScore == max(regressionTable$regressionScore))])
  proq3Index <- which(modelScores$ProQ3D == max(modelScores$ProQ3D))

  # ## Correlation plots and linear regressions
  #
  # require(ggplot2)
  #
  # ### Max TM-Score
  # summary(lm(modelScores$`TM-Score`[-maxIndex] ~ optimumSimilarityTable[maxIndex, ][-maxIndex]))
  # ggplot(data.frame(NativeTMScore = modelScores$`TM-Score`[-maxIndex],
  #                   ModelTMScore = optimumSimilarityTable[maxIndex, ][-maxIndex]), aes(x=NativeTMScore, y=ModelTMScore)) +
  #   geom_point(size = 2) + geom_smooth(method="lm")
  #
  # ### Max Davies Score
  # summary(lm(modelScores$`TM-Score`[-daviesIndex] ~ optimumSimilarityTable[daviesIndex, ][-daviesIndex]))
  # ggplot(data.frame(NativeTMScore = modelScores$`TM-Score`[-daviesIndex],
  #                   ModelTMScore = optimumSimilarityTable[daviesIndex, ][-daviesIndex]), aes(x=NativeTMScore, y=ModelTMScore)) +
  #   geom_point(size = 2) + geom_smooth(method="lm")
  #
  # ### Max Regression Score
  # summary(lm(modelScores$`TM-Score`[-regressionIndex] ~ optimumSimilarityTable[regressionIndex, ][-regressionIndex]))
  # ggplot(data.frame(NativeTMScore = modelScores$`TM-Score`[-regressionIndex],
  #                   ModelTMScore = optimumSimilarityTable[regressionIndex, ][-regressionIndex]), aes(x=NativeTMScore, y=ModelTMScore)) +
  #   geom_point(size = 2) + geom_smooth(method="lm")
  #
  # ### Max Proq3D Score
  # summary(lm(modelScores$`TM-Score`[-proq3Index] ~ optimumSimilarityTable[proq3Index, ][-proq3Index]))
  # ggplot(data.frame(NativeTMScore = modelScores$`TM-Score`[-proq3Index],
  #                   ModelTMScore = optimumSimilarityTable[proq3Index, ][-proq3Index]), aes(x=NativeTMScore, y=ModelTMScore)) +
  #   geom_point(size = 2) + geom_smooth(method="lm")
  #
  # ## ProQ3 Plots
  #
  # ### ProQ2D
  # ggplot(data.frame(NativeTMScore = modelScores$`TM-Score`, ProQ2D = modelScores$ProQ2D), aes(x=NativeTMScore, y=ProQ2D)) +
  #   geom_point(size = 2)
  #
  # ### ProQ3D
  # ggplot(data.frame(NativeTMScore = modelScores$`TM-Score`, ProQ3D = modelScores$ProQ3D), aes(x=NativeTMScore, y=ProQ3D)) +
  #   geom_point(size = 2)
  #
  # ### ProQRosetta Centroid
  # ggplot(data.frame(NativeTMScore = modelScores$`TM-Score`, ProQRosCenD = modelScores$ProQRosCenD), aes(x=NativeTMScore, y=ProQRosCenD)) +
  #   geom_point(size = 2)
  #
  # ### ProQRosetts Full-Atom
  # ggplot(data.frame(NativeTMScore = modelScores$`TM-Score`, ProQRosFAD = modelScores$ProQRosFAD), aes(x=NativeTMScore, y=ProQRosFAD)) +
  #   geom_point(size = 2)

  # Protocol Termination

  ## Building the constraint lists

  bislist <- rownames(restrictionScores)[order(-restrictionScores$bis)][1:nconst]
  rscorelist <- rownames(restrictionScores)[order(-restrictionScores$rscore)][1:nconst]
  biscorelist <- rownames(restrictionScores)[order(-restrictionScores$biscore)][1:nconst]
  biscore_nativelist <- rownames(restrictionScores)[order(-restrictionScores$biscore_native)][1:nconst]
  biscore_bestlist <- rownames(restrictionScores)[order(-restrictionScores$biscore_best)][1:nconst]
  biscore_regressionlist <- rownames(restrictionScores)[order(-restrictionScores$biscore_regression)][1:nconst]

  ## Writing the constraint files

  appendname <- unlist(strsplit(simulationName, split = "[.]"))
  appendname <- as.numeric(appendname[length(appendname)])

  if ("freq" %in% indicator) {

    write.table(x = write.rosetta.constraints(freqlist),
                file = paste0(outputfile, (appendname+1)),
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)

  }

  if ("bis" %in% indicator) {

    write.table(x = write.rosetta.constraints(bislist),
                file = paste0(outputfile, (appendname+1)),
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)

  }

  if ("rscore" %in% indicator) {

    write.table(x = write.rosetta.constraints(rscorelist),
                file = paste0(outputfile, (appendname+1)),
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)

  }

  if ("biscore" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscorelist),
                file = paste0(outputfile, (appendname+1)),
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  if ("biscore_native" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscore_nativelist),
                file = paste0(outputfile, (appendname+1)),
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  if ("biscore_best" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscore_bestlist),
                file = paste0(outputfile, (appendname+1)),
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  if ("biscore_regression" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscore_regressionlist),
                file = paste0(outputfile, (appendname+1)),
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

}
