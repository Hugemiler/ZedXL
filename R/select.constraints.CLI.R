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

  bestcstcol <- optimumXlinkMirttable[which(modelScores$`TM-Score` == max(modelScores$`TM-Score`)), ]
  bestcstlist <- colnames(bestcstcol)[bestcstcol == 1]

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

  if(computeproq3 == T) {

    restrictionScores$biscore_proq3 <- -apply(optimumXlinkMirttable, 2, function(x) {
      ltm::biserial.cor(optimumSimilarityTable[which(modelScores$ProQ3D == max(modelScores$ProQ3D)), ], x)})

    restrictionScores$biscore_proq3.TM <- -apply(optimumXlinkMirttable, 2, function(x) {
      ltm::biserial.cor(optimumSimilarityTable[which(modelScores$ProQ3D.TM == max(modelScores$ProQ3D.TM)), ], x)})

  }

  restrictionScores <- attribute.crys.and.opt(restrictionScores, cryslistLocation, optlistLocation)

  # Correlations and Charts

  ## Index Assignment

  maxIndex <- which(modelScores$`TM-Score` == max(modelScores$`TM-Score`))
  daviesIndex <- which(modelScores$`TM-Score` ==
                         modelScores$`TM-Score`[which(modelScores$davisconsensus == max(modelScores$davisconsensus))])

  if(computeproq3 == T) {
    proq3Index <- which(modelScores$ProQ3D == max(modelScores$ProQ3D))
    proq3_TMIndex <- which(modelScores$ProQ3D.TM == max(modelScores$ProQ3D.TM))
  }
  # Protocol Termination

  ## Building the constraint lists

  bislist <- rownames(restrictionScores)[order(-restrictionScores$bis)][1:nconst]
  rscorelist <- rownames(restrictionScores)[order(-restrictionScores$rscore)][1:nconst]
  biscorelist <- rownames(restrictionScores)[order(-restrictionScores$biscore)][1:nconst]
  biscore_nativelist <- rownames(restrictionScores)[order(-restrictionScores$biscore_native)][1:nconst]
  biscore_bestlist <- rownames(restrictionScores)[order(-restrictionScores$biscore_best)][1:nconst]
  if(computeproq3 == T) {
    biscore_proq3list <- rownames(restrictionScores)[order(-restrictionScores$biscore_proq3)][1:nconst]
    biscore_proq3_TMlist <- rownames(restrictionScores)[order(-restrictionScores$biscore_proq3.TM)][1:nconst]
  }

  ## Writing the constraint files

  if ("freq" %in% indicator) {

    write.table(x = write.rosetta.constraints(freqlist, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)

  }

  if ("bestcst" %in% indicator) {

    write.table(x = write.rosetta.constraints(bestcstlist, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  if ("bis" %in% indicator) {

    write.table(x = write.rosetta.constraints(bislist, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)

  }

  if ("rscore" %in% indicator) {

    write.table(x = write.rosetta.constraints(rscorelist, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)

  }

  if ("biscore" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscorelist, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  if ("biscore_native" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscore_nativelist, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  if ("biscore_best" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscore_bestlist, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  if ("biscore_regression" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscore_regressionlist, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  if ("biscore_proq3" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscore_proq3list, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  if ("biscore_proq3.TM" %in% indicator) {

    write.table(x = write.rosetta.constraints(biscore_proq3_TMlist, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

}
