# quickbuild

  require(ZedXL)
  set.seed(102528)

  modelScores <- compute.model.scores(type = "gscore",
                                      gscoreLogPath,
                                      lovoalignLogPath,
                                      proq3listLocation,
                                      computeproq3)

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

    restrictionScores$biscore_proq3.TM <- -apply(optimumXlinkMirttable, 2, function(x) {
      ltm::biserial.cor(optimumSimilarityTable[which(modelScores$ProQ3D.TM == max(modelScores$ProQ3D.TM)), ], x)})

  }

  restrictionScores <- attribute.crys.and.opt(restrictionScores, cryslistLocation, cryslistLocation)

  # Correlations and Charts

  ## Index Assignment

  maxIndex <- which(modelScores$`TM-Score` == max(modelScores$`TM-Score`))
  daviesIndex <- which(modelScores$`TM-Score` ==
                         modelScores$`TM-Score`[which(modelScores$davisconsensus == max(modelScores$davisconsensus))])

  if(computeproq3 == T) {
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
    biscore_proq3_TMlist <- rownames(restrictionScores)[order(-restrictionScores$biscore_proq3.TM)][1:nconst]
  }

  ## Writing the constraint files

    write.table(x = write.rosetta.constraints(freqlist, table.location = distanceTableLocation),
                file = 'xl_freq', quote = FALSE, col.names = FALSE, row.names = FALSE)

    write.table(x = write.rosetta.constraints(bislist, table.location = distanceTableLocation),
                file = 'xl_bis', quote = FALSE, col.names = FALSE, row.names = FALSE)

    write.table(x = write.rosetta.constraints(biscorelist, table.location = distanceTableLocation),
                file = 'xl_biscre_consensus', quote = FALSE, col.names = FALSE, row.names = FALSE)

    write.table(x = write.rosetta.constraints(biscore_nativelist, table.location = distanceTableLocation),
                file = 'xl_biscore_native', quote = FALSE, col.names = FALSE, row.names = FALSE)

    write.table(x = write.rosetta.constraints(biscore_bestlist, table.location = distanceTableLocation),
                file = 'xl_biscore_best', quote = FALSE, col.names = FALSE, row.names = FALSE)

    write.table(x = write.rosetta.constraints(biscore_proq3_TMlist, table.location = distanceTableLocation),
                file = 'xl_biscore_proq3', quote = FALSE, col.names = FALSE, row.names = FALSE)
