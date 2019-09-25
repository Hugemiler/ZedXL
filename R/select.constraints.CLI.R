select.constraints.CLI <- function(inputfile, outputfile = "xl") {

  source(inputfile)

  ###############
  # Script for the command-line selection methodology
  # Sep 2019
  #
  # Guilherme Fahur Bottino
  ###############

  #####
  #  0. Initial Package Load and Environment Set
  #####

  require(ZedXL)
  set.seed(102528)

  #####
  #  1. Parsing input data
  #####

  ## 1.1 Model scores (alignment against crystallographix structure, G-score (optional) and ProQ3D-score (optional)).

  modelScores <- compute.model.scores(lovoalignLogPath,
                                      computegscore = F,
                                      computeproq3 = ("BISCORE_PROQ3D_TM" %in% indicator | "BISCORE_PROQ3D" %in% indicator))

  ## 1.2 Computing the binary N x M (Models x Constraints) Matrix.

  optimumXlinkMirttable <- create.xlink.mirttable(
    prepare.topolink.logs(
      read.topolink.output(loglistLocation,
                           topolinkLogsDirectory)
    )
  )

  ### 1.2.1. Extracting the most frequent XLs before digestion.

  if ("FREQ" %in% indicator) {

    FREQ_constraintList <- colnames(optimumXlinkMirttable)[order(-colSums(optimumXlinkMirttable))][1:nconst]

  }

  ### 1.2.2. Extracting the list of constraints from the best.

  if ("BEST" %in% indicator) {

    bestcstcol <- optimumXlinkMirttable[which(modelScores$`TM-Score` == max(modelScores$`TM-Score`)), ]
    BEST_constraintList <- colnames(bestcstcol)[bestcstcol == 1]

  }

  ## 1.3 Computing the similarity Matrix from the nxn alignments.

  if ("BISCORE_BEST" %in% indicator | "BISCORE_CONSENSUS" %in% indicator) {
    optimumSimilarityTable <- create.dissimilarity.matrix(mode = "similarity",
                                                          diagonal = 1,
                                                          alignlistPath,
                                                          compactlogPath,
                                                          nmodels = nrow(optimumXlinkMirttable))

    ### 1.3.1 Appending Davis-QA-Consensus to modelScores table

    modelScores$davisconsensus <- imperialist.score(mode = "consensus",
                                                    models = rownames(optimumSimilarityTable),
                                                    similarityTable = optimumSimilarityTable)

  }

  #####
  #  2. Analyzing Data
  #####

  ###
  ## 2.1 Digesting the optimumXlinkMirttable (removing 0% and 100% constraints because they don't have variance).

  optimumXlinkMirttable <- digest.xlink.mirttable(optimumXlinkMirttable)

  ###
  ## 2.2 Descriptive Statistics (employing LTM package's descript class object).

  optimumDescript <- ltm::descript(optimumXlinkMirttable, chi.squared = F)

  ###
  ## 2.3 Initializing the restrictionScores data frame.

  restrictionScores <- data.frame(bis = optimumDescript$bisCorr,
                                  freq = optimumDescript$perc[,2],
                                  logit = optimumDescript$perc[,3])

  ### 2.3.1. Listing the top BIS constraints
  BIS_constraintList <- rownames(restrictionScores)[order(-restrictionScores$bis)][1:nconst]

  ###
  ## 2.4 Calculating additional scores, according to inputfile.

  ### 2.4.1.

  if ("BISCORE_CONSENSUS" %in% indicator) {
    restrictionScores$biscore_consensus <- -apply(optimumXlinkMirttable, 2, function(x) {
      ltm::biserial.cor(optimumSimilarityTable[which(modelScores$davisconsensus == max(modelScores$davisconsensus)), ], x)})

    BISCORE_CONSENSUS_constraintList <- rownames(restrictionScores)[order(-restrictionScores$biscore_consensus)][1:nconst]
  }

  ### 2.4.2.

  if ("BISCORE_BEST" %in% indicator) {
    restrictionScores$biscore_best <- -apply(optimumXlinkMirttable, 2, function(x) {
      ltm::biserial.cor(optimumSimilarityTable[which(modelScores$`TM-Score` == max(modelScores$`TM-Score`)), ], x)})

    BISCORE_BEST_constraintList <- rownames(restrictionScores)[order(-restrictionScores$biscore_best)][1:nconst]
  }

  ### 2.4.3.

  if ("BISCORE_CRYS" %in% indicator) {
    restrictionScores$biscore_crys <- -apply(optimumXlinkMirttable, 2, function(x) {
      ltm::biserial.cor(modelScores$`TM-Score`, x)})

    BISCORE_CRYS_constraintList <- rownames(restrictionScores)[order(-restrictionScores$biscore_crys)][1:nconst]
  }

  ### 2.4.4.

  if ("BISCORE_PROQ3D" %in% indicator) {
    restrictionScores$biscore_proq3 <- -apply(optimumXlinkMirttable, 2, function(x) {
      ltm::biserial.cor(optimumSimilarityTable[which(modelScores$ProQ3D == max(modelScores$ProQ3D)), ], x)})

    BISCORE_PROQ3D_constraintList <- rownames(restrictionScores)[order(-restrictionScores$biscore_proq3)][1:nconst]
  }

  ### 2.4.5.

  if ("BISCORE_PROQ3D_TM" %in% indicator) {
    restrictionScores$biscore_proq3.TM <- -apply(optimumXlinkMirttable, 2, function(x) {
      ltm::biserial.cor(optimumSimilarityTable[which(modelScores$ProQ3D.TM == max(modelScores$ProQ3D.TM)), ], x)})

    BISCORE_PROQ3D_TM_constraintList <- rownames(restrictionScores)[order(-restrictionScores$biscore_proq3.TM)][1:nconst]
  }

  ###
  ## 2.5 Appending crystallographic constraint list, if available.

  restrictionScores <- attribute.crys.and.opt(restrictionScores, cryslistLocation, optlistLocation)

  ###
  ## 2.6 Exporting constraint scores table (uncomment to export).
  #
  # write.csv(restrictionScores, file = './constraintScores.csv')

  #####
  #  3. Outputting constraint files
  #####

  ##
  # 3.1 Frequency

  if ("FREQ" %in% indicator) {

    write.table(x = write.rosetta.constraints(FREQ_constraintList, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)

  }

  ##
  # 3.2 Constraints from the best model

  if ("BESTCST" %in% indicator) {

    write.table(x = write.rosetta.constraints(BEST_constraintL, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)

  }

  ##
  # 3.3 Point-biserial correlation of each constraint versus the total number of constraints per model

  if ("BISCORE_TOTALCST" %in% indicator) {

    write.table(x = write.rosetta.constraints(BIS_constraintList, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)

  }

  ##
  # 3.4 Point-biserial correlation of each constraint versus the alignment against the largest Davis Consensus score model

  if ("BISCORE_CONSENSUS" %in% indicator) {

    write.table(x = write.rosetta.constraints(BISCORE_CONSENSUS_constraintList, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  ##
  # 3.5 Point-biserial correlation of each constraint versus the alignment against the crystallographic structure

  if ("BISCORE_CRYS" %in% indicator) {

    write.table(x = write.rosetta.constraints(BISCORE_CRYS_constraintList, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  ##
  # 3.6 Point-biserial correlation of each constraint versus the alignment against the model most similar to the crystallographic structure

  if ("BISCORE_BEST" %in% indicator) {

    write.table(x = write.rosetta.constraints(BISCORE_BEST_constraintList, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  ##
  # 3.7 Point-biserial correlation of each constraint versus the alignment against the largest Proq3-sscore model

  if ("BISCORE_PROQ3" %in% indicator) {

    write.table(x = write.rosetta.constraints(BISCORE_PROQ3_constraintList, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

  ##
  # 3.8 Point-biserial correlation of each constraint versus the alignment against the largest Proq3-TMscore model

  if ("BISCORE_PROQ3_TM" %in% indicator) {

    write.table(x = write.rosetta.constraints(BISCORE_PROQ3_TM_constraintList, table.location = distanceTableLocation),
                file = outputfile,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }

}
