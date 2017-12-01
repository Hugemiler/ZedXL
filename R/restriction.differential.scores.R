restriction.differential.scores <- function(type = "TM-Score",
                                            xlinkMirtTable,
                                            modelScores) {

  averageTrueScore <- apply(xlinkMirtTable, 2, function(x) {
    if (type == "TM-Score") {mean(modelScores[which(x == 1),][,type])}
    else {1/mean(modelScores[which(x == 1),][,type])}})

  deviationTrueScore <- apply(xlinkMirtTable, 2, function(x) {
                                sd(modelScores[which(x == 1),][,type])})

  countTrueScore <- apply(xlinkMirtTable, 2, function(x) {
                            length(modelScores[which(x == 1),][,type])})

  averageFalseScore <- apply(xlinkMirtTable, 2, function(x) {
                              if (type == "TM-Score") {mean(modelScores[which(x == 0),][,type])}
                              else {1/mean(modelScores[which(x == 0),][,type])}})

  deviationFalseScore <- apply(xlinkMirtTable, 2, function(x) {
                                sd(modelScores[which(x == 0),][,type])})

  countFalseScore <- apply(xlinkMirtTable, 2, function(x) {
                            length(modelScores[which(x == 0),][,type])})

  ## Uncomment columns as needed for Exploratory Data Analysis.

  differentialScores <- data.frame(
    "trueScore" = averageTrueScore,
    # "trueDev" = deviationTrueScore,
    # "trueCount" = countTrueScore,
    "falseScore" = averageFalseScore,
    # "falseDev" = deviationFalseScore,
    # "falseCount" = countFalseScore,
    # "rscore" = (averageTrueScore/averageFalseScore),
    "adjustedRscore" = (averageTrueScore*countTrueScore)/(averageFalseScore*countFalseScore),
    "logRscore" = log((averageTrueScore*countTrueScore)/(averageFalseScore*countFalseScore)))

  return(differentialScores)
}
