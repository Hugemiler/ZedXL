restriction.differential.scores <- function(type = "TM-Score",
                                            xlinkMirtTable,
                                            modelScores) {

  scoreTrue <- apply(xlinkMirtTable, 2, function(x) {
    if (type == "TM-Score") {mean(modelScores[which(x == 1),][,type])}
    else if (type == "G-Score") {mean(1/modelScores[which(x == 1),][,type])}
    else if (type == "Wdegree") {mean(modelScores[which(x == 1),][,type])}})

  devTrue <- apply(xlinkMirtTable, 2, function(x) {
                                sd(modelScores[which(x == 1),][,type])})

  freqTrue <- apply(xlinkMirtTable, 2, function(x) {
    length(modelScores[which(x == 1),][,type])/length(x)})

  scoreFalse <- apply(xlinkMirtTable, 2, function(x) {
    if (type == "TM-Score") {mean(modelScores[which(x == 0),][,type])}
    else if (type == "G-Score") {mean(1/modelScores[which(x == 0),][,type])}
    else if (type == "Wdegree") {mean(modelScores[which(x == 0),][,type])}})

  devFalse <- apply(xlinkMirtTable, 2, function(x) {
    sd(modelScores[which(x == 0),][,type])})

  freqFalse <- apply(xlinkMirtTable, 2, function(x) {
    length(modelScores[which(x == 0),][,type])/length(x)})

  ## Uncomment columns as needed for Exploratory Data Analysis.

  logit <- log(freqTrue/(1 - freqTrue))

  rscorescale <- function(x){(x-min(x))/(max(x)-min(x))}

  scaledLogit <- rscorescale(logit)

  differentialScores <- data.frame(
    "trueScore" = scoreTrue,
    # "trueDev" = deviationTrueScore,
    # "trueCount" = countTrueScore,
    "falseScore" = scoreFalse,
    # "falseDev" = deviationFalseScore,
    # "falseCount" = countFalseScore,
    # "rscore" = (averageTrueScore/averageFalseScore),
    "logscore" = logit,
    "rscore" = scaledLogit + rscorescale(scoreTrue - scoreFalse))

  return(differentialScores)
}
