attribute.crys.and.opt <- function(restrictionScores,
                                   cryslistLocation,
                                   optListLocation){

  restriction.list <- rownames(restrictionScores)

  #CRYSTALLOGRAPHIC list
  crystallographic.list <- readLines(cryslistLocation)
  crystallographyc.factor <- restriction.list %in% crystallographic.list
  restrictionScores$crys <- as.logical(crystallographyc.factor)

  #OPTIMAL list (Crystallographic + 5A)
  optimum.list <- readLines(optListLocation)
  optimum.factor <- restriction.list %in% optimum.list
  restrictionScores$opt <- as.logical(optimum.factor)

  return(restrictionScores)
}
