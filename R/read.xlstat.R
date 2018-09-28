read.xlstat <- function(xlstatLogLocation = "/home/guilherme/software/xl_stat/salb3output.csv") {

  xldf <- read.table(xlstatLogLocation)
  colnames(xldf) <- c("restriction",
                      "crys",
                      "avgscore1",
                      "avgscore2",
                      "maxscore1",
                      "maxscore2",
                      "sumscore1",
                      "sumscore2",
                      "nscans",
                      "nspecies")

}
