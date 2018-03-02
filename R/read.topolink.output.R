read.topolink.output <- function(mode = "observed",
                                 simulationName,
                                 loglistLocation = paste0("~/datasets/",simulationName,"/loglist.txt"),
                                 logsLocation = paste0("~/datasets/",simulationName,"/topolink_",mode,"/")){

  ######
  #
  # This script relies on logs obtained by TopoLink Software. TopoLink is open-source software.
  # Whenever you use TopoLink software, please cite:
  #
  # Ferrari, Allan & Martínez, Leandro & Cesar Gozzo, Fabio. (2017).
  # TopoLink: A software to validate structural models using chemical crosslinking constraints.
  # Nature Protocol Exchange. 10.1038/protex.2017.035.
  #
  # Obtain TopoLink software at http://leandro.iqm.unicamp.br/topolink/home.shtml
  #
  ######

  logList <- readLines(loglistLocation)
  logTableList <- list()

  for (i in 1:length(logList)){

    logPath <- paste0(logsLocation, logList[i])
    logTableList[[logList[i]]] <- readLines(logPath)
    filterVector <- grep(pattern = 'LINK:', x = logTableList[[logList[i]]])

    logTableList[[logList[i]]] <- read.table(file = logPath,
                                             header = FALSE,
                                             skip = (min(filterVector)-1),
                                             nrows = length(filterVector))

    # logTableList[[logList[i]]] <- as.data.frame(logTableList[[logList[i]]][filterVector])

    }
  names(logTableList) = logList
  return(logTableList)
}
