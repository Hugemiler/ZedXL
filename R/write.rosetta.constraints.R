write.rosetta.constraints <- function(restrictionVector,
                                      table.location = NULL,
                                      mode = "xlff"){

  distanceTable <- read.table(table.location)

  rownames(distanceTable) <- paste0(distanceTable$V4,
                                    distanceTable$V5,
                                     "-",
                                    distanceTable$V8,
                                    distanceTable$V9)

  if (mode == "flat_linear") {
    distanceTable$input <- apply(distanceTable, 1, function(line) {
      paste("AtomPair",
            line[5],
            line[4],
            line[9],
            line[8],
            "LINEAR_PENALTY",
            as.numeric(line[14])/2,
            0,
            as.numeric(line[14])/2,
            1,
            collapse = " ")})

  } else if (mode == "xlff") {

    distanceTable$input <- apply(distanceTable, 1, function(line) {
      paste("observed",
            line[2],
            line[3],
            line[4],
            line[6],
            line[7],
            line[8],
            collapse = " ")})

  }

  distanceTable <- subset(distanceTable,
                          rownames(distanceTable) %in% restrictionVector,
                          select = "input")

  return(distanceTable)

}
