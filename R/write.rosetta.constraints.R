write.rosetta.constraints <- function(restrictionVector,
                                      table.location = "~/datasets/distance_table.log"){

  distanceTable <- read.table(table.location)

  rownames(distanceTable) <- paste0(distanceTable$V4,
                                    distanceTable$V5,
                                     "-",
                                    distanceTable$V8,
                                    distanceTable$V9)

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

  distanceTable <- subset(distanceTable,
                          rownames(distanceTable) %in% restrictionVector,
                          select = "input")

  return(distanceTable)

}
