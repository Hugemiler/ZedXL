write.rosetta.constraints <- function(restriction.vector,
                                      table.location = "~/Documents/datasets/distance_table.log"){

  distance.table <- read.table(table.location)

  restriction.vector <- gsub("CB", "", restriction.vector)
  restriction.vector <- strsplit(restriction.vector, '-')
  rosetta.df <- list()

  for (i in 1:length(restriction.vector)){

    temp.row <- subset(distance.table, (V4 == restriction.vector[[i]][1] & V8 == restriction.vector[[i]][2]))
    temp.vector <- data.frame("#ATOMPAIR" = "AtomPair",
                              "ATOM1"= "CB",
                              "RESIDUE1" = temp.row[1,4],
                              "ATOM2" = "CB",
                              "RESIDUE2" = temp.row[1,8],
                              "FUNCTION" = "LINEAR_PENALTY",
                              "DIST1" =  temp.row[1,14]/2,
                              "CENTER" = 0,
                              "DIST2" = temp.row[1,14]/2,
                              "TYPE" = 1)

    rosetta.df <- rbind(rosetta.df, temp.vector)
  }

  return(rosetta.df)

}
