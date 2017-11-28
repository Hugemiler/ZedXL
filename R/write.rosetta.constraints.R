write.rosetta.constraints <- function(restriction.vector,
                                      table.location = "~/datasets/distance_table.log"){

  distance.table <- read.table(table.location)

  rownames(distance.table) <- paste0(distance.table$V4,
                                     distance.table$V5,
                                     "-",
                                     distance.table$V8,
                                     distance.table$V9)

    for (i in 1:length(restriction.vector)){

      distance.table[which(rownames(distance.table) == restriction.vector[i])]

    temp.row <- subset(distance.table, (V4 == restriction.vector[[i]][1] & V8 == restriction.vector[[i]][2]))
    relevant <- c(temp.row[,4], temp.row[,8], temp.row[,13])
    temp.vector <- data.frame("#ATOMPAIR" = "AtomPair",
                              "ATOM1"= "CB",
                              "RESIDUE1" = relevant[1],
                              "ATOM2" = "CB",
                              "RESIDUE2" = relevant[2],
                              "FUNCTION" = "LINEAR_PENALTY",
                              "DIST1" =  as.numeric(relevant[3])/2,
                              "CENTER" = 0,
                              "DIST2" = as.numeric(relevant[3])/2,
                              "TYPE" = 1)

    rosetta.df <- rbind(rosetta.df, temp.vector)
  }

  return(rosetta.df)

}
