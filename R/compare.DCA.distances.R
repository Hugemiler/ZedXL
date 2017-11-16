compare.DCA.distances <- function(DCA.distances.list,
                                  utility.directory = "~/Documents/datasets/ABCt.DCA1.001/utility/"){

  distance.comparison.table <- data.frame()
  dca.truedistances <- read.table(paste0(utility.directory, "DCA302_truedistances"))$V3
  dca.preddistances <- read.table(paste0(utility.directory, "DCA302_preddist"))$V3

  for (i in 1:length(DCA.distances.list)){

    grep.line <- DCA.distances.list[[i]]
    distance.comparison.table <- rbind(distance.comparison.table, (grep.line <= dca.preddistances))
    }

  colnames(distance.comparison.table) <- names(DCA.distances.list[[1]])
  rownames(distance.comparison.table) <- names(DCA.distances.list)

  return(distance.comparison.table)
}
