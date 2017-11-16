prepare.topolink.logs <- function(logTableList){

  # optimizedTableNames <- c("Restriction","Eucldist", "Topodist", "Observed", "Result")
  # optimizedListOfTables <- list()

  for (i in 1:length(logTableList)){

    logTableList[[i]]$Restriction <- paste0(logTableList[[i]]$V4,
                                            logTableList[[i]]$V5,
                                            "-",
                                            logTableList[[i]]$V8,
                                            logTableList[[i]]$V9)
    logTableList[[i]]$Result <- paste0(logTableList[[i]]$V15,
                                       logTableList[[i]]$V16)

    logTableList[[i]] <- subset(logTableList[[i]], select = c("Restriction", "Result"))

     # optimized.list.of.tables[[i]] <- cbind(paste0(topolink.output.list[[i]]$V4,
     #                                              topolink.output.list[[i]]$V5,
     #                                              "-",
     #                                              topolink.output.list[[i]]$V8,
     #                                              topolink.output.list[[i]]$V9),
     #                                       topolink.output.list[[i]]$V10,
     #                                       topolink.output.list[[i]]$V11,
     #                                       topolink.output.list[[i]]$V12,
     #                                       paste0(topolink.output.list[[i]]$V15,
     #                                              topolink.output.list[[i]]$V16))
     # colnames(optimized.list.of.tables[[i]]) <- optimized.table.names
  }
  # names(optimized.list.of.tables) = names(topolink.output.list)
  return(logTableList)
}
