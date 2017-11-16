create.xlink.mirttable <- function(optimizedLogs){
  xlink.mirttable <- data.frame()
  for (i in 1:length(optimizedLogs)){
    grep.line <- optimizedLogs[[i]]$Result

    ## Restriction-type-wise, according to Ferrari et al.

    grep.line <- gsub("OK:FOUND", TRUE, grep.line) # A valid topological distance was found, and it is consistent with all observations.
    grep.line <- gsub("OK:LONG", FALSE, grep.line) # A valid topological distance was found, is it longer than the linker length, and the link was NOT observed, thus the result is consistent with observations.
    grep.line <- gsub("OK:EUCL", FALSE, grep.line) # The euclidean distance is already too long, and the link was NOT observed, thus the result is consistent with observations.
    grep.line <- gsub("OK:NOTFOUND", FALSE, grep.line) # No valid topological distance was found, but the euclidean distance is shorter than the linker length. The link was also NOT observed, such that the result is consistent with observations.
    grep.line <- gsub("BAD:SHORT", TRUE, grep.line) # A topological distance was found, which is shorter than a linker length for which the link was NOT observed. Since it was NOT observed, that linker should be a lower bound to the distance, which is violated.
    grep.line <- gsub("BAD:LONG", FALSE, grep.line) # A topological distance was found, and the link was observed. However, the distance is too long for the linker length. The link, therefore, is not consistent with observations.
    grep.line <- gsub("BAD:EUCL", FALSE, grep.line) # The link was observed, but the euclidean distance in the model is already too long to be consistent with that observation.
    grep.line <- gsub("BAD:NOTFOUND", FALSE, grep.line) # The link was observed, the euclidean distance is fine, but no valid topological distance was found. Thus, the link is not consistent with the observation.
    grep.line <- gsub("BAD:MISSING", TRUE, grep.line) # The link was NOT observed, and the topological distance is such that it should have been.

    grep.line <- as.logical(grep.line)
    grep.line <- as.numeric(grep.line)

    xlink.mirttable <- rbind(xlink.mirttable, grep.line)
  }

  colnames(xlink.mirttable) <- optimizedLogs[[1]]$Restriction
  rownames(xlink.mirttable) <- names(optimizedLogs)

  return(xlink.mirttable)
}
