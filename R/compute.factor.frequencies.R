compute.factor.frequencies <- function(optimized.list.of.tables){

  factor.list <- vector()

  for (i in 1:length(optimized.list.of.tables)){

    factor.list <- c(factor.list, as.character(as.data.frame(optimized.list.of.tables[[i]])$Result))

  }

  return(factor.list)

}

# ggplot(data.frame(factor.list),aes(factor.list)) + geom_bar()
