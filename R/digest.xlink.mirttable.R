digest.xlink.mirttable <- function(xlink.mirttable){

  all.zeroes <- subset(xlink.mirttable, select = (colSums(xlink.mirttable) == 0))
  all.ones <- subset(xlink.mirttable, select = (colSums(xlink.mirttable) == nrow(xlink.mirttable)))
  digest.list <- c(colnames(all.ones), colnames(all.zeroes))
  xlink.mirtdigested <- subset(xlink.mirttable, select = !(colnames(xlink.mirttable) %in% digest.list))

  return(xlink.mirtdigested)
}
