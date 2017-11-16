compute.restriction.parameters <- function(xlink.mirttable){

  digest.xlink.mirttable <- function(){
    all.zeroes <- subset(xlink.mirttable, select = (colSums(xlink.mirttable) == 0))
    all.ones <- subset(xlink.mirttable, select = (colSums(xlink.mirttable) == nrow(xlink.mirttable)))
    digest.list <- c(colnames(all.ones), colnames(all.zeroes))
    xlink.mirtdigested <- subset(xlink.mirttable, select = !(colnames(xlink.mirttable) %in% digest.list))
    return(xlink.mirtdigested)
  }

  xlink.mirtdigested <- digest.xlink.mirttable()

  ##mirt and frequentist analysis

  xlink.mirtObject <- mirt(xlink.mirtdigested, 1, itemtype = '3PL')
  descript.xlink.mirttable <- descript(xlink.mirtdigested)

  #data bind

  restriction.parameters <- as.data.frame(coef(xlink.mirtObject, simplify = TRUE)$items)
  restriction.parameters <- cbind(restriction.parameters, bis = descript.xlink.mirttable$bisCorr)
  restriction.parameters <- cbind(restriction.parameters, descript.xlink.mirttable$perc)

  # Factor Attibution (curated and recovered)

  attribute.cur.and.rec <- function(){
    restriction.list <- rownames(restriction.parameters)
    curated.list <- readLines("~/Documents/datasets/SALB3.22CUR.001/utility/curlist.txt")
    curated.factor <- restriction.list %in% curated.list
    restriction.parameters <- cbind(restriction.parameters, cur = as.factor(curated.factor))
    recovered.list <- readLines("~/Documents/datasets/SALB3.22CUR.001/utility/reclist.txt")
    recovered.factor <- restriction.list %in% recovered.list
    restriction.parameters <- cbind(restriction.parameters, rec = as.factor(recovered.factor))

    colnames(restriction.parameters)[7] <- "right"
    colnames(restriction.parameters)[6] <- "wrong"

    return(restriction.parameters)
  }

  restriction.parameters <- attribute.cur.and.rec()

  #Exploratory Data Analysis

  return(restriction.parameters)
}
