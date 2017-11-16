direct.mirttable <- function(x = xlink.mirttable){
  x <- create.xlink.mirttable(prepare.topolink.logs(read.topolink.output()))
  return(x)
}
