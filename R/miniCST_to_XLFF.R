######
# For debug purposes:
# fastaObject <- read.fasta(file = "~/[Backup]SAFE/papers/2019_DCA_Ricardo/FASTA/1amm.fasta", seqtype = "AA")
# cst_table <- read.table("~/[Backup]SAFE/papers/2019_DCA_Ricardo/CST/1AMM_XL.dat")
#
# observed LYS A 123 SER A 14 short
# observed SER A 15 ASP A 24 zl
# observed GLU A 17 ASP A 44 long
######

miniCST_to_XLFF <- function(fastaObject,
                            cst_table) {

  buildList <- apply(cst_table, 1, function(ll) {

    paste0(

      "observed ",
      fastaDictionary(fastaObject[[1]][ll[1]]),
      " A ",
      ll[1],
      " ",
      fastaDictionary(fastaObject[[1]][ll[2]]),
      " A ",
      ll[2]
      )
  }
  )

    return(buildList)

}
