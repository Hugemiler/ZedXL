fastaDictionary <- function(query,
                            type = "1to3") {

  lettercode <- c(

    "ALA" = "A", # Alanine
    "ASX" = "B", # Aspartic acid (D) or Asparagine (N)
    "CYS" = "C", # Cysteine
    "ASP" = "D", # Aspartic Acid
    "GLU" = "E", # Glutamic Acid
    "PHE" = "F", # Phenylalanine
    "GLY" = "G", # Glycine
    "HIS" = "H", # Histidine
    "ILE" = "I", # Isoleucine
    "LEX" = "J", # Leucine (L) or Isoleucine (I)
    "LYS" = "K", # Lysine
    "LEU" = "L", # Leucine
    "MET" = "M", # Methionine
    "ASN" = "N", # Asparagine
    "PYL" = "O", # Pyrrolysine
    "PRO" = "P", # Proline
    "GLN" = "Q", # Glutamine
    "ARG" = "R", # Arginine
    "SER" = "S", # Serine
    "THR" = "T", # Threonine
    "SEC" = "U", # Selenocysteine
    "VAL" = "V", # Valine
    "TRP" = "W", # Tryptophan
    "TYR" = "Y", # Tyrosine
    "CLX" = "Z", # Glutamic acid (E) or Glutamine (Q)
    "ANY" = "X", # any
    "STP" = "*", # translation stop
    "GAP" = "-" # gap of indeterminate length
    )

  return(names(lettercode)[which(lettercode == query)])

}
