#' Creates complement of DNA sequence (as vector)
#'
#' @param chrstring character vector describing DNA sequence
#'
#' @return Returns character vector of DNA sequence complement
#' @export
#'
#' @examples
#' dna_complement(c("g", "a", "t", "t", "a", "c", "a"))
dna_complement <- function(chrstring){
  newstr <- chrstring
  for (n in 1:length(chrstring)){
    chr <- chrstring[n]
    if (chr == "t" | chr == "T") newstr[n] <- "a"
    if (chr == "a" | chr == "A") newstr[n] <- "t"
    if (chr == "g" | chr == "G") newstr[n] <- "c"
    if (chr == "c" | chr == "C") newstr[n] <- "g"
    if (chr == "n" | chr == "N") newstr[n] <- "n"
    if (chr == "r" | chr == "R") newstr[n] <- "y"
    if (chr == "y" | chr == "Y") newstr[n] <- "r"
    if (chr == "s" | chr == "S") newstr[n] <- "s"
    if (chr == "w" | chr == "W") newstr[n] <- "w"
    if (chr == "k" | chr == "K") newstr[n] <- "m"
    if (chr == "m" | chr == "M") newstr[n] <- "k"
    if (chr == "b" | chr == "B") newstr[n] <- "v"
    if (chr == "v" | chr == "V") newstr[n] <- "b"
    if (chr == "d" | chr == "D") newstr[n] <- "h"
    if (chr == "h" | chr == "H") newstr[n] <- "d"
  }
  return(newstr)
}