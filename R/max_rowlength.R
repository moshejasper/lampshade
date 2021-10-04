#' Max rowlength
#'
#' @param n Number of bases in DNA sequence
#'
#' @return Returns integer number of bases in a DNA row (currently under A5 page assumptions)
#' @export
#'
#' @examples
#' max_rowlength(100)
max_rowlength <- function(n){
  rnum <- sqrt(n / 12.5) * 12.5
  return(rnum %/% 1)
}