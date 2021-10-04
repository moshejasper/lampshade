#' Create DNA frame list object
#'
#' @param charstring concatenated string of DNA letters
#'
#' @return DNA frame list object (not currently defined)
#' @export
#'
#' @examples
#' dna_frame("gattaca")
dna_frame <- function(charstring){
  charstring <- stringr::str_to_lower(charstring)
  charstring <- strsplit(charstring, split="")[[1]]
  dnatib <- tibble::tibble(position = 0, fwd = "N", .rows = 0)
  for (n in 1:length(charstring)){
    dnatib <- tibble::add_row(dnatib, position = n, fwd = charstring[[n]])
  }
  return(list(seqdata = dplyr::mutate(dnatib, cmp = dna_complement(.data$fwd)), 
              fprimers = NULL, rcprimers = NULL, rprimers = NULL, cprimers = NULL))
}