#' Diagnose standard size requirements for A5 page from sequence
#'
#' @param dnaobj DNA or LAMP frame list object
#' @param rowcount Number of bases per row
#' @param size Font size to print bases
#' @param sepr Vertical separation between bases
#' @param family Font family; default is the custom 'DNAMosaic'
#'
#' @return returns list of key graph parameter settings under 'default' values
#' @export
#'
#' @examples
#' lframe <- dna_frame(random_sequence_generator(500))
#' dnagraph_diagnostic(lframe)
dnagraph_diagnostic <- function(dnaobj, rowcount=max_rowlength(nrow(dnaobj$seqdata)), size=5*50/rowcount, 
                           sepr=.007*size, family = "DNAMosaic"){
  dnaframe <- dnaobj$seqdata
  sepr <- sepr/2
  if (nrow(dnaframe) %% rowcount == 0) {maxrow <- nrow(dnaframe) %/% rowcount}
  else maxrow <- nrow(dnaframe) %/% rowcount + 1
  dna_mat <- mutate(dnaframe, xgrid = (.data$position-1) %% rowcount + 1, 
                                 ygrid = ((.data$position-1) %/% rowcount + 1))
  dna_mat <- mutate(dna_mat, ygrid = 1 - ((.data$ygrid-0.5) / (maxrow)))
  return(list(rowcount = rowcount, size = size, sepr = sepr * 2, maxrow = maxrow))
}
