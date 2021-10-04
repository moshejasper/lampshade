#' Create LAMP Frame DNA object list
#'
#' @param dnastring concatenated string of DNA letters
#' @param f3 charstring - F3 primer
#' @param f2 charstring - F2 primer
#' @param f1c charstring - F1c primer
#' @param b1c charstring - B1c primer
#' @param b2 charstring - B2 primer
#' @param b3 charstring - B3 primer
#' @param lf charstring - LF primer
#' @param lb charstring - LB primer
#'
#' @return LAMP Frame list object
#' @export
#'
#' @examples
#' lamp_frame("gattaca", f3 = "gatt", b3 = "tgta")
lamp_frame <- function(dnastring, f3=NULL, f2=NULL, f1c=NULL, b1c=NULL, b2=NULL, b3=NULL, lf=NULL, lb=NULL){
  dnastring <- stringr::str_to_lower(dnastring)
  dnavect <- strsplit(dnastring, split="")[[1]]
  dnatib <- tibble::tibble(position = 0, fwd = "N", .rows = 0)
  for (n in 1:length(dnavect)){
    dnatib <- tibble::add_row(dnatib, position = n, fwd = dnavect[[n]])
  }
  lframe <- list(seqdata = dplyr::mutate(dnatib, cmp = dna_complement(.data$fwd)), 
                 fprimers = NULL, rcprimers = NULL, rprimers = NULL, cprimers = NULL,
                 f3 = f3, f2 = f2, f1c = f1c, b1c = b1c, b2 = b2, b3 = b3, lf = lf, lb = lb,
                 f3pos = NULL, f2pos = NULL, f1cpos = NULL, b1cpos = NULL, b2pos = NULL, b3pos = NULL,
                 lfpos = NULL, lbpos = NULL)
  if (! is.null(f3)) {
    lframe <- primer_add(lframe, f3, "f")
    lframe$f3pos <- primer_match(lframe, f3, "f")
  }
  if (! is.null(f2)) {
    lframe <- primer_add(lframe, f2, "f")
    lframe$f2pos <- primer_match(lframe, f2, "f")
  }
  if (! is.null(b3)) {
    lframe <- primer_add(lframe, b3, "rc")
    lframe$b3pos <- primer_match(lframe, b3, "rc")
  }
  if (! is.null(b2)) {
    lframe <- primer_add(lframe, b2, "rc")
    lframe$b2pos <- primer_match(lframe, b2, "rc")
  }
  if (! is.null(f1c)) {
    lframe <- primer_add(lframe, f1c, "rc")
    lframe$f1cpos <- primer_match(lframe, f1c, "rc")
  }
  if (! is.null(b1c)) {
    lframe <- primer_add(lframe, b1c, "f")
    lframe$b1cpos <- primer_match(lframe, b1c, "f")
  }
  if (! is.null(lf)) {
    lframe <- primer_add(lframe, lf, "rc")
    lframe$lfpos <- primer_match(lframe, lf, "rc")
  }
  if (! is.null(lb)) {
    lframe <- primer_add(lframe, lb, "f")
    lframe$lbpos <- primer_match(lframe, lb, "f")
  }
  return(lframe)
}