#' Match primer with DNA frame object
#'
#' @param seqframe DNA frame object
#' @param primr sequence to match char
#' @param matchtype "f", "r", "rc" or "c"
#'
#' @return Returns integer vector of sequence match locations
#' @export
#'
#' @examples
#' dframe <- dna_frame("gattaca")
#' primer_match(dframe, "atta", "f")
primer_match <- function(seqframe, primr, matchtype = "f"){
  primr <- stringr::str_to_lower(primr)
  if (matchtype == "f"){
    fseq <- paste0(seqframe$seqdata$fwd, collapse = "")
    matches <- stringr::str_locate_all(fseq, primr)
    for (n in 1:nrow(matches[[1]])){
      match <- matches[[1]][n,]
      match <- match[1]:match[2]
      return(match)
    }
  }
  if (matchtype == "rc"){
    primr <- stringi::stri_reverse(primr)
    rcseq <- paste0(seqframe$seqdata$cmp, collapse = "")
    matches <- stringr::str_locate_all(rcseq, primr)
    for (n in 1:nrow(matches[[1]])){
      match <- matches[[1]][n,]
      match <- match[1]:match[2]
      return(match)
    }
  }
  if (matchtype == "r"){
    primr <- stringi::stri_reverse(primr)
    rseq <- paste0(seqframe$seqdata$fwd, collapse = "") 
    matches <- stringr::str_locate_all(rseq, primr)
    for (n in 1:nrow(matches[[1]])){
      match <- matches[[1]][n,]
      match <- match[1]:match[2]
      return(match)
    }
  }
  if (matchtype == "c"){
    cseq <- paste0(seqframe$seqdata$cmp, collapse = "")
    matches <- stringr::str_locate_all(cseq, primr)
    for (n in 1:nrow(matches[[1]])){
      match <- matches[[1]][n,]
      match <- match[1]:match[2]
      return(match)
    }
  }
  return(NA)
}