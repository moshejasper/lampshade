#' Add primer to DNA frame object
#'
#' @param seqframe DNA frame list object
#' @param primr String sequence of Primer
#' @param matchtype Matchtype 'f', 'r', 'c', or 'rc'. Must be exact. 
#'
#' @return Returns seqframe list with added primers. 
#' @export
#'
#' @examples
#' dframe <- dna_frame("gattaca")
#' primer_add(dframe, "atta", "f")
primer_add <- function(seqframe, primr, matchtype = "f"){
  primr <- stringr::str_to_lower(primr)
  
  if (matchtype == "f"){
    fseq <- paste0(seqframe$seqdata$fwd, collapse = "")
    matches <- stringr::str_locate_all(fseq, primr)
    for (n in 1:nrow(matches[[1]])){
      match <- matches[[1]][n,]
      match <- match[1]:match[2]
      seqframe$fprimers <- unique(c(seqframe$fprimers, match))
    }
  }
  if (matchtype == "rc"){
    primr <- stringi::stri_reverse(primr)
    rcseq <- paste0(seqframe$seqdata$cmp, collapse = "")
    matches <- stringr::str_locate_all(rcseq, primr)
    for (n in 1:nrow(matches[[1]])){
      match <- matches[[1]][n,]
      match <- match[1]:match[2]
      seqframe$rcprimers <- unique(c(seqframe$rcprimers, match))
    }
  }
  if (matchtype == "r"){
    primr <- stringi::stri_reverse(primr)
    rseq <- paste0(seqframe$seqdata$fwd, collapse = "") 
    matches <- stringr::str_locate_all(rseq, primr)
    for (n in 1:nrow(matches[[1]])){
      match <- matches[[1]][n,]
      match <- match[1]:match[2]
      seqframe$rprimers <- unique(c(seqframe$rprimers, match))
    }
  }
  if (matchtype == "c"){
    cseq <- paste0(seqframe$seqdata$cmp, collapse = "")
    matches <- stringr::str_locate_all(cseq, primr)
    for (n in 1:nrow(matches[[1]])){
      match <- matches[[1]][n,]
      match <- match[1]:match[2]
      seqframe$cprimers <- unique(c(seqframe$cprimers, match))
    }
  }
  return(seqframe)
}