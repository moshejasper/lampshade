#' Setup system for custom DNA graphing
#'
#' @param addfont T/F. Choose to add DNAMosaic custom font to system
#' @param showtext F/F. Choose to deploy showtext 
#'
#' @return invisible
#' @export
#'
dnagraph_setup <- function(addfont = TRUE, showtext = TRUE){
  if (addfont == TRUE){
    sysfonts::font_add("DNAMosaic", paste0(find.package("lampshade"), "/extdata/Dnamosaic-Regular-2.ttf"))
  }
  if (showtext == TRUE){
    showtext::showtext_auto()
  }
}

#' Stop giving showtext control of graphics
#'
#' @return invisible
#' @export
#'
dnagraph_detach <- function(){
  showtext::showtext_auto(enable = FALSE)
}