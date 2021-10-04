#' Draw dsDNA to fill A5 scale from DNA frame object
#'
#' @param dnaobj DNA frame list object
#' @param rowcount Number of bases per row
#' @param size Font size to print bases
#' @param sepr Vertical separation between bases
#' @param family Font family; default is the custom 'DNAMosaic'
#' @param upper Whether to print DNA in UPPERCASE (default FALSE)
#'
#' @return returns ggplot object scaled to A5 dimensions with DNA sequence drawn. 
#' @export
#'
#' @examples
#' dframe <- dna_frame("gattaca")
#' simple_dnagraph_double(dframe, family = "serif")
simple_dnagraph_double <- function(dnaobj, rowcount=max_rowlength(nrow(dnaobj$seqdata)), size=5*50/rowcount, 
                            sepr=0.035 * 50 / rowcount, family = "DNAMosaic", upper = FALSE){
  
  
  dnaframe <- dnaobj$seqdata
  if (upper){
    dnaframe <- dplyr::mutate(dnaframe, fwd = stringr::str_to_upper(.data$fwd), 
                              cmp = stringr::str_to_upper(.data$cmp))
    size <- size * 4 / 5
  }
  sepr <- sepr/2
  sepr2 <- sepr*2
  if (! family == "DNAMosaic") {
    sepr2 <- sepr
    sepr <- sepr * 2/3
  }
  
  if (nrow(dnaframe) %% rowcount == 0) {maxrow <- nrow(dnaframe) %/% rowcount}
  else maxrow <- nrow(dnaframe) %/% rowcount + 1
  dna_mat <- dplyr::mutate(dnaframe, xgrid = (.data$position-1) %% rowcount + 1, 
                                 ygrid = ((.data$position-1) %/% rowcount + 1))
  dna_mat <- dplyr::mutate(dna_mat, ygrid = 1 - ((.data$ygrid-0.5) / (maxrow)))
  gg <- ggplot(dna_mat) + aes(x = .data$xgrid, y = .data$ygrid+sepr, label = .data$fwd) + 
    geom_text(size = size, family = family) + 
    geom_text(mapping = aes(label = .data$cmp, y = .data$ygrid -sepr), size = size, family = family)+
    coord_cartesian(ylim = c(0, 1))+
    theme_void()
  
  if (! is.null(dnaobj$fprimers)){
    ndna_mat <- dplyr::filter(dna_mat, .data$position %in% dnaobj$fprimers)
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 3), label = .data$fwd), 
                size = size, family = family, colour = "blue")
  }
  if (! is.null(dnaobj$rcprimers)){
    ndna_mat <- dplyr::filter(dna_mat, .data$position %in% dnaobj$rcprimers)
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 3), label = .data$cmp), 
                size = size, family = family, colour = "red")
  }
  if (! is.null(dnaobj$rprimers)){
    ndna_mat <- dplyr::filter(dna_mat, .data$position %in% dnaobj$rprimers)
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 3 + sepr2), label = .data$fwd), 
                size = size, family = family, colour = "orange")
  }
  if (! is.null(dnaobj$cprimers)){
    ndna_mat <- dplyr::filter(dna_mat, .data$position %in% dnaobj$cprimers)
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 3 + sepr2), label = .data$cmp), 
                size = size, family = family, colour = "green")
  }
  return(gg)
}