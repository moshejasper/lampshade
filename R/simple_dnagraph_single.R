#' Draw ssDNA to fill A5 scale from DNA frame object
#'
#' @param dnaobj DNA frame list object
#' @param rowcount Number of bases per row
#' @param size Font size to print bases
#' @param sepr Vertical separation between bases
#' @param family Font family; default is the custom 'DNAMosaic'
#'
#' @return returns ggplot object scaled to A5 dimensions with DNA sequence drawn.
#' @export
#'
#' @examples
#' dframe <- dna_frame("gattaca")
#' simple_dnagraph_single(dframe, family = "serif")
simple_dnagraph_single <- function(dnaobj, rowcount=max_rowlength(nrow(dnaobj$seqdata)), size=5*50/rowcount, 
                                   sepr=0.035 * 50 / rowcount, family = "DNAMosaic"){
  dnaframe <- dnaobj$seqdata
  sepr <- sepr/2
  if (nrow(dnaframe) %% rowcount == 0) {maxrow <- nrow(dnaframe) %/% rowcount}
  else maxrow <- nrow(dnaframe) %/% rowcount + 1
  dna_mat <- dplyr::mutate(dnaframe, xgrid = (.data$position-1) %% rowcount + 1, 
                           ygrid = ((.data$position-1) %/% rowcount + 1))
  dna_mat <- dplyr::mutate(dna_mat, ygrid = 1 - ((.data$ygrid-0.5) / (maxrow)))
  gg <- ggplot(dna_mat) + aes(x = .data$xgrid, y = .data$ygrid, label = .data$fwd) + 
    geom_text(size = size, family = family) + 
    coord_cartesian(ylim = c(0, 1))+
    theme_void()
  
  if (! is.null(dnaobj$fprimers)){
    ndna_mat <- dplyr::filter(dna_mat, .data$position %in% dnaobj$fprimers)
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 2), label = .data$fwd), 
                size = size, family = family, colour = "blue")
  }
  if (! is.null(dnaobj$rcprimers)){
    ndna_mat <- dplyr::filter(dna_mat, .data$position %in% dnaobj$rcprimers)
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 2), label = .data$fwd), 
                size = size, family = family, colour = "red")
  }
  if (! is.null(dnaobj$rprimers)){
    ndna_mat <- dplyr::filter(dna_mat, .data$position %in% dnaobj$rprimers)
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 4), label = .data$fwd), 
                size = size, family = family, colour = "orange")
  }
  if (! is.null(dnaobj$cprimers)){
    ndna_mat <- dplyr::filter(dna_mat, .data$position %in% dnaobj$cprimers)
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 4), label = .data$fwd), 
                size = size, family = family, colour = "green")
  }
  return(gg)
}