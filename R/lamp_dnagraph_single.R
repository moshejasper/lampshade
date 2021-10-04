#' Draw ssDNA of LAMP sequence and primers to fill A5 scale from LAMP frame object
#'
#' @param dnaobj LAMP frame list object
#' @param rowcount Number of bases per row
#' @param size Font size to print bases
#' @param sepr Vertical separation between bases
#' @param family Font family; default is the custom 'DNAMosaic'
#'
#' @return returns ggplot object scaled to A5 dimensions with DNA sequence and LAMP primers
#' @export
#'
#' @examples
#' lframe <- lamp_frame("gattaca", f3 = "gatt", b3 = "tgta")
#' lamp_dnagraph_single(lframe, family = "serif")
lamp_dnagraph_single <- function(dnaobj, rowcount=max_rowlength(nrow(dnaobj$seqdata)), size=5*50/rowcount, 
                                 sepr=.007*size, family = "DNAMosaic"){
  dnaframe <- dnaobj$seqdata
  if (family == "DNAMosaic"){
    dnaframe <- mutate(dnaframe, cmp = .data$fwd)
  }
  sepr <- sepr/2
  if (nrow(dnaframe) %% rowcount == 0) {maxrow <- nrow(dnaframe) %/% rowcount}
  else maxrow <- nrow(dnaframe) %/% rowcount + 1
  dna_mat <- mutate(dnaframe, xgrid = (.data$position-1) %% rowcount + 1, 
                                 ygrid = ((.data$position-1) %/% rowcount + 1))
  dna_mat <- mutate(dna_mat, ygrid = 1 - ((.data$ygrid-0.5) / (maxrow)))
  gg <- ggplot(dna_mat) + aes(x = .data$xgrid, y = .data$ygrid, label = .data$fwd) + 
    geom_text(size = size, family = family) + 
    coord_cartesian(ylim = c(0, 1))+
    theme_void()
  
  if (! is.null(dnaobj$f3)){
    ndna_mat <- filter(dna_mat, .data$position %in% dnaobj$f3pos)
    if (length(unique(ndna_mat$ygrid)) ==2 ){
      yvals <- ndna_mat$ygrid
      ndna_mat1 <- filter(ndna_mat, .data$ygrid == max(yvals))
      gg <- gg + 
        geom_text(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 2), label = .data$fwd), 
                  size = size, family = family, colour = "blue")+
        geom_line(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 3)), 
                  colour = "blue", size = size/10)
      ndna_mat <- filter(ndna_mat, .data$ygrid == min(yvals))
    }
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 2), label = .data$fwd), 
                size = size, family = family, colour = "blue")+
      geom_line(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 3)), 
                colour = "blue", size = size/10, arrow = arrow(length = unit(size, "points")))+
      annotate("text", x = stats::median(ndna_mat$xgrid), y = stats::median(ndna_mat$ygrid) + (sepr * 4),
               label = "F3", size = size, colour = "blue")
  }
  if (! is.null(dnaobj$f2)){
    ndna_mat <- filter(dna_mat, .data$position %in% dnaobj$f2pos)
    if (length(unique(ndna_mat$ygrid)) ==2 ){
      yvals <- ndna_mat$ygrid
      ndna_mat1 <- filter(ndna_mat, .data$ygrid == max(yvals))
      gg <- gg + 
        geom_text(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 2), label = .data$fwd), 
                  size = size, family = family, colour = "blue")+
        geom_line(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 3)), 
                  colour = "blue", size = size/10)
      ndna_mat <- filter(ndna_mat, .data$ygrid == min(yvals))
    }
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 2), label = .data$fwd), 
                size = size, family = family, colour = "blue")+
      geom_line(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 3)), 
                colour = "blue", size = size/10, arrow = arrow(length = unit(size, "points")))+
      annotate("text", x = stats::median(ndna_mat$xgrid), y = stats::median(ndna_mat$ygrid) + (sepr * 4),
               label = "F2", size = size, colour = "blue")
  }
  if (! is.null(dnaobj$lb)){
    ndna_mat <- filter(dna_mat, .data$position %in% dnaobj$lbpos)
    if (length(unique(ndna_mat$ygrid)) ==2 ){
      yvals <- ndna_mat$ygrid
      ndna_mat1 <- filter(ndna_mat, .data$ygrid == max(yvals))
      gg <- gg + 
        geom_text(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 2), label = .data$fwd), 
                  size = size, family = family, colour = "blue")+
        geom_line(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 3)), 
                  colour = "blue", size = size/10)
      ndna_mat <- filter(ndna_mat, .data$ygrid == min(yvals))
    }
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 2), label = .data$fwd), 
                size = size, family = family, colour = "blue")+
      geom_line(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 3)), 
                colour = "blue", size = size/10, arrow = arrow(length = unit(size, "points")))+
      annotate("text", x = stats::median(ndna_mat$xgrid), y = stats::median(ndna_mat$ygrid) + (sepr * 4),
               label = "LB", size = size, colour = "blue")
  }
  if (! is.null(dnaobj$b3)){
    ndna_mat <- filter(dna_mat, .data$position %in% dnaobj$b3pos)
    if (length(unique(ndna_mat$ygrid)) ==2 ){
      yvals <- ndna_mat$ygrid
      ndna_mat1 <- filter(ndna_mat, .data$ygrid == min(yvals))
      gg <- gg + 
        geom_text(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 2), label = .data$cmp), 
                  size = size, family = family, colour = "red")+
        geom_line(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 3.5)), colour = "red", 
                  size = size/10)
      ndna_mat <- filter(ndna_mat, .data$ygrid == max(yvals))
    }
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 2), label = .data$cmp), 
                size = size, family = family, colour = "red")+
      geom_line(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 3.5)), colour = "red", 
                size = size/10, arrow = arrow(length = unit(size, "points"), ends = "first"))+
      annotate("text", x = stats::median(ndna_mat$xgrid), y = stats::median(ndna_mat$ygrid) - (sepr * 4.5),
               label = "B3", size = size, colour = "red")
  }
  if (! is.null(dnaobj$b2)){
    ndna_mat <- filter(dna_mat, .data$position %in% dnaobj$b2pos)
    if (length(unique(ndna_mat$ygrid)) ==2 ){
      yvals <- ndna_mat$ygrid
      ndna_mat1 <- filter(ndna_mat, .data$ygrid == min(yvals))
      gg <- gg + 
        geom_text(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 2), label = .data$cmp), 
                  size = size, family = family, colour = "red")+
        geom_line(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 3.5)), colour = "red", 
                  size = size/10)
      ndna_mat <- filter(ndna_mat, .data$ygrid == max(yvals))
    }
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 2), label = .data$cmp), 
                size = size, family = family, colour = "red")+
      geom_line(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 3.5)), colour = "red", 
                size = size/10, arrow = arrow(length = unit(size, "points"), ends = "first"))+
      annotate("text", x = stats::median(ndna_mat$xgrid), y = stats::median(ndna_mat$ygrid) - (sepr * 4.5),
               label = "B2", size = size, colour = "red")
  }
  if (! is.null(dnaobj$lf)){
    ndna_mat <- filter(dna_mat, .data$position %in% dnaobj$lfpos)
    if (length(unique(ndna_mat$ygrid)) ==2 ){
      yvals <- ndna_mat$ygrid
      ndna_mat1 <- filter(ndna_mat, .data$ygrid == min(yvals))
      gg <- gg + 
        geom_text(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 2), label = .data$cmp), 
                  size = size, family = family, colour = "red")+
        geom_line(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 3.5)), colour = "red", 
                  size = size/10)
      ndna_mat <- filter(ndna_mat, .data$ygrid == max(yvals))
    }
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 2), label = .data$cmp), 
                size = size, family = family, colour = "red")+
      geom_line(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 3.5)), colour = "red", 
                size = size/10, arrow = arrow(length = unit(size, "points"), ends = "first"))+
      annotate("text", x = stats::median(ndna_mat$xgrid), y = stats::median(ndna_mat$ygrid) - (sepr * 4.5),
               label = "LF", size = size, colour = "red")
  }
  if (! is.null(dnaobj$f1c)){
    ndna_mat <- filter(dna_mat, .data$position %in% dnaobj$f1cpos)
    if (length(unique(ndna_mat$ygrid)) ==2 ){
      yvals <- ndna_mat$ygrid
      ndna_mat1 <- filter(ndna_mat, .data$ygrid == max(yvals))
      gg <- gg + 
        geom_text(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 4), label = .data$fwd), 
                  size = size, family = family, colour = "green")+
        geom_line(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 5.5)), colour = "green", 
                  size = size/10)
      ndna_mat <- filter(ndna_mat, .data$ygrid == min(yvals))
    }
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 4), label = .data$fwd), 
                size = size, family = family, colour = "green")+
      geom_line(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid - (sepr * 5.5)), colour = "green", 
                size = size/10)+
      annotate("text", x = stats::median(ndna_mat$xgrid), y = stats::median(ndna_mat$ygrid) - (sepr * 2),
               label = "F1c", size = size, colour = "green")
  }
  if (! is.null(dnaobj$b1c)){
    ndna_mat <- filter(dna_mat, .data$position %in% dnaobj$b1cpos)
    if (length(unique(ndna_mat$ygrid)) ==2 ){
      yvals <- ndna_mat$ygrid
      ndna_mat1 <- filter(ndna_mat, .data$ygrid == min(yvals))
      gg <- gg + 
        geom_text(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 4), label = .data$cmp), 
                  size = size, family = family, colour = "orange")+
        geom_line(data = ndna_mat1, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 5)), colour = "orange", 
                  size = size/10)
      ndna_mat <- filter(ndna_mat, .data$ygrid == max(yvals))
    }
    gg <- gg + 
      geom_text(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 4), label = .data$cmp), 
                size = size, family = family, colour = "orange")+
      geom_line(data = ndna_mat, mapping = aes(x = .data$xgrid, y = .data$ygrid + (sepr * 5)), colour = "orange", 
                size = size/10)+
      annotate("text", x = stats::median(ndna_mat$xgrid), y = stats::median(ndna_mat$ygrid) + (sepr * 2),
               label = "B1c", size = size, colour = "orange")
  }
  if (! is.null(dnaobj$f1c) & ! is.null(dnaobj$f2)){
    f2startpos <- min(dnaobj$f2pos)
    f1cendpos <- min(dnaobj$f1cpos)
    
    f2data <- dna_mat[dna_mat$position == f2startpos,]
    f1cdata <- dna_mat[dna_mat$position == f1cendpos,]
    startx <- f1cdata$xgrid
    starty <- f1cdata$ygrid - (sepr * 5.5)
    endx <- f2data$xgrid
    endy <- f2data$ygrid + (sepr * 3)
    print(endy)
    gg <- gg + 
      annotate("curve", x = startx, y = starty, xend = endx, yend = endy, 
               linetype = 2, curvature = -0.3, colour = "blue", 
               alpha = 0.5, size = size/10)
  }
  if (! is.null(dnaobj$b1c) & ! is.null(dnaobj$b2)){
    b2startpos <- max(dnaobj$b2pos)
    b1cendpos <- max(dnaobj$b1cpos)
    
    b2data <- dna_mat[dna_mat$position == b2startpos,]
    b1cdata <- dna_mat[dna_mat$position == b1cendpos,]
    startx <- b1cdata$xgrid
    starty <- b1cdata$ygrid + (sepr * 5)
    endx <- b2data$xgrid
    endy <- b2data$ygrid - (sepr * 3.5)
    print(endy)
    gg <- gg + 
      annotate("curve", x = startx, y = starty, xend = endx, yend = endy, 
               linetype = 2, curvature = -0.3, colour = "red",  
               alpha = 0.5, size = size/10)
  }
  return(gg)
}