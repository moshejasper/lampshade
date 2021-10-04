#' Load PrimerInfo file from Primer Explorer v5
#'
#' @param filename file location of  PrimerInfo file
#'
#' @return returns LAMP frame object with 'core' primers attached
#' @export
#'
pe_load_primerinfo <- function(filename){
  # assume the file is standard
  container = list()
  
  lines <- readr::read_lines(filename)
  
  dnafwd <- strsplit(lines[3], "=")[[1]][2]
  fwd <- strsplit(dnafwd, "")[[1]]
  dnacmp <- paste0(dna_complement(strsplit(dnafwd, "")[[1]]), collapse = "")
  cmp <- strsplit(dnacmp, "")[[1]]
  f3s <- as.integer(strsplit(lines[10], "=")[[1]][2])
  f3e <- as.integer(strsplit(lines[11], "=")[[1]][2])
  f2s <- as.integer(strsplit(lines[12], "=")[[1]][2])
  f2e <- as.integer(strsplit(lines[13], "=")[[1]][2])
  f1cs <- as.integer(strsplit(lines[14], "=")[[1]][2])
  f1ce <- as.integer(strsplit(lines[15], "=")[[1]][2])
  b1cs <- as.integer(strsplit(lines[16], "=")[[1]][2])
  b1ce <- as.integer(strsplit(lines[17], "=")[[1]][2])
  b2s <- as.integer(strsplit(lines[18], "=")[[1]][2])
  b2e <- as.integer(strsplit(lines[19], "=")[[1]][2])
  b3s <- as.integer(strsplit(lines[20], "=")[[1]][2])
  b3e <- as.integer(strsplit(lines[21], "=")[[1]][2])
  
  # now build
  
  f3pos <- f3s:f3e
  f3 <- paste0(fwd[f3pos], collapse = "")
  f2pos <- f2s:f2e
  f2 <- paste0(fwd[f2pos], collapse = "")
  f1cpos <- f1cs:f1ce
  f1c <- stringi::stri_reverse(paste0(cmp[f1cpos], collapse = ""))
  b1cpos <- b1cs:b1ce
  b1c <- paste0(fwd[b1cpos], collapse = "")
  b2pos <- b2s:b2e
  b2 <- stringi::stri_reverse(paste0(cmp[b2pos], collapse = ""))
  b3pos <- b3s:b3e
  b3 <- stringi::stri_reverse(paste0(cmp[b3pos], collapse = ""))
  
  # now assign
  
  
  return(lamp_frame(dnafwd, f3 = f3, f2 = f2, f1c = f1c, b1c = b1c, b2 = b2, b3 = b3))
}