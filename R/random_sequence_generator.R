#' Generate random DNA sequence of length n
#'
#' @param n Length of sequence to generate (integer)
#'
#' @return Charstring of length n with random DNA letters (including nonstandard)
#' @export
#'
#' @examples
#' random_sequence_generator(10)
random_sequence_generator <- function(n){
  nb <- n%/% 5
  startval <- c(replicate(nb, "a"), replicate(nb, "t"), replicate(nb, "g"), 
                replicate(nb, "c"), replicate(nb%/% 6, c("s","w","y","r","k","m")))
  startval <- c(startval, replicate(n - length(startval), "n"))
  return(paste0(sample(startval, n, TRUE), collapse = ""))
}