
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lampshade

<!-- badges: start -->
<!-- badges: end -->

The goal of lampshade is to manipulate and visualise short DNA sequences
and associated LAMP primers

## Installation

You can install the released version of lampshade from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("lampshade") # doesn't currently exist
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("moshejasper/lampshade")
```

## Example

Here we create a simple DNA container object

``` r
library(lampshade)
## basic example code

dna_frame(random_sequence_generator(10))
#> $seqdata
#> # A tibble: 10 x 3
#>    position fwd   cmp  
#>       <dbl> <chr> <chr>
#>  1        1 n     n    
#>  2        2 a     t    
#>  3        3 t     a    
#>  4        4 g     c    
#>  5        5 a     t    
#>  6        6 a     t    
#>  7        7 a     t    
#>  8        8 c     g    
#>  9        9 a     t    
#> 10       10 n     n    
#> 
#> $fprimers
#> NULL
#> 
#> $rcprimers
#> NULL
#> 
#> $rprimers
#> NULL
#> 
#> $cprimers
#> NULL
```

Moshe Jasper October 4, 2021
