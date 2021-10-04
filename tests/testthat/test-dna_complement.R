test_that("DNA complement works as intended",{
          expect_equal(dna_complement(c("g", "a", "t", "t", "a", "c", "a")), c("c","t","a","a","t","g","t"))
})