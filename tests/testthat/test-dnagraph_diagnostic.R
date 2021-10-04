set.seed(123)


test_that("Diagnostics run correctly", {
  lframe <- dna_frame(random_sequence_generator(500))
  diagn <- dnagraph_diagnostic(lframe)
  expect_equal(length(diagn), 4)
})