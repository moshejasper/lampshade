test_that("Random sequence generator works", {
  set.seed(125)
  rseq <- random_sequence_generator(100)
  expect_snapshot_output(rseq)
  expect_equal(length(strsplit(rseq, "")[[1]]), 100)
})