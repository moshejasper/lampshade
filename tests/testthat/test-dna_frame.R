
test_that("DNA frame created properly", {
  frametest <- dna_frame("gattaca")
  expect_equal(frametest$seqdata$cmp[7], "t")
})