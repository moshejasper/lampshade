test_that("LAMP frame created properly", {
  frametest <- lamp_frame("gattaca")
  expect_equal(frametest$seqdata$cmp[7], "t")
})