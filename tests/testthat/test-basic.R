test_that("LRDE basic workflow runs", {
  set.seed(123)

  mat <- matrix(rnbinom(600, size = 5, mu = 5), nrow = 100)
  grp <- factor(rep(c("A", "B"), each = 3))

  y <- prepareDGE(mat, grp)
  y <- sizeFactorsEst(y)
  y <- tagwiseEst(y)

  expect_true(!is.null(y$tagwise.disp))
  expect_equal(length(y$tagwise.disp), nrow(mat))

  y <- hurdle.LRT(y)
  y <- hurdle.Wald.Test(y)

  expect_true(!is.null(y$p.values))
  expect_equal(length(y$p.values), nrow(mat))
})
