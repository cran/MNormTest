test_that("single mean test works", {
  data(iris)
  X <- iris[, 1:4]
  mu0 <- c(5.8, 3.0, 4.3, 1.3)

  expect_path <- normalizePath("./Expect/meanTestsingleExpected.RData")
  load(expect_path)

  test1 <- meanTest.single(X, mu0)
  test2 <- meanTest.single(X, mu0, Sigma0 = diag(1, 4))

  expect_equal(test1, output1)
  expect_equal(test2, output2)
})
