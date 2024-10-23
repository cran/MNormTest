test_that("Single Covariance Test works", {
  data(iris)
  X <- iris[, 1:4]

  expect_path <- normalizePath("./Expect/covTestsingleExpected.RData")
  load(expect_path)
  
  test1 <- covTest.single(X, diag(1, 4))
  test2 <- covTest.single(X, diag(1, 4), ball = TRUE)
  test3 <- covTest.single(X, diag(2, 4), ball = TRUE)

  expect_equal(test1, output1, tolerance = 1e-10)
  expect_equal(test2, output2, tolerance = 1e-10)
  expect_equal(test3, output3, tolerance = 1e-10)
})
