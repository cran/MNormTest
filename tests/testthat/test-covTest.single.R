test_that("Single Covariance Test works", {
  data(iris)
  X <- iris[, 1:4]

  expect_path <- normalizePath("./Expect/covTestsingleExpected.RData")
  load(expect_path)
  
  test1 <- covTest.single(X, diag(1, 4))
  test2 <- covTest.single(X, diag(1, 4), ball = TRUE)
  test3 <- covTest.single(X, diag(2, 4), ball = TRUE)

  expect_identical(test1, output1)
  expect_identical(test2, output2)
  expect_identical(test3, output3)
})
