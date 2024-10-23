test_that("Mean Covariance Test works", {
  data(iris)
  chart <- iris[, 1:4]
  species <- iris[, 5]

  expect_path <- normalizePath("./Expect/meancovTestExpected.RData")
  load(expect_path)

  test <- meancov.Test(chart, species)

  expect_equal(test, output)
})
