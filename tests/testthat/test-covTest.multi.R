test_that("Multi Covariance Test works", {
  data(iris)
  chart <- iris[, 1:4]
  species <- iris[, 5]
  
  expect_path <- normalizePath("./Expect/covTestmultiExpected.RData")
  load(expect_path)

  test <- covTest.multi(chart, species)

  expect_identical(test, output)
})
