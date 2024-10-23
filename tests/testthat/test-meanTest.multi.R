test_that("Multi Mean Test works", {
  data(iris)
  chart <- iris[, 1:4]
  species <- iris[, 5]

  expect_path <- normalizePath("./Expect/meanTestmultiExpected.RData")
  load(expect_path)

  test <- meanTest.multi(chart, species)

  expect_equal(test, output)
})
