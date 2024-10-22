test_that("independent Test works", {
  data(iris)
  chart <- iris[, 1:4]

  expect_path <- normalizePath("./Expect/indTestmultiExpected.RData")
  load(expect_path)

  test1 <- indTest.multi(chart)
  test2 <- indTest.multi(chart, subdim = c(2, 1, 1))

  expect_identical(output1, test1)
  expect_identical(output2, test2)
})
