test_that("Two Mean Test works", {
  data(iris)
  X <- iris[1:50, 1:4]
  Y <- iris[51:100, 1:4]

  expect_path <- normalizePath("./Expect/meanTesttwoExpected.RData")
  load(expect_path)

  test1 <- meanTest.two(X, Y)
  test2 <- meanTest.two(X, Y, equal = FALSE, method = "Coupled")
  test3 <- meanTest.two(X, Y, equal = FALSE, method = "Transformed")

  expect_identical(test1, output1)
  expect_identical(test2, output2)
  expect_identical(test3, output3)
})
