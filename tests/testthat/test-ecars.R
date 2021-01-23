test_that("ecars works", {
  X <- matrix(1, 4, 4)
  expect_equal(ecars(X), 1)
})
