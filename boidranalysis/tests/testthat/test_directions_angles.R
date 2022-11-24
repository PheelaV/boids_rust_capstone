# to execute:
# > devtools::test()


test_that("Acquiring headings", {
  x <- c(0, 2, 3, 5, 6)
  y <- c(0, 2, 5, 6, 3)
  expected_results <- c(45.00000, 18.43495, 63.43495, 161.56505)

  acquired_results <- round(get_headings(x, y), 13)
  expect_equal(acquired_results, expected_results, tolerance = 1e-05)
})

test_that("Acquiring bearings", {
  x <- c(0, 2, 3, 5, 6)
  y <- c(0, 2, 5, 6, 3)
  expected_results <- c(26.56505, -45.00000, -98.13010)

  acquired_results <- round(get_bearings(get_headings(x, y)), 13)
  expect_equal(acquired_results, expected_results, tolerance = 1e-05)
})
