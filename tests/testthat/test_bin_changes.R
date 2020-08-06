library(testthat)
library(Claddis)

test_that("bin_changes returns correct bin counts", {
  expect_equal(bin_changes(change_times = c(1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 3.5, 3.5, 3.5), time_bins = c(4, 3, 2, 1)), c("4-3" = 3, "3-2" = 3, "2-1" = 3))
  expect_equal(bin_changes(change_times = c(1, 1.5, 2, 2.5, 3, 3.5, 4), time_bins = c(4, 3, 2, 1)), c("4-3" = 2, "3-2" = 2, "2-1" = 2))
})

# Undo any set.seed usage:
set.seed(Sys.time())
