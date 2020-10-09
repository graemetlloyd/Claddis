library(testthat)
library(Claddis)

time_bins <- matrix(data = c(4, 3, 3, 2, 2, 1), ncol = 2, byrow = TRUE, dimnames = list(LETTERS[1:3], c("fad", "lad")))
class(time_bins) <- "timeBins"

test_that("bin_changes returns correct bin counts", {
  expect_equal(bin_changes(change_times = c(1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 3.5, 3.5, 3.5), time_bins = time_bins), c("A" = 3, "B" = 3, "C" = 3))
  expect_equal(bin_changes(change_times = c(1, 1.5, 2, 2.5, 3, 3.5, 4), time_bins = time_bins), c("A" = 2, "B" = 2, "C" = 2))
})

# Undo any set.seed usage:
set.seed(Sys.time())
