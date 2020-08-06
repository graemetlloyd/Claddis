library(testthat)
library(Claddis)

test_that("find_linked_edges returns correct matrix", {
  expect_equal(find_linked_edges(tree = ape::read.tree(text = "(A,(B,(C,D)));")), matrix(data = c(0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0), ncol = 6, byrow = TRUE, dimnames = list(1:6, 1:6)))
})

# Undo any set.seed usage:
set.seed(Sys.time())
