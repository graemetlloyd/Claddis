library(testthat)
library(Claddis)

test_that("find_mrca returns correct node numbers", {
  expect_equal(find_mrca(descendant_names = c("B", "C", "D"), tree = ape::read.tree(text = "(A,(B,(C,D)));")), 6)
  expect_equal(find_mrca(descendant_names = c("A", "D"), tree = ape::read.tree(text = "(A,(B,(C,D)));")), 5)
  expect_equal(find_mrca(descendant_names = c("C", "D"), tree = ape::read.tree(text = "(A,(B,(C,D)));")), 7)
})

# Undo any set.seed usage:
set.seed(Sys.time())
