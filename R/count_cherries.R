#' Counts the number of cherries in a tree
#'
#' @description
#'
#' Given a set of phylogenetic tree(s) returns the number of cherries in each one.
#'
#' @param tree A tree (phylo or multiPhylo object).
#'
#' @details
#'
#' Cherries are components of a phylogenetic tree defined as internal nodes with exactly two terminal descendants.
#'
#' This function simply counts the number present in a given tree.
#'
#' Note that any fully dichotomous phylogenetic tree must have at least one cherry.
#'
#' @return Returns a vector of cherry counts for each tree retaining the order in which they were supplied.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create simple two-cherry tree:
#' tree <- ape::read.tree(text = "((A,B),(C,D));")
#'
#' # Show count of cherries is two:
#' count_cherries(tree = tree)
#'
#' # Create a star tree:
#' tree <- ape::read.tree(text = "(A,B,C,D);")
#'
#' # Show count of cherries is zero:
#' count_cherries(tree = tree)
#'
#' @export count_cherries
count_cherries <- function(tree) {
  
  # If a single tree return scalar of pre-terminal nodes with two descendants (i.e., cherries):
  if (inherits(x = tree, what = "phylo")) return(value = sum(x = rle(x = sort(x = tree$edge[match(x = 1:ape::Ntip(tree), table = tree$edge[, 2]), 1]))$lengths == 2))
  
  # If multiple trees return vector of pre-terminal nodes with two descendants (i.e., cherries):
  if (inherits(x = tree, what = "multiPhylo")) return(unlist(lapply(X = tree, FUN = function(x) value = sum(x = rle(x = sort(x = x$edge[match(x = 1:ape::Ntip(x), table = x$edge[, 2]), 1]))$lengths == 2))))
}
