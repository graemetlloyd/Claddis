#' Find ancestor
#'
#' @description
#'
#' Finds the last common ancestor (node) of a set of two or more descendant tips.
#'
#' @param descendant_names A vector of mode character representing the tip names for which an ancestor is sought.
#' @param tree The tree as a phylo object.
#'
#' @details
#'
#' Intended for use as an internal function for \link{trim_matrix}, but potentially of more general use.
#'
#' @return \item{ancestor_node}{The ancestral node number.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a simple four-taxon tree:
#' tree <- ape::read.tree(text = "(A,(B,(C,D)));")
#'
#' # Plot the tree:
#' ape::plot.phylo(tree)
#'
#' # Add nodelabels and show that the most recent common
#' # ancestor of B, C, and D is node 6:
#' ape::nodelabels()
#'
#' # Use find_mrca to show that the most recent common
#' # ancestor of B, C, and D is node 6:
#' find_mrca(
#'   descendant_names = c("B", "C", "D"),
#'   tree = tree
#' )
#' @export find_mrca
find_mrca <- function(descendant_names, tree) {

  # Get tip numbers:
  tip_numbers <- match(descendant_names, tree$tip.label)

  # Get ancestral nodes in order:
  ancestor_node <- sort(x = unique(x = tree$edge[, 1][match(tip_numbers, tree$edge[, 2])]))

  # Keep going until a single ancestral node is converged upon:
  while (length(x = ancestor_node) > 1) {

    # Get node with highest number (definitely not ancestor):
    highest_node <- ancestor_node[length(x = ancestor_node)]

    # Remove this node from the list:
    ancestor_node <- ancestor_node[-length(x = ancestor_node)]

    # Find its ancestor and add to unique list:
    ancestor_node <- sort(x = unique(x = c(ancestor_node, tree$edge[match(highest_node, tree$edge[, 2]), 1])))
  }

  # Return ancestral node:
  ancestor_node
}
