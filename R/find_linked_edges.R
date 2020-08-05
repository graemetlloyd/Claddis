#' Find linked edges for a tree
#'
#' @description
#'
#' Given a tree finds edges that are linked to each other.
#'
#' @param tree A tree (phylo object).
#'
#' @details
#'
#' Finds all edges that link (share a node) with each edge of a tree.
#'
#' This is intended as an internal function, but may be of use to someone else.
#'
#' @return Returns a matrix where links are scored 1 and everything else 0. The diagonal is left as zero.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a simple four-taxon tree:
#' tree <- ape::read.tree(text = "(A,(B,(C,D)));")
#'
#' # Find linked (1) edges matrix for tree:
#' find_linked_edges(tree)
#' @export find_linked_edges
find_linked_edges <- function(tree) {

  # Create matrix of linked edges (all starting as unlinked, 0):
  edge_link_matrix <- matrix(0, nrow = nrow(tree$edge), ncol = nrow(tree$edge), dimnames = list(1:nrow(tree$edge), 1:nrow(tree$edge)))

  # For each edge:
  for (i in 1:nrow(tree$edge)) {

    # Find linked edges:
    links <- setdiff(x = union(which(x = apply(tree$edge == tree$edge[i, 1], 1, sum) == 1), which(x = apply(tree$edge == tree$edge[i, 2], 1, sum) == 1)), y = i)

    # Code 1 (linked) for linked edges in matrix:
    edge_link_matrix[i, links] <- edge_link_matrix[links, i] <- 1
  }

  # Return edge links matrix:
  edge_link_matrix
}
