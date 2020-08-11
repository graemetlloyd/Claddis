#' Gets descendant edges of an internal node
#'
#' @description
#'
#' Returns all descendant edges of an internal node for a phylo object.
#'
#' @param n An integer corresponding to the internal node for which the descendant edges are sought.
#' @param tree A tree as a phylo object.
#'
#' @details
#'
#' Returns a vector of integers corresponding to row numbers in \code{$edge} or cells in \code{$edge.length} of the descendant edges of the internal node supplied.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create simple four-taxon tree:
#' tree <- ape::read.tree(text = "(A,(B,(C,D)));")
#'
#' # Plot tree:
#' plot(tree)
#'
#' # Show nodelabels:
#' nodelabels()
#'
#' # Show edgelabels (note that edges 5 and 6
#' # are descendants of node 7):
#' edgelabels()
#'
#' # Use find_descendant_edges to show that edges
#' # 5 and 6 are descendants of node 7:
#' find_descendant_edges(n = 7, tree = tree)
#' @export find_descendant_edges
find_descendant_edges <- function(n, tree) {

  # Find number of tips:
  n_tips <- ape::Ntip(phy = tree)

  # Find number of terminals (i.e. stopping point):
  n_terminals <- length(x = strap::FindDescendants(n = n, tree = tree))

  # Create vector to store internal nodes:
  nodes <- n

  # Create vector to store edge numbers (i.e. row numbers for tree$edge):
  edges <- grep(n, tree$edge[, 1])

  # Keep going until all descendant edges are found:
  while (length(x = which(x = tree$edge[edges, 2] <= n_tips)) < n_terminals) {

    # Get internal nodes found so far:
    nodes <- tree$edge[edges, 2][which(x = tree$edge[edges, 2] > n_tips)]

    # For each node add any new descendant edges:
    for (i in nodes) edges <- sort(x = unique(x = c(edges, which(x = tree$edge[, 1] == i))))
  }

  # Return edges vector:
  edges
}
