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
#' tree <- read.tree(text = "(A,(B,(C,D)));")
#' 
#' # Plot tree:
#' plot(tree)
#' 
#' # Add nodelabels:
#' nodelabels()
#' 
#' # Add edgelabels (note that edges 5 and 6
#' # are descendants of node 7):
#' edgelabels()
#' 
#' # Use GdtDescendantEdges to show that edges
#' # 5 and 6 are descendants of node 7:
#' GetDescendantEdges(7, tree)
#' 
#' @export GetDescendantEdges
GetDescendantEdges <- function(n, tree) {
  
  # Find number of terminals (i.e. stopping point):
  n.terminals <- length(strap::FindDescendants(n, tree))
  
  # Create vector to store internal nodes:
  nodes <- n
  
    # Create vector to store edge numbers (i.e. row numbers for tree$edge):
    edges <- grep(n, tree$edge[, 1])
    
  # Keep going until all descendant edges are found:
  while(length(which(tree$edge[edges, 2] <= ape::Ntip(tree))) < n.terminals) {
    
    # Get internal nodes found so far:
    nodes <- tree$edge[edges, 2][which(tree$edge[edges, 2] > ape::Ntip(tree))]
        
    # For each node add any new descendant edges:
    for(i in nodes) edges <- sort(unique(c(edges, which(tree$edge[, 1] == i))))
    
  }
  
  # Return edges vector:
  return(edges)

}
