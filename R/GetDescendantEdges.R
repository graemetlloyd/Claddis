#' Gets descendant edges of an internal node
#' 
#' Returns all descendant edges of an internal node for a phylo object.
#' 
#' Returns a vector of integers corresponding to row numbers in \code{$edge} or
#' cells in \code{$edge.length} of the descendant edges of the internal node
#' supplied.
#' 
#' @param n An integer corresponding to the internal node for which the
#' descendant edges are sought.
#' @param tree A tree as a phylo object.
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#' @examples
#' 
#' # Create simple four-taxon tree:
#' tree <- read.tree(text="(A,(B,(C,D)));")
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
    n.terminals <- length(FindDescendants(n, tree))
  
    # Cretae vector to store internal nodes:
    nodes <- n
  
    # Create vector to store edge numbers (i.e. row numbers for tree$edge):
    edges <- grep(n, tree$edge[, 1])
    
  # Keep going until all descendant edges are found:
  while(length(grep(TRUE, tree$edge[edges, 2] <= Ntip(tree))) < n.terminals) {
    
    # Get internal nodes found so far:
    nodes <- tree$edge[edges, 2][grep(TRUE, tree$edge[edges, 2] > Ntip(tree))]
        
    # For each node add any new descendant edges:
    for(i in nodes) edges <- sort(unique(c(edges, grep(TRUE, tree$edge[, 1] == i))))
    
  }
  
  # Return edges vector:
  return(edges)

}
