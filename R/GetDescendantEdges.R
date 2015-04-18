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
