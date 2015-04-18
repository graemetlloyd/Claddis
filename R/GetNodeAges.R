GetNodeAges <- function(tree) {
  
  # Store root node number:
  rootnode <- Ntip(tree) + 1
  
  # Create initial paths list with end nodes (terminal and internal, excluding the root):
  paths <- split(c(1:Ntip(tree), (Ntip(tree) + 2):(Ntip(tree) + Nnode(tree))), f=1:(Ntip(tree) + Nnode(tree) - 1))
  
  # Strip names:
  names(paths) <- NULL

  # For each path:
  for (i in 1:length(paths)) {
    
    # Set counter as 1:
    j <- 1
    
    # Identify current node:
    currentnode <- paths[[i]][j]
    
    # While current node is not the root (path has not terminated):
    while (currentnode != rootnode) {
      
      # Update current node and add to path:
      currentnode <- paths[[i]][j + 1] <- tree$edge[match(currentnode, tree$edge[, 2]), 1]
      
      # Update counter:
      j <- j + 1
      
    }
    
  }
  
  # Create vector to store node ages:
  nodeages <- vector(mode="numeric", length=Ntip(tree) + Nnode(tree))
  
  # For each path:
  for (i in 1:length(paths)) {
    
    # Store path lengths from root:
    nodeages[paths[[i]][1]] <- sum(tree$edge.length[match(paths[[i]][1:(length(paths[[i]]) - 1)], tree$edge[, 2])])
  
  }
  
  # Subtract path lengths from root time:
  nodeages <- tree$root.time - nodeages
  
  # Add node numbers:
  names(nodeages) <- 1:(Ntip(tree) + Nnode(tree))
  
  # Return node ages:
  return(nodeages)
  
}
