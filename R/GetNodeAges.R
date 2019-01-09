#' Returns node ages for a time-scaled tree
#'
#' @description
#'
#' Given a tree with branch-lengths scaled to time and a value for \code{$root.time} will return a vector of node ages.
#'
#' @param tree A tree (phylo object) with branch lengths representing time and a value for \code{$root.time}.
#'
#' @details
#'
#' Returns a vector of node ages (terminal and internal) labelled by their node number.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#' 
#' # Create simple four-taxon tree with edge lengths all
#' # set to 1 Ma:
#' tree <- read.tree(text="(A:1,(B:1,(C:1,D:1):1):1);")
#' 
#' # Set root.time as 10 Ma:
#' tree$root.time <- 10
#' 
#' # Get node ages:
#' GetNodeAges(tree)
#' 
#' @export GetNodeAges
GetNodeAges <- function(tree) {
  
  # Store root node number:
  rootnode <- ape::Ntip(tree) + 1
  
  # If tree is a complete polytomy:
  if(tree$Nnode == 1) {
      
    # Create paths for just tips:
    paths <- as.list(1:ape::Ntip(tree))
    
    # Add root to each path:
    for(i in 1:length(paths)) paths[[i]] <- c(paths[[i]], ape::Ntip(tree) + 1)
    
  # If tree is not a complete polytomy:
  } else {
      
    # Create initial paths list with end nodes (terminal and internal, excluding the root):
    paths <- split(c(1:ape::Ntip(tree), (ape::Ntip(tree) + 2):(ape::Ntip(tree) + ape::Nnode(tree))), f = 1:(ape::Ntip(tree) + ape::Nnode(tree) - 1))
      
    # Strip names:
    names(paths) <- NULL
      
    # For each path:
    for(i in 1:length(paths)) {
          
      # Set counter as 1:
      j <- 1
          
      # Identify current node:
      currentnode <- paths[[i]][j]
          
      # While current node is not the root (path has not terminated):
      while(currentnode != rootnode) {
              
        # Update current node and add to path:
        currentnode <- paths[[i]][j + 1] <- tree$edge[match(currentnode, tree$edge[, 2]), 1]
              
        # Update counter:
        j <- j + 1
              
      }
          
    }
      
  }

  # Create vector to store node ages:
  nodeages <- vector(mode = "numeric", length = ape::Ntip(tree) + ape::Nnode(tree))
  
  # For each path:
  for(i in 1:length(paths)) {
    
    # Store path lengths from root:
    nodeages[paths[[i]][1]] <- sum(tree$edge.length[match(paths[[i]][1:(length(paths[[i]]) - 1)], tree$edge[, 2])])
  
  }
  
  # Subtract path lengths from root time:
  nodeages <- tree$root.time - nodeages
  
  # Add node numbers:
  names(nodeages) <- 1:(ape::Ntip(tree) + ape::Nnode(tree))
  
  # Return node ages:
  return(nodeages)
  
}
