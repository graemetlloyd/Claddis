#' Find ancestor
#'
#' @description
#'
#' Finds the last common ancestor (node) of a set of two or more descendant tips.
#'
#' @param descs A vector of mode character representing the tip names for which an ancestor is sought.
#' @param tree The tree as a phylo object.
#'
#' @details
#'
#' Intended for use as an internal function for \link{TrimMorphDistMatrix}, but potentially of more general use.
#'
#' @return \item{anc.node}{The ancestral node number.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a simple four-taxon tree:
#' tree <- read.tree(text = "(A,(B,(C,D)));")
#' 
#' # Plot the tree:
#' plot(tree)
#' 
#' # Add nodelabels and show that the most recent common
#' # ancestor of B, C, and D is node 6:
#' nodelabels()
#' 
#' # Use FindAncestor to show that the most recent common
#' # ancestor of B, C, and D is node 6:
#' FindAncestor(c("B", "C", "D"), tree)
#' 
#' @export FindAncestor
FindAncestor <- function(descs, tree) {

  # Get tip numbers:
  tipnos <- match(descs, tree$tip.label)
  
  # Get ancestral nodes in order:
  anc.node <- sort(unique(tree$edge[, 1][match(tipnos, tree$edge[, 2])]))
  
  # Keep going until a single ancestral node is converged upon:
  while(length(anc.node) > 1) {
    
    # Get node with highest number (definitely not ancestor):
    highestnode <- anc.node[length(anc.node)]
    
    # Remove this node from the list:
    anc.node <- anc.node[-length(anc.node)]
    
    # Find its ancestor and add to unique list:
    anc.node <- sort(unique(c(anc.node, tree$edge[match(highestnode, tree$edge[, 2]), 1])))

  }
  
  # Return ancestral node:
  return(anc.node)

}
