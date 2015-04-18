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
