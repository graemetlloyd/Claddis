#' Trims a morphological distance matrix
#'
#' @description
#'
#' Trims a morphological distance matrix by removing objects that cause empty cells.
#'
#' @param dist.matrix A distance matrix in the format created by \link{MorphDistMatrix}.
#' @param Tree If the distance matrix includes ancestors this should be the tree (phylo object) used to reconstruct them.
#'
#' @details
#'
#' Trims a morphological distance matrix by removing nodes (terminal or internal) that cause empty cells allowing it to be passed to an ordination function such as \link{cmdscale}.
#'
#' Some distances are not calculable from cladistic matrices if there are taxa that have no coded characters in common. This algorithm iteratively removes the taxa responsible for the most empty cells until the matrix is complete (no empty cells).
#'
#' If the matrix includes reconstructed ancestors the user should also provide the tree used (as the \code{tree} argument). The function will then also remove the tips from the tree and where reconstructed ancestors also cause empty cells will prune the minimum number of descendants of that node. The function will then renumber the nodes in the distance matrix so they match the pruned tree.
#'
#' @return
#'
#' \item{DistMatrix}{A complete distance matrix with all cells filled. If there were no empty cells will return original.}
#' \item{Tree}{A tree (if supplied) with the removed taxa (see below) pruned. If no taxa are dropped will return the same tree as inputted. If no tree is supplied this is set to NULL.}
#' \item{RemovedTaxa}{A character vector listing the taxa removed. If none are removed this will be set to NULL.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get morphological distances for Michaux (1989) data set:
#' distances <- MorphDistMatrix(Michaux1989)
#'
#' # Attempt to trim max.dist.matrix:
#' TrimMorphDistMatrix(distances$DistanceMatrix)
#'
#' @export TrimMorphDistMatrix
TrimMorphDistMatrix <- function(dist.matrix, Tree = NULL) {
  
  # If input is of class "dist" first convert to a regular matrix:
  if(class(dist.matrix) == "dist") dist.matrix <- as.matrix(dist.matrix)
  
  # Check the input is a distance matrix:
  if(!is.matrix(dist.matrix)) stop("ERROR: Input must be a distance matrix (i.e., either an object of class \"dist\" or a square matrix).")
  
  # Case if there is no tree:
  if(is.null(Tree)) {
    
    # Case if distance matrix is already complete:
    if(length(which(is.na(dist.matrix))) == 0) {
      
      # Warn user:
      print("There are no gaps in the distance matrix")
      
      # There are no taxa to be removed:
      removed.taxa <- NULL
      
      # Compile data in single variable:
      out <- list(dist.matrix, Tree, removed.taxa)
      
      # Add names:
      names(out) <- c("DistMatrix", "Tree", "RemovedTaxa")
      
      # Output:
      return(out)
      
    # Case if distance matrix has gaps:
    } else {
      
      # Vector to store taxa that need removing:
      removes <- vector(mode = "character")
      
      # Whilst there are still gaps in the matrix:
      while(length(which(is.na(dist.matrix))) > 0) {
        
        # Vector to store number of NAs for each row of the matrix:
        na.lengths <- vector(mode="numeric")
        
        # For each row of the matrix get the number of NAs:
        for(i in 1:length(dist.matrix[, 1])) na.lengths[i] <- length(which(is.na(dist.matrix[i, ])))
        
        # Take row with most NAs and make it the taxon to delete:
        taxon.to.delete <- rownames(dist.matrix)[which(na.lengths == max(na.lengths))[1]]
        
        # Find the taxons row:
        delete.row <- which(rownames(dist.matrix) == taxon.to.delete)
        
        # Remove it from the distance matrix:
        dist.matrix <- dist.matrix[-delete.row, -delete.row]
        
        # Add taxon to removes list:
        removes <- c(removes, taxon.to.delete)
        
      }
      
      # Compile data in single variable:
      out <- list(dist.matrix, Tree, removes)
      
      # Add names:
      names(out) <- c("DistMatrix", "Tree", "RemovedTaxa")
      
      # Output:
      return(out)
      
    }
    
  # Case if there is a tree:
  } else {
    
    # Case if distance matrix is already complete:
    if(length(which(is.na(dist.matrix))) == 0) {
      
      # Warn user:
      print("There are no gaps in the distance matrix")
      
      # Create NULL vector for removed taxa:
      removed.taxa <- NULL
      
      # Compile data in single variable:
      out <- list(dist.matrix, Tree, removed.taxa)
      
      # Add names:
      names(out) <- c("DistMatrix", "Tree", "RemovedTaxa")
      
      # Return unmodified matrix and tree:
      return(out)
      
    # Case if distance matrix has gaps:
    } else {
      
      # Get list of node numbers as text:
      node.nos <- as.character((ape::Ntip(Tree) + 1):(ape::Ntip(Tree) + ape::Nnode(Tree)))
      
      # Rename dist.matrix rownames by descendant taxa that define them:
      for(i in match(node.nos, rownames(dist.matrix))) colnames(dist.matrix)[i] <- rownames(dist.matrix)[i] <- paste(sort(Tree$tip.label[FindDescendants(rownames(dist.matrix)[i], Tree)]), collapse="%%")
      
      # Vector to store taxa and nodes that need removing:
      removes <- vector(mode="character")
      
      # Temporarily store complete distance matrix:
      temp.dist.matrix <- dist.matrix
      
      # Whilst there are still gaps in the matrix:
      while(length(which(is.na(dist.matrix))) > 0) {
        
        # Vector to store number of NAs for each row of the matrix:
        na.lengths <- vector(mode="numeric")
        
        # For each row of the matrix get the number of NAs:
        for(i in 1:length(dist.matrix[, 1])) na.lengths[i] <- length(which(is.na(dist.matrix[i, ])))
        
        # Take row with most NAs and make it the taxon to delete:
        taxon.to.delete <- rownames(dist.matrix)[which(na.lengths == max(na.lengths))[1]]
        
        # Find the taxons row:
        delete.row <- which(rownames(dist.matrix) == taxon.to.delete)
        
        # Remove it from the distance matrix:
        dist.matrix <- dist.matrix[-delete.row, -delete.row]
        
        # Add taxon to removes list:
        removes <- c(removes, taxon.to.delete)
        
      }
      
      # Restore complete distance matrix:
      dist.matrix <- temp.dist.matrix
      
      # Find tips to remove:
      tips.to.remove <- removes[sort(match(Tree$tip.label, removes))]
      
      # Find nodes to remove:
      nodes.to.remove <- removes[grep("%%", removes)]
      
      # Vector to store nodes to delete (i.e. from where tips to delete originate):
      tip.name.nodes.to.remove <- vector(mode="numeric")
      
      # For each tip to delete:
      for(i in tips.to.remove) {
        
        # Get originating node for tip:
        originating.node <- Tree$edge[match(which(Tree$tip.label == i), Tree$edge[, 2]), 1]
        
        # Get name of originating node:
        originating.node.name <- paste(sort(Tree$tip.label[FindDescendants(originating.node, Tree)]), collapse = "%%")
        
        # Add to nodes to delete vector:
        tip.name.nodes.to.remove <- unique(c(tip.name.nodes.to.remove, originating.node.name))
        
      }
      
      # Check if nodes need to be deleted that are not related to tips to be deleted:
      if(length(setdiff(nodes.to.remove, tip.name.nodes.to.remove)) > 0) {
        
        # Identify nodes left to remove:
        nodes.left.to.remove <- setdiff(nodes.to.remove, tip.name.nodes.to.remove)
        
        # Go through each node:
        for(i in nodes.left.to.remove) {
          
          # Find node number:
          originating.node <- FindAncestor(strsplit(i, "%%")[[1]], Tree)
          
          # Find descendants of that node:
          descendant.nodes <- Tree$edge[which(Tree$edge[, 1] == originating.node), 2]
          
          # Case if one of the descendants is a tip:
          if(length(which(descendant.nodes <= ape::Ntip(Tree))) > 0) {
            
            # Get taxon to exclude:
            taxon.to.exclude <- Tree$tip.label[min(descendant.nodes)]
            
          # Case if all descendants are internal nodes:
          } else {
            
            # If first descendant node has more taxa:
            if(length(FindDescendants(descendant.nodes[1], Tree)) >= length(FindDescendants(descendant.nodes[2], Tree))) {
              
              # Get taxon to exclude:
              taxon.to.exclude <- Tree$tip.label[FindDescendants(descendant.nodes[2], Tree)]
              
            # If second descendant node has more taxa:
            } else {
              
              # Get taxon to exclude:
              taxon.to.exclude <- Tree$tip.label[FindDescendants(descendant.nodes[1], Tree)]
              
            }
            
          }
          
          # Add to removes list:
          removes <- c(removes, taxon.to.exclude)
          
          # Add to tips to remove list:
          tips.to.remove <- c(tips.to.remove, taxon.to.exclude)
          
        }
        
      }
      
      # Prune matrix until complete:
      dist.matrix <- dist.matrix[-match(removes, rownames(dist.matrix)), -match(removes, rownames(dist.matrix))]
      
      # Loop to prune removed taxa from node names:
      for(i in 1:length(dist.matrix[, 1])) {
        
        # Split name into constituent taxa:
        nms <- strsplit(rownames(dist.matrix)[i], "%%")[[1]]
        
        # Get only still present taxa and rename rows and columns:
        rownames(dist.matrix)[i] <- colnames(dist.matrix)[i] <- paste(sort(nms[which(is.na(match(nms, tips.to.remove)))]), collapse="%%")
        
      }
      
      # Establish if there are redundant rows (nodes defined by a single tip or by no tips at all):
      redundant.rows <- c(which(duplicated(rownames(dist.matrix))), which(rownames(dist.matrix) == ""))
      
      # If there are any redundant rows:
      if(length(redundant.rows) > 0) {
        
        # Remove them from the distance matrix:
        dist.matrix <- dist.matrix[-redundant.rows, -redundant.rows]
        
      }
      
      # Remove pruned taxa from tree:
      Tree <- drop.tip(Tree, tips.to.remove)
      
      # Find node names:
      node.names <- rownames(dist.matrix)[grep("%%", rownames(dist.matrix))]
      
      # Find node numbers:
      for(j in 1:length(node.names)) names(node.names)[j] <- FindAncestor(strsplit(node.names[j], "%%")[[1]], Tree)
      
      # Replace node names with numbers:
      colnames(dist.matrix)[match(node.names, colnames(dist.matrix))] <- rownames(dist.matrix)[match(node.names, rownames(dist.matrix))] <- names(node.names)
      
      # Compile data in single variable:
      out <- list(dist.matrix, Tree, removes)
      
      # Add names:
      names(out) <- c("DistMatrix", "Tree", "RemovedTaxa")
      
      # Return answer:
      return(out)
      
    }
    
  }
  
}
