TrimMorphDistMatrix <- function(dist.matrix, tree=NULL) {

  # If input is iof class "dist" first convert to a regular matrix:
  if(class(dist.matrix) == "dist") dist.matrix <- as.matrix(dist.matrix)

  # Check the input is a distance matrix:
  if(!is.matrix(dist.matrix)) stop("ERROR: Input must be a distance matrix (i.e., either an object of class \"dist\" or a square matrix).")

  # Case if there is no tree:
  if(is.null(tree)) {
    
    # Case if distance matrix is already complete:
    if(length(grep(TRUE, is.na(dist.matrix))) == 0) {

      # Warn user:
      print("There are no gaps in the distance matrix")
    
      # There are no taxa to be removed:
      removed.taxa <- NULL
      
      # Compile data in single variable:
      out <- list(dist.matrix, tree, removed.taxa)
      
      # Add names:
      names(out) <- c("dist.matrix", "tree", "removed.taxa")
      
      # Output:
      return(out)
      
    # Case if distance matrix has gaps:
    } else {
      
      # Vector to store taxa that need removing:
      removes <- vector(mode="character")
      
      # Whilst there are still gaps in the matrix:
      while(length(grep(TRUE, is.na(dist.matrix))) > 0) {
        
        # Vector to store number of NAs for each row of the matrix:
        na.lengths <- vector(mode="numeric")
        
        # For each row of the matrix get the number of NAs:
        for(i in 1:length(dist.matrix[, 1])) na.lengths[i] <- length(grep(TRUE, is.na(dist.matrix[i, ])))
        
        # Take row with most NAs and make it the taxon to delete:
        taxon.to.delete <- rownames(dist.matrix)[grep(TRUE, na.lengths == max(na.lengths))[1]]
        
        # Find the taxons row:
        delete.row <- grep(TRUE, rownames(dist.matrix) == taxon.to.delete)
        
        # Remove it from the distance matrix:
        dist.matrix <- dist.matrix[-delete.row, -delete.row]
        
        # Add taxon to removes list:
        removes <- c(removes, taxon.to.delete)

      }

      # Compile data in single variable:
      out <- list(dist.matrix, tree, removes)
      
      # Add names:
      names(out) <- c("dist.matrix", "tree", "removed.taxa")
      
      # Output:
      return(out)

    }
    
  # Case if there is a tree:
  } else {
    
    # Case if distance matrix is already complete:
    if(length(grep(TRUE, is.na(dist.matrix))) == 0) {
      
      # Warn user:
      print("There are no gaps in the distance matrix")
      
      # Create NULL vector for removed taxa:
      removed.taxa <- NULL
      
      # Compile data in single variable:
      out <- list(dist.matrix, tree, removed.taxa)
      
      # Add names:
      names(out) <- c("dist.matrix", "tree", "removed.taxa")
      
      # Return unmodified matrix and tree:
      return(out)
      
    # Case if distance matrix has gaps:
    } else {
      
      # Get list of node numbers as text:
      node.nos <- as.character((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
      
      # Rename dist.matrix rownames by descendant taxa that define them:
      for(i in match(node.nos, rownames(dist.matrix))) colnames(dist.matrix)[i] <- rownames(dist.matrix)[i] <- paste(sort(tree$tip.label[FindDescendants(rownames(dist.matrix)[i], tree)]), collapse="%%")
      
      # Vector to store taxa and nodes that need removing:
      removes <- vector(mode="character")
      
      # Temporarily store complete distance matrix:
      temp.dist.matrix <- dist.matrix
      
      # Whilst there are still gaps in the matrix:
      while(length(grep(TRUE, is.na(dist.matrix))) > 0) {
        
        # Vector to store number of NAs for each row of the matrix:
        na.lengths <- vector(mode="numeric")
        
        # For each row of the matrix get the number of NAs:
        for(i in 1:length(dist.matrix[, 1])) na.lengths[i] <- length(grep(TRUE, is.na(dist.matrix[i, ])))
        
        # Take row with most NAs and make it the taxon to delete:
        taxon.to.delete <- rownames(dist.matrix)[grep(TRUE, na.lengths == max(na.lengths))[1]]
        
        # Find the taxons row:
        delete.row <- grep(TRUE, rownames(dist.matrix) == taxon.to.delete)
        
        # Remove it from the distance matrix:
        dist.matrix <- dist.matrix[-delete.row, -delete.row]
        
        # Add taxon to removes list:
        removes <- c(removes, taxon.to.delete)

      }
      
      # Restore complete distance matrix:
      dist.matrix <- temp.dist.matrix
      
      # Find tips to remove:
      tips.to.remove <- removes[sort(match(tree$tip.label, removes))]
      
      # Find nodes to remove:
      nodes.to.remove <- removes[grep("%%", removes)]
      
      # Vector to store nodes to delete (i.e. from where tips to delete originate):
      tip.name.nodes.to.remove <- vector(mode="numeric")
      
      # For each tip to delete:
      for(i in tips.to.remove) {
        
        # Get originating node for tip:
        originating.node <- tree$edge[match(grep(TRUE, tree$tip.label == i), tree$edge[, 2]), 1]
        
        # Get name of originating node:
        originating.node.name <- paste(sort(tree$tip.label[FindDescendants(originating.node, tree)]), collapse="%%")
        
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
          originating.node <- FindAncestor(strsplit(i, "%%")[[1]], tree)
          
          # Find descendants of that node:
          descendant.nodes <- tree$edge[grep(TRUE, tree$edge[, 1] == originating.node), 2]
          
          # Case if one of the descendants is a tip:
          if(length(grep(TRUE,descendant.nodes <= Ntip(tree))) > 0) {
            
            # Get taxon to exclude:
            taxon.to.exclude <- tree$tip.label[min(descendant.nodes)]
            
          # Case if all descendants are internal nodes:
          } else {
            
            # If first descendant node has more taxa:
            if(length(FindDescendants(descendant.nodes[1], tree)) >= length(FindDescendants(descendant.nodes[2], tree))) {
              
              # Get taxon to exclude:
              taxon.to.exclude <- tree$tip.label[FindDescendants(descendant.nodes[2], tree)]
              
            # If second descendant node has more taxa:
            } else {
              
              # Get taxon to exclude:
              taxon.to.exclude <- tree$tip.label[FindDescendants(descendant.nodes[1], tree)]

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
        rownames(dist.matrix)[i] <- colnames(dist.matrix)[i] <- paste(sort(nms[grep(TRUE, is.na(match(nms, tips.to.remove)))]), collapse="%%")
        
      }

      # Estbalish if there are redundant rows (nodes defined by a single tip or by no tips at all):
      redundant.rows <- c(grep(TRUE, duplicated(rownames(dist.matrix))), grep(TRUE, rownames(dist.matrix) == ""))
      
      # If there are any redundant rows:
      if(length(redundant.rows) > 0) {
        
        # Remove them from the distance matrix:
        dist.matrix <- dist.matrix[-redundant.rows, -redundant.rows]

      }
      
      # Remove pruned taxa from tree:
      tree <- drop.tip(tree, tips.to.remove)
      
      # Find node names:
      node.names <- rownames(dist.matrix)[grep("%%", rownames(dist.matrix))]
      
      # Find node numbers:
      for(j in 1:length(node.names)) names(node.names)[j] <- FindAncestor(strsplit(node.names[j], "%%")[[1]], tree)
      
      # Replace node names with numbers:
      colnames(dist.matrix)[match(node.names, colnames(dist.matrix))] <- rownames(dist.matrix)[match(node.names, rownames(dist.matrix))] <- names(node.names)
      
      # Compile data in single variable:
      out <- list(dist.matrix, tree, removes)
      
      # Add names:
      names(out) <- c("dist.matrix", "tree", "removed.taxa")
      
      # Return answer:
      return(out)

    }

  }

}
