AncStateEstMatrix <- function(morph.matrix, tree, estimate.allchars=FALSE, estimate.tips=FALSE) {

  # Catch problem with trees with no branch lengths:
  if(is.null(tree$edge.length)) stop("ERROR:\n Tree must have branch lengths.")

  # Catch problem with polytomies:
  if(tree$Nnode < (Ntip(tree) - 1)) stop("ERROR:\n Tree must be fully bifurcating.")

  # Catch problem with zero-length branches:
  if(any(tree$edge.length == 0)) stop("ERROR:\n Tree must not have zero-length branches.")

  # Remove node labels from tree (causes bug downstream):
  tree$node.label <- NULL

  # Collapse matrix to vectors for each character (state and ordering combination):
  collapse.matrix <- apply(rbind(morph.matrix$matrix, morph.matrix$ordering), 2, paste, collapse="")

  # Find just unique characters (no point repeating ancestral state reconstruction if codings and ordering are identical):
  unique.characters <- match(unique(collapse.matrix), collapse.matrix)

  # Case if only estimating states for ancestral nodes:
  if(estimate.tips == FALSE) {
    
    # Create ancestral storage matrix:
    anc.lik.matrix <- matrix(nrow=Nnode(tree), ncol=length(morph.matrix$matrix[1, ]))
    
    # Label matrix to record ancestral state estimates:
    rownames(anc.lik.matrix) <- c((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
    
  # Case if also estimating states for tip nodes:
  } else {

    # Create ancestral storage matrix (including tips):
    anc.lik.matrix <- matrix(nrow=Ntip(tree) + Nnode(tree), ncol=length(morph.matrix$matrix[1, ]))

    # Label matrix to record ancestral state estimates:
    rownames(anc.lik.matrix) <- c(tree$tip.label, (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)))
    
  }
  
  # For each unique character:
  for(i in unique.characters) {
        
    # Get minimum value for character:
    minval <- morph.matrix$min.vals[i]

    # Get maximum value for character:
    maxval <- morph.matrix$max.vals[i]

    # Case that not all characters are NA:
    if(!is.na(minval) && !is.na(maxval)) {
    
      # Only proceed if character is variable (non-constant):
      if(maxval != minval) {
            
        # If estimating states for all taxa then treat missing values as all possible states:
        if(estimate.allchars) morph.matrix$matrix[grep(TRUE, is.na(morph.matrix$matrix[, i])), i] <- paste(minval:maxval, collapse="&")
            
        # Find tips which cannot be used due to missing data:
        tipstogo <- rownames(morph.matrix$matrix)[grep(TRUE, is.na(morph.matrix$matrix[, i]))]
            
        # Only continue if at least three tips in pruned tree:
        if(length(tipstogo) < (Ntip(tree) - 2)) {
      
          # If there are tips to remove create a pruned tree:
          if(length(tipstogo) > 0) chartree <- drop.tip(tree, tipstogo)
        
          # If there are no tips to remove just use full tree:
          if(length(tipstogo) == 0) chartree <- tree
            
          # Get tip values for the pruned tree:
          tipvals <- morph.matrix$matrix[chartree$tip.label, i]
                
          # Continue if not gap-coded (or apparently gap-coded) or is unordered:
          if(!any(diff(sort(as.numeric(unlist(strsplit(tipvals, "&"))))) > 1) || morph.matrix$ordering[i] == "unord") {
                    
            # Set discrete character estimation model if unordered and or binary character:
            if(maxval - minval == 1 || maxval - minval > 1 && morph.matrix$ordering[i] == "unord") mymodel <- "ER"
                    
            # Set discrete character estimation model if ordered multistate character:
            if(maxval - minval > 1 && morph.matrix$ordering[i] == "ord") {
                        
              # Create all zero matrix:
              mymodel <- matrix(0, nrow=(maxval - minval) + 1, ncol=(maxval - minval) + 1)
                        
              # Name rows and columns as states:
              rownames(mymodel) <- colnames(mymodel) <- minval:maxval
                        
              # Enter one for all the off-diagonal diagonals (an ordered change model):
              for(j in 1:(length(mymodel[1, ]) - 1)) mymodel[j + 1, j] <- mymodel[j, j + 1] <- 1

            }
                    
            # Create matrix to store probabilities of tip values:
            tipvals.mat <- matrix(0, nrow=length(tipvals), ncol=maxval - minval + 1)
      
            # Add rownames (tip labels):
            rownames(tipvals.mat) <- names(tipvals)

            # Add colunames (state values):
            colnames(tipvals.mat) <- minval:maxval
                    
            # Fill all probabilities equal to one (non-polymorphisms):
            for(j in colnames(tipvals.mat)) tipvals.mat[grep(TRUE, tipvals == j), j] <- 1
                    
            # If there are polymorphisms make all observed states equally probable:
            if(any(apply(tipvals.mat, 1, sum) == 0)) {
                        
              # Get list of tip values with polymorphisms:
              polymorphism.values <- grep(TRUE, apply(tipvals.mat, 1, sum) == 0)
                        
              # Go through each polymorphism:
              for(j in polymorphism.values) {
                            
                # Get list of each state:
                states <- strsplit(tipvals[j], "&")[[1]]
                            
                # Make each state equally probable in tip values matrix:
                tipvals.mat[j, states] <- 1 / length(states)

              }

            }

            # Get likelihoods for each state in taxa and ancestors:
            state.likelihoods <- rerootingMethod(chartree, tipvals.mat, model=mymodel)$marginal.anc

            # Remove any potential node labels on the character tree to avoid an error from rerootingMethod():
            chartree$node.label<-NULL 

            # Get likelihoods for each state in taxa and ancestors:
            state.likelihoods <- rerootingMethod(chartree, tipvals.mat, model=mymodel)$marginal.anc
                    
            # Get maximum likelihood:
            max.lik <- apply(state.likelihoods, 1, max)
                    
            # Create vector to store maximum likelihood states:
            max.lik.state <- vector(mode="character", length=length(rownames(state.likelihoods)))
      
            # Add names:
            names(max.lik.state) <- rownames(state.likelihoods)
                    
            # For each tip and node find most likely state(s):
            for(j in 1:length(max.lik)) max.lik.state[j] <- paste(colnames(state.likelihoods)[grep(TRUE, state.likelihoods[j, ] == max.lik[j])], collapse="&")

          # Case if gap-coded or apparently gap-coded:
          } else {
                    
            # If there are polymorphisms present:
            if(length(grep("&", tipvals)) > 0) {
                        
              # List rows with polymorphisms:
              rows <- grep("&", tipvals)
                        
              # For each polymorphism:
              for(j in 1:length(rows)) {
                            
                # Replace polymorphisms with mean values:
                tipvals[rows[j]] <- mean(as.numeric(strsplit(tipvals[rows[j]], "&")[[1]]))
      
              }

            }
                    
            # Treat character as continuous to get ancestral state estimates:
            max.lik.state <- round(ace(as.numeric(tipvals), chartree)$ace)

          }
                
          # Get node numbers for pruned tree:
          chartreenodes <- (Ntip(chartree) + 1):(Ntip(chartree) + Nnode(chartree))
                
          # For each node on pruned tree:
          for(j in chartreenodes) {
                    
            # Find descendants in pruned tree:
            descs <- sort(chartree$tip.label[FindDescendants(j, chartree)])
                    
            # Store ancestral values in state estimation matrix:
            anc.lik.matrix[as.character(FindAncestor(descs, tree)), i] <- max.lik.state[as.character(j)]
 
          }

          # If also estimating tip states:
          if(estimate.tips == TRUE) {
          
            # Store tip states for taxa:
            anc.lik.matrix[chartree$tip.label, i] <- max.lik.state[match(chartree$tip.label, names(max.lik.state))]

          }

        # Case if pruned tree too small to inform (leaves matrix as NAs)
        }
            
      # Case if constant character:
      } else {
            
        # Enter constant value for all nodes:
        anc.lik.matrix[, i] <- minval

      }

    # Case if all characters are NA:
    } else {

      # Enter NA for all nodes:
      anc.lik.matrix[, i] <- NA

    }

  }

  # Fill in repeating characters to give full matrix:
  anc.lik.matrix <- anc.lik.matrix[, unique.characters[match(collapse.matrix, unique(collapse.matrix))]]
    
  # Return completed discrete ancestral character estimation matrix:
  return(anc.lik.matrix)

}
