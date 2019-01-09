#' Stochastic Character Map For Dollo Character
#'
#' @description
#'
#' Given a tree with binary tip states produces a stochastic Dollo character map.
#'
#' @param tree A tree in phylo format with barnch lengths and a value for \code{$root.time}.
#' @param tip.states A vector of tip states (must be 0 or 1) with names matching \code{tree$tip.label}.
#'
#' @details
#'
#' A non-ideal solution to the problem of generating a stochastic character map for a Dollo character (i.e., a single gain of the derived state (1) with any number of losses).
#'
#' The function operates as follows:
#'
#' 1) Establishes least inclusive clade exhibiting the derived state (1).
#' 2) Assumes single gain occurred on branch subtending this clade and with equal probability of occurring at any point along the branch.
#' 3) Treats inclusive clade as a subtree and places a strong prior on the root of the derived state (1).
#' 4) Calls \code{make.simmap} from the \code{phytools} package to generate a stochastic character map using a model where only losses are possible.
#' 5) Outputs both the stochastic character map (time spent in each state on each branch) and a matrix of state changes.
#'
#' (NB: As the map is stochastic the answer will be different each time the function is run and multiple replicates are strongly advised in order to ascertain uncertainty.)
#'
#' This was the method used in Tarver et al. (2018).
#'
#' @return
#'
#' \item{Changes}{A matrix of all changes (gains and losses).}
#' \item{SCM}{The stochastic character map.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Tarver, J. E., Taylor, R. S., Puttick, M. N., Lloyd, G. T., Pett, W., Fromm, B., Schirrmeister, B. E., Pisani, D., Peterson, K. J. and Donoghue, P. C. J., 2018. Well-annotated microRNAomes do not evidence pervasive miRNA loss. Genome Biology and Evolution, 6, 1457-1470.
#'
#' @examples
#' 
#' # Create a random 10-taxon tree:
#' tree <- rtree(10)
#' 
#' # Arbitrarily add a root.time value of 100 Ma:
#' tree$root.time <- 100
#'
#' # Generate random tip states (0s and 1s):
#' tip.states <- sample(c(0, 1), 10, replace = TRUE)
#'
#' # Add labels to tip states:
#' names(tip.states) <- tree$tip.label
#'
#' # Get a single stochastic character map:
#' out <- DolloSCM(tree, tip.states)
#'
#' # View matrix of changes:
#' out$Changes
#'
#' # View stochastic character map (time spent in each state on each branch):
#' out$SCM
#'
#' @export DolloSCM
DolloSCM <- function(tree, tip.states) {
  
  # Check tree has branch lengths:
  if(is.null(tree$edge.length)) stop("Tree must have branch lengths.")
  
  # Check tree has a root time:
  if(is.null(tree$root.time)) stop("Tree must have a $root.time value.")
  
  # Record non-matching names:
  non.matches <- c(setdiff(names(tip.states), tree$tip.label), setdiff(tree$tip.label, names(tip.states)))
  
  # If there are any non-matching names stop and notify user:
  if(length(non.matches) > 0) stop(paste("Non-matching names between tree and tips states found:", paste(non.matches, collapse = ", ")))
  
  # Check for states other than zero or one:
  non.binary.states <- setdiff(unique(tip.states), c(0, 1))
  
  # If states other than zero or one found warn user:
  if(length(non.binary.states) > 0) stop("States other than 0 or 1 found. All states must be 0 or 1.")
  
  # Make Dollo-like model where only losses (1 -> 0) are allowed (will later force root state to be one):
  Dollo.model <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2, dimnames = list(c("0", "1"), c("0", "1")))
  
  # Find new root (indicating least inclusive clade of
  new.root <- FindAncestor(names(which(tip.states == 1)), tree)
  
  # If new root is really the old root:
  if((ape::Ntip(tree) + 1) == new.root && sum(tip.states == 1) > 1) {
    
    # Set acqusition branch to zero as precedes root:
    acquisition.branch <- 0
    
    # Set acquisition time to root time:
    acquisition.time <- tree$root.time
    
    # If tip states vary:
    if(length(unique(tip.states)) > 1) {
      
      # Get stochastic character map for full tree using Dollo model and a strong root prior of one:
      SCM <- make.simmap(tree, tip.states[tree$tip.label], model = Dollo.model, pi = c(0, 1), message = FALSE)$maps
      
    # If tip states are constant:
    } else {
      
      # Create SCM with no losses:
      SCM <- as.list(tree$edge.length)
      
      # Label all states as derived:
      for(i in 1:length(SCM)) names(SCM[[i]]) <- as.character(unique(tip.states))
      
    }
    
  # If new root reflects a subtree:
  } else {
    
    # Case if character is an autapomorphy:
    if(sum(tip.states == 1) == 1) {
      
      # Create base stochastic character map (SCM:
      SCM <- as.list(c(1:nrow(tree$edge)))
      
      # For each branch in the SCM:
      for(i in 1:length(SCM)) {
        
        # Create a null value (vector equal to branch length):
        x <- tree$edge.length[i]
        
        # Label null vector with null (root) state of zero:
        names(x) <- "0"
        
        # Update base SCM with null value:
        SCM[[i]] <- x
        
      }
      
      # Get branch along which (single) acqusition of derived state occurs:
      acquisition.branch <- match(which(tip.states == 1), tree$edge[, 2])
      
      # Get bounds for acquisition times:
      acquisition.bounds <- GetNodeAges(tree)[tree$edge[acquisition.branch, ]]
      
      # Draw an acusition time from a uniform distribution between the bounds:
      acquisition.time <- runif(1, min = acquisition.bounds[2], max = acquisition.bounds[1])
      
      # Create an SCM branch for the acquisition:
      acquisition.branch.SCM <- c(acquisition.bounds[1] - acquisition.time, acquisition.time - acquisition.bounds[2])
      
      # Add labels to the acquisition branch:
      names(acquisition.branch.SCM) <- c("0", "1")
      
      # Add acquisition branch to the SCM:
      SCM[[acquisition.branch]] <- acquisition.branch.SCM
      
    # Case if at least two taxa exhibit the derived state:
    } else {
      
      # Find members of the least inclusive clade:
      clade.members <- tree$tip.label[FindDescendants(new.root, tree)]
      
      # Find non-members of the least inclusive clade:
      nonclade.members <- setdiff(tree$tip.label, clade.members)
      
      # Case if only two members in clade:
      if(length(clade.members) == 2) {
        
        # Get branch along which (single) acqusition of derived state occurs:
        acquisition.branch <- match(new.root, tree$edge[, 2])
        
        # Get bounds for acquisition times:
        acquisition.bounds <- GetNodeAges(tree)[tree$edge[acquisition.branch, ]]
        
        # Draw an acusition time from a uniform distribution between the bounds:
        acquisition.time <- runif(1, min = acquisition.bounds[2], max = acquisition.bounds[1])
        
        # Create base stochastic character map (SCM):
        SCM <- as.list(c(1:nrow(tree$edge)))
        
        # For each branch in the SCM:
        for(i in 1:length(SCM)) {
          
          # Create a null value (vector equal to branch length):
          x <- tree$edge.length[i]
          
          # Label null vector with null (root) state of zero:
          names(x) <- "0"
          
          # Update base SCM with null value:
          SCM[[i]] <- x
          
        }
        
        # Create an SCM branch for the acquisition:
        acquisition.branch.SCM <- c(acquisition.bounds[1] - acquisition.time, acquisition.time - acquisition.bounds[2])
        
        # Add labels to the acqusition branch:
        names(acquisition.branch.SCM) <- c("0", "1")
        
        # Add acquisition branch to the SCM:
        SCM[[acquisition.branch]] <- acquisition.branch.SCM
        
        # Get two descending edges:
        descendant.edges <- GetDescendantEdges(new.root, tree)
        
        # Update state of descending edges to derived (1):
        for(i in descendant.edges) names(SCM[[i]]) <- "1"
        
      # Case if at least three members in clade:
      } else {
        
        # Prune taxa external to least inclusive clade to create a pruned tree:
        new.tree <- drop.tip(tree, nonclade.members)
        
        # Ensure root time is correct:
        new.tree <- CorrectRootTime(tree, new.tree)
        
        # Update tip states for new pruned tree:
        new.tips <- tip.states[clade.members]
        
        # Get branch along which (single) acqusition of derived state occurs:
        acquisition.branch <- match(new.root, tree$edge[, 2])
        
        # Get bounds for acquisition times:
        acquisition.bounds <- GetNodeAges(tree)[tree$edge[acquisition.branch, ]]
        
        # Draw an acusition time from a uniform distribution between the bounds:
        acquisition.time <- runif(1, min = acquisition.bounds[2], max = acquisition.bounds[1])
        
        # Create base stochastic character map (SCM):
        SCM <- as.list(c(1:nrow(tree$edge)))
        
        # For each branch in the SCM:
        for(i in 1:length(SCM)) {
          
          # Create a null value (vector equal to branch length):
          x <- tree$edge.length[i]
          
          # Label null vector with null (root) state of zero:
          names(x) <- "0"
          
          # Update base SCM with null value:
          SCM[[i]] <- x
          
        }
        
        # Create an SCM branch for the acquisition:
        acquisition.branch.SCM <- c(acquisition.bounds[1] - acquisition.time, acquisition.time - acquisition.bounds[2])
        
        # Add labels to the acqusition branch:
        names(acquisition.branch.SCM) <- c("0", "1")
        
        # Add acquisition branch to the SCM:
        SCM[[acquisition.branch]] <- acquisition.branch.SCM
        
        # If tip states vary:
        if(length(unique(new.tips)) > 1) {
          
          # Now do real SCM on pruned tree using Dollo model and a strong root prior of one:
          SCM_real <- make.simmap(new.tree, new.tips[new.tree$tip.label], model = Dollo.model, pi = c(0, 1), message = FALSE)$maps
          
        # If tip state is constant:
        } else {
          
          # Create SCM with no losses:
          SCM_real <- as.list(new.tree$edge.length)
          
          # Label all states as derived:
          for(i in 1:length(SCM_real)) names(SCM_real[[i]]) <- "1"
          
        }
        
        # Create edge matrix for original tree:
        orig.edges <- tree$edge
        
        # Create edge matrix for pruned tree:
        new.edges <- new.tree$edge
        
        # Update tip names for original tree edge matrix:
        for(i in 1:ape::Ntip(tree)) orig.edges[which(orig.edges[, 2] == i), 2] <- tree$tip.label[i]
        
        # Update tip names for pruned tree edge matrix:
        for(i in 1:ape::Ntip(new.tree)) new.edges[which(new.edges[, 2] == i), 2] <- new.tree$tip.label[i]
        
        # Update node names for original tree edge matrix:
        for(i in (ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree))) orig.edges[which(orig.edges == i)] <- paste(sort(tree$tip.label[FindDescendants(i, tree)]), collapse="")
        
        # Update node names for pruned tree edge matrix:
        for(i in (ape::Ntip(new.tree) + 1):(ape::Ntip(new.tree) + ape::Nnode(new.tree))) new.edges[which(new.edges == i)] <- paste(sort(new.tree$tip.label[FindDescendants(i, new.tree)]), collapse="")
        
        # Collapse original edge matrix to from-to straings for matching:
        orig.edges <- apply(orig.edges, 1, paste, collapse="%%TO%%")
        
        # Collapse pruned edge matrix to from-to straings for matching:
        new.edges <- apply(new.edges, 1, paste, collapse="%%TO%%")
        
        # Match edges between pruned and original trees:
        edge.matches <- match(new.edges, orig.edges)
        
        # Store real SCM in SCM for full tree:
        SCM[edge.matches] <- SCM_real
        
      }
      
    }
    
  }
  
  # Get node ages for tree in advance of establishing character change times:
  all.node.ages <- GetNodeAges(tree)
  
  # Create changes matrix:
  changes.matrix <- matrix(nrow = 0, ncol = 4, dimnames = list(c(), c("Branch", "From", "To", "Time")))
  
  # Get branches that record changes:
  change.branches <- which(lapply(SCM, length) > 1)
  
  # Special case of acquisition occurring prior to root (add extra line to changes matrix):
  if(acquisition.branch == 0) changes.matrix <- rbind(changes.matrix, c(acquisition.branch, 0, 1, acquisition.time))
  
  # As long as there are changes to record (i.e., it is not a constant character):
  if(length(change.branches) > 0) {
    
    # For each change branch:
    for(i in change.branches) {
      
      # Record branch number:
      change.branch <- i
      
      # Record from and to states:
      change.states <- as.numeric(names(SCM[[i]]))
      
      # Record change time:
      change.time <- all.node.ages[tree$edge[i, 1]] - SCM[[i]][1]
      
      # Make changes line ready for inserting into matrix:
      change.line <- c(change.branch, change.states, change.time)
      
      # Add changes to matrix:
      changes.matrix <- rbind(changes.matrix, change.line)
      
    }
    
  }
  
  # Clean up row names:
  rownames(changes.matrix) <- NULL
  
  # Build output:
  output <- list(changes.matrix, SCM)
  
  # Add names to output:
  names(output) <- c("Changes", "SCM")
  
  # Return output:
  return(output)
  
}
