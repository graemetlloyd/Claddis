#' Reconstruct ancestral states using parsimony
#'
#' @description
#'
#' Given a tree, discrete states for each tip, and either a specific stepmatrix or character type, returns the most parsimonious ancestral reconstruction(s).
#'
#' @param tree A tree (phylo object).
#' @param tip_states A labelled vector of tip states. These should be discrete and mathc the row and column headings in \code{stepmatrix}, with labels matching the tip labels of \code{tree}.
#' @param stepmatrix Either the character type (one of \code{"ordered"} or \code{"unordered"}) or a custom-specified step matrix. If the latter this must be square with identical row and column names correspnding to the values in \code{tip_states}. The diagonal must be zero, but off-diagonal values can be any non-negative real number. I.e., they needn't be integers, and the matrix does not need to be symmetric. Note that for transitions the rows are considered the "from" values, and the columns the "to" value. Thus a cost for the transition from state "0" to state "1" would be specified by \code{stepmatrix["0", "1"]}. (This is only relevant where matrices are asymmetric in tranistion costs.)
#'
#' @details
#'
#' Text.
#'
#' # Uses Swofford and Maddison 1992 general solution as covers more base's but is also slower.
#' # Rooted tres only as direction fucking matters and roots have meaning (stepmatrices that are asymmetric)
#' # Option to fix root state (polarize) somehow
#' # From is row, to is column

#' # MANUAL: generalized solution as accounts for more possibilities
#' # MANUAL: really only for rooted trees (unrooted possible using same methods (see Swofford and Maddison(1992), but unclear why you would want to do this with trees without direction of time (i.e., a root).
#' # MANUAL: works with multifirctaions and doesn't reuire a fully bifurcating tree, but does assume polytomies are therefore "hard".
#' # MANUAL: Swofford and Maddison (1992; p210): "Both [ACCTRAN and DELTRAN] start by either assuming that the ancestral state is known or by choosing a state from the set of optimal assignments to the root node" - i.e., always an arbitrainess problem at root!
#' MISSING DATA (NA) IS DEALT WITH BY SETTING ALL TIP VALUES FOR A MISSING STATE TO ZERO (AS PER SWOFFORD AND MADDISON 1992).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. \emph{In}  R. L. Mayden (ed.) Systematics, Historical Ecology, and North American Freshwater Fishes. Stanford University Press, Stanford, p187-223.
#'
#' @return
#'
#' A matrix where rows correspond to ALL nodes (i.e., terminal and internal) in the order numbered by \code{ape}, and columns correspond to unique most parsimnious reconstruction(s). I.e., if there is only one most parsimnious reconstruction there will be only one column.
#'
#' Users may also wish to refer to the more complex whole-matrix, likelihood-based function \link{estimate_ancestral_states}. Although, note that eventually parsimony will simply be an option to that function.
#'
#' @seealso
#'
#' \link{estimate_ancestral_states}
#'
#' @examples
#'
#' # Set up the example from Swofford and Maddison 1992:
#' tree <- ape::read.tree(text = "((C,(A,B)),(D,E));")
#' tip_states <- c(A = 1, B = 2, C = 2, D = 1, E = 0)
#' stepmatrix <- matrix(
#'   data = c(0, 1, 2, 1, 0, 1, 2, 1, 0),
#'   nrow = 3,
#'   ncol = 3,
#'   dimnames = list(0:2, 0:2)
#' )
#'
#' # Look at stepmatrix to confirm it is symmetric and ordered:
#' stepmatrix
#'
#' # Reconstruct ancestral states:
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )
#'
#' # Plot all most parsimonious reconstructions on the tree:
#' par(mfrow = c(1, ncol(x = node_estimates)))
#' for(i in 1:ncol(x = node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Repeat using the simpler "ordered" option:
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(x = node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Repeat but considering the character "unordered":
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "unordered"
#' )
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(x = node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Repeat using a stepmatrix that penalises gains:
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = matrix(
#'     data = c(0, 1, 2, 999, 0, 1, 999, 999, 0),
#'     nrow = 3,
#'     ncol = 3,
#'     dimnames = list(0:2, 0:2)
#'   )
#' )
#' par(mfrow = c(1, ncol(x = node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Repeat using a stepmatrix that penalises losses:
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = matrix(
#'     data = c(0, 999, 999, 1, 0, 999, 2, 1, 0),
#'     nrow = 3,
#'     ncol = 3,
#'     dimnames = list(0:2, 0:2)
#'   )
#' )
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(x = tree, show.tip.label = FALSE)
#'   tiplabels(text = tip_states[tree$tip.label], cex = 2)
#'   nodelabels(text = node_estimates[6:9, i], cex = 2)
#' }
#'
#' # Ordered with polytomies:
#' tree <- ape::read.tree(text = "((C,A,B),(D,E));") # Two non-root internal
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(tree, show.tip.label = FALSE)
#'   tiplabels(tip_states[tree$tip.label], cex = 2)
#'   nodelabels(node_estimates[6:8, i], cex = 2)
#' }
#' tree <- ape::read.tree(text = "((C,A,B),D,E);") # One non-root internal
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(tree, show.tip.label = FALSE)
#'   tiplabels(tip_states[tree$tip.label], cex = 2)
#'   nodelabels(node_estimates[6:7, i], cex = 2)
#' }
#' tree <- ape::read.tree(text = "(C,A,B,D,E);") # Star tree
#' node_estimates <- reconstruct_ancestral_states(
#'   tree = tree,
#'   tip_states = tip_states,
#'   stepmatrix = "ordered"
#' )
#' par(mfrow = c(1, ncol(node_estimates)))
#' for(i in 1:ncol(node_estimates)) {
#'   plot(tree, show.tip.label = FALSE)
#'   tiplabels(tip_states[tree$tip.label], cex = 2)
#'   nodelabels(node_estimates[6, i], cex = 2)
#' }
#'
#' @export reconstruct_ancestral_states
reconstruct_ancestral_states <- function(tree, tip_states, stepmatrix) {
  
  # TO DO
  # - Allow for polymorphisms (maybe conditionals for counting changes?) This is a hard problem!
  # - How to allow for missing or uncertainty? This bit should be easy as set all states to zero, or all uncertain states to zero at tips.
  # - ACCTRAN/DELTRAN/(Random?)/Branch lengths as weights to help collapse?/Uncertainity (e.g., 0/1)
  # - Character weights? (Just multiply through stepmatrix by weight?)
  
  # CHECKS TO WRITE
  # - states are discrete
  # - states in tip states are all available in stepmatrix and vice versa
  # - stepmatrix values are zero or positive (don't have to be integers?)
  # - Check there are enough states to do something with!
  # - Put checks at higher level to avoid repeating for every character in a matrix (slow?)
  
  # If stepmatrix is not already specified as a matrix (i.e., it is simply "ordered", "unordrered" etc.):
  if(!is.matrix(x = stepmatrix)) {
    
    # NEED TO SET MORE OPTIONS LIKE DOLLO, CAMIN-SOKAL AND....?
    
    # First check stepmatrix is of a vaid type and stop wand warn user if not:
    if(length(x = setdiff(x = stepmatrix, y = c("ordered", "unordered"))) > 0) stop("If not a specific matrix, then stepmatrix must be one of \"ordered\" or \"unordered\".")
    
    # Get tip states as numerics:
    tip_state_numbers <- as.numeric(x = tip_states)
    
    # Get the range (min-max) of tip values:
    tip_state_range <- range(x = tip_state_numbers[!is.na(x = tip_state_numbers)])
    
    # Get the size (n rows and columns) required for the stepmatrix:
    stepmatrix_size <- diff(x = tip_state_range) + 1
    
    # Create base stepmatrix (to be altered further):
    base_stepmatrix <- matrix(data = NA, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(tip_state_range[1]:tip_state_range[2], tip_state_range[1]:tip_state_range[2]))
    
    # Add diagonal of zero as this should always be true regardless of character type:
    diag(x = base_stepmatrix) <- 0
    
    # Case if character is ordered:
    if(stepmatrix == "ordered") {
      
      # Populate upper triangle (cost is difference between states):
      base_stepmatrix[upper.tri(x = base_stepmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(stepmatrix_size - 1)), FUN = function(steps) steps:1))
      
      # Populate lower triangle (cost is difference between states):
      base_stepmatrix[lower.tri(x = base_stepmatrix)] <- unlist(x = lapply(X = as.list(x = (stepmatrix_size - 1):1), FUN = function(steps) 1:steps))
      
    }
    
    # Case if character is unordered all off-diagonal values cost one step):
    if(stepmatrix == "unordered") base_stepmatrix[is.na(x = base_stepmatrix)] <- 1
    
    # Overwrite stepmatrix with the fully formed stepmatrix:
    stepmatrix <- base_stepmatrix
    
  }
  
  # Reorder tip states (1 to N):
  tip_states <- tip_states[tree$tip.label]
  
  # Establish number of tips:
  tip_count <- ape::Ntip(phy = tree)
  
  # Establish total number of nodes (terminal and internal):
  node_count <- tip_count + tree$Nnode
  
  # Initialise node_values matrix:
  node_values <- matrix(data = NA, nrow = node_count, ncol = ncol(x = stepmatrix), dimnames = list(c(), colnames(stepmatrix)))
  
  # Begin by inserting Inf as default tip value:
  node_values[1:tip_count, ] <- Inf
  
  # Populate tip values (zeroes for assigned state(s):
  for(i in 1:tip_count) node_values[i, colnames(x = node_values) == tip_states[i]] <- 0
  
  # If any tips are missing then set all states for those as zero:
  if(any(x = is.na(x = tip_states))) node_values[which(x = is.na(x = tip_states)), ] <- 0
  
  # First pass (traverse tree from tips to root):
  for(needle in (tip_count + node_count - tip_count):(tip_count + 1)) {
    
    # Find decsendants of current node:
    descendants <- tree$edge[tree$edge[, 1] == needle, 2]
    
    # Calculate and store new node values:
    node_values[needle, ] <- unlist(x = lapply(X = as.list(x = colnames(x = node_values)), FUN = function(fromstate) sum(x = unlist(x = lapply(X = as.list(x = descendants), FUN = function(descendant) min(x = node_values[descendant, ] + stepmatrix[fromstate, ]))))))
    
  }
  
  # ADD OPTION TO FIX ROOT HERE
  # CHECK OPTION DOESN'T FORCE INFINITY!
  
  tree_length <- min(x = node_values[tip_count + 1, ])
  
  # STOP AFTER HERE IF ONLY WANT TREE LENGTH?
  
  # WILL NEED TO MODIFY BELOW TO DEAL WITH UNCERTAINTIES AND POLYMORPHISMS AT TIPS
  
  # Build node estimates from single state results:
  node_estimates <- matrix(data = apply(X = node_values, MARGIN = 1, FUN = function(x) {x <- names(x = x[x == min(x = x)]); ifelse(test = length(x = x) == 1, x, NA)}), ncol = 1, nrow = node_count)
  
  # Only continue if there is any ambiguity:
  if(any(x = is.na(x = node_estimates))) {
    
    # Find all possible root states:
    possible_root_states <- colnames(x = node_values)[node_values[tip_count + 1, ] == min(x = node_values[tip_count + 1, ])]
    
    # Make new node estimates with every possible root state:
    node_estimates <- do.call(what = cbind, args = lapply(X = as.list(x = possible_root_states), FUN = function(x) {y <- node_estimates; y[tip_count + 1, ] <- x; y}))
      
    # For each internal node above the root to the tips:
    for(needle in (tip_count + 2):(tip_count + node_count - tip_count)) {
      
      # Establish ancestor of current node:
      ancestor_node <- tree$edge[tree$edge[, 2] == needle, 1]
      
      # Reformat node_estimates as a list to enable lapply:
      node_estimates <- split(x = node_estimates, f = col(x = node_estimates))
      
      # Permute all possible values and reformat as matrix:
      node_estimates <- do.call(what = cbind, args = lapply(X = node_estimates, FUN = function(x) {
        
        # Get updated tree lengths for current node as per Swofford and Maddison 1992 second pass:
        updated_tree_lengths <- stepmatrix[x[ancestor_node], ] + node_values[needle, ]
        
        # Store all possible most parsimonious state(s) for current node:
        possible_states <- names(x = updated_tree_lengths[updated_tree_lengths == min(x = updated_tree_lengths)])
        
        # Store total number of possible states:
        n_states <- length(x = possible_states)
        
        # Create new node estimates to store possible state(s):
        new_estimates <- matrix(data = rep(x = x, times = n_states), ncol = n_states)
        
        # Add possible state(s)
        new_estimates[needle, ] <- possible_states
        
        # Return new estimates only:
        new_estimates
        
      }))
      
    }
    
  }

  # Return node estimates:
  node_estimates

  # OUTPUT:
  # - Ancestral states
  # - Ambiguities?
  # - Stepmatrix used
  # - Tree used
  # - Tip states used
  # - Tree length
  # - How root value was chosen (arbitraty forced, unambiguous?)
  # - Algorithm (parsimony, ML, etc.)
  # - CI?
  # - RI?
  # - N reversals
  # - N parallelisms
  # - MAKE THIS A CLASS TOO!
  
}
