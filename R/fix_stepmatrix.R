#' Fixes a stepmatrix that has inconsistent costs
#'
#' @description
#'
#' Given a stepmatrix where transition costs are not self-consistent finds and returns a stepmatrix that does.
#'
#' @param stepmatrix A stepMatrix object.
#' @param message A logical indicating whether messages should be output (defaults to \code{TRUE}).
#'
#' @details
#'
#' A user may wish to consider a complicated set of transition costs between states when modelling discrete character evolution. This can be achieved with a custom stepmatrix in Claddis (and elsewhere). However, some caution is urged when using such matrices to ensure that these costs are \emph{self-consistent}. More specifically, no direct state-to-state transition cost (excluding infinite costs) should be greater than is possible with an indirect path via one or more intermediate states.
#'
#' This problem was first recognised (to the best of my knowledge) by Maddison and Maddison (2003), although they offered no direct solution.
#'
#' This function offers a solution through an algorithm that will iteratively alter a stepmatrix until all direct transition costs are self-consistent. It does so by finding the shortest state-to-state path for all possible transitions using the \link{find_shortest_stepmatrix_path} function. Because the first solution may itself be inconsistent (as it relied on costs that have since updated) the algorithm is repeated until an equilbirum is reached. (This scenario is unlikely in most real world cases, but may be possible with very large matrices representing many states so was implemented for safety.)
#'
#' Note: that TNT (Goloboff et al. 2008; Goloboff and Catalano, 2016) offers a correction based on the triangle inequality which appears to be an attempt to solve the same problem.
#'
#' Note: other issues with a stepmatrix may arise that are better revealed by using the \link{check_stepMatrix} function which has informative error messages for how to fix any errors in format.
#'
#' @return
#'
#' A stepMatrix object with self-consistent transition costs.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics. \emph{Cladistics}, \bold{32}, 221-238.
#'
#' Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. \emph{Cladistics}, \bold{24}, 774-786.
#'
#' Maddison, D. R. and Maddison, W. P., 2003. \emph{MacClade 4: Analysis of phylogeny and character evolution}. Version 4.06. Sinauer Associates, Sunderland, Massachusetts.
#'
#' @examples
#'
#' # Build a custom stepmatrix with non-self consistent path lengths:
#' stepmatrix <- make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "irreversible",
#'   polymorphism_shape = "simplex",
#'   polymorphism_distance = "euclidean"
#' )
#' stepmatrix$stepmatrix[1:9] <- c(0, 2, 4, 1, 0, 3, 5, 3, 0)
#' stepmatrix$symmetry <- "Asymmetric"
#' stepmatrix$type <- "custom"
#' stepmatrix
#'
#' # Fix stepmatrix:
#' fixed_stepmatrix <- fix_stepmatrix(stepmatrix = stepmatrix)
#'
#' # Compare transition costs:
#' stepmatrix$stepmatrix
#' fixed_stepmatrix$stepmatrix
#' @export fix_stepmatrix
fix_stepmatrix <- function(stepmatrix, message = TRUE) {
  
  # Check object is at least a stepmatrix and stop and warn user if not:
  if (class(stepmatrix) != "stepMatrix") stop("stepmatrix must be an object of class \"stepMatrix\".")
  
  # Sanity check (stop and warn user if all values are zero):
  if (all(stepmatrix$stepmatrix == 0)) stop("Stepmatrix is all zero values! (A fix is impossible.)")
  
  # Ensure stepmatrix is a custom type (only type that possibly needs fixing):
  if (all(stepmatrix$type != "custom")) stop("Stepmatrix is not of type \"custom\" - a fix is illogical.")
  
  # Store starting stepmatrix (for later comparisons if changes):
  starting_stepmatrix <- stepmatrix
  
  # Subfunction to make shortest path stepmatrix:
  make_shortest_path_stepmatrix <- function(stepmatrix) {
  
    # Get stepmatrix states:
    states <- rownames(x = stepmatrix$stepmatrix)
  
    # Get length of each shortest path:
    path_lengths <- apply(
      X = expand.grid(start = states, end = states),
      MARGIN = 1,
      FUN = function(x) {
        path <- find_shortest_stepmatrix_path(
          stepmatrix = stepmatrix,
          start = as.character(x = x[1]),
          end = as.character(x = x[2])
        )
        lengths <- lapply(
          X = as.list(x = 2:length(x = path)),
          function(i) stepmatrix$stepmatrix[path[(i - 1)], path[i]]
        )
        sum(x = unlist(x = lengths))
      }
    )
    
    # Make new shortest path stepmatrix:
    stepmatrix$stepmatrix <- matrix(data = path_lengths, nrow = stepmatrix$size, dimnames = list(states, states))
    
    # Return stepmatrix:
    stepmatrix
  }
  
  # Generate a (first) shortest path stepmatrix:
  shortest_path_stepmatrix <- make_shortest_path_stepmatrix(stepmatrix = stepmatrix)
  
  # If at least one shorter path (smaller transition cost) exists:
  if (!all(stepmatrix$stepmatrix == shortest_path_stepmatrix$stepmatrix)) {
  
    # As long as the matrix has not reached an equilbrium:
    while(!all(stepmatrix$stepmatrix == shortest_path_stepmatrix$stepmatrix)) {
      
      # Reset stepmatrix as current shortest path stepmatrix:
      stepmatrix <- shortest_path_stepmatrix
      
      # Generate a new shortest path stepmatrix from current stepmatrix:
      shortest_path_stepmatrix <- make_shortest_path_stepmatrix(stepmatrix = stepmatrix)
    }
    
    # Sanity check for new stepmatrix symmetry:
    new_symmetry <- ifelse(test = isSymmetric(object = stepmatrix$stepmatrix), yes = "Asymmetric", no = "Symmetric")
    
    # If new shortest stepmatrix changes symmetry of input matrix:
    if (starting_stepmatrix$symmetry != new_symmetry) {
      
      # Return message to user regarding symmetry change (if messages are requested):
      if (message) cat(paste("stepmatrix symmetry has changed from ", starting_stepmatrix$symmetry, " to ", new_symmetry, ".\n", sep = ""))
      
      # Update symmetry value before outputting:
      stepmatrix$symmetry <- new_symmetry
    }
    
    # Return new shortest path stepmatrix:
    return(stepmatrix)
    
  # Case if stepmatrix:
  } else {
    
    # Return message to user (if messages are requested):
    if (message) cat("stepmatrix doesn't need fixing - returning original stepmatrix.\n")
    
    # Return stepmatrix unaltered:
    return(stepmatrix)
  }
}
