#' Fixes a costmatrix that has inconsistent costs
#'
#' @description
#'
#' Given a costmatrix where transition costs are not self-consistent finds and returns a costmatrix that does.
#'
#' @param costmatrix A costMatrix object.
#' @param message A logical indicating whether messages should be output (defaults to \code{TRUE}).
#'
#' @details
#'
#' A user may wish to consider a complicated set of transition costs between states when modelling discrete character evolution. This can be achieved with a custom costmatrix in Claddis (and elsewhere). However, some caution is urged when using such matrices to ensure that these costs are \emph{self-consistent} (Maddison and Maddison 2003). More specifically, no direct state-to-state transition cost should be greater than is possible with an indirect path via one or more intermediate states.
#'
#' This function offers a solution through an algorithm that will iteratively alter a costmatrix until all direct transition costs are self-consistent. It does so by finding the shortest state-to-state path for all possible transitions using the \link{find_shortest_costmatrix_path} function. Because the first solution may itself be inconsistent (as it relied on costs that have since updated) the algorithm is repeated until an equilibrium is reached. (This scenario is unlikely in most real world cases, but may be possible with very large matrices representing many states so was implemented here for safety.)
#'
#' Note: infinite costs are allowed in costmatrices but unless they fill entire rows or columns (excluding the diagonal) they will not be self-consistent as there will always be a cheaper indirect cost.
#'
#' Note: that both PAUP* (Swofford 2003) TNT (Goloboff et al. 2008; Goloboff and Catalano, 2016) offerthe same correction using the triangle inequality.
#'
#' Note: other issues with a costmatrix may arise that are better revealed by using the \link{check_costMatrix} function, which returns informative error messages (with fix suggestions) where issues are found.
#'
#' @return
#'
#' A costMatrix object with self-consistent transition costs.
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
#' Swofford, D. L., 2003. \emph{PAUP*. Phylogenetic Analysis Using Parsimony (*and Other Methods). Version 4}. Sinauer Associates, Sunderland, Massachusetts.
#'
#' @examples
#'
#' # Build a custom costmatrix with non-self consistent path lengths:
#' costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "irreversible"
#' )
#' costmatrix$costmatrix[1:9] <- c(0, 2, 4, 1, 0, 3, 5, 3, 0)
#' costmatrix$symmetry <- "Asymmetric"
#' costmatrix$type <- "custom"
#'
#' # Fix costmatrix:
#' fixed_costmatrix <- fix_costmatrix(costmatrix = costmatrix)
#'
#' # Compare transition costs:
#' costmatrix$costmatrix
#' fixed_costmatrix$costmatrix
#' @export fix_costmatrix
fix_costmatrix <- function(costmatrix, message = TRUE) {
  
  # Check object is at least a costmatrix and stop and warn user if not:
  if (!inherits(x = costmatrix, what = "costMatrix")) stop("costmatrix must be an object of class \"costMatrix\".")
  
  # Sanity check (stop and warn user if all values are zero):
  if (all(costmatrix$costmatrix == 0)) stop("Costmatrix is all zero values! (A fix is impossible.)")
  
  # Ensure costmatrix is a custom type (only type that possibly needs fixing):
  if (all(costmatrix$type != "custom")) stop("Costmatrix is not of type \"custom\" - a fix is illogical.")
  
  # Store starting costmatrix (for later comparisons if changes):
  starting_costmatrix <- costmatrix
  
  # Subfunction to make shortest path costmatrix:
  make_shortest_path_costmatrix <- function(costmatrix) {
  
    # Get costmatrix states:
    states <- rownames(x = costmatrix$costmatrix)
  
    # Get length of each shortest path:
    path_lengths <- apply(
      X = expand.grid(start = states, end = states),
      MARGIN = 1,
      FUN = function(x) {
        path <- find_shortest_costmatrix_path(
          costmatrix = costmatrix,
          start = as.character(x = x[1]),
          end = as.character(x = x[2])
        )
        lengths <- lapply(
          X = as.list(x = 2:length(x = path)),
          function(i) costmatrix$costmatrix[path[(i - 1)], path[i]]
        )
        sum(x = unlist(x = lengths))
      }
    )
    
    # Make new shortest path costmatrix:
    costmatrix$costmatrix <- matrix(data = path_lengths, nrow = costmatrix$size, dimnames = list(states, states))
    
    # Return costmatrix:
    costmatrix
  }
  
  # Generate a (first) shortest path costmatrix:
  shortest_path_costmatrix <- make_shortest_path_costmatrix(costmatrix = costmatrix)
  
  # If at least one shorter path (smaller transition cost) exists:
  if (!all(costmatrix$costmatrix == shortest_path_costmatrix$costmatrix)) {
  
    # As long as the matrix has not reached an equilbrium:
    while(!all(costmatrix$costmatrix == shortest_path_costmatrix$costmatrix)) {
      
      # Reset costmatrix as current shortest path costmatrix:
      costmatrix <- shortest_path_costmatrix
      
      # Generate a new shortest path costmatrix from current costmatrix:
      shortest_path_costmatrix <- make_shortest_path_costmatrix(costmatrix = costmatrix)
    }
    
    # Sanity check for new costmatrix symmetry:
    new_symmetry <- ifelse(test = isSymmetric(object = costmatrix$costmatrix), yes = "Symmetric", no = "Asymmetric")
    
    # If new shortest costmatrix changes symmetry of input matrix:
    if (starting_costmatrix$symmetry != new_symmetry) {
      
      # Return message to user regarding symmetry change (if messages are requested):
      if (message) cat(paste("costmatrix symmetry has changed from ", starting_costmatrix$symmetry, " to ", new_symmetry, ".\n", sep = ""))
      
      # Update symmetry value before outputting:
      costmatrix$symmetry <- new_symmetry
    }
    
    # Return new shortest path costmatrix:
    return(costmatrix)
    
  # Case if costmatrix:
  } else {
    
    # Return message to user (if messages are requested):
    if (message) cat("costmatrix doesn't need fixing - returning original costmatrix.\n")
    
    # Return costmatrix unaltered:
    return(costmatrix)
  }
}
