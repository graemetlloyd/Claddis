#' Calculate the maximum possible tree length, gmax, under parsimony
#'
#' @description
#'
#' Given a costmatrix and number of tips, calculates the longest possible tree length under maximum parsimony.
#'
#' @param costmatrix An object of class \code{costMatrix}.
#' @param n_taxa The number of taxa (ree tips) to consider.
#' @param allow_zeroes Logical indicating whether some states are allowed to be assigned zero tips. Defaults to \code{FALSE} (recommended). Note that if the number of states exceeds the number of tips then \code{allow_zeroes = TRUE} is forced.
#'
#' @details
#'
#' The concept of \emph{gmax} was introduced by Hoyal Cuthill et al. (2010) and in the strictest sense represents the integer state frequency that maximizes \emph{g} with the restriction that every state is sampled at least once. This function is intended to relax that slightly using the \code{allow_zeroes} option and indeed this is not possible in rare cases such as either high levels of missing data or gap-weighted (Thiele 1993) characters where it is possible for the number of states to exceed the (effective) number of tips. For a more detailed description of the problem and implemented solution(s) see Hoyal Cuthill and Lloyd (in press).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Jen Hoyal Cuthill \email{j.hoyal-cuthill@@essex.ac.uk}
#'
#' @references
#'
#' Hoyal Cuthill, J. F. and Lloyd, G. T., in press. Measuring homoplasy I: comprehensive measures of maximum and minimum cost under parsimony across discrete cost matrix character types. \emph{Cladistics}, bold{}, .
#'
#' Hoyal Cuthill, J. F., Braddy, S. J. and Donoghue, P. C. J., 2010. A formula for maximum possible steps in multistate characters: isolating matrix parameter effects on measures of evolutionary convergence. \emph{Cladistics}, \bold{26}, 98-102.
#'
#' Thiele, K.. 1993. The Holy Grail of the perfect character: the cladistic treatment of morphometric data. \emph{Cladistics}, \bold{9}, 275-304.
#'
#' @return
#'
#' A numeric indicating the maximum possible cost, \emph{gmax}.
#'
#' @seealso
#'
#' \link{calculate_g}
#'
#' @examples
#'
#' # Create a Type I character costmatrix:
#' constant_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 0,
#'   character_type = "unordered"
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = constant_costmatrix,
#'   n_taxa = 10
#' )
#'
#' # Create a Type II character costmatrix:
#' binary_symmetric_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 1,
#'   character_type = "unordered"
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = binary_symmetric_costmatrix,
#'   n_taxa = 10
#' )
#'
#' # Create a Type III character costmatrix:
#' unordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered"
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = unordered_costmatrix,
#'   n_taxa = 10
#' )
#'
#' # Create a Type IV character costmatrix:
#' linear_ordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "ordered"
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = linear_ordered_costmatrix,
#'   n_taxa = 10
#' )
#'
#' # Create a Type V character costmatrix:
#' nonlinear_ordered_costmatrix <- convert_adjacency_matrix_to_costmatrix(
#'   adjacency_matrix = matrix(
#'     data = c(
#'       0, 1, 0, 0,
#'       1, 0, 1, 1,
#'       0, 1, 0, 0,
#'       0, 1, 0, 0
#'     ),
#'     nrow = 4,
#'     dimnames = list(0:3, 0:3)
#'   )
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = nonlinear_ordered_costmatrix,
#'   n_taxa = 10
#' )
#'
#'
#' # Create a Type VI character costmatrix:
#' binary_irreversible_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 1,
#'   character_type = "irreversible"
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = binary_irreversible_costmatrix,
#'   n_taxa = 10
#' )
#'
#' # Create a Type VII character costmatrix:
#' multistate_irreversible_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "irreversible"
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = multistate_irreversible_costmatrix,
#'   n_taxa = 10
#' )
#'
#' #' # Create a Type VIII character costmatrix:
#' binary_dollo_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 1,
#'   character_type = "dollo"
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = binary_dollo_costmatrix,
#'   n_taxa = 10
#' )
#'
#' # Create a Type IX character costmatrix:
#' multistate_dollo_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "dollo"
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = multistate_dollo_costmatrix,
#'   n_taxa = 10
#' )
#'
#' # Create a Type X character costmatrix:
#' multistate_symmetric_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 5,
#'   character_type = "ordered"
#' )
#' multistate_symmetric_costmatrix$type <- "custom"
#' multistate_symmetric_costmatrix$costmatrix <- matrix(
#'   data = c(
#'     0, 1, 2, 3, 2, 3,
#'     1, 0, 3, 2, 1, 2,
#'     2, 3, 0, 3, 2, 1,
#'     3, 2, 3, 0, 1, 2,
#'     2, 1, 2, 1, 0, 1,
#'     3, 2, 1, 2, 1, 0
#'   ),
#'   nrow = multistate_symmetric_costmatrix$size,
#'   ncol = multistate_symmetric_costmatrix$size,
#'   byrow = TRUE,
#'   dimnames = list(
#'     multistate_symmetric_costmatrix$single_states,
#'     multistate_symmetric_costmatrix$single_states
#'   )
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = multistate_symmetric_costmatrix,
#'   n_taxa = 10
#' )
#'
#' # Create a Type XI character costmatrix:
#' binary_asymmetric_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 1,
#'   character_type = "ordered"
#' )
#' binary_asymmetric_costmatrix$type <- "custom"
#' binary_asymmetric_costmatrix$costmatrix <- matrix(
#'   data = c(
#'     0, 1,
#'     10, 0
#'   ),
#'   nrow = binary_asymmetric_costmatrix$size,
#'   ncol = binary_asymmetric_costmatrix$size,
#'   byrow = TRUE,
#'   dimnames = list(
#'     binary_asymmetric_costmatrix$single_states,
#'     binary_asymmetric_costmatrix$single_states
#'   )
#' )
#' binary_asymmetric_costmatrix$symmetry <- "Asymmetric"
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = binary_asymmetric_costmatrix,
#'   n_taxa = 10
#' )
#'
#' # Create a Type XII character costmatrix:
#' multistate_asymmetric_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "ordered"
#' )
#' multistate_asymmetric_costmatrix$type <- "custom"
#' multistate_asymmetric_costmatrix$costmatrix <- matrix(
#'   data = c(
#'     0, 1, 1,
#'     1, 0, 1,
#'     10, 10, 0
#'   ),
#'   nrow = multistate_asymmetric_costmatrix$size,
#'   ncol = multistate_asymmetric_costmatrix$size,
#'   byrow = TRUE,
#'   dimnames = list(
#'     multistate_asymmetric_costmatrix$single_states,
#'     multistate_asymmetric_costmatrix$single_states
#'   )
#' )
#' multistate_asymmetric_costmatrix$symmetry <- "Asymmetric"
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = multistate_asymmetric_costmatrix,
#'   n_taxa = 10
#' )
#'
#' @export calculate_gmax
calculate_gmax <- function(costmatrix, n_taxa, allow_zeroes = FALSE) {
  
  ### IGNORES POLYMORPHISMS AND UNCERTAINTIES FOR TIPS, BUT COULD PERHAPS ALLOW FORMER FOR ANCESTRAL STATE?
  ### CHECK N TAXA IS A SENSIBLE NUMBER! (Below assumes at least 1 but zero, negatives or decimals are bad)
  
  # Check if allow_zeroes is consistent with tip and state counts:
  if (allow_zeroes && costmatrix$n_states > n_taxa) {
    
    # Reset allow_zeroes:
    allow_zeroes <- TRUE
    
    # Warn user:
    print("Zeroes must be allowed as costmatrix$n_states exceeds n_taxa, allow_zeroes is reset to TRUE.")
  }
  
  # Special case of a Type I constant character (return zero):
  if (costmatrix$n_states == 1 || n_taxa <= 1) return(0)
  
  # Special case of a binary character (will never want to allow zeroes for these as would then be a constant character and gmax > 0 as long as there are at least two states and two tips):
  if (costmatrix$n_states == 2) {
    
    # Special case of a Type II standard binary character (return gmax according to formula from Hoyal Cuthill and Lloyd):
    if (length(x = setdiff(x = costmatrix$type, y = c("ordered", "unordered"))) == 0) return(costmatrix$costmatrix[1, 2] * (n_taxa - ceiling(x = n_taxa / 2)))
    
    # Special case of a Type VI standard binary irreversible character (return gmax according to formula from Hoyal Cuthill and Lloyd):
    if (length(x = setdiff(x = costmatrix$type, y = c("irreversible", "stratigraphy"))) == 0) return(costmatrix$costmatrix[1, 2] * (n_taxa - 1))
    
    # Special case of a Type XI custom binary character (or non-standard Type II or Type VI character) (return gmax according to formula from Hoyal Cuthill and Lloyd):
    if (length(x = setdiff(x = costmatrix$type, y = c("custom"))) == 0) return(max(x = min(x = (n_taxa - floor(x = (n_taxa * (costmatrix$costmatrix[1, 2] / (costmatrix$costmatrix[1, 2] + costmatrix$costmatrix[2, 1]))))) * costmatrix$costmatrix[1, 2], floor(x = (n_taxa * (costmatrix$costmatrix[1, 2] / (costmatrix$costmatrix[1, 2] + costmatrix$costmatrix[2, 1])))) * costmatrix$costmatrix[2, 1]), min(x = (n_taxa - ceiling((n_taxa * (costmatrix$costmatrix[1, 2] / (costmatrix$costmatrix[1, 2] + costmatrix$costmatrix[2, 1]))))) * costmatrix$costmatrix[1, 2], ceiling(x = (n_taxa * (costmatrix$costmatrix[1, 2] / (costmatrix$costmatrix[1, 2] + costmatrix$costmatrix[2, 1])))) * costmatrix$costmatrix[2, 1])))
    
    # Special case of a Type VIII binary Dollo character (return gmax according to formula from Hoyal Cuthill and Lloyd):
    if (length(x = setdiff(x = costmatrix$type, y = c("dollo"))) == 0) {
      
      # If n tips and n states are same (i.e., two) return gmax according to formula from Hoyal Cuthill and Lloyd:
      if (n_taxa == 2) return(costmatrix$costmatrix[2, 1] * (n_taxa - 1))
    
      # If more tips than states return gmax according to formula from Hoyal Cuthill and Lloyd:
      if (n_taxa >= 2) return(costmatrix$costmatrix[2, 1] * (n_taxa - 2))
    }
  }
  
  # Special case of a multistate character:
  if (costmatrix$n_states > 2) {
    
    # Special case of a Type III multistate unordered character (return gmax according to formula from Hoyal Cuthill and Lloyd):
    if (length(x = setdiff(x = costmatrix$type, y = c("unordered"))) == 0) return(costmatrix$costmatrix[1, 2] * (n_taxa - ceiling(x = n_taxa / costmatrix$n_states)))
    
    # If zeroes are not permitted:
    if (!allow_zeroes) {
    
      # Special case of a Type IV multistate linear ordered character (return gmax according to formula from Hoyal Cuthill and Lloyd):
      if (length(x = setdiff(x = costmatrix$type, y = c("ordered"))) == 0) return(costmatrix$costmatrix[1, 2] * floor(x = (costmatrix$n_states - 3) + (n_taxa - (costmatrix$n_states - 2)) * ((costmatrix$n_states - 1) / 2)))
    
      # Special case of a Type VII multistate irreversible character (return gmax according to formula from Hoyal Cuthill and Lloyd):
      if (length(x = setdiff(x = costmatrix$type, y = c("irreversible"))) == 0) return(costmatrix$costmatrix[1, 2] * ((n_taxa - costmatrix$n_states + 1) * costmatrix$costmatrix[1, costmatrix$n_states] + sum(costmatrix$costmatrix[1, (1:costmatrix$n_states - 1)])))
    
      # Special case of a Type IX multistate Dollo character:
      if (length(x = setdiff(x = costmatrix$type, y = c("dollo"))) == 0) {
      
        # If n tips and n states are same return gmax according to formula from Hoyal Cuthill and Lloyd:
        if (n_taxa < costmatrix$n_states + 1) return(costmatrix$costmatrix[2, 1] * ((costmatrix$n_states - 2) * (n_taxa - costmatrix$n_states + 1) + (1 / 2) * (costmatrix$n_states - 2) * (costmatrix$n_states - 3) + 1))
      
        # If more tips than states return gmax according to formula from Hoyal Cuthill and Lloyd:
        if (n_taxa >= costmatrix$n_states + 1) return(costmatrix$costmatrix[2, 1] * ((costmatrix$n_states - 1) * (n_taxa - costmatrix$n_states) + (1 / 2) * (costmatrix$n_states - 1) * (costmatrix$n_states - 2)))
      }
      
      # If still going must be one of remaining types:
      #
      # Type V - Multistate non-linear ordered
      # Type X - Multisatte custom symmetric
      # Type XII - Multistate custom asymmetric

      # Before adopting exhaustive solution check it is computationally feasible:
      exhaustive_problem_size <- choose(n = n_taxa, k = costmatrix$n_states)
      if (exhaustive_problem_size > 40000000) stop("Cost matrix is of a type without a direct solution and combination of state and tip counts is too large for an exhaustive solution.")
      
      # Permute all possible state frequencies:
      permuted_state_frequencies <- permute_restricted_compositions(n = n_taxa, m_labels = costmatrix$single_states, allow_zero = FALSE)
      
      # Calculate gmax via an exhaustive search of possible state frequencies:
      exhaustive_gmax <- max(
        x = apply(
          X = permuted_state_frequencies,
          MARGIN = 1,
          FUN = function(i) {
            cost_for_each_ancestral_state <- apply(
              X = matrix(data = rep(x = i, each = costmatrix$n_states), nrow = costmatrix$n_states) * costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states],
              MARGIN = 1,
              FUN = sum
            )
            # Maximum parsimony cost:
            min(x = cost_for_each_ancestral_state)
          }
        )
      )
      
      # Return exhaustive gmax:
      return(exhaustive_gmax)
    }
    
    # If zeroes are permitted:
    if (allow_zeroes) {
      
      # Special case of a Type IV multistate linear ordered character:
      if (length(x = setdiff(x = costmatrix$type, y = c("ordered"))) == 0) return(costmatrix$costmatrix[1, 2] * floor(x = (n_taxa / 2) * (costmatrix$n_states - 1)))
      
      # Special case of a Type VII multistate irreversible character:
      if (length(x = setdiff(x = costmatrix$type, y = c("irreversible"))) == 0) return(costmatrix$costmatrix[1, 2] * ((n_taxa - 1) * (costmatrix$n_states - 1)))
      
      # Special case of a Type IX multistate Dollo character:
      if (length(x = setdiff(x = costmatrix$type, y = c("dollo"))) == 0) {
        
        # If n tips and n states are same:
        if (n_taxa == 2) return(costmatrix$costmatrix[2, 1] * (costmatrix$n_states - 1))
        
        # If more tips than states:
        if (n_taxa >= 2) return(costmatrix$costmatrix[2, 1] * ((costmatrix$n_states - 1) * (n_taxa - 2)))
      }
      
      # If still going must be one of remaining types:
      #
      # Type V - Multistate non-linear ordered
      # Type X - Multisatte custom symmetric
      # Type XII - Multistate custom asymmetric
      
      # Before adopting exhaustive solution check it is computationally feasible:
      exhaustive_problem_size <- choose(n = n_taxa + costmatrix$n_states - 1, k = costmatrix$n_states - 1)
      if (exhaustive_problem_size > 40000000) stop("Cost matrix is of a type without a direct solution and combination of state and tip counts is too large for an exhaustive solution.")
      
      # Permute all possible state frequencies:
      permuted_state_frequencies <- permute_restricted_compositions(n = n_taxa, m_labels = costmatrix$single_states, allow_zero = TRUE)

      # Calculate gmax via an exhaustive search of possible state frequencies:
      exhaustive_gmax <- max(
        x = apply(
          X = permuted_state_frequencies,
          MARGIN = 1,
          FUN = function(i) {
            cost_for_each_ancestral_state <- apply(
              X = matrix(data = rep(x = i, each = costmatrix$n_states), nrow = costmatrix$n_states) * costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states],
              MARGIN = 1,
              FUN = sum
            )
            # Maximum parsimony cost:
            min(x = cost_for_each_ancestral_state)
          }
        )
      )
      
      # Return exhaustive gmax:
      return(exhaustive_gmax)
    }
  }
  
  # Error check if not already returned out of function:
  stop("Something has gone horribly wrong if you see this error as every possibility should have already returned a value from calculate_gmax!")
}

#THIS IS NOT A PROPER TYPE V CHARACTER. MUST HAVE A VERTEX OF DEGREE 3!
#graph_size <- 6
#tip_count <- 10
#adjacency_matrix <- matrix(
#  data = 0,
#  nrow = graph_size,
#  ncol = graph_size,
#  dimnames = list(0:(graph_size - 1), 0:(graph_size - 1))
#)
# Add minimum links (creating a subgraph that is a path graph:
#for(i in 1:(graph_size - 1)) adjacency_matrix[i, (i + 1)] <- adjacency_matrix[(i + 1), i] <- 1
#min_links <- graph_size - 1
#max_links <- (((graph_size ^ 2) - graph_size) / 2) - 1
#adjacency_matrices <- list(adjacency_matrix)
#if (max_links > min_links) {
#  for(i in 1:(max_links - min_links)) {
#    sampled_coordinates <- sample(x = 1:graph_size, size = 2, replace = TRUE)
#    while(length(x = unique(x = sampled_coordinates)) == 1 || adjacency_matrices[[length(adjacency_matrices)]][sampled_coordinates[1], sampled_coordinates[2]] == 1) sampled_coordinates <- sample(x = 1:graph_size, size = 2, replace = TRUE)
#    adjacency_matrix <- adjacency_matrices[[length(adjacency_matrices)]]
#    adjacency_matrix[sampled_coordinates[1], sampled_coordinates[2]] <- adjacency_matrix[sampled_coordinates[2], sampled_coordinates[1]] <- 1
#    adjacency_matrices[[(length(x = adjacency_matrices) + 1)]] <- adjacency_matrix
#  }
#}
#type_v_costmatrices <- lapply(X = adjacency_matrices, FUN = Claddis::convert_adjacency_matrix_to_costmatrix)
#permuted_tipstates <- permute_tipstates(
#  tree = ape::read.tree(text = paste("(", paste(make_labels(tip_count), collapse = ","), ");", sep = "")),
#  states = as.character(x = 0:(graph_size - 1)),
#  all_states_present = TRUE
#)
#exhaustive_gmax <- apply(
#  X = apply(
#    X = permuted_tipstates,
#    MARGIN = 1,
#    FUN = function(i) {
#      unlist(
#        x = lapply(
#          X = type_v_costmatrices,
#          FUN = function(j) {
#            calculate_g(
#              costmatrix = j,
#              tip_states = i
#            )
#          }
#        )
#      )
#    }
#  ),
#  MARGIN = 1,
#  FUN = max
#)
# Assumes starting path from 0 to n is longest and assigns states equally to ends (DOES NOT WORK):
#integer_gmax_1 <- unlist(
#  x = lapply(
#    X = type_v_costmatrices,
#    FUN = function(i) calculate_g(costmatrix = i, tip_states = c(rep("0", length.out = floor(x = (tip_count - graph_size) / 2)), as.character(x = 0:(graph_size - 1)), rep(as.character(x = (graph_size - 1)), length.out = ceiling(x = (tip_count - graph_size) / 2))))
#  )
#)
# Assumes can dsitribute spare states to peripheral vertices (DOES NOT WORK):
#integer_gmax_2 <- unlist(
#  x = lapply(
#    X = type_v_costmatrices,
#    FUN = function(i) {
#      vertices <- convert_costmatrix_to_stategraph(i)$vertices
#      peripheral_vertices <- vertices[vertices[, "periphery"] == 1, "label"]
#      calculate_g(costmatrix = i, tip_states = tip_states <- sort(x = c(as.character(x = 0:(graph_size - 1)), peripheral_vertices[((0:(tip_count - graph_size - 1)) %% length(x = peripheral_vertices)) + 1])))
#    }
#  )
#)
# TRY ACTUALLY FINDING A LONGEST PATH AND THEN ASSIGNING SPARE STATES EQUALLY TO ENDS
#cbind(exhaustive_gmax, integer_gmax_1, integer_gmax_2)
