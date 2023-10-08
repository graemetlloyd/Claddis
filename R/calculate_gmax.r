#' Calculate the maximum possible tree length, gmax, under parsimony
#'
#' @description
#'
#' Given a costmatrix and number of tips, calculates the longest possible tree length under maximum parsimony.
#'
#' @param costmatrix An object of class \code{costMatrix}.
#' @param n_taxa The number of taxa (ree tips) to consider.
#' @param allow_zeroes Logical indicating whether some states are allowed to be assigned zero tips. Defaults to \code{TRUE}. Note that if the number of states exceeds the number of tips then \code{allow_zeroes = TRUE} is forced. Additionally, this rstriction only applies to integer gmax.
#'
#' @details
#'
#' The concept of \emph{gmax} was introduced by Hoyal Cuthill et al. (2010) and in the strictest sense represents the integer state frequency that maximises \emph{g} given every state is sampled at least once. This function is intended to relax that slightly using the \code{allow_zeroes} option and indeed this is not possible in rare cases such as either high levels of missing data or gap-weighted (Thiele 1993) characters where it is possible for the number of states to exceed the (effective) number of tips. In practice the function also estimates an algebraic version of \emph{gmax} which allows for fractional or negative state frequencies (Hoyal Cuthill and Lloyd in revision).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Jen Hoyal Cuthill \email{j.hoyal-cuthill@@essex.ac.uk}
#'
#' @references
#'
#' Hoyal Cuthill, J. F., Braddy, S. J. and Donoghue, P. C. J., 2010. A formula for maximum possible steps in multistate characters: isolating matrix parameter effects on measures of evolutionary convergence. \emph{Cladistics}, \bold{26}, 98-102.
#'
#' Thiele, K.. 1993. The Holy Grail of the perfect character: the cladistic treatment of morphometric data. \emph{Cladistics}, \bold{9}, 275-304.
#'
#' @return
#'
#' A list containings two items:
#'
#' \item{algebraic_gmax}{The absolute limit gmax can reach if not restricted to positive or integer state frequencies.}
#' \item{integer_gmax}{The value gmax can reach when restfricted to integer values and whatever was chosen for the \code{allow_zeroes} option.}
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
#'   max_state= 2,
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
#'   max_state= 2,
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
#'   max_state= 1,
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
#'   max_state= 2,
#'   character_type = "irreversible"
#' )
#'
#' # Calculate gmax for ten tips:
#' calculate_gmax(
#'   costmatrix = multistate_irreversible_costmatrix,
#'   n_taxa = 10
#' )
#'
# Create a Type VIII character costmatrix:
#' binary_dollo_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 1,
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
#'   max_state= 2,
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
#'   max_state= 5,
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
#'   max_state= 1,
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
#'   max_state= 2,
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
calculate_gmax <- function(costmatrix, n_taxa, allow_zeroes = TRUE) {
  
  ### IGNORES POLYMORPHISMS AND UNCERTAINTIES FOR TIPS, BUT COULD PERHAPS ALLOW FORMER FOR ANCESTRAL STATE?
  ### CHECK N TAXA IS A SENSIBLE NUMBER!
  ### CHECK FOR COSTMATRIX OF TYPE I (OUTPUT IMMEDIATELY IS ZERO!)
  ### CURRENTLY CUSTOM INTEGER GMAX DOES NOT DEAL WITH ALLOW ZEROES = FALSE CASE
  
  # Set logical value for Dollo status (used later):
  costmatrix_is_dollo <- ifelse(
    test = costmatrix$type == "dollo",
    yes = TRUE,
    no = FALSE
  )
  
  # Special case of an ordered character:
  if (costmatrix$type == "ordered") {
      
    # If can safely assign at least one taxon to each state:
    if (!allow_zeroes && n_taxa >= costmatrix$n_states) {
        
      # Assign one taxon to each state to initialise state frequencies:
      state_frequencies <- rep(x = 1, times = costmatrix$n_states)
        
    # If cannot assign a taxon to each state:
    } else {
        
      # Initialise state frequencies as all zeros:
      state_frequencies <- rep(x = 0, times = costmatrix$n_states)
    }
      
    # Find number of tips left to assign to a state:
    tips_to_assign <- n_taxa - sum(x = state_frequencies)
      
    # As long as there are still tips to assign:
    if (tips_to_assign > 0) {
        
      # Assign half (rounding up) of taxa to first state:
      state_frequencies[1] <- state_frequencies[1] + ceiling(x = tips_to_assign / 2)
        
      # Assign half (rounding down) of taxa to last state:
      state_frequencies[costmatrix$n_states] <- state_frequencies[costmatrix$n_states]  + floor(x = tips_to_assign / 2)
    }
    
    # Special case of ordered algebraic gmax:
    algebraic_gmax <- (max(x = costmatrix$costmatrix[1, costmatrix$single_states]) / 2) * n_taxa
  }
    
  # Special case of an unordered character:
  if (costmatrix$type == "unordered") {
      
    # Initialise state frequencies with largest integer that can be assigned to each state:
    state_frequencies <- rep(x = floor(x = n_taxa / costmatrix$n_states), times = costmatrix$n_states)
      
    # If there are still some taxa left without a state assignment:
    if (sum(x = state_frequencies) < n_taxa) {
        
      # Add 1 to each of the first N state frequencies where N is the number of taxa without a state assignment:
      state_frequencies[1:(n_taxa - sum(x = state_frequencies))] <- state_frequencies[1:(n_taxa - sum(x = state_frequencies))] + 1
    }
    
    # Special case of algenraic gmax for unordered characters:
    algebraic_gmax <- sum(x = costmatrix$costmatrix[1, costmatrix$single_states] * rep(x = n_taxa / costmatrix$n_states, times = costmatrix$n_states))
  }
    
  # Special case of a Dollo character:
  if (costmatrix_is_dollo) {
    
    # Convert costmatrix into a linear ordered character to maintain correct costs below:
    new_costmatrix <- make_costmatrix(
      min_state = 0,
      max_state = costmatrix$n_states - 1,
      character_type = "ordered",
      include_polymorphisms = costmatrix$includes_polymorphisms,
      include_uncertainties = costmatrix$includes_uncertainties,
      polymorphism_costs = costmatrix$polymorphism_costs,
      polymorphism_geometry = costmatrix$polymorphism_geometry,
      polymorphism_distance = costmatrix$polymorphism_distance,
      dollo_penalty = costmatrix$dollo_penalty
    )
    
    # Update weight term:
    new_costmatrix$weight <- costmatrix$weight
    
    # Overwrite Dollo costmatrix with new costmatrix:
    costmatrix <- new_costmatrix

    # If can safely assign at least one taxon to each state:
    if (!allow_zeroes && n_taxa >= costmatrix$n_states) {
        
      # Assign one taxon to each state to initialise state frequencies:
      state_frequencies <- rep(x = 1, times = costmatrix$n_states)
        
    # If cannot assign a taxon to each state:
    } else {
        
      # Initialise state frequencies by assigning one taxon to the two terminal states:
      state_frequencies <- c(1, rep(x = 0, times = costmatrix$n_states - 2), 1)
    }
      
    # Find number of tips left to assign to a state:
    tips_to_assign <- n_taxa - sum(x = state_frequencies)
    
    # If allowing zeroes and still tips to assign:
    if (allow_zeroes && tips_to_assign > 0) {
      
      # Make most derived ttate 2 to force as root (maximising length on star tre):
      state_frequencies[costmatrix$n_states] <- 2
      
      # Lower tips to assign by one:
      tips_to_assign <- tips_to_assign - 1
    }
    
    # As long as there are still tips to assign assign them to the lowest state (maximising losses):
    if (tips_to_assign > 0) state_frequencies[1] <- state_frequencies[1] + tips_to_assign
    
    # Set Dollo star root as second most derived individual tip state (forces no multiple acquisitions):
    dollo_star_root <- unlist(
      x = lapply(
        X = as.list(x = 1:costmatrix$n_states),
        FUN = function(state) rep(x = costmatrix$single_states[state], times = state_frequencies[state]))
    )[(n_taxa - 1)]
    
    # Special case of algebraic gmax for Dollo (at limit there are 1 plus 1/Inf of derived state and hence n - 1 losses to least derived state - with any intermediate states not sampled):
    algebraic_gmax <- (n_taxa - 1) * costmatrix$costmatrix["0", costmatrix$n_states]
  }
    
  # Special case of an irreversible or stratigraphic character:
  if (costmatrix$type == "irreversible" || costmatrix$type == "stratigraphy") {
      
    # If can safely assign at least one taxon to each state:
    if (!allow_zeroes && n_taxa >= costmatrix$n_states) {
        
      # Assign one taxon to each state to initialise state frequencies:
      state_frequencies <- rep(x = 1, times = costmatrix$n_states)
        
    # If cannot assign a taxon to each state:
    } else {
        
      # Initialise state frequencies by assigning one taxon to the two terminal states:
      state_frequencies <- c(1, rep(x = 0, times = costmatrix$n_states - 2), 1)
    }

    # Find number of tips left to assign to a state:
    tips_to_assign <- n_taxa - sum(x = state_frequencies)
      
    # As long as there are still tips to assign assign them to the smallest state (most reversals):
    if (tips_to_assign > 0) state_frequencies[costmatrix$n_states] <- state_frequencies[1] + tips_to_assign
    
    # Special case of irreversible or stratigraphic algebraic gmax (frequency of most derived state tends to n_taxa but least derived state still considered "present"):
    algebraic_gmax <- costmatrix$costmatrix[1, costmatrix$n_states] * n_taxa
  }
  
  # If the general case:
  if (costmatrix$type == "custom") {
    
    # Subfunction to traverse simplex of all possible state frequencies to find integer gmax:
    traverse_simplex_gmax <- function(costmatrix, n_taxa) {
      
      ### HOW TO MAKE THIS WORK FOR THE ALL STATES MUST BE PRESENT AT LEAST ONCE CASE? CAN WE JUST EXCLUDE THE "PERIMETER" LAYER SOMEHOW?
      
      # Little fix for infinite characters (just replace with a very large value as this will never be part of g anyway):
      if (any(costmatrix$costmatrix == Inf)) costmatrix$costmatrix[costmatrix$costmatrix == Inf] <- max(x = costmatrix$costmatrix[costmatrix$costmatrix != Inf]) * 9999999
      
      # Starting state frequency:
      state_frequency <- c(n_taxa, rep(x = 0, times = costmatrix$n_states - 1))
      names(x = state_frequency) <- costmatrix$single_states
      
      # Build state frequency matrix:
      state_frequency_matrix <- matrix(
        data = rep(x = state_frequency, times = costmatrix$n_states),
        nrow = costmatrix$n_states,
        ncol = costmatrix$n_states,
        byrow = TRUE,
        dimnames = list(costmatrix$single_states, costmatrix$single_states)
      )
      
      # Get costs of each possible root:
      root_costs <- apply(
        X = state_frequency_matrix * costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states],
        MARGIN = 1,
        FUN = sum
      )

      # Get minimum star tree cost (g):
      g = min(x = root_costs)
      
      # Populate start point (all tips assigned to first state):
      start_point <- list(
        step = 0,
        state_frequency = state_frequency,
        g = g,
        root_state = names(x = which(x = root_costs == g))
      )
      
      # Create empty list to store path taken in gmax search:
      path <- list()
      
      # Add start point as first item in path:
      path[[1]] <- start_point
      
      # Create empty vector to store previous roots:
      previous_roots <- vector(mode = "character")
      
      # Set current root as first state:
      current_root <- start_point$root_state
      
      # Set starting next root as second state:
      next_root <- costmatrix$single_states[2]
      
      # Until all roots have been discovered:
      while(!is.na(x = next_root)) {
        
        # Generate candidate state frequencies:
        candidate_state_frequencies <- do.call(
          what = rbind,
          args = lapply(
            X = as.list(x = c(previous_roots, current_root)),
            FUN = function(old_state) {
              new_state_frequency <- state_frequency
              if (state_frequency[old_state] > 0) {
                new_state_frequency[old_state] <- new_state_frequency[old_state] - 1
                new_state_frequency[next_root] <- new_state_frequency[next_root] + 1
                return(value = new_state_frequency)
              } else {
                return(value = matrix(ncol = costmatrix$n_states, nrow = 0, dimnames = list(c(), costmatrix$single_states)))
              }
            }
          )
        )

        # If there is more than one candidate state freqeuncy (first look for maximum roots):
        if (nrow(x = candidate_state_frequencies) > 1) {
          
          # Get some stats on each candidate state frequency:
          candidate_stats <- apply(
            X = candidate_state_frequencies,
            MARGIN = 1,
            FUN = function(candidate_state_frequency) {
              
              # Build state frequency matrix:
              state_frequency_matrix <- matrix(
                data = rep(x = candidate_state_frequency, times = costmatrix$n_states),
                nrow = costmatrix$n_states,
                ncol = costmatrix$n_states,
                byrow = TRUE,
                dimnames = list(costmatrix$single_states, costmatrix$single_states)
              )
              
              # Get costs of each possible root:
              root_costs <- apply(
                X = state_frequency_matrix * costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states],
                MARGIN = 1,
                FUN = sum
              )
              
              # Get g:
              g <- min(x = root_costs)
              
              # Get root state(s):
              root_state <- names(x = which(x = root_costs == g))
              
              # Return g, root states:
              list(
                state_frequency = candidate_state_frequency,
                g = g,
                root_state = root_state,
                n_roots = length(root_state)
              )
            }
          )
          
          # Get total number of roots for each candidate state frequency:
          n_roots <- unlist(x = lapply(X = candidate_stats, FUN = function(x) x$n_roots))
          
          # Get maximum number of roots:
          max_n_roots <- max(x = n_roots)
          
          # Prune out lower root counts:
          candidate_state_frequencies <- candidate_state_frequencies[n_roots == max_n_roots, , drop = FALSE]
          candidate_stats <- candidate_stats[n_roots == max_n_roots]
        }

        # If there is still more than one candidate state freqeuncy (choose most different root(s)):
        if (nrow(x = candidate_state_frequencies) > 1) {

          # Check last root in path:
          last_root <- path[[length(x = path)]]$root_state
            
          # Identify number of diffreent roots to last root:
          different_roots <- unlist(x = lapply(X = candidate_stats, FUN = function(x) length(x = setdiff(x = x$root_state, y = last_root))))
          
          # Get maximum different roots value:
          max_different_roots <- max(x = different_roots)
          
          # Prune out lower different root counts:
          candidate_state_frequencies <- candidate_state_frequencies[different_roots == max_different_roots, , drop = FALSE]
          candidate_stats <- candidate_stats[different_roots == max_different_roots]
        }

        # If there is still more than one candidate state frequency (choose largest g value):
        if (nrow(x = candidate_state_frequencies) > 1) {
          
          # Get all g values:
          g_values <- unlist(x = lapply(X = candidate_stats, FUN = function(x) x$g))

          # Get maximum g value:
          max_g_value <- max(x = g_values)
          
          # Prune out lower g values:
          candidate_state_frequencies <- candidate_state_frequencies[g_values == max_g_value, , drop = FALSE]
          candidate_stats <- candidate_stats[g_values == max_g_value]
        }
        
        # Should now either be a single caniddtae state frequency (or no further reason to choose one over another) so take first:
        state_frequency <- candidate_state_frequencies[1, ]

        # Build state frequency matrix:
        state_frequency_matrix <- matrix(
          data = rep(x = state_frequency, times = costmatrix$n_states),
          nrow = costmatrix$n_states,
          ncol = costmatrix$n_states,
          byrow = TRUE,
          dimnames = list(costmatrix$single_states, costmatrix$single_states)
        )
        
        # Get costs of each possible root:
        root_costs <- apply(
          X = state_frequency_matrix * costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states],
          MARGIN = 1,
          FUN = sum
        )

        # Get minimum star tree cost (g):
        g = min(x = root_costs)
        
        # Add new state frequency to path:
        path[[(length(x = path) + 1)]] <- list(
          step = length(x = path),
          state_frequency = state_frequency,
          g = g,
          root_state = names(x = which(x = root_costs == g))
        )
        
        # Check for next_root and update if found:
        if (any(path[[length(x = path)]]$root_state == next_root)) {
          previous_roots <- c(previous_roots, current_root)
          current_root <- next_root
          next_root <- costmatrix$single_states[which(x = costmatrix$single_states == current_root) + 1]
        }
      }
      
      # Get final state frequency (last step in path):
      final_state_frequency <- path[[length(x = path)]]$state_frequency
      
      # Make final state frequencies to check:
      final_state_frequencies <- do.call(
        what = rbind,
        args = lapply(
          X = as.list(x = names(x = which(x = final_state_frequency > 0))),
          FUN = function(state) {
            do.call(
              what = rbind,
              args = lapply(
                X = as.list(x = setdiff(x = costmatrix$single_states, y = state)),
                FUN = function(other_state) {
                  final_state_frequency[state] <- final_state_frequency[state] - 1
                  final_state_frequency[other_state] <- final_state_frequency[other_state] + 1
                  final_state_frequency
                }
              )
            )
          }
        )
      )
      
      # Make final path additions by checking area around final step of path:
      final_path_additions <- apply(
        X = final_state_frequencies,
        MARGIN = 1,
        FUN = function(state_frequency) {
          
          # Build state frequency matrix:
          state_frequency_matrix <- matrix(
            data = rep(x = state_frequency, times = costmatrix$n_states),
            nrow = costmatrix$n_states,
            ncol = costmatrix$n_states,
            byrow = TRUE,
            dimnames = list(costmatrix$single_states, costmatrix$single_states)
          )

          # Get costs of each possible root:
          root_costs <- apply(
            X = state_frequency_matrix * costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states],
            MARGIN = 1,
            FUN = sum
          )
          
          # Get g:
          g <- min(x = root_costs)
          
          # Output data in path format:
          list(
            step = NA,
            state_frequency = state_frequency,
            g = g,
            root_state = costmatrix$single_states[root_costs == g]
          )
        },
        simplify = FALSE
      )
      
      # Add to path object:
      path <- c(path, final_path_additions)
      
      # Isolate g values:
      g_values <- unlist(x = lapply(X = path, FUN = function(x) x$g))
      
      # Calculate gmax:
      gmax <- max(x = g_values)
      
      # return first state frequecny with g maximized:
      path[[which(x = g_values == gmax)[1]]]$state_frequency
    }
    
    # Define a subfunction to get state frequencies for algebraic gmax:
    algebraic_gmax_state_frequencies <- function(costmatrix, n_taxa) {
      
      # Little fix for infinite characters (just replace with a very large value as this will never be part of g anyway):
      if (any(costmatrix$costmatrix == Inf)) costmatrix$costmatrix[costmatrix$costmatrix == Inf] <- max(x = costmatrix$costmatrix[costmatrix$costmatrix != Inf]) * 9999999
      
      # Prepare zero value matrices to fill:
      coefficient_matrix <- matrix(
        data = 0,
        nrow = costmatrix$n_states,
        ncol = costmatrix$n_states
      )
      constant_matrix <- matrix(
        data = 0,
        nrow = costmatrix$n_states,
        ncol = 1
      )
      
      # Enter frequency coefficients of 1 given that state frequencies sum to the number of taxa:
      coefficient_matrix[1, ] <- 1

      # Enter constant t given that state frequencies sum to the number of taxa:
      constant_matrix[1] <- n_taxa

      # Enter remaining n - 1 rows of frequency coefficients based on subtraction of rows in the cost matrix:
      coefficient_matrix[2:costmatrix$n_states, ] <- do.call(
        what = rbind,
        args = apply(
          X = costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states][2:costmatrix$n_states, , drop = FALSE],
          MARGIN = 1,
          FUN = function(row) costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states][1, ] - row,
          simplify = FALSE
        )
      )
        
      # Get frequency of states:
      state_frequency_solution <- base::solve(a = coefficient_matrix, b = constant_matrix, tol = 1.85037e-18)
      rownames(x = state_frequency_solution) <- costmatrix$single_states
      
      # Return algebraic gmax state frequencies:
      state_frequency_solution[, 1]
    }
    
    # Set state frequencies for integer gmax:
    state_frequencies <- traverse_simplex_gmax(costmatrix = costmatrix, n_taxa = n_taxa)
    
    # Set satte frequencies for algebraic gmax:
    algebraic_state_frequencies <- algebraic_gmax_state_frequencies(costmatrix = costmatrix, n_taxa = n_taxa)

    # Calculate algebraic gmax:
    algebraic_gmax <- min(
      x = apply(
        X = matrix(
          data = rep(x = algebraic_state_frequencies, times = costmatrix$n_states),
          nrow = costmatrix$n_states,
          byrow = TRUE
        ) * costmatrix$costmatrix,
        MARGIN = 1,
        FUN = sum
      )
    )
  }
  
  # Add state names to state frequencies:
  names(x = state_frequencies) <- costmatrix$single_states
    
  # Build state frequencies into a matrix ready to calculate g_max:
  frequency_matrix <- matrix(
    data = rep(x = state_frequencies, times = costmatrix$n_states),
    nrow = costmatrix$n_states,
    ncol = costmatrix$n_states,
    byrow = TRUE,
    dimnames = list(names(x = state_frequencies), names(x = state_frequencies))
  )
    
  # Get costs for each possible ancestral state:
  costs_for_ancestral_state <- apply(
    X = frequency_matrix * costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states],
    MARGIN = 1,
    FUN = function(x) sum(x = as.numeric(x = gsub(pattern = NaN, replacement = 0, x = x)))
  )
  
  # Special case of a Dollo character:
  if (costmatrix_is_dollo) {

    # Calculate integer gmax for Dollo character (uses forced root and state frequencies):
    integer_gmax <- sum(x = costmatrix$costmatrix[dollo_star_root, ] * state_frequencies)
    
  # If any other type of character:
  } else {
    
    # Calculate integer g_max for non-Dollo character:
    integer_gmax <- min(x = costs_for_ancestral_state)
  }
  
  # Return gmax values:
  list(
    algebraic_gmax = algebraic_gmax,
    integer_gmax = integer_gmax
  )
}
