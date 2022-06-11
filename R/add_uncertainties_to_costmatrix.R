#' Adds uncertainties to a costmatrix
#'
#' @description
#'
#' Adds uncertainty transitions to a costmatrix.
#'
#' @param costmatrix An object of class "costMatrix".
#' @param message Logical indicating whether (\code{TRUE}, the default) or not (\code{FALSE}) to provide messages about the process of the function.
#'
#' @details
#'
#' Under a parsimony algorithm uncertainties indicate that the observed state cannot be restricted to a single value but some finite set of values. In #NEXUS format these are indicated with curly braces, {}, and within Claddis by a slash. For example, state 0 or 1 would be 0/1.
#'
#' For most operations (calculating the length of a tree, estimating ancestral states) there is no need to formally add uncertainties to a costmatrix. Instead, the uncertainty can be established via setting up the tip states in the Swofford and Maddison (1992) approach. However, there are some usages within Claddis where it is useful to include them in the costmatrix itself, for example calculating minimum step counts for homoplasy metrics.
#'
#' Here, uncertainties are added to a costmatrix in a similar way to polymorphisms. I.e., all possible uncertainty combinations are calculated and then the associated transition cost calculated and inserted. Note that here only one approach for calculating this cost is applied - the minimu rule of Hoyal Cuthill and Lloyd (in prep.).
#'
#' Importantly costs of transitioning \emph{from} uncertainties are assigned infinite cost. There aree two main reasons for this. Firstly, under parsimony and the minimum rule allowing uncertainties as ancestral states would lead to a tree length of zero and an ancestral state for every internal node that includes every possible state being favoured, regardless of the states at the tips. In other words, it would have no analytical value. Secondly, assigning uncertainty of ancestral states is best done by generating all (single or polymorphic) most parsimonious reconstructions and summarising these.
#'
#' Note that in many practical cases an uncertainty involving all possible states - e.g., the uncertianty 0/1/2 for a three-state character - operates the same as a missing value (NA). However, it may still be useful in some circumstances. For example, to denote that the coding is definitively applicable (as opposed to inapplicable) or to properly record uncertainty rather than missing values for summary statistical purposes.
#'
#' Finally, it should be noted that uncertainties are the only allowable breaking of the costmatrix rule that all off-diagonal transition costs must be positive. I.e., costs of zero from a single or polymorphic state to an uncertainty state are permitted.
#'
#' @return
#'
#' An object of class "costMatrix".
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Hoyal Cuthill and Lloyd, in prep.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. In R. L. Mayden (ed.), \emph{Systematics, Historical Ecology, and North American Freshwater Fishes}. Stanford University Press, Stanford. pp187-223.
#'
#' @examples
#'
#' # Generate an example three-state unordered character costmatrix:
#' # Generate an example three-state unordered character costmatrix:
#' unordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered"
#' )
#'
#' # Generate an example three-state ordered character costmatrix:
#' ordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "ordered"
#' )
#'
#' # Generate an example three-state ordered character costmatrix with polymorphism included:
#' unordered_polymorphism_costmatrix <- list(
#'   size = 7,
#'   n_states = 3,
#'   single_states = c("0", "1", "2"),
#'   type = "unordered",
#'   costmatrix = matrix(
#'     data = c(
#'       0, 2, 2, 1, 1, 3, 2,
#'       2, 0, 2, 1, 3, 1, 2,
#'       2, 2, 0, 3, 1, 1, 2,
#'       1, 1, 3, 0, 2, 2, 1,
#'       1, 3, 1, 2, 0, 2, 1,
#'       3, 1, 1, 2, 2, 0, 1,
#'       2, 2, 2, 1, 1, 1, 0
#'     ),
#'     nrow = 7,
#'     byrow = TRUE,
#'     dimnames = list(
#'       c(as.character(x = 0:2), "0&1", "0&2", "1&2", "0&1&2"),
#'       c(as.character(x = 0:2), "0&1", "0&2", "1&2", "0&1&2")
#'     )
#'   ),
#'   symmetry = "Symmetric",
#'   includes_polymorphisms = TRUE,
#'   polymorphism_costs = "geometric",
#'   polymorphism_geometry = "hypercube",
#'   polymorphism_distance = "manhattan",
#'   includes_uncertainties = FALSE,
#'   pruned = FALSE,
#'   dollo_penalty = 999,
#'   base_age = 1,
#'   weight = 0.5
#' )
#' class(unordered_polymorphism_costmatrix) <- "costMatrix"
#'
#' # Add uncertainties to unordered costmatrix:
#' unordered_costmatrix_uncertainties <- add_uncertainties_to_costmatrix(
#'   costmatrix = unordered_costmatrix
#' )
#'
#' # Show unordered costmatrix with uncertainties included:
#' unordered_costmatrix_uncertainties$costmatrix
#'
#' # Add uncertainties to ordered costmatrix:
#' ordered_costmatrix_uncertainties <- add_uncertainties_to_costmatrix(
#'   costmatrix = ordered_costmatrix
#' )
#'
#' # Show ordered costmatrix with uncertainties included:
#' ordered_costmatrix_uncertainties$costmatrix
#'
#' # Add uncertainties to unordered costmatrix with polymorphisms:
#' unordered_polymorphism_costmatrix_uncertainties <- add_uncertainties_to_costmatrix(
#'   costmatrix = unordered_polymorphism_costmatrix
#' )
#'
#' # Show unordered polymorphism costmatrix with uncertainties included:
#' unordered_polymorphism_costmatrix_uncertainties$costmatrix
#'
#' @export add_uncertainties_to_costmatrix
add_uncertainties_to_costmatrix <- function(
  costmatrix,
  message = TRUE
) {
  
  ### NEED TO DECIDE WHAT CHARACTER TYPE IS FOR POLYMORPHIC/UNCERTAINTY MATRICES PLUS WHAT MATRIX SYMMETRY IS (MAYBE NEED TO CHECK WHAT THIS IS USED FOR ACROSS PACKAGE?)
  ### CHECK costmatrix is of class costmatrix! (Need to rewrite that function for new structure)
  
  # Check if costmatrix is size 1 (no possibility for uncertainties):
  if (costmatrix$size == 1) {
    
    # Message user about this if requested:
    if (message) print("Costmatrix is size 1 meaning no uncertianties are possible. Returning original costmatrix.")
    
    # Return original costmatrix:
    return(value = costmatrix)
  }

  # Check if uncertainties are already present:
  if (costmatrix$includes_uncertainties) {
    
    # Message user about this if requested:
    if (message) print("costmatrix already includes uncertainties. Returning original costmatrix.")
    
    # Return original costmatrix:
    return(value = costmatrix)
  }
  
  # Store original state labels (before adding polymorphisms to them):
  original_state_labels <- colnames(x = costmatrix$costmatrix)
  
  # Permute new uncertainty states:
  uncertainty_states <- permute_all_uncertainties(single_states = costmatrix$single_states)
  
  # Calculate number of new states:
  n_new_states <- length(x = uncertainty_states)
  
  # Check data are not too big (>= 2^14 states) and stop and warn user if so:
  if (n_new_states >= 16384) stop("Costmatrix would be too large. Use fewer states.")

  # Store original costmatrix for later use:
  original_costmatrix <- costmatrix$costmatrix
  
  # Update costmatrix size with new states:
  costmatrix$size <- costmatrix$size + n_new_states
  
  # Initialise new costmatrix with infinities:
  costmatrix$costmatrix <- matrix(
    data = Inf,
    nrow = costmatrix$size,
    ncol = costmatrix$size,
    dimnames = list(
      c(original_state_labels, uncertainty_states),
      c(original_state_labels, uncertainty_states)
    )
  )

  # Populate portion representing original costmatrix with original values:
  costmatrix$costmatrix[original_state_labels, original_state_labels] <- original_costmatrix
  
  # Set diagonal to zero (always true):
  diag(x = costmatrix$costmatrix) <- 0
  
  # Add uncertianty costs using minimum rule:
  costmatrix$costmatrix[costmatrix$single_states, uncertainty_states] <- do.call(
    what = cbind,
    args = lapply(
      X = as.list(x = uncertainty_states),
      FUN = function(uncertainty) {
        apply(
          X = costmatrix$costmatrix[costmatrix$single_states, strsplit(x = uncertainty, split = "/")[[1]]],
          MARGIN = 1,
          FUN = min
        )
      }
    )
  )
  
  # If polymorphisms were included (will need to calculate costs for polymorphism to uncertainty transitions):
  if (costmatrix$includes_polymorphisms) {
    
    # First establish which polymorphisms are present:
    polymorphic_states <- colnames(costmatrix$costmatrix)[grep(pattern = "&", x = colnames(costmatrix$costmatrix))]
    
    # Insert cost for minimum polymorphic to single state transition for each uncertainty:
    costmatrix$costmatrix[polymorphic_states, uncertainty_states] <- do.call(
      what = cbind,
      args = lapply(
        X = as.list(x = uncertainty_states),
        FUN = function(uncertainty) {
          uncertainty_components <- strsplit(x = uncertainty, split = "/")[[1]]
          apply(
            X = costmatrix$costmatrix[polymorphic_states, uncertainty_components],
            MARGIN = 1,
            FUN = min
          )
        }
      )
    )
    
    # Sort costmatrix so that polymorphic states come last (will match output from add polymorphisms function if uncertainties were included prior):
    costmatrix$costmatrix <- costmatrix$costmatrix[
      c(costmatrix$single_states, uncertainty_states, polymorphic_states),
      c(costmatrix$single_states, uncertainty_states, polymorphic_states)
    ]
  }
  
  # Update includes_uncertainties to TRUE:
  costmatrix$includes_uncertainties <- TRUE
  
  # Return costmatrix with uncertainties added:
  return(value = costmatrix)
}
