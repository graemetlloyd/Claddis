#' XXX
#'
#' @description
#'
#' XXX.
#'
#' @param N The number of labels required,
#'
#' @details
#'
#' XXX.
#'
#' @return XXX.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # XXX:
#' equalize_weights()
#'
#' @export equalize_weights
equalize_weights <- function(cladistic_matrix, method) {

  # method being for disparity or phylogenetics
  # Need to calculate what weights would work TNT and not tht not all of them need to be integers if they would not cause FPEs, e.g., 0.5 is fine
  
  if (method == "disparity") {
    
    # Pull this bit from read_nexus_matrix and then re-add this function there with method call
  }
  if (method == "cladistic") {
    
    # Take cladistic_matrix and find m and g for each character
    # Then calculate ranges as g - m for each character
    
    ### EXAMPLE
    ranges <- c(3, 1, 1, 2, 8, 1, 1, 0, 0, 5, 9, 2, 5, 1, 1, 3, 4, 1, 3, 4, 5, 2, 8, 6, 1, 0, 8, 2, 4, 2, 1, 6, 3, 3, 1, 4, 0, 1, 7, 5, 9, 3, 3, 0, 5, 6, 4, 4, 2, 9, 7, 10, 5, 7, 10, 3, 0, 10, 4, 12, 4, 2, 2, 3, 2, 2, 5, 4, 5, 0, 3)
    
    # Get just unique ranges:
    unique_ranges <- sort(x = unique(x = ranges))
    
    # If there are any zero ranges then remove these from calculation:
    if (unique_ranges[1] == 0) unique_ranges<- unique_ranges[-1]
    
    # Set common factor not found to TRUE:
    common_factor_not_found <- TRUE
    
    # Initialise base weight (will not actually use zero!):
    base_weight <- 0
    
    # While a common factor has not been found:
    while(common_factor_not_found) {
      
      # Increment base weight:
      base_weight <- base_weight + 1
      
      # Check if current base_weight is a factor of every range:
      if (all(x = base_weight %% unique_ranges == 0)) {
        
        # Common factor found! We all know what that means:
        common_factor_not_found <- FALSE
      }
    }
    
    # Can now define unique (integer) weights:
    unique_weights <- base_weight / unique_ranges
    
    unique_weights
    
    # Return cladistic_matrix object with updated character_weights
    
    # Return updated cladistic matrix:
    return(cladistic_matrix)
  }
}
