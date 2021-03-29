#' Make a stepmatrix for a given set of states
#'
#' @description
#'
#' Given a set of discrete states and a character type will make the approriate stepmatrix for them.
#'
#' @param min_state The minimum character state (defaults to \code{0}).
#' @param max_state The maximum character state. Must be \code{1} or greater.
#' @param character_type The type of character desired. Must be one of: \code{"ordered"}, \code{"unordered"}, \code{"dollo"}, \code{"irreversible"}, or \code{"stratigraphy"}.
#' @param include_polymorphisms Logical indicating whether or not to include polymorphic state combinations (defaults to \code{FALSE}).
#' @param polymorphism_shape The shape to use for assigning polymorphism coordinates. Must be one of: \code{"hypercube"}, \code{"hypersphere"}, or \code{"simplex"}. (Only relevant if \code{character_type = "unordered"} and \code{include_polymorphisms = TRUE}.)
#' @param polymorphism_distance The distance to use to set transformation costs between polymorphic states. Must be one of: \code{"euclidean"}, \code{"great_circle"}, or \code{"manhattan"}.
#' @param state_ages A vector of ages assigned to each state (only relevant if \code{character_type = stratigraphy}).
#' @param dollo_penalty The size of the cost penalty for the acquisition of a Dollo character (defaults to \code{999}). Note: this should always be a positive real value greater than one, and never infinity (\code{Inf}), as at least one acquisition is expected.
#'
#' @details
#'
#' Text.
#'
#' \bold{Character types}
#' \emph{ordered} Wagner States must be in order of transistion, e.g., 0, 1, 2 means 0<->1<->2
#' \emph{unordered} Fitch
#' \emph{dollo} States must be in order.
#' \emph{irreversible} States must be in order. Aka. Camin-Sokal.
#' \emph{stratigraphy} States must be in order from oldest to youngest. In usage zero must also be root state or else it makes no damn sense. Special case of irreversible and needs additional input of ages for differences!
#'
#' \bold{Polymorphic characters}
#'
#' Depends on character type too. E.g., for stratigraphy polymorphisms uses oldest state - leads to scenario where off-diagonal values can also be zero which may be important if looking for most parsimonious nacestral state reconstructions later. Not sure about Dollo and Irreversible?
#'
#' \emph{Shape of polymorphic state space}
#'
#' "hypercube", "hypersphere", "simplex"
#'
#' \emph{Distances in polymorphic state space}
#'
#' "manhattan", "euclidean", "great_circle" (assumes radius of one as does coordinates of shape)
#'
#' \bold{Character weights}
#'
#' Tempting to multiply space by weight to get weighted form, but no need as dealt with in parsimony function and can be dangerous with some values (e.g., weights of zero - may still want relative weights for ancestral state estimation even if not "charging" cost of transitions). Tranistions are mostly recslaed to a weight of one (except stratigraphy?).
#'
#' \bold{More complex stepmatrices}
#'
#' Can customise by mkaing one here as a starting point then manually modifying it. [CITE HOOKER]
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' REF
#'
#' @return
#'
#' A square stepmatrix with rows representing "from" states, columns "to" states, and individual cells the cost in steps of that transition.
#'
#' @seealso
#'
#' \link{make_all_polymorphisms}
#'
#' @examples
#'
#' # Make an unordered stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make an unordered stepmatrix including polymorphisms:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered",
#'   include_polymorphisms = TRUE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make an ordered stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "ordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make an ordered stepmatrix including polymorphisms:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "ordered",
#'   include_polymorphisms = TRUE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make a Dollo stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "dollo",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle",
#'   dollo_penalty = 100
#' )
#'
#' # Make an irreversible stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "irreversible",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Make a stratigraphic stepmatrix:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "stratigraphy",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle",
#'   state_ages = c(52, 34, 12)
#' )
#'
#' # Make a stratigraphic stepmatrix with polymorphisms:
#' make_stepmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "stratigraphy",
#'   include_polymorphisms = TRUE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle",
#'   state_ages = c(52, 34, 12)
#' )
#'
#' @export make_stepmatrix
make_stepmatrix <- function(min_state = 0, max_state, character_type, include_polymorphisms = FALSE, polymorphism_shape, polymorphism_distance, state_ages, dollo_penalty = 999) {
  
  # - Need more incoming checks probably?
  # - Option to disallow polymorphisms as ancestors (still allows costs for transitions to polymorphic states to be taken into account)
  
  # Check character_type is a valid value and stop and warn user if not:
  if(length(x = setdiff(x = character_type, y = c("ordered", "unordered", "dollo", "irreversible", "stratigraphy"))) > 0) stop("character_type must be one of \"ordered\", \"unordered\", \"dollo\", \"irreversible\", or \"stratigraphy\".")
  
  # Check polymorphism_shape is a valid value and stop and warn user if not:
  if(length(x = setdiff(x = polymorphism_shape, y = c("hypercube", "hypersphere", "simplex"))) > 0) stop("character_type must be one of \"hypercube\", \"hypersphere\", \"simplex\".")
  
  # Check polymorphism_distance is a valid value and stop and warn user if not:
  if(length(x = setdiff(x = polymorphism_distance, y = c("manhattan", "euclidean", "great_circle"))) > 0) stop("character_type must be one of \"manhattan\", \"euclidean\", \"great_circle\".")
  
  # Get single states:
  single_states <- min_state:max_state
  
  # Case if character is ordered:
  if(character_type == "ordered") {
    
    # If including polymorphisms:
    if(include_polymorphisms) {
      
      # Generate and store all possible states:
      all_states <- make_all_polymorphisms(single_states = single_states)
      
      # Check data are not too big (>= 2^14 states) and stop and warn user if so:
      if(length(x = all_states) >= 16384) stop("Stepmatrix would be too large. Use fewer states.")
      
      # Calculate stepmatrix size:
      stepmatrix_size <- length(x = all_states)
      
      # Initialise adjacency atrix with zeroes:
      adjacency_matrix <- matrix(data = 0, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(all_states, all_states))
      
      # For each state save the last one:
      for(state in all_states[-length(x = all_states)]) {
        
        # Find matching states for current state:
        state_matches <- apply(X = do.call(what = rbind, args = lapply(X = as.list(x = strsplit(x = state, split = "&")[[1]]), FUN = function(x) single_states == x)), MARGIN = 2, FUN = any)
        
        # Find states that are addable (reachable) from current state(s):
        addable_states <- single_states[apply(X = rbind(!state_matches, unlist(lapply(X = as.list(x = 1:length(x = state_matches)), FUN = function(x) any(x = c(state_matches[x - 1], state_matches[x + 1]), na.rm = TRUE)))), MARGIN = 2, FUN = all)]
        
        # Compile state(s) adjacent to current state:
        adjacent_states <- unlist(x = lapply(X = as.list(x = addable_states), FUN = function(new_state) paste(x = sort(x = c(new_state, strsplit(x = state, split = "&")[[1]])), collapse = "&")))
        
        # Populate adjaceney matrix with these states:
        adjacency_matrix[state, adjacent_states] <- adjacency_matrix[adjacent_states, state] <- 1
        
      }
      
      # Convert adjacency matrix to step matrix:
      stepmatrix <- convert_adjacency_matrix_to_stepmatrix(adjacency_matrix = adjacency_matrix)
      
      # Return stepmatrix rescaled such that adjacent singe states are one step:
      return(stepmatrix / stepmatrix[1, 2])
    
    # If excluding polymorphisms:
    } else {
      
      # Calculate stepmatrix size:
      stepmatrix_size <- length(x = single_states)
      
      # Initialise stepmatrix with all values set to zero:
      stepmatrix <- matrix(data = 0, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(single_states, single_states))
      
      # Store ordered values in upper and lower triangles:
      stepmatrix[upper.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(stepmatrix_size - 1)), FUN = function(steps) steps:1))
      stepmatrix[lower.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = (stepmatrix_size - 1):1), FUN = function(steps) 1:steps))
      
      # Return stepmatrix:
      return(stepmatrix)

    }
  
  }

  # Case if character is unordered:
  if(character_type == "unordered") {
    
    # If including polymorphisms:
    if(include_polymorphisms) {
      
      # Generate and store all possible states:
      all_states <- make_all_polymorphisms(single_states = single_states)
      
      # Check data are not too big (>= 2^14 states) and stop and warn user if so:
      if(length(x = all_states) >= 16384) stop("Stepmatrix would be too large. Use fewer states.")
      
      # Create coordinate matrix and initialise with zeroes:
      state_presence_matrix <- matrix(data = 0, nrow = length(x = single_states), ncol = length(x = all_states), dimnames = list(single_states, all_states))
      
      # If using the hypercube coordinate space assign coordinates accordingly:
      if(polymorphism_shape == "hypercube") for(i in 1:ncol(x = state_presence_matrix)) state_presence_matrix[strsplit(x = colnames(x = state_presence_matrix)[i], split = "&")[[1]], i] <- 1
      
      # If using the hypersphere or simplex coordinate space:
      if(polymorphism_shape == "hypersphere" || polymorphism_shape == "simplex") {
        
        # For each coding:
        for(i in 1:ncol(x = state_presence_matrix)) {
          
          # Isolate components of polymorphism:
          components <- strsplit(x = colnames(x = state_presence_matrix)[i], split = "&")[[1]]
          
          # If using the hypersphere coordinate space then apply coordinates accordingly (the square root of 1/N states in polymorphism on each axis state is present):
          if(polymorphism_shape == "hypersphere") state_presence_matrix[components, i] <- sqrt(x = 1 / length(x = components))
          
          # If using the simplex coordinate space then apply coordinates accordingly (1/N states in polymorphism on each axis state is present):
          if(polymorphism_shape == "simplex") state_presence_matrix[components, i] <- 1 / length(x = components)
          
        }
        
      }
      
      # If using a manhattan or euclidean distance, calculate distance directly from coordinate-space:
      if(polymorphism_distance == "euclidean" || polymorphism_distance == "manhattan") stepmatrix <- as.matrix(x = dist(x = t(x = state_presence_matrix), method = polymorphism_distance, diag = TRUE, upper = TRUE))
      
      # If using a great circle distance:
      if(polymorphism_distance == "great_circle") {
        
        # Start by calculating the euclidean distances between each point:
        stepmatrix <- as.matrix(x = dist(x = t(x = state_presence_matrix), method = "euclidean", diag = TRUE, upper = TRUE))
        
        # Treating distances as the chord length of a circle of radius one transform them to the arc length for the same:
        stepmatrix <- 2 * asin(x = stepmatrix / 2)
        
      }
      
      # Return a stepmatrix rescaled such that single state distances (e.g., 0 to 1) are one (i.e., a normal unordered character):
      return(stepmatrix / stepmatrix[1, 2])
      
    # If excluding polymorphisms:
    } else {
      
      # Initialise stepmatrix with all values set to one:
      stepmatrix <- matrix(data = 1, nrow = length(x = single_states), ncol = length(x = single_states), dimnames = list(single_states, single_states))
      
      # Make diagional of stepmatrix zero:
      diag(x = stepmatrix) <- 0
      
      # Return stepmatrix:
      return(stepmatrix)
      
    }
    
  }
  
  # Case if character is Dollo:
  if(character_type == "dollo") {
    
    # If including polymorphisms:
    if(include_polymorphisms) {
      
      # Generate and store all possible states:
      all_states <- make_all_polymorphisms(single_states = single_states)
      
      # Check data are not too big (>= 2^14 states) and stop and warn user if so:
      if(length(x = all_states) >= 16384) stop("Stepmatrix would be too large. Use fewer states.")

      ### TO FINISH
      
    # If excluding polymorphisms:
    } else {
      
      # Calculate stepmatrix size:
      stepmatrix_size <- length(x = single_states)
      
      # Initialise stepmatrix with all values set to zero:
      stepmatrix <- matrix(data = 0, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(single_states, single_states))
      
      # Set upper triangle as an ordered character times the dollo_penalty
      stepmatrix[upper.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(stepmatrix_size - 1)), FUN = function(steps) steps:1)) * dollo_penalty
      
      # Store regular ordered values in lower triangles::
      stepmatrix[lower.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = (stepmatrix_size - 1):1), FUN = function(steps) 1:steps))
      
      # Return stepmatrix:
      return(stepmatrix)
      
    }
    
  }
  
  # Case if character is irreversible:
  if(character_type == "irreversible") {
    
    # If including polymorphisms:
    if(include_polymorphisms) {
      
      # Generate and store all possible states:
      all_states <- make_all_polymorphisms(single_states = single_states)
      
      # Check data are not too big (>= 2^14 states) and stop and warn user if so:
      if(length(x = all_states) >= 16384) stop("Stepmatrix would be too large. Use fewer states.")

      ### TO FINISH
      
    # If excluding polymorphisms:
    } else {
      
      # Calculate stepmatrix size:
      stepmatrix_size <- length(x = single_states)
      
      # Initialise stepmatrix with all values set to zero:
      stepmatrix <- matrix(data = 0, nrow = stepmatrix_size, ncol = stepmatrix_size, dimnames = list(single_states, single_states))
      
      # Exclude losses by setting lower triangle values to infinity:
      stepmatrix[lower.tri(x = stepmatrix)] <- Inf
      
      # Set upper triangle as an ordered character:
      stepmatrix[upper.tri(x = stepmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(stepmatrix_size - 1)), FUN = function(steps) steps:1))
      
      # Return stepmatrix:
      return(stepmatrix)
      
    }
    
  }
  
  # Case if character is stratigraphy:
  if(character_type == "stratigraphy") {
    
    # Add state names to ages:
    names(x = state_ages) <- min_state:max_state
    
    # If including polymorphisms:
    if(include_polymorphisms) stop("If character_type is \"stratigraphy\" then include_polymorphisms cannot be TRUE. If the age of the OTU is uncertain then code it as such (use / instead of & between state). If the OTU is truly present in multiple units then only code it as the oldest one or break it up into multiple OTUs.")
    
    # Generate initial stepmatrix using temporal distances:
    stepmatrix <- as.matrix(x = dist(x = state_ages, diag = TRUE, upper = TRUE))
    
    # Set infinite cost to all reversals:
    stepmatrix[lower.tri(x = stepmatrix)] <- Inf
    
    # Return stepmatrix:
    return(stepmatrix)
    
  }
  
}
