#' Convert a costmatrix to a minimal state graph
#'
#' @description
#'
#' Given a costmatrix, returns the smallest possible state graph (fewest arcs).
#'
#' @param costmatrix An object of class \code{costMatrix}.
#'
#' @details
#'
#' A costmatrix summarises all possible state-to-state transition costs and hence each entry can also be considered as an arc of a directed state graph. However, many of these arcs could be removed and a complete description of the graph still be provided. For example, the diagonal (any transition from a state to itself - a loop) can be removed, as can any arc with infinite cost (as this arc would never be traversed in practice). Finally, some arcs are redundant as indirect paths already represent the same cost.
#'
#' As an example, we can consider the linear ordered costmatrix:
#'
#' \preformatted{    -------------
#'     | 0 | 1 | 2 |
#' -----------------
#' | 0 | 0 | 1 | 2 |
#' -----------------
#' | 1 | 1 | 0 | 1 |
#' -----------------
#' | 2 | 2 | 1 | 0 |
#' -----------------}
#'
#' A maximum directed graph representation would thus be:
#'
#' \preformatted{----------------------
#' | from | to | weight |
#' ----------------------
#' |   0  | 0  |   0    |
#' ----------------------
#' |   0  | 1  |   1    |
#' ----------------------
#' |   0  | 2  |   2    |
#' ----------------------
#' |   1  | 0  |   1    |
#' ----------------------
#' |   1  | 1  |   0    |
#' ----------------------
#' |   1  | 2  |   1    |
#' ----------------------
#' |   2  | 0  |   2    |
#' ----------------------
#' |   2  | 1  |   1    |
#' ----------------------
#' |   2  | 2  |   0    |
#' ----------------------}
#'
#' But the following description is still complete, and minimal:
#'
#' \preformatted{----------------------
#' | from | to | weight |
#' ----------------------
#' |   0  | 1  |   1    |
#' ----------------------
#' |   1  | 0  |   1    |
#' ----------------------
#' |   1  | 2  |   1    |
#' ----------------------
#' |   2  | 1  |   1    |
#' ----------------------}
#'
#' This function effectively generates the latter (the minimal directed graph representation as a matrix of arcs).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' An object of class \code{stateGraph}.
#'
#' @seealso
#'
#' \link{convert_adjacency_matrix_to_costmatrix}
#'
#' @examples
#'
#' # Make a six-state unordered character costmatrix:
#' unordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 5,
#'   character_type = "unordered"
#' )
#'
#' # Find the minimal directed graph representation:
#' convert_costmatrix_to_stategraph(costmatrix = unordered_costmatrix)
#'
#' # Make a six-state ordered character costmatrix:
#' ordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 5,
#'   character_type = "ordered"
#' )
#'
#' # Find the minimal directed graph representation:
#' convert_costmatrix_to_stategraph(costmatrix = ordered_costmatrix)
#'
#' # Make a six-state stratigraphic character costmatrix:
#' stratigraphic_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 5,
#'   character_type = "stratigraphy",
#'   state_ages = c(103, 91.4, 78.2, 73.4, 66.0, 59.7)
#' )
#'
#' # Find the minimal directed graph representation:
#' convert_costmatrix_to_stategraph(costmatrix = stratigraphic_costmatrix)
#'
#' @export convert_costmatrix_to_stategraph
convert_costmatrix_to_stategraph <- function(costmatrix) {
  
  ### ADD SOME INPUT CHECKS!
  
  # Set sampled states:
  sampled_states <- rownames(x = costmatrix$costmatrix)
  
  # Create table of all possible edges:
  all_edges <- expand.grid(from = sampled_states, to = sampled_states)
  
  # Add weight to edges:
  all_edges <- data.frame(
    from = all_edges[, "from"],
    to = all_edges[, "to"],
    weight = apply(
      X = all_edges,
      MARGIN = 1,
      FUN = function(x) costmatrix$costmatrix[x[1], x[2]]
    )
  )
  
  # Prune loops (i.e., diagonal):
  all_edges <- all_edges[-which(x = all_edges[, "from"] == all_edges[, "to"]), , drop = FALSE]
  
  # Prune any infinite weight edges:
  if (any(all_edges[, "weight"] == Inf)) all_edges <- all_edges[-which(x = all_edges[, "weight"] == Inf), , drop = FALSE]
  
  # Make vector of redundant edges to remove (may be empty):
  to_remove <- unlist(
    x = lapply(
      X = as.list(x = 1:nrow(x = all_edges)),
      FUN = function(i) {
    
        # Define from and to states:
        start_state <- all_edges[i, "from"]
        end_state <- all_edges[i, "to"]
    
        # Make from and to edges that match start and end state:
        from_edges <- all_edges[all_edges[, "from"] == start_state, , drop = FALSE]
        to_edges <- all_edges[all_edges[, "to"] == end_state, , drop = FALSE]
    
        # Merge edges on from-to to find linked edges (may be none):
        merged_edges <- merge(
          x = from_edges,
          y = to_edges,
          by.x = "to",
          by.y = "from"
        )
    
        # If there are indirect routes from start_state to end_state:
        if (nrow(x = merged_edges)) {
      
          # Find minimum indirect cost:
          minimum_cost <- min(
            x = apply(
              X = merged_edges[, c("weight.x", "weight.y")],
              MARGIN = 1,
              FUN = sum
            )
          )
      
          # Find direct cost:
          direct_cost <- all_edges[i, "weight"]
      
          # Important check that matrix is self-consistent (there is no indirect cost shorter than a direct cost) and stop and warn user if not:
          if (minimum_cost < direct_cost) stop("Matrix is not self-consistent!")
      
          # If minimum cost is same as direct cost:
          if (minimum_cost == direct_cost) {
        
            # Can remove ith edge so output this:
            i
        
          # If minimum cost s greater than direct cost:
          } else {
        
            # Return empty vector:
            c()
          }
      
        # If there are no indirect routes:
        } else {
      
          # Return empty vector:
          c()
        }
      }
    )
  )
  
  # If there are redundant edges to remove then remove them:
  if (length(x = to_remove) > 0) all_edges <- all_edges[-to_remove, ]
  
  # Create empty adjacency matrix:
  adjacency_matrix <- matrix(
    data = 0,
    nrow = nrow(x = costmatrix$costmatrix),
    ncol = ncol(x = costmatrix$costmatrix),
    dimnames = list(
      rownames(x = costmatrix$costmatrix),
      colnames(x = costmatrix$costmatrix)
    )
  )
  
  # Populate adjacency matrix:
  for (i in 1:nrow(x = all_edges)) {
    adjacency_matrix[all_edges[i, "from"], all_edges[i, "to"]] <- 1
    adjacency_matrix[all_edges[i, "to"], all_edges[i, "from"]] <- 1
  }
  
  # Create vertices data frame:
  vertices <- data.frame(
    label = sampled_states,
    in_degree = unlist(
      x = lapply(
        X = as.list(x = sampled_states),
        FUN = function(i) sum(x = all_edges[, "to"] == i)
      )
    ),
    out_degree = unlist(
      x = lapply(
        X = as.list(x = sampled_states),
        FUN = function(i) sum(x = all_edges[, "from"] == i)
      )
    ),
    eccentricity = unlist(
      x = lapply(
        X = as.list(x = sampled_states),
        FUN = function(i) {
          visited_states <- vector(mode = "character")
          current_states <- i
          path_length <- 0
          while(
            length(x =
              intersect(
                x = all_edges[, "from"],
                y = setdiff(
                  x = current_states,
                  y = visited_states
                )
              )
            ) > 0
          ) {
            visited_states <- c(visited_states, current_states)
            current_states <- setdiff(
              x = unique(
                x = unlist(
                  x = lapply(
                    X = as.list(x = current_states),
                    FUN = function(j) all_edges[all_edges[, "from"] == j, "to"]
                  )
                )
              ),
              y = visited_states
            )
            if (length(x = current_states) > 0) path_length <- path_length + 1
          }
          path_length
        }
      )
    ),
    periphery = 0,
    centre = 0
  )
  vertices[, "periphery"] <- as.numeric(x = vertices[, "eccentricity"] == max(x = vertices[, "eccentricity"]))
  vertices[, "centre"] <- as.numeric(x = vertices[, "eccentricity"] == min(x = vertices[, "eccentricity"]))
  
  # Compile output as stategraph object:
  stategraph <- list(
    n_vertices = costmatrix$size,
    n_arcs = nrow(x = all_edges),
    n_states = costmatrix$n_states,
    single_states = costmatrix$single_states,
    type = costmatrix$type,
    arcs = all_edges,
    vertices = vertices,
    radius = min(x = vertices[, "eccentricity"]),
    diameter = max(x = vertices[, "eccentricity"]),
    adjacency_matrix = adjacency_matrix,
    directed = ifelse(
      test = costmatrix$symmetry == "Symmetric",
      yes = FALSE,
      no = TRUE
    ),
    includes_polymorphisms = costmatrix$includes_polymorphisms,
    polymorphism_costs = costmatrix$polymorphism_costs,
    polymorphism_geometry = costmatrix$polymorphism_geometry,
    polymorphism_distance = costmatrix$polymorphism_distance,
    includes_uncertainties = costmatrix$includes_uncertainties,
    pruned = costmatrix$pruned,
    dollo_penalty = costmatrix$dollo_penalty,
    base_age = costmatrix$base_age,
    weight = costmatrix$weight
  )
  
  # Set object class as stateGraph:
  class(x = stategraph) <- "stateGraph"
  
  # Return stategraph object:
  stategraph
}
