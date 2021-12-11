#' Convert a costmatrix to a minimal graph
#'
#' @description
#'
#' Given a costmatrix, returns the smallest possible graph (fewest edges).
#'
#' @param costmatrix An object of class \code{costMatrix}.
#'
#' @details
#'
#' A costmatrix summarises all possible state-to-state transition costs and hence each entry could also be considered an edge of a directed state graph. However, many of these edges could be removed and a complete description of the graph still be provided. For example, the diagonal (any transition from a state to itself) can be removed, as can any edge with infinite cost (as this edge would never be traversed in practice). Finally, some edges are redundant as indirect paths already represent the same cost.
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
#' This function effectively generates the latter (the minimal directed graph representation as a matrix of edges).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' A matrix representing the minimum graph as directed edges (from, to, and edge weight). For an undirected graph both ways an edge can be defined are included.
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
#'   character_type = "unordered",
#'   polymorphism_shape = "simplex",
#'   polymorphism_distance = "euclidean"
#' )
#'
#' # Find the minimal directed graph representation:
#' convert_costmatrix_to_graph(costmatrix = unordered_costmatrix)
#'
#' # Make a six-state ordered character costmatrix:
#' ordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 5,
#'   character_type = "ordered",
#'   polymorphism_shape = "simplex",
#'   polymorphism_distance = "euclidean"
#' )
#'
#' # Find the minimal directed graph representation:
#' convert_costmatrix_to_graph(costmatrix = ordered_costmatrix)
#'
#' # Make a six-state stratigraphic character costmatrix:
#' stratigraphic_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 5,
#'   character_type = "stratigraphy",
#'   polymorphism_shape = "simplex",
#'   polymorphism_distance = "euclidean",
#'   state_ages = c(103, 91.4, 78.2, 73.4, 66.0, 59.7)
#' )
#'
#' # Find the minimal directed graph representation:
#' convert_costmatrix_to_graph(costmatrix = stratigraphic_costmatrix)
#'
#' @export convert_costmatrix_to_graph
convert_costmatrix_to_graph <- function(costmatrix) {
  
  # Store original states for using as output:
  original_states <- colnames(x = costmatrix$costmatrix)
  
  # Set sampled states and costmatrix rows and columns to simple numbers to work with optrees:
  sampled_states <- rownames(x = costmatrix$costmatrix) <- colnames(x = costmatrix$costmatrix) <- 1:costmatrix$size
  
  # Create table of all possible edges:
  all_edges <- expand.grid(from = sampled_states, to = sampled_states)
  
  # Add weight to edges:
  all_edges <- cbind(from = all_edges[, "from"], to = all_edges[, "to"], weight = apply(X = all_edges, MARGIN = 1, function(x) costmatrix$costmatrix[x[1], x[2]]))
  
  # Prune all zero weight edges (i.e., diagonal):
  all_edges <- all_edges[-which(x = all_edges[, "weight"] == 0), , drop = FALSE]
  
  # Prune any Infinite weight edges:
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
  
  ###
  # CHECK GRAPH IS *CONNECTED*
  #optrees::checkArbor
  #optrees::checkGraph
  
  # Return graph as edges:
  cbind(
    from = as.numeric(x = original_states[all_edges[, "from"]]),
    to = as.numeric(x = original_states[all_edges[, "to"]]),
    weight = as.numeric(x = all_edges[, "weight"])
  )
}
