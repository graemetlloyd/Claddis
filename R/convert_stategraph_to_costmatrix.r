#' Convert a minimal state graph to a costmatrix
#'
#' @description
#'
#' Given a smallest possible state graph (fewest edges), returns the corresponding costmatrix.
#'
#' @param state_graph A matrix describing a state graph with columns indicating from states, to states and weights (costs).
#' @param fix_ratio A logical indicating whether (\cost{TRUE}) or not (\code{FALSE}, the default) to fix costs such that the minimum cost is one.
#'
#' @details
#'
#' A state graph describes the relationship between character states in terms of a graph of vertices (character states) and edges (weights). Edges can be "symmetric" (the graph is undirected) or "asymmetric" (the graph is directed, a digraph).
#'
#' For example, a simple symmetric binary state graph might look ike this:
#'
#' \preformatted{  1
#' A---B}
#'
#' Here the two states are A and B and the weight (cost of transitioning from A to B or B to A) is one.
#'
#' In Claddis this graph can be represented using a data.frame of edges:
#'
#' \preformatted{from to weight
#'    0  1      1
#'    1  0      1}
#'
#' Each row represents an edge, from one vertex to another, and the corresponding weight (transition cost). Note that for symmetric graphs the edge is stated twice, once for each direction of the transition.
#'
#' This function converts these state transitions to costmatrices, including interpolating any missing transitions. Here the costmatrix would look quite simple:
#'
#' \preformatted{    ---------
#'     | 0 | 1 |
#' -------------
#' | 0 | 0 | 1 |
#' -------------
#' | 1 | 1 | 0 |
#' -------------}
#'
#' The function also offers the option to "fix" the ratio of costs such that the minimum edge weight (minimum direct state-to-state transition) is one. Note: this is unlikely to be something most users will need and so defaults to false.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' An object of class \code{costMatrix}.
#'
#' @seealso
#'
#' \link{convert_adjacency_matrix_to_costmatrix}
#' \link{convert_costmatrix_to_stategraph}
#'
#' @examples
#'
#' # Make state graph for a four-state linear ordered character:
#' state_graph <- data.frame(
#'   from = c("0", "1", "2", "3", "2", "1"),
#'   to = c("1", "2", "3", "2", "1", "0"),
#'   weight = c(1, 1, 1, 1, 1, 1),
#'   stringsAsFactors = FALSE
#' )
#'
#' # View state graph:
#' state_graph
#'
#' # Convert state graph to a costmatrix:
#' costmatrix <- convert_stategraph_to_costmatrix(state_graph = state_graph)
#'
#' # Show costmatrix reflects linear ordered costs:
#' costmatrix$costmatrix
#'
#' @export convert_stategraph_to_costmatrix
convert_stategraph_to_costmatrix <- function(state_graph, fix_ratio = FALSE) {
  
  ### ADD CHECKS AT END FOR COMMON COSTMATRIX TYPES?
  
  # Identify state labels from state graph:
  states <- unique(x = sort(x = c(state_graph[, "from"], state_graph[, "to"])))
  
  # Set number of states:
  n_states <- length(x = states)
  
  # Create base costmatrix object (will need populating):
  costmatrix <- matrix(
    data = Inf,
    nrow = n_states,
    ncol = n_states,
    dimnames = list(states, states)
  )
  
  # Set diagonal as all zeroes:
  diag(x = costmatrix) <- 0
  
  # Populate costmatrix with edge weights:
  for(i in 1:nrow(x = state_graph)) costmatrix[state_graph[i, "from"], state_graph[i, "to"]] <- state_graph[i, "weight"]
  
  # Format as a costmatrix object:
  costmatrix <- list(
    size = n_states,
    type = "custom",
    costmatrix = costmatrix,
    symmetry = ifelse(
      test = isSymmetric(object = costmatrix),
      yes = "Symmetric",
      no = "Asymmetric"
    ),
    includes_polymorphisms = FALSE
  )
  
  # Set class to costMatrix:
  class(x = costmatrix) <- "costMatrix"
  
  # Fix any missing costs by interpolation:
  costmatrix <- fix_costmatrix(costmatrix = costmatrix, message = FALSE)
  
  # If fixing ratio (i.e., minimum cost as one):
  if (fix_ratio) {
  
    # Identify minimum cost:
    minimum_cost <- min(
      x = c(
        costmatrix$costmatrix[upper.tri(x = costmatrix$costmatrix)],
        costmatrix$costmatrix[lower.tri(x = costmatrix$costmatrix)]
      )
    )
    
    # Divide through by minimum cost:
    costmatrix$costmatrix <- costmatrix$costmatrix / minimum_cost
  }
  
  # Return costmatrix object:
  costmatrix
}
