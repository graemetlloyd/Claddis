#' Convert a minimal state graph to a costmatrix
#'
#' @description
#'
#' Given a state graph returns the corresponding costmatrix.
#'
#' @param stategraph An object of class \code{stateGraph}.
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
#' In Claddis this graph can be represented using a \code{stateGraph} object including the arcs that describe the (di)graph:
#'
#' \preformatted{from to weight
#'    0  1      1
#'    1  0      1}
#'
#' Each row represents an arc from one vertex to another, and the corresponding weight (transition cost). Note that for symmetric graphs the edge is stated using two arcs (from \emph{i} to \emph{j} and from \emph{j} to \emph{i}).
#'
#' This function converts these state transitions to costmatrices, including interpolating any missing transitions by using the shortest indirect cost. In the example used here the costmatrix would look quite simple:
#'
#' \preformatted{    ---------
#'     | 0 | 1 |
#' -------------
#' | 0 | 0 | 1 |
#' -------------
#' | 1 | 1 | 0 |
#' -------------}
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
#' # Make state graph for a six-state linear ordered character:
#' stategraph <- list(
#'   n_vertices = 6,
#'   n_arcs = 10,
#'   n_states = 6,
#'   single_states = c("0", "1", "2", "3", "4", "5"),
#'   type = "ordered",
#'   arcs = data.frame(
#'     from = c("1", "0", "2", "1", "3", "2", "4", "3", "5", "4"),
#'     to = c("0", "1", "1", "2", "2", "3", "3", "4", "4", "5"),
#'     weight = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#'   ),
#'   vertices = data.frame(
#'     label = c("0", "1", "2", "3", "4", "5"),
#'     in_degree = c(1, 2, 2, 2, 2, 1),
#'     out_degree = c(1, 2, 2, 2, 2, 1),
#'     eccentricity = c(5, 4, 3, 3, 4, 5),
#'     periphery = c(1, 0, 0, 0, 0, 1),
#'     centre = c(0, 0, 1, 1, 0, 0)
#'   ),
#'   radius = 3,
#'   diameter = 5,
#'   adjacency_matrix = matrix(
#'     data = c(
#'       0, 1, 0, 0, 0, 0,
#'       1, 0, 1, 0, 0, 0,
#'       0, 1, 0, 1, 0, 0,
#'       0, 0, 1, 0, 1, 0,
#'       0, 0, 0, 1, 0, 1,
#'       0, 0, 0, 0, 1, 0
#'     ),
#'     ncol = 6,
#'     byrow = TRUE,
#'     dimnames = list(
#'       c("0", "1", "2", "3", "4", "5"),
#'       c("0", "1", "2", "3", "4", "5")
#'     )
#'   ),
#'   directed = FALSE,
#'   includes_polymorphisms = FALSE,
#'   polymorphism_costs = "additive",
#'   polymorphism_geometry = "simplex",
#'   polymorphism_distance = "euclidean",
#'   includes_uncertainties = FALSE,
#'   pruned = FALSE,
#'   dollo_penalty = 999,
#'   base_age = 1,
#'   weight = 1
#' )
#'
#' # Set calss as stateGraph:
#' class(x = stategraph) <- "stateGraph"
#'
#' # View state graph:
#' stategraph
#'
#' # Convert state graph to a costmatrix:
#' costmatrix <- convert_stategraph_to_costmatrix(stategraph = stategraph)
#'
#' # Show costmatrix reflects linear ordered costs:
#' costmatrix$costmatrix
#'
#' @export convert_stategraph_to_costmatrix
convert_stategraph_to_costmatrix <- function(stategraph) {
  
  # Identify state labels from stategraph:
  states <- stategraph$vertices$label
  
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
  for(i in 1:nrow(x = stategraph$arcs)) costmatrix[stategraph$arcs[i, "from"], stategraph$arcs[i, "to"]] <- stategraph$arcs[i, "weight"]
  
  # Format as a costmatrix object:
  costmatrix <- list(
    size = stategraph$n_vertices,
    n_states = stategraph$n_states,
    single_states = stategraph$single_states,
    type = stategraph$type,
    costmatrix = costmatrix,
    symmetry = ifelse(
      test = stategraph$directed,
      yes = "Asymmetric",
      no = "Symmetric"
    ),
    includes_polymorphisms = stategraph$includes_polymorphisms,
    polymorphism_costs = stategraph$polymorphism_costs,
    polymorphism_geometry = stategraph$polymorphism_geometry,
    polymorphism_distance = stategraph$polymorphism_distance,
    includes_uncertainties = stategraph$includes_uncertainties,
    pruned = stategraph$pruned,
    dollo_penalty = stategraph$dollo_penalty,
    base_age = stategraph$base_age,
    weight = stategraph$weight
  )

  # Set class to costMatrix:
  class(x = costmatrix) <- "costMatrix"
  
  # Fill in Inf states with shortest paths:
  costmatrix$costmatrix[states, states] <- do.call(
    what = rbind,
    args = lapply(
      X = as.list(x = states),
      FUN = function(state) {
        apply(
          X = cbind(state, states),
          MARGIN = 1,
          FUN = function(fromtopair) {
            shortest_path <- find_shortest_costmatrix_path(
              costmatrix = costmatrix,
              start = fromtopair[1],
              end = fromtopair[2]
            )
            path_length <- 0
            for(i in 2:length(shortest_path)) path_length <- path_length + costmatrix$costmatrix[shortest_path[(i - 1)], shortest_path[i]]
            path_length
          }
        )
      }
    )
  )
  
  # Return costmatrix object:
  costmatrix
}
