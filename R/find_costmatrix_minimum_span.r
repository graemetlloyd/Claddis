#' Finds a minimum spanning tree of a costmatrix
#'
#' @description
#'
#' Given a costmatrix, returns a shortest tree connecting every state.
#'
#' @param costmatrix An object of class \code{costMatrix}.
#'
#' @details
#'
#' The minimum parsimony length a phylogenetic hypothesis can have depends on both the costmatrix of transition costs and the states actually sampled. If the costmatrix rows and columns already represent sampled states (the assumption here) then this minimum length is reduced to a graph theory problem - the minimum spanning tree or, if a directed graph, the minimum weight spanning arboresence. (NB: if there are unsampled states then the \code{find_steiner_tree_of_stategraph} function should be used instead.) This function returns one such shortest tree (although others may exist). The sum of the weights of the edges or arcs returned is the minimum cost.
#'
#' As the algorithms used are graph theory based the function operates by simply calling \link{convert_costmatrix_to_stategraph} and \link{find_stategraph_minimum_span}. In practice, if the costmatrix represents a graph (transition costs are all symmetric) then Kruskal's algorithm is applied (Kruskal 1956). If costs are asymmetric, however, then the graph representation is a directed graph (or digraph) and so a version of Edmonds' algorithm is applied (Edmonds 1967).
#'
#' Note that Dollo characters represent a special case solution as although a penalty weight is applied to the edges intended to only ever be traversed once this weight should not be used when calculating tree lengths. The function catches this and returns the edges with the weight that would actually be counted for a minimum weight spanning tree.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Edmonds, J., 1967. Optimum branchings. \emph{Journal of Research of the National Bureau of Standards Section B}, \bold{71B}, 233-240.
#'
#' Kruskal, J. B., 1956. On the shortest spanning subtree of a graph and the traveling salesman problem. \emph{Proceedings of the American Mathematical Society}, \bold{7}, 48-50.
#'
#' @return
#'
#' A \code{data.frame} object describing a minimum spanning tree or minimum weight arboresence as a series of edges or arcs.
#'
#' @seealso
#'
#' \link{find_shortest_costmatrix_path}
#'
#' @examples
#'
#' # Make a four-state ordered character costmatrix:
#' ordered_costmatrix <- make_costmatrix(
#'   min_state = "0",
#'   max_state = "3",
#'   character_type = "ordered"
#' )
#'
#' # Find length of shortest spanning tree of costmatrix:
#' find_costmatrix_minimum_span(costmatrix = ordered_costmatrix)
#'
#' # Make a four-state unordered character costmatrix:
#' unordered_costmatrix <- make_costmatrix(
#'   min_state = "0",
#'   max_state = "3",
#'   character_type = "unordered"
#' )
#'
#' # Find length of shortest spanning tree of costmatrix:
#' find_costmatrix_minimum_span(costmatrix = unordered_costmatrix)
#'
#' # Make a four-state irreversible character costmatrix:
#' irreversible_costmatrix <- make_costmatrix(
#'   min_state = "0",
#'   max_state = "3",
#'   character_type = "irreversible"
#' )
#'
#' # Find length of shortest spanning tree of costmatrix:
#' find_costmatrix_minimum_span(costmatrix = irreversible_costmatrix)
#'
#' # Make a four-state Dollo character costmatrix:
#' dollo_costmatrix <- make_costmatrix(
#'   min_state = "0",
#'   max_state = "3",
#'   character_type = "dollo"
#' )
#'
#' # Find length of shortest spanning tree of costmatrix:
#' find_costmatrix_minimum_span(costmatrix = dollo_costmatrix)
#'
#' @export find_costmatrix_minimum_span
find_costmatrix_minimum_span <- function(costmatrix) {
  
  ### ADD DATA CHECK
  ### CANNOT DEAL WTH POLYMORPHISMS OR UNCERTAINTIES WITHOUT MODIFICATION!
  
  # Special case of no edges (return empty matrix of edges):
  if (costmatrix$size == 1) return(data.frame(from = "0", to = "1", weight = 1)[-1, ])
  
  # Begin by converting costmatrix to state graph ready for using graph theory algorithms:
  stategraph <- convert_costmatrix_to_stategraph(costmatrix = costmatrix)

  # Find and return minimum span:
  find_stategraph_minimum_span(stategraph = stategraph)
}
