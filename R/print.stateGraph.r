#' Compact display of a stategraph
#'
#' @description
#'
#' Displays a compact summary of a stateGraph object.
#'
#' @param x An object of class \code{"stateGraph"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#'
#' Displays some basic summary information on a stateGraph object.
#'
#' @return
#'
#' Nothing is directly returned, instead a text summary describing a \code{"stateGraph"} object is printed to the console.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Make an unordered costmatrix:
#' example_stategraph <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered"
#' )
#'
#' # Show print.stateGraph version:
#' print.stateGraph(x = example_stategraph)
#'
#' @export print.stateGraph
print.stateGraph <- function(x, ...) {
  

#  "n_vertices" The total number of vertices
#  "n_edges" The total number of edges
#  "n_states" The number of unique single states
#  "single_states" The labels of the unique single states
#  "type" The type of stateGraph (ordered, unordered, dollo, irreversibe, stratigraphy, custom)
#  "edges" Edges matrix (from, to, weight)
#  "vertices" Vertices matrix (label, in degree, out degree, eccentricity)
#  "radius" The radius of the graph (minimum vertex eccenticity)
#  "diameter" The diameter of the graph (maximum vertex eccentricity)
#  "adjacency_matrix" Adjacency matrix representation of the state (di)graph, with rows as "from" states and columns as "to" states.
#  "directed" If TRUE is a digraph, if FALSE is a graph
#  "includes_polymorphisms" Whether (TRUE) or not (FALSE) polymorphic state vertices are included
#  "polymorphism_costs" The means by which costs (edge weights) were assigned to polymorphisms. Must be one of: \code{"additive"}, \code{"geometric"}, \code{"maddison"}, or \code{"stratigraphic"}
#  "polymorphism_geometry" The geometry by which costs (edge weights) were assigned to polymorphisms. Must be one of: \code{"hypercube"}, \code{"hypersphere"}, or \code{"simplex"}
#  "polymorphism_distance" The distance metric by which costs were assigned to polymorphisms. Must be one of: \code{"euclidean"}, \code{"great_circle"}, or \code{"manhattan"}
#  "includes_uncertainties" Whether (TRUE) or not (FALSE) uncertainty state vertices are included
#  "pruned" Whether (TRUE) or not (FALSE) the state graph represents a pruned version
#  "dollo_penalty" A numeric value indicating the penalty used for a Dollo character
#  "base_age" A numeric value indicating the base (oldest) age used for a stratigraphic character
#  "weight" A numeric value indicating the character weight
  
  
  # Check x has class stateGraph and stop and warn user if not:
  if (!inherits(x = x, what = "stateGraph")) stop("x must be an object of class \"stateGraph\".")
  
  # If not a valid stateGraph object then stop and provide feedback to user on what is wrong:
  if (!is.stateGraph(x = x)) stop(check_stateGraph(stategraph = x)[1])
  
  # Return summary information about object:
  cat(
    paste0(
      x$symmetry, ### REPLACE WITH GRAPH/DIGRAPH BASED ON $directed
      " ",
      x$type,
      " stateGraph object containing ",
      x$n_states,
      " unique states",
      ifelse(
        test = all(c(x$includes_polymorphisms, x$includes_uncertainties)),
        yes = paste0(
          " (plus ",
          length(x = grep(pattern = "&", x = colnames(x = x$costmatrix))),
          " polymorphic and ",
          length(x = grep(pattern = "/", x = colnames(x = x$costmatrix))),
          " uncertain states)"
        ),
        no = ""
      ),
      ifelse(
        test = all(c(x$includes_polymorphisms, !x$includes_uncertainties)),
        yes = paste0(
          " (plus ",
          length(x = grep(pattern = "&", x = colnames(x = x$costmatrix))),
          " polymorphic states)"
        ),
        no = ""
      ),
      ifelse(
        test = all(c(!x$includes_polymorphisms, x$includes_uncertainties)),
        yes = paste0(
          " (plus ",
          length(x = grep(pattern = "/", x = colnames(x = x$costmatrix))),
          " uncertain states)"
        ),
        no = ""
      ),
      "."
    )
  )
}
