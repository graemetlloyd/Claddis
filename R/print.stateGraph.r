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
#' # Make an example stategraph:
#' example_stategraph <- list(
#'   n_vertices = 6,
#'   n_arcs = 12,
#'   n_states = 6,
#'   single_states = as.character(x = 0:5),
#'   type = "custom",
#'   arcs = data.frame(
#'     from = as.character(x = c(0, 1, 0, 2, 2, 5, 1, 4, 5, 4, 3, 4)),
#'     to = as.character(x = c(1, 0, 2, 0, 5, 2, 4, 1, 4, 5, 4, 3)),
#'     weight = c(1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1)
#'   ),
#'   vertices = data.frame(
#'     label = as.character(x = 0:5),
#'     in_degree = c(2, 2, 2, 1, 3, 2),
#'     out_degree = c(2, 2, 2, 1, 3, 2),
#'     eccentricity = c(3, 2, 3, 3, 2, 2),
#'     periphery = c(1, 0, 1, 1, 0, 0),
#'     centre = c(0, 1, 0, 0, 1, 1)
#'   ),
#'   radius = 2,
#'   diameter = 3,
#'   adjacency_matrix = matrix(
#'     data = c(
#'       0, 1, 1, 0, 0, 0,
#'       1, 0, 0, 0, 1, 0,
#'       1, 0, 0, 0, 0, 1,
#'       0, 0, 0, 0, 1, 0,
#'       0, 1, 0, 1, 0, 1,
#'       0, 0, 1, 0, 1, 0
#'     ),
#'     nrow = 6,
#'     byrow = TRUE,
#'     dimnames = list(0:5, 0:5)
#'   ),
#'   directed = FALSE,
#'   includes_polymorphisms = FALSE,
#'   polymorphism_costs = "additive",
#'   polymorphism_geometry = "simplex",
#'   polymorphism_distance = "euclidean",
#'   includes_uncertainties = FALSE,
#'   pruned = FALSE,
#'   dollo_penalty = 999,
#'   base_age = 100,
#'   weight = 1
#' )
#'
#' # Set class as stateGraph:
#' class(x = example_stategraph) <- "stateGraph"
#'
#' # Show print.stateGraph version:
#' print.stateGraph(x = example_stategraph)
#'
#' @export print.stateGraph
print.stateGraph <- function(x, ...) {
  
#  "n_vertices" The total number of vertices
#  "n_arcs" The total number of arcs (2 arcs for an edge)
#  "n_states" The number of unique single states
#  "single_states" The labels of the unique single states
#  "type" The type of stateGraph (ordered, unordered, dollo, irreversibe, stratigraphy, custom)
#  "arcs" Arcs matrix (from, to, weight) - an edge is composed of two symmetric arcs
#  "vertices" Vertices matrix (label, in degree, out degree, eccentricity, periphery (1/0), centre (0/1))
#  "radius" The radius of the graph (minimum vertex eccenticity)
#  "diameter" The diameter of the graph (maximum vertex eccentricity)
#  "adjacency_matrix" Adjacency matrix representation of the state graph, with rows as "from" states and columns as "to" states.
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
      paste0(
        toupper(x = strsplit(x = x$type, split = "")[[1]][1]),
        paste(strsplit(x = x$type, split = "")[[1]][-1], collapse = "")
      ),
      " state ",
      ifelse(test = x$directed, yes = "digraph", no = "graph"),
      " composed of ",
      x$n_arcs,
      " arcs connecting ",
      x$n_states,
      " unique states",
      ifelse(
        test = all(c(x$includes_polymorphisms, x$includes_uncertainties)),
        yes = paste0(
          " (plus ",
          length(x = grep(pattern = "&", x = colnames(x = x$adjacency_matrix))),
          " polymorphic and ",
          length(x = grep(pattern = "/", x = colnames(x = x$adjacency_matrix))),
          " uncertain states)"
        ),
        no = ""
      ),
      ifelse(
        test = all(c(x$includes_polymorphisms, !x$includes_uncertainties)),
        yes = paste0(
          " (plus ",
          length(x = grep(pattern = "&", x = colnames(x = x$adjacency_matrix))),
          " polymorphic states)"
        ),
        no = ""
      ),
      ifelse(
        test = all(c(!x$includes_polymorphisms, x$includes_uncertainties)),
        yes = paste0(
          " (plus ",
          length(x = grep(pattern = "/", x = colnames(x = x$adjacency_matrix))),
          " uncertain states)"
        ),
        no = ""
      ),
      "."
    )
  )
}
