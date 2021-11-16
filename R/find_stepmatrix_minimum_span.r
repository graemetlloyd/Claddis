#' Finds the minimum spanning tree length for a stepmatrix
#'
#' @description
#'
#' Given a stepmatrix, returns the length of the shortest tree that connects every state.
#'
#' @param stepmatrix An object of class \code{stepMatrix}.
#'
#' @details
#'
#' The minimum parsimony length a phylogenetic hypothesis can have depends on the stepmatrix of transition costs and the states actually sampled. If the stepmatrix rows and columns already represent sampled states (the assumption here) then this minimum length is reduced to a graph theory problem - the minimum spanning tree. This function returns the length of that tree and hence the minimum cost (as used in stratocladistics to help define "parsimony debt" Fisher 1994, or in homoplasy indices, see Hoyal Cuthill 2015).
#'
#' \bold{Special case solutions}
#'
#' In order for the function to run as fast as possible the following special case solutions are used.
#'
#' \emph{Invariant character solution}
#'
#' If the character is invariant then the stepmatrix will have size one (one-by-one), and more simply the minimum spanning tree is length zero as there re no transitions to make.
#'
#' \emph{Unordered character solution}
#'
#' If the character is unordered all off-diagonal values of the stepmatrix are one and the minimum steps is the size of the stepmatrix minus one (as there is no cost to visit the root state).
#'
#' \emph{Ordered/irreversible/stratigraphic character solution}
#'
#' For an ordered, irreversible, or stratigraphic character the solution is also simple, and identical. This is the cost to go between the furthest two states which is simply the value in the upper right of the stepmatrix.
#'
#' \emph{Dollo character solution}
#'
#' Despite other complexities Dollo characters also have a simple solution that mirrors that for ordered, irreversible, or stratigraphic characters. I.e., the derived states must all be acquired once which is simply the value in the lower left of the stepmatrix.
#'
#' \bold{General case solution}
#'
#' For characters that do not constitute any of the above special case solutions a slower running general case solution is required. This is effectively limited to so-called "custom" characters - \code{stepmatrix$type = "custom"} that can only be represented by either a state tree or stepmatrix.
#'
#' Strictly speaking the solution used depends on whether the stepmatrix is symmetric (representing an undirected graph, or simply a graph) or asymmetric (representing a directed graph, or digraph). In both cases the edge weights of the graph/digraph are the transition costs from the stepmatrix. For a graph the problem is the "regular" minimum spanning tree, and if for a digraph then it is the arboresence. For regular minimum spanning trees Kruskal's (1956) algorithm is used and for aboresence's Edmonds' (1967) algorithm is used. (The function automatically makes this choice using the value of \code{stepmatrix$symmetry}.) In practice the implementations used come from the optrees package (Fontenla 2014), and specifically the functions \code{getMinimumSpanningTree} and \code{msArborEdmonds}.
#'
#' Note that regardless of method the function only returns the \emph{length} of the minimum spanning tree or arboresence (i.e., the total minimum cost), not the actual path (pattern of states traversed). Similarly, multiple paths (trees or arboresences) may exist that have this minimum length.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Edmonds, J., 1967. Optimum branchings. \emph{Journal of Research of the National Bureau of Standards Section B}, \bold{71B}, 233-240.
#'
#' Fisher, D. C., 1994. Stratocladistics: morphological and temporal patterns and their relation to phylogenetic process. In L. Grande and O. Rieppel (eds.), \emph{Interpreting the Hierarchy of Nature}. Academic Press, San Diego. pp133–171.
#'
#' Fontenla, M., 2014. optrees: Optimal Trees in Weighted Graphs. R package version 1.0. https://CRAN.R-project.org/package=optrees
#'
#' Hoyal Cuthill, J., 2015. The size of the character state space affects the occurrence and detection of homoplasy: modelling the probability of incompatibility for unordered phylogenetic characters. \emph{Journal ofTheoretical Biology}, \bold{366}, 24-32.
#'
#' Kruskal, J. B., 1956. On the shortest spanning subtree of a graph and the traveling salesman problem. \emph{Proceedings of the American Mathematical Society}, \bold{7}, 48-50.
#'
#' @return
#'
#' A scalar indicating the length in steps of the minimum spanning tree or arboresence.
#'
#' @seealso
#'
#' \link{find_shortest_stepmatrix_path}
#'
#' @examples
#'
#' # Make a six-state ordered character stepmatrix:
#' stepmatrix <- make_stepmatrix(
#'   min_state = "0",
#'   max_state = "5",
#'   character_type = "ordered",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle"
#' )
#'
#' # Find length of shortest spanning tree of stepmatrix:
#' find_stepmatrix_minimum_span(stepmatrix = stepmatrix)
#'
#' @export find_stepmatrix_minimum_span
find_stepmatrix_minimum_span <- function(stepmatrix) {
  
  # ADD SOME DATA CHECKS?
  # CANNOT DEAL WTH POLYMORPHISMS OR UNCERTAINTIES WITHOUT MODIFICATION!

  # Things are simple if character is invariant (one-by-one matrix) - there are no changes:
  if (stepmatrix$size == 1) return(0)

  # Things are simple if is an unordered character (there must minimally be N - 1 links):
  if (stepmatrix$type == "unordered") return(stepmatrix$size - 1)
  
  # Things are also simple if it is an ordered, irreversible or stratigraphic character (top right value is minimum steps):
  if (stepmatrix$type == "ordered" || stepmatrix$type == "irreversible" || stepmatrix$type == "stratigraphy") return(stepmatrix$stepmatrix[1, stepmatrix$size])
  
  # Things are simple if is an Dollo character (lower left value is minimum steps):
  if (stepmatrix$type == "dollo") return(stepmatrix$stepmatrix[stepmatrix$size, 1])
  
  # Set sampled states and stepmatrix rows and columns to simple numbers to work with optrees:
  sampled_states <- rownames(x = stepmatrix$stepmatrix) <- colnames(x = stepmatrix$stepmatrix) <- 1:stepmatrix$size
  
  # Build matrix of shortest incoming edge(s) (those available for the aboresence):
  available_edges <- do.call(
    what = rbind,
    args = lapply(
      X = as.list(x = sampled_states),
      FUN = function(state) {
        cbind(
          from = names(x = which(x = stepmatrix$stepmatrix[, state] == min(x = stepmatrix$stepmatrix[-which(x = rownames(x = stepmatrix$stepmatrix) == state), state]))),
          to = state
        )
      }
    )
  )
  
  # Add weights of edges:
  available_edges <- cbind(from = as.numeric(x = available_edges[, "from"]), to = as.numeric(x = available_edges[, "to"]), weight = apply(X = available_edges, MARGIN = 1, function(x) stepmatrix$stepmatrix[x[1], x[2]]))
  
  # If graph (stepmatrix) is symmetric then can use Kruskal:
  if (stepmatrix$symmetry == "Symmetric") {
    
    # First remove duplicate edges (due to symmetries):
    available_edges <- matrix(
      data = as.numeric(x = unlist(x = strsplit(x = unique(x = apply(
        X = available_edges,
        MARGIN = 1,
        FUN = function(x) {
          x[1:2] <- sort(x = x[1:2])
          paste(x = x[1:3], collapse = "%")
        }
      )), split = "%"))),
      ncol = 3,
      byrow = TRUE,
      dimnames = list(c(), c("from", "to", "weight"))
    )
    
    # Return length of minimum spanning tree:
    return(sum(x = optrees::getMinimumSpanningTree(
      nodes = sampled_states,
      arcs = available_edges,
      algorithm = "Kruskal",
      start.node = 1,
      show.data = FALSE,
      show.graph = FALSE,
      check.graph = FALSE
    )$tree.arcs[, 3]))
  }
  
  # If graph (stepmatrix) is asymmetric then use Edmonds' algorithm:
  if (stepmatrix$symmetry == "Asymmetric") {
    
    # Return length of minimum arboresence:
    return(sum(x = optrees::msArborEdmonds(
      nodes = as.numeric(x = sampled_states),
      arcs = available_edges,
      source.node = 1,
      stages.data = FALSE
    )$tree.arcs[, 3]))
  }
}
