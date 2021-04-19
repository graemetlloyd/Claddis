#' Finds the shortest path between two states in a stepmatrix
#'
#' @description
#'
#' Given a start and end state, returns the shortest path through a stepmatrix.
#'
#' @param stepmatrix An object of class \code{stepMatrix}.
#' @param start The start state for the requested path.
#' @param end The end state of the requested path.
#'
#' @details
#'
#' A common problem in graph theory is identifying the shortest path to take between two vertices of a connected graph. A stepmatrix also describes a graph, with transition costs representing edge weights that can be asymmetric (e.g., going from 0 to 1 can have a different weight than going from 1 to 0). Finding the shortest path between states - i.e., the path that minimises total weiight (cost) - from a stepmatrix has two main applications in Claddis: 1. to check a stepmatrix is internally consistent (no cost is misidentified due to a "cheaper" route being available - solving a problem identified in MADDISON), and 2. to identify the minimum number of steps a character could have on a tree (an essential aspect of various homoplasy metrics, Hoyal Cuthill et al. 2010).
#'
#' This function implements a slightly modified version of Dijkstra's algorithm (Dijkstra 1959) to deal with the special case of a cladistic-style stepmatrix. The core code is taken from \link[https://www.r-bloggers.com/2020/10/finding-the-shortest-path-with-dijkstras-algorithm/]{this} R-bloggers post.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Dijkstra, E. W., 1959. A note on two problems in connexion with graphs. \emph{Numerische Mathematik}, \bold{1}, 269–271.
#'
#' Hoyal Cuthill, J. F., Braddy, S. J. and Donoghue, P. C. J., 2010. A formula for maximum possible steps in multistate characters: isolating matrix parameter effects on measures of evolutionary convergence. \emph{Cladistics}, \bold{26}, 98-102.
#'
#' @return
#'
#' A vector of states describing a path in order from start to end.
#'
#' @seealso
#'
#' \link{make_all_polymorphisms}
#'
#' @examples
#'
#' # Make a four-state Dollo stepmatrix:
#' stepmatrix <- Claddis::make_stepmatrix(
#'   min_state = 0,
#'   max_state = 3,
#'   character_type = "dollo",
#'   include_polymorphisms = FALSE,
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle",
#'   dollo_penalty = 100
#' )
#'
#' # Find the shortest path from state 0 to state 3:
#' find_shortest_path(
#'   stepmatrix = stepmatrix,
#'   start = "0",
#'   end = "3"
#' )
#'
#' # Show that this is directional by asking for the reverse path:
#' find_shortest_path(
#'   stepmatrix = stepmatrix,
#'   start = "3",
#'   end = "0"
#' )
#'
#' @export find_shortest_path
find_shortest_path <- function(stepmatrix, start, end) {
  
  # Check stepmatrix has class stepMatrix and stop and warn user if not:
  if (!inherits(x = stepmatrix, what = "stepMatrix")) stop("stepmatrix must be an object of class \"stepMatrix\".")
  
  ### RETURN ZERO IF START AND END ARE SAME
  ### NEED TO CHECK IF Inf COST FALLS IN ALL Inf TRIANGLE OF STEPMATRIX AND RETURN JUST DIRECT ROUTE

  # Isolate states:
  states <- rownames(x = stepmatrix$stepmatrix)
  
  # Restate stepmatrix as graph:
  graph <- lapply(
    X = as.list(x = states),
    FUN = function(state) setdiff(x = states, y = state)
  )
  
  # Use transition costs as weights of graph:
  weights <- lapply(
    X = as.list(x = states),
    FUN = function(state) unname(obj = stepmatrix$stepmatrix[state, setdiff(x = states, y = state)])
  )
  
  # Add names to graph:
  names(graph) <- names(weights) <- states
  
  # NEED TO DO BOTH WAYS AND SKIP CALCULATING IF A WHOLE TRIANGLE IS Inf VALUES
  
  
  ### BELOW IS CODE TAKEN FROM https://www.r-bloggers.com/2020/10/finding-the-shortest-path-with-dijkstras-algorithm/ WITH LITTLE MODIFICATION
  
  # Create edgelist with weights
  G <- data.frame(stack(graph), weights = stack(weights)[[1]])
  
  # Sub function to give path length:
  path_length <- function(path) {
    
    # If path is NULL return infinite length:
    if (is.null(path)) return(Inf)
    
    # Get all consecutive nodes:
    pairs <- cbind(values = path[-length(path)], ind = path[-1])
    
    # Join with G and sum over weights
    sum(merge(x = pairs, y = G)[ , "weights"])
  }
  
  # Recursive subfunction to find shortest path:
  find_the_shortest_path <- function(graph, start, end, path = c()) {
    
    # if there are no nodes linked from current node (= dead end) return NULL
    if (is.null(graph[[start]])) return(NULL)
    
    # add next node to path so far
    path <- c(path, start)
    
    # base case of recursion: if end is reached return path
    if (start == end) return(path)
    
    # initialize shortest path as NULL
    shortest <- NULL
    
    # loop through all nodes linked from the current node (given in start)
    for (node in graph[[start]]) {
      
      # proceed only if linked node is not already in path
      if (!(node %in% path)) {
        
        # recursively call function for finding shortest path with node as start and assign it to newpath
        newpath <- find_the_shortest_path(graph, node, end, path)
        
        # if newpath is shorter than shortest so far assign newpath to shortest
        if (path_length(newpath) < path_length(shortest))
        shortest <- newpath
      }
    }
    
    # return shortest path
    shortest
  }
  
  # Return shortest path:
  find_the_shortest_path(graph = graph, start = start, end = end)
}

# HOMOPLASY FUNCTION
# - Option to define outgroup (root state?) as will affect min steps value
# - Otherwise need to do path both ways, e.g., 0->1->2 and 2->1->0 in case of asymmetric characters
# - Keep check for variance as no changes at all if invariant.
# - How to handle uncertainties?
