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
#' A common problem in graph theory is identifying the shortest path to take between two vertices of a connected graph. A stepmatrix also describes a graph, with transition costs representing edge weights that can be asymmetric (e.g., going from 0 to 1 can have a different weight than going from 1 to 0). Finding the shortest path between states - i.e., the path that minimises total weight (cost) - from a stepmatrix has two main applications in Claddis: 1. to check a stepmatrix is internally consistent (no cost is misidentified due to a "cheaper" route being available - solving a problem identified in Maddison and Maddison 2003), and 2. to identify the minimum number of steps a character could have on a tree (an essential aspect of various homoplasy metrics, see Hoyal Cuthill et al. 2010).
#'
#' The function returns a vector describing (one) shortest (i.e., lowest cost) path between \code{start} and \code{end}. If the direct path is shortest this will be simply \code{start} and \code{end}, but if an indirect route is cheaper then other node(s) will appear between these values.
#'
#' In operation the function is inspired by Dijkstra's algorithm (Dijkstra 1959) but differs in some aspects to deal with the special case of a cladistic-style stepmatrix. Essentially multiple paths are considered with the algorithm continuing until either the destination node (\code{end}) is reached or the accumulated cost (path length) exceeds the direct cost (meaning the path cannot be more optimal, lower cost, than the direct one).
#'
#' Note: that because infinite costs are allowed in step matrices to convey that a particular transition is not allowed these are not considered by the function and by default the direct path (\code{start} -> \code{end}) will be returned. (If this is the user's aim they should replaced any \\code{Inf} value with a high number instead.)
#'
#' Note: negative costs are not allowed in stepmatrices and will confound the algorithm.
#'
#' Note: if multiple equally optimal solutions are possible, the function will only return one of them. I.e., just because a solution is not presented it cannot be assumed it is suboptimal. For example, for any ordered character (of three or more states) there will always be multiple equally optimal solutions.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Dijkstra, E. W., 1959. A note on two problems in connexion with graphs. \emph{Numerische Mathematik}, \bold{1}, 269–271.
#'
#' Hoyal Cuthill, J. F., Braddy, S. J. and Donoghue, P. C. J., 2010. A formula for maximum possible steps in multistate characters: isolating matrix parameter effects on measures of evolutionary convergence. \emph{Cladistics}, \bold{26}, 98-102.
#'
#' Maddison, D. R. and Maddison, W. P., 2003. \emph{MacClade 4: Analysis of phylogeny and character evolution}. Version 4.06. Sinauer Associates, Sunderland, Massachusetts.
#'
#' @return
#'
#' A vector of states describing (one) of the optimal path(s) in the order \code{start} to \code{end}.
#'
#' @seealso
#'
#' \link{convert_adjacency_matrix_to_stepmatrix}, \link{make_stepmatrix}
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
#' find_shortest_stepmatrix_path(
#'   stepmatrix = stepmatrix,
#'   start = "0",
#'   end = "3"
#' )
#'
#' # Show that this is directional by asking for the reverse path:
#' find_shortest_stepmatrix_path(
#'   stepmatrix = stepmatrix,
#'   start = "3",
#'   end = "0"
#' )
#'
#' @export find_shortest_stepmatrix_path
find_shortest_stepmatrix_path <- function(stepmatrix, start, end) {

  # Check stepmatrix is correct class and stop and warn user if not:
  if (class(stepmatrix) != "stepMatrix") stop("stepmatrix must be an object o class \"stepMatrix\".")
  
  # Check start and end are in character format so that the correct stepmatrix rows and columns are called:
  if (any(c(!is.character(x = start), !is.character(x = end)))) stop("start and end must be characters (i.e., \"0\" not 0).")
  
  # If start and end are the same cost must be zero (and hence unbeatable) so return start and end:
  if (start == end) return(c(start, end))
  
  # If an infinite cost then stop and return direct path as this indicates the path is not accessible:
  if (stepmatrix$stepmatrix[start, end] == Inf) return(c(start, end))
  
  # Get all character states:
  states <- rownames(x = stepmatrix$stepmatrix)

  # Get direct cost of transition (default option unless a lower cost path is found):
  direct_cost <- stepmatrix$stepmatrix[start, end]

  # Set up a startng path format:
  path <- list(
    start_node = start,
    destination_node = end,
    current_node = start,
    unvisited_nodes = setdiff(x = states, y = start),
    visited_nodes = start,
    length = 0,
    reached_destination = FALSE
  )
  
  # Set up a list of possile path(s):
  paths <- list(path)
  
  # Start a while loop that will run until either destination node is always reached and/or costs are always higher than direct cost:
  while(!all(unlist(x = lapply(X = paths, FUN = function(x) x$reached_destination))) || any(unlist(x = lapply(X = paths, FUN = function(x) ifelse(x$reached_destination, Inf, x$length))) < direct_cost)) {
  
    # Find which paths need to continue (have not reached destination or exceeded drect cost):
    continuing_paths <- which(x = !unlist(x = lapply(X = paths, FUN = function(x) x$reached_destination)))
    
    # For each continuing path (in reverse order so i doesn't get confused):
    for(i in rev(x = continuing_paths)) {
    
      # Isolate unvisited nodes:
      unvisited_nodes <- paths[[i]]$unvisited_nodes
    
      # Store previous node (about to update current one):
      previous_node <- paths[[i]]$current_node
      
      # Build new paths list:
      new_paths <- rep(paths[i], length(x = unvisited_nodes))
    
      # Set current node as each unviisted node:
      for(j in 1:length(x = new_paths)) new_paths[[j]]$current_node <- unvisited_nodes[j]
    
      # Update new paths:
      new_paths <- lapply(
        X = new_paths,
        FUN = function(x) {
          x$unvisited_nodes <- setdiff(x = x$unvisited_nodes, y = x$current_node)
          x$length <- x$length + stepmatrix$stepmatrix[previous_node, x$current_node]
          x$visited_nodes <- c(x$visited_nodes, x$current_node)
          if (x$current_node == x$destination_node) x$reached_destination <- TRUE
          x
        }
      )
    
      # Add new paths to list:
      paths <- c(paths, new_paths)
    }
  
    # Remove (now) old paths:
    paths <- paths[-continuing_paths]
  }
  
  # Make vector of discovered path lengths:
  path_lengths <- unlist(x = lapply(X = paths, FUN = function(x) x$length))

  # If any path lengths are lower than direct cost:
  if (any(path_lengths < direct_cost)) {
  
    # Return first shortest path found:
    return(paths[[which(x = path_lengths == min(x = path_lengths))[1]]]$visited_nodes)
  
  # If direct cost is best:
  } else {
  
    # Return direct path:
    return(c(start, end))
  }
}
