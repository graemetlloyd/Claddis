#' Finds the shortest path between two states in a costmatrix
#'
#' @description
#'
#' Given a start and end state, returns the shortest path through a costmatrix.
#'
#' @param costmatrix An object of class \code{costMatrix}.
#' @param start The start state for the requested path.
#' @param end The end state of the requested path.
#'
#' @details
#'
#' A common problem in graph theory is identifying the shortest path to take between two vertices of a connected graph. A costmatrix also describes a graph, with transition costs representing edge weights that can be asymmetric (e.g., going from 0 to 1 can have a different weight than going from 1 to 0). Finding the shortest path between states - i.e., the path that minimises total weight (cost) - from a costmatrix has two main applications in Claddis: 1. to check a costmatrix is internally consistent (no cost is misidentified due to a "cheaper" route being available - solving a problem identified in Maddison and Maddison 2003), and 2. to identify the minimum cost a character could have on a tree (an essential aspect of various homoplasy metrics, see Hoyal Cuthill et al. 2010).
#'
#' The function returns a vector describing (one) shortest (i.e., lowest cost) path between \code{start} and \code{end}. If the direct path is shortest this will be simply \code{start} and \code{end}, but if an indirect route is cheaper then other node(s) will appear between these values.
#'
#' In operation the function is inspired by Dijkstra's algorithm (Dijkstra 1959) but differs in some aspects to deal with the special case of a cladistic-style costmatrix. Essentially multiple paths are considered with the algorithm continuing until either the destination node (\code{end}) is reached or the accumulated cost (path length) exceeds the direct cost (meaning the path cannot be more optimal, lower cost, than the direct one).
#'
#' Note: that because infinite costs are allowed in costmatrices to convey that a particular transition is not allowed these are always likely to be replaced (meanng the path does become possible) unless they apply to entire rows or columns of a costmatrix.
#'
#' Note: negative costs are not allowed in costmatrices as they will confound the halting criteria of the algorithm. (They also do not make logical sense in a costmatrix anyway!)
#'
#' Note: if multiple equally optimal solutions are possible, the function will only return one of them. I.e., just because a solution is not presented it cannot be assumed it is suboptimal.
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
#' \link{convert_adjacency_matrix_to_costmatrix}, \link{make_costmatrix}
#'
#' @examples
#'
#' # Make a four-state Dollo costmatrix:
#' costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 3,
#'   character_type = "dollo",
#'   dollo_penalty = 100
#' )
#'
#' # Find the shortest path from state 0 to state 3:
#' find_shortest_costmatrix_path(
#'   costmatrix = costmatrix,
#'   start = "0",
#'   end = "3"
#' )
#'
#' # Show that this is directional by asking for the reverse path:
#' find_shortest_costmatrix_path(
#'   costmatrix = costmatrix,
#'   start = "3",
#'   end = "0"
#' )
#'
#' @export find_shortest_costmatrix_path
find_shortest_costmatrix_path <- function(costmatrix, start, end) {

  # Check costmatrix is correct class and stop and warn user if not:
  if (!inherits(x = costmatrix, what = "costMatrix")) stop("costmatrix must be an object of class \"costMatrix\".")
  
  # Check start and end are in character format so that the correct costmatrix rows and columns are called:
  if (any(c(!is.character(x = start), !is.character(x = end)))) stop("start and end must be characters (i.e., \"0\" not 0).")
  
  # If start and end are the same cost must be zero (and hence unbeatable) so return start and end:
  if (start == end) return(c(start, end))
  
  # Get all character states:
  states <- rownames(x = costmatrix$costmatrix)

  # Get direct cost of transition (default option unless a lower cost path is found):
  direct_cost <- costmatrix$costmatrix[start, end]

  # Set up a startng path format:
  path <- list(
    start_node = start,
    destination_node = end,
    current_node = start,
    unvisited_nodes = setdiff(x = states, y = start),
    visited_nodes = start,
    length = 0,
    reached_destination = FALSE,
    exceeded_direct_cost = FALSE
  )
  
  # Set up a list of possible path(s):
  paths <- list(path)
  
  # Start a while loop that will run until either destination node is always reached and/or costs are always higher than direct cost:
  while(
    any(
      unlist(
        x = lapply(
          X = paths,
          FUN = function(path) all(x = c(!path$reached_destination, !path$exceeded_direct_cost))
        )
      )
    )
  ) {
  
    # Find which paths need to continue (have not reached destination or exceeded direct cost):
    continuing_paths <- which(
      x = unlist(
        x = lapply(
          X = paths,
          FUN = function(path) all(x = !c(path$reached_destination, path$exceeded_direct_cost))
        )
      )
    )
    
    # For each continuing path (in reverse order so i doesn't get confused):
    for(i in rev(x = continuing_paths)) {
    
      # Isolate unvisited nodes:
      unvisited_nodes <- paths[[i]]$unvisited_nodes
    
      # Store previous node (about to update current one):
      previous_node <- paths[[i]]$current_node
      
      # Build new paths list:
      new_paths <- rep(paths[i], length(x = unvisited_nodes))
    
      # Set current node as each unvisited node:
      for(j in 1:length(x = new_paths)) new_paths[[j]]$current_node <- unvisited_nodes[j]
    
      # Update new paths:
      new_paths <- lapply(
        X = new_paths,
        FUN = function(x) {
          x$unvisited_nodes <- setdiff(x = x$unvisited_nodes, y = x$current_node)
          x$length <- x$length + costmatrix$costmatrix[previous_node, x$current_node]
          x$visited_nodes <- c(x$visited_nodes, x$current_node)
          if (x$current_node == x$destination_node) x$reached_destination <- TRUE
          if (x$length >= direct_cost) x$exceeded_direct_cost <- TRUE
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
