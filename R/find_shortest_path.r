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
#' This function implements a slightly modified version of Dijkstra's algorithm (Dijkstra 1959) to deal with the special case of a cladistic-style stepmatrix. The core code is taken from \link[https://www.r-bloggers.com/2020/10/finding-the-shortest-path-with-dijkstras-algorithm/]{this} R-bloggers post.
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
#' A vector of states describing one optimal path in the order \code{start} to \code{end}.
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
  
  # THIS WOULD IDEALLY BE REFACTORED TO JUST USE STEPMATRIX DIRECTLY (I.E., WITHOUT CONVERTING TO THE GRAPH FORMAT.)
  
  # If start and end are identical (on diagonal of matrix):
  if (start == end) {
    
    # Just return start and end (cost is zero):
    return(value = c(start, end))
    
  # If start and end differ (not on diagonal of matrix):
  } else {
    
    # Check stepmatrix has class stepMatrix and stop and warn user if not:
    if (!inherits(x = stepmatrix, what = "stepMatrix")) stop("stepmatrix must be an object of class \"stepMatrix\".")
    
    # Create position version of matrix to check for all Inf triangles:
    position_matrix <- matrix(data = 1:stepmatrix$size ^ 2, ncol = stepmatrix$size, dimnames = list(rownames(x = stepmatrix$stepmatrix), colnames(x = stepmatrix$stepmatrix)))
    
    # Find out whether the transition in upper triangle of stepmatrix:
    is_in_upper_triangle <- ifelse(test = any(x = which(x = upper.tri(x = position_matrix)) == position_matrix[start, end]), yes = TRUE, no = FALSE)
    
    # If in upper triangle and upper triangle is all Inf values:
    if (is_in_upper_triangle && all(x = stepmatrix$stepmatrix[upper.tri(x = stepmatrix$stepmatrix)] == Inf)) {
      
      # Return start and end directly (cost is Inf):
      return(value = c(start, end))
      
    # If not upper triangle or upper triangle is not all Inf values:
    } else {
      
      # If in lower triangle and lower triangle is all Inf values:
      if (!is_in_upper_triangle && all(x = stepmatrix$stepmatrix[lower.tri(x = stepmatrix$stepmatrix)] == Inf)) {
        
        # Return start and end directly (cost is Inf):
        return(value = c(start, end))
        
      # If matrix triangle present in has finite values (can proceed to Dijkstra's algorithm):
      } else {
        
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
        
        # Case if Inf edges are found (will need to be removed):
        if (any(x = unlist(x = weights) == Inf)) {
          
          # Remove Inf edges from graph and weights:
          graph <- mapply(FUN = function(graph, weights) graph[weights != Inf], graph = graph, weights = weights, SIMPLIFY = FALSE)
          weights <- mapply(FUN = function(graph, weights) weights[weights != Inf], graph = graph, weights = weights, SIMPLIFY = FALSE)
        }
        
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
          
          # Loop through all nodes linked from the current node (given in start):
          for (node in graph[[start]]) {
            
            # Proceed only if linked node is not already in path:
            if (!(node %in% path)) {
              
              # Recursively call function for finding shortest path with node as start and assign it to newpath:
              newpath <- find_the_shortest_path(graph, node, end, path)
              
              # If newpath is shorter than shortest so far assign newpath to shortest:
              if (path_length(newpath) < path_length(shortest))
              shortest <- newpath
            }
          }
          
          # Return shortest path:
          shortest
        }
        
        # Return shortest path:
        return(value = find_the_shortest_path(graph = graph, start = start, end = end))
      }
    }
  }
}
