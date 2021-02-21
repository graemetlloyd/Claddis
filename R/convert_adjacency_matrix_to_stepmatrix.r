#' Converts an adjacency matrix to a stepmatrix
#'
#' @description
#'
#' Takes an adjacency matrix as input and returns the corresponding stepmatrix.
#'
#' @param adjacency_matrix A labelled square matrix with zeroes denoting non-adjacencies and ones denoting adjacencies.
#'
#' @details
#'
#' This function is intended for internal use by \link{convert_state_tree_to_stepmatrix}, but as it also generalizes to solving any graph theory problem - generating a distance matrix corresponding to each shortest path between vertices of a connected graph represented as an adjacency matrix - it is made available explicitly here.
#'
#' The process is best understood with an example. Imagine we have a graph like this:
#'
#' \preformatted{  0-1-2
#'   |
#'   3}
#'
#' I.e., we have four labelled vertices, 0-3, and three edges (connections) between them:
#'
#' \preformatted{  0-1
#'   1-2
#'   0-3}
#'
#' Note: here we assume symmetry, 0-1 = 1-0.
#'
#' Graphs like this can be explicitly captured as adjacency matrices, where a one denotes two vertices are "adjacent" (connected by an edge) and a zero that they are not.
#'
#' \preformatted{ |0|1|2|3|
#' 0|0|1|0|1|
#' 1|1|0|1|0|
#' 2|0|1|0|0|
#' 3|1|0|0|0|}
#'
#' But what such matrices do not tell us is how far every vertex-to-vertex path is in terms of edge counts. E.g., the path length from vertex 3 to vertex 2.
#'
#' This function simply takes the adjacency matrix and returns the corresponding stepmatrix, corresponding to every minimum vertex-to-vertex path length.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' A stepmatrix where values indicate the minimum number of edges (shortest path) between each pair of vertices.
#'
#' @seealso
#'
#' \link{convert_state_tree_to_stepmatrix}, link{locate_bracket_positions}
#'
#' @examples
#'
#' # Build the example adjacency matrix for the graph above:
#' adjacency_matrix <- matrix(
#'   data = c(0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0),
#'   nrow = 4,
#'   ncol = 4,
#'   dimnames = list(0:3, 0:3)
#' )
#'
#' # Convert this to a step matrix:
#' convert_adjacency_matrix_to_stepmatrix(
#'   adjacency_matrix = adjacency_matrix
#' )
#'
#' @export convert_adjacency_matrix_to_stepmatrix
convert_adjacency_matrix_to_stepmatrix <- function(adjacency_matrix) {
  
  # TO DO:
  #
  # - Check that matrix is connected somehow (cannot solve all distances if not)
  
  # Set initial stepmatrix as adjacency matrix:
  stepmatrix <- adjacency_matrix
  
  # Set all zero values to NA (steps that need to be calculated):
  stepmatrix[stepmatrix == 0] <- NA
  
  # Set diagional as zero (the distance from any value to itself):
  diag(stepmatrix) <- 0
  
  # Only continue as long as there are unknown path lengths:
  while(any(x = is.na(x = stepmatrix))) {
    
    # Start with first path of unknown length:
    current_path <- which(x = is.na(x = stepmatrix), arr.ind = TRUE)[1, ]
    
    # Isolate path start:
    path_start <- colnames(x = adjacency_matrix)[current_path][1]
    
    # Isolate path end:
    path_end <- colnames(x = adjacency_matrix)[current_path][2]
    
    # Little subfunction to find which verti(ces) current vertex is connected to:
    linked_to <- function(node, adjacency_matrix) names(x = which(x = adjacency_matrix[node, ] == 1))
    
    # Populate path list with connections to current verti(ces):
    path_list <- list(linked_to(node = colnames(adjacency_matrix)[current_path[1]], adjacency_matrix = adjacency_matrix))
    
    # As long as the path end has not been reached:
    while(!any(x = path_list[[length(x = path_list)]] == path_end)) {
      
      # Add connections from current verti(ces):
      path_list[[length(x = path_list) + 1]] <- unique(x = unlist(x = lapply(X = as.list(x = path_list[[length(x = path_list)]]), FUN = linked_to, adjacency_matrix = adjacency_matrix)))
      
    }
    
    # Update stepmatrix with path length (minimum length to connect path_start to path_end):
    stepmatrix[current_path[1], current_path[2]] <- stepmatrix[current_path[2], current_path[1]] <- length(path_list)
    
  }
  
  # Return complete stepmatrix:
  stepmatrix
  
}
