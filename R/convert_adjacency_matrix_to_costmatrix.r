#' Converts an adjacency matrix to a costmatrix
#'
#' @description
#'
#' Takes an adjacency matrix as input and returns the corresponding costmatrix.
#'
#' @param adjacency_matrix A labelled square matrix with zeroes denoting non-adjacencies and ones denoting adjacencies.
#'
#' @details
#'
#' This function is intended for internal use, but as it also generalizes to solving a general graph theory problem - generating a distance matrix corresponding to each shortest path between vertices of a connected graph represented as an adjacency matrix - it is made available explicitly here.
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
#' \preformatted{    _________________
#'     | 0 | 1 | 2 | 3 |
#' ---------------------
#' | 0 | 0 | 1 | 0 | 1 |
#' ---------------------
#' | 1 | 1 | 0 | 1 | 0 |
#' ---------------------
#' | 2 | 0 | 1 | 0 | 0 |
#' ---------------------
#' | 3 | 1 | 0 | 0 | 0 |
#' ---------------------}
#'
#' But what such matrices do not tell us is how far every vertex-to-vertex path is in terms of edge counts. E.g., the path length from vertex 3 to vertex 2.
#'
#' This function simply takes the adjacency matrix and returns the corresponding costmatrix, corresponding to every minimum vertex-to-vertex path length.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' An object of class \code{costMatrix}.
#'
#' @seealso
#'
#' \link{convert_state_tree_to_adjacency_matrix}, \link{locate_bracket_positions}
#'
#' @examples
#'
#' # Build the example adjacency matrix for the graph above:
#' adjacency_matrix <- matrix(
#'   data = c(
#'     0, 1, 0, 1,
#'     1, 0, 1, 0,
#'     0, 1, 0, 0,
#'     1, 0, 0, 0
#'   ),
#'   nrow = 4,
#'   ncol = 4,
#'   dimnames = list(0:3, 0:3)
#' )
#'
#' # Convert this to a costmatrix:
#' convert_adjacency_matrix_to_costmatrix(
#'   adjacency_matrix = adjacency_matrix
#' )
#'
#' @export convert_adjacency_matrix_to_costmatrix
convert_adjacency_matrix_to_costmatrix <- function(adjacency_matrix) {
  
  # Check adjacency_matrix is symmetric and stop and warn user if not:
  if (!isSymmetric(object = adjacency_matrix)) stop("adjacency_matrix must be symmetric. Fix and try again.")
  
  # Check adjacency_matrix is for a connected graph:
  if (!is_graph_connected(adjacency_matrix = adjacency_matrix)) stop("adjacency_matrix must represent a connected graph. Fix and try again")
  
  # Check for polymorphisms and uncertainties and remove if found:
  if (length(x = grep(pattern = "&", x = colnames(x = adjacency_matrix)))) stop("Polymorphisms (&) are not permitted by this function.")
  if (length(x = grep(pattern = "/", x = colnames(x = adjacency_matrix)))) stop("Uncertainties (/) are not permitted by this function.")

  # Calculate size of adjacency matrix:
  matrix_size <- nrow(x = adjacency_matrix)
  
  # Set initial costmatrix as adjacency matrix:
  costmatrix <- adjacency_matrix
  
  # Set all zero values to NA (costs that need to be calculated):
  costmatrix[costmatrix == 0] <- NA
  
  # Set diagonal as zero (the distance from any value to itself):
  diag(costmatrix) <- 0
  
  # Only continue as long as there are unknown path lengths:
  while(any(x = is.na(x = costmatrix))) {
    
    # Start with first path of unknown length:
    current_path <- which(x = is.na(x = costmatrix), arr.ind = TRUE)[1, ]
    
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
    
    # Update costmatrix with path length (minimum length to connect path_start to path_end):
    costmatrix[current_path[1], current_path[2]] <- costmatrix[current_path[2], current_path[1]] <- length(path_list)
  }
  
  # Initialise character type as custom:
  character_type <- "custom"
  
  # If adjacency matrix matches an unordered character (everything is adjacent) then store character type as such:
  if (all(x = as.dist(m = adjacency_matrix) == 1)) character_type <- "unordered"
  
  # If adjacency matrix matches an ordered character (maximum vertex degree is 2) then store character type as such:
  if (all(adjacency_matrix[c(2, cumsum(x = rep(x = matrix_size + 1, length.out = matrix_size - 2)) + 2)] == 1) && sum(adjacency_matrix) == ((2 * matrix_size) - 2)) character_type <- "ordered"
  
  # Format as a costmatrix object:
  costmatrix <- list(
    size = ncol(x = costmatrix),
    n_states = ncol(x = costmatrix),
    single_states = colnames(x = costmatrix),
    type = character_type,
    costmatrix = costmatrix,
    symmetry = "Symmetric",
    includes_polymorphisms = FALSE,
    polymorphism_costs = "additive",
    polymorphism_geometry = "simplex",
    polymorphism_distance = "euclidean",
    includes_uncertainties = FALSE,
    pruned = FALSE,
    dollo_penalty = 99,
    base_age = 100,
    weight = 1
  )
  
  # Set class of output as costMatrix:
  class(costmatrix) <- "costMatrix"
  
  # Return complete costmatrix:
  costmatrix
}
