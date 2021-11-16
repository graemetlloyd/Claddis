#' Split adjacency matrix into subgraphs
#'
#' @description
#'
#' Given a graph represented by an adjacency matrix splits into all connected subgraphs.
#'
#' @param adjacency_matrix An adjacency matrix where the diagonal is zeroes and the off-diagonal either ones (if the two vertices are directly connected) or zeroes (if not directly connected).
#'
#' @details
#'
#' This functions take any undirected graph (connected or unconnected) represented as an adjacency matrix and identifies all connected subgraphs and returns these as a list of adjacency matr(ices).
#'
#' @return A list of all connected subgraphs represented as adjacency matri(ces).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create an adjacency matrix representing an unconnected graph:
#' adjacency_matrix <- matrix(
#'   data = c(
#'     0, 0, 0, 1, 1, 0,
#'     0, 0, 1, 0, 0, 1,
#'     0, 1, 0, 0, 0, 1,
#'     1, 0, 0, 0, 0, 0,
#'     1, 0, 0, 0, 0, 0,
#'     0, 1, 1, 0, 0, 0
#'   ),
#'   ncol = 6,
#'   byrow = TRUE,
#'   dimnames = list(LETTERS[1:6], LETTERS[1:6])
#' )
#'
#' # Check graph is connected:
#' split_out_subgraphs(adjacency_matrix = adjacency_matrix)
#'
#' @export split_out_subgraphs
split_out_subgraphs <- function(adjacency_matrix) {
  
  # Create empty list to store subgraphs:
  subgraphs <- list()
  
  # Identify any isolated vert(ices) (single vertex subgraphs):
  isolated_vertices <- names(x = which(x = apply(X = adjacency_matrix, MARGIN = 2, FUN = sum) == 0))
  
  # If isolated vert(ices) exist:
  if (length(x = isolated_vertices) > 0) {
    
    # Add these as sepaarte subgrahs to the su graphs output:
    subgraphs <- c(subgraphs, lapply(X = as.list(x = isolated_vertices), FUN = function(i) adjacency_matrix[i, i, drop = FALSE]))
    
    # Identify remaining vert(ices):
    remaining_vertices <- setdiff(x = rownames(x = adjacency_matrix), y = isolated_vertices)
    
    # If there are remaining vert(ices):
    if (length(x = remaining_vertices) > 0) {
      
      # Collapse adjacency_matrix to just remaining vert(ices):
      adjacency_matrix <- adjacency_matrix[remaining_vertices, remaining_vertices, drop = FALSE]
      
      # If there are no remaning vert(ices):
    } else {
      
      # Can return subgraphs as already done:
      return(subgraphs)
    }
  }
  
  # If remaining adjacency_matrix is a connected graph:
  if (is_graph_connected(adjacency_matrix = adjacency_matrix)) {
    
    # Add adjacency_matrix to sugraphs:
    subgraphs[[(length(x = subgraphs) + 1)]] <- adjacency_matrix
    
    # Return subgraphs as already done:
    return(subgraphs)
  }
  
  # Create power_matrix from adjacency_matrix:
  power_matrix <- adjacency_matrix
  
  # Set power_matrix diagonal as one:
  diag(x = power_matrix) <- 1
  
  # Create clean power_matrix for exponenting:
  clean_power_matrix <- power_matrix
  
  # Get kth power of power_matrix (any zeroes left indicate unconnected subgraphs):
  for(i in 1:(dim(x = power_matrix)[1] - 2)) power_matrix <- power_matrix %*% clean_power_matrix
  
  # As long as there are still subgraphs to split off:
  while(length(x = power_matrix) > 0) {
    
    # Get vertices of first remaining subgraph:
    first_subgraph_vertices <- names(x = which(x = power_matrix[1, ] > 0))
    
    # Store this subgraph in subgraph:
    subgraphs[[(length(x = subgraphs) + 1)]] <- adjacency_matrix[first_subgraph_vertices, first_subgraph_vertices, drop = FALSE]
    
    # Identify remaining vertices:
    remaining_vertices <- setdiff(x = rownames(x = power_matrix), y = first_subgraph_vertices)
    
    # Remove this sbgraph from the power matrix:
    power_matrix <- power_matrix[remaining_vertices, remaining_vertices, drop = FALSE]
  }
  
  # Return complete list of subgraphs:
  return(subgraphs)
}
