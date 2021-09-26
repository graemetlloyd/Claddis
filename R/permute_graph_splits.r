#' Permute all ways to split a graph
#'
#' @description
#'
#' Given a graph represented by an adjacency matrix, permutes all ways this graph could be split apart.
#'
#' @param adjacency_matrix A labelled adjacency matrix where the diagonal is zeroes and the off-diagonal either ones (if the two vertices are directly connected) or zeroes (if not directly connected). Labels must match across row names and column names,
#'
#' @details
#'
#' This function takes a connected graph and considers all the ways it \emph{could} be split by removing every possible combination of edges (inclusive of none and all).
#'
#' @return A vector of splits where connected vert(ices) are listed with "+" joining them and disconnected vert(ices) are separated by "|".
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a connected graph matrix:
#' adjacency_matrix <- matrix(
#'   data = c(
#'     0, 1, 0, 0, 1, 0,
#'     1, 0, 1, 0, 1, 0,
#'     0, 1, 0, 1, 0, 0,
#'     0, 0, 1, 0, 1, 1,
#'     1, 1, 0, 1, 0, 0,
#'     0, 0, 0, 1, 0, 0
#'   ),
#'   ncol = 6,
#'   byrow = TRUE,
#'   dimnames = list(
#'     LETTERS[1:6],
#'     LETTERS[1:6]
#'   )
#' )
#'
#' # Check graph is connectedPermute all ways to split graph (remove edges):
#' permute_graph_splits(adjacency_matrix = adjacency_matrix)
#'
#' @export permute_graph_splits
permute_graph_splits <- function(adjacency_matrix) {
  
  # CHECK LABELS ARE PRESENT AND MATCH AND MATRIX IS TRULY ADJACANECY MATRIX (SQUARE, DIAGONAL OF ZERO, MIRRORED ONES).
  
  # Check graph is connected and stop and warn user if not:
  if (!is_graph_connected(adjacency_matrix = adjacency_matrix)) stop("adjacency_matrix must be a connected graph.")
  
  # Create coordinate matrix from adjacency matrix:
  coordinate_matrix <- adjacency_matrix
  
  # Remove half of matrix (as redundant with other half):
  coordinate_matrix[upper.tri(x = coordinate_matrix)] <- 0
  
  # Subfunction to convert positions in a matrix into coordinates:
  pos2coord <- function(pos = NULL, dim.mat = NULL){
    coord <- matrix(NA, nrow = length(pos), ncol = 2)
    coord[,1] <- ((pos - 1)  %% dim.mat[1]) + 1
    coord[,2] <- ((pos - 1) %/% dim.mat[1]) + 1
    coord
  }
  
  # Find edge coordinates:
  edge_coordinates <- pos2coord(pos = which(x = coordinate_matrix == 1), dim.mat = dim(x = coordinate_matrix))
  
  # Make matrix of all "switch" combinations (i.e., whether an edge is included or not):
  switches_matrix <- expand.grid(lapply(X = 1:sum(x = coordinate_matrix), FUN = function(x) c(0, 1)))
  
  # Split into subgraphs and identfy vertice components of each split:
  permuted_splits <- unlist(
    x = lapply(
      X = 1:nrow(x = switches_matrix),
      FUN = function(i) {
        switches_i <- switches_matrix[i, ]
        adjacency_i <- adjacency_matrix
        edges_to_remove <- which(x = switches_i == 0)
        if (length(x = edges_to_remove) > 0) {
          for(j in edges_to_remove) adjacency_i[edge_coordinates[j, 1], edge_coordinates[j, 2]] <- adjacency_i[edge_coordinates[j, 2], edge_coordinates[j, 1]] <- 0
        }
        paste(
          x = sort(
            x = unlist(
              x = lapply(
                X = split_out_subgraphs(adjacency_matrix = adjacency_i),
                FUN = function(j) paste(x = sort(x = colnames(x = j)), collapse = " + ")
              )
            )
          ),
          collapse = " | "
        )
      }
    )
  )
  
  # Sort and remove duplicates:
  permuted_splits <- sort(x = unique(x = permuted_splits))
  
  # Return permuted splits:
  permuted_splits
}
