#' Get edges of minimum spanning tree
#'
#' @description
#'
#' Returns edges of a minimum spanning tree given a distance matrix.
#'
#' @param distance_matrix A square matrix of distances between objects.
#'
#' @details
#'
#' This function is a wrapper for \link{mst} in the \link{ape} package, but returns a vector of edges rather than a square matrix of links.
#'
#' @return A vector of named edges (X->Y) with their distances. The sum of this vector is the length of the minimum spanning tree.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a simple square matrix of distances:
#' distance_matrix <- matrix(c(0, 1, 2, 3, 1, 0, 1, 2, 2, 1, 0, 1, 3, 2, 1, 0),
#'   nrow = 4,
#'   dimnames = list(LETTERS[1:4], LETTERS[1:4])
#' )
#'
#' # Show matrix to confirm that the off diagonal has the shortest
#' # distances:
#' distance_matrix
#'
#' # Use find_minimum_spanning_edges to get the edges for the minimum spanning
#' # tree:
#' find_minimum_spanning_edges(distance_matrix)
#'
#' # Use sum of find_minimum_spanning_edges to get the length of the minimum
#' # spanning tree:
#' sum(find_minimum_spanning_edges(distance_matrix))
#' @export find_minimum_spanning_edges
find_minimum_spanning_edges <- function(distance_matrix) {

  # Convert to matrix and set up links matrix:
  distance_matrix <- as.matrix(distance_matrix)

  # Get links matrix for minimum spanning tree:
  links_matrix <- ape::mst(distance_matrix)

  # Create empty matrix to store edges for minimum spanning tree:
  minimum_spanning_tree_edges <- matrix(nrow = 0, ncol = 2, dimnames = list(c(), c("From", "To")))

  # For each row:
  for (i in 1:(nrow(links_matrix) - 1)) {

    # For each column:
    for (j in (i + 1):ncol(links_matrix)) {

      # If there is a link then record it:
      if (links_matrix[i, j] == 1) {
        minimum_spanning_tree_edges <- rbind(
          minimum_spanning_tree_edges,
          c(
            rownames(x = links_matrix)[i],
            colnames(x = links_matrix)[j]
          )
        )
      }
    }
  }

  # Get distances:
  distances <- diag(x = distance_matrix[minimum_spanning_tree_edges[, "From"], minimum_spanning_tree_edges[, "To"]])

  # Add names to distances:
  names(distances) <- apply(minimum_spanning_tree_edges, 1, paste, collapse = "->")

  # Return distances for minimum spanning tree:
  distances
}
