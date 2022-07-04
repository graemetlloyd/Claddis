#' Is a graph connected?
#'
#' @description
#'
#' Is a graph represented by an adjacenecy matrix connected?
#'
#' @param adjacency_matrix An adjacency matrix where the diagonal is zeroes and the off-diagonal either ones (if the two vertices are directly connected) or zeroes (if not directly connected).
#'
#' @details
#'
#' Any undirected graph can be represented as an adjacency matrix and the properties of this matrix can be used to determine whether or not the graph is connected (i.e., a path between any two vertices exists) or not.
#'
#' For example, the following graph:
#'
#' \preformatted{6---4---5
#'     |   |\
#'     |   | 1
#'     |   |/
#'     3---2}
#'
#' Has the adjacency matrix:
#'
#' \preformatted{    _________________________
#'     | 1 | 2 | 3 | 4 | 5 | 6 |
#' -----------------------------
#' | 1 | 0 | 1 | 0 | 0 | 1 | 0 |
#' -----------------------------
#' | 2 | 1 | 0 | 1 | 0 | 1 | 0 |
#' -----------------------------
#' | 3 | 0 | 1 | 0 | 1 | 0 | 0 |
#' -----------------------------
#' | 4 | 0 | 0 | 1 | 0 | 1 | 1 |
#' -----------------------------
#' | 5 | 1 | 1 | 0 | 1 | 0 | 0 |
#' -----------------------------
#' | 6 | 0 | 0 | 0 | 1 | 0 | 0 |
#' -----------------------------}
#'
#' This functions run through the following checks in order to confirm the connectivity of the graph:
#'
#' \enumerate{
#'   \item As the graph has more than one vertex then further checks are required (a single vertex graph is considered connected).
#'   \item As no vertice has degree zero (no links) then the graph \emph{may} be connected (further checks are required).
#'   \item As there are more than two vertices the graph may be disconnected (further checks are required).
#'   \item As no missing paths are found (that would separate the graph into one or more disconnected subgraphs) then the graph must be connected.
#' }
#'
#' This ordering means more complex queries are not triggered unless simpler tests do not provide a definitive answer.
#'
#' @return A logical (TRUE or FALSE).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create the connected graph matrix:
#' x <- matrix(
#'   data = c(
#'     0, 1, 0, 0, 1, 0,
#'     1, 0, 1, 0, 1, 0,
#'     0, 1, 0, 1, 0, 0,
#'     0, 0, 1, 0, 1, 1,
#'     1, 1, 0, 1, 0, 0,
#'     0, 0, 0, 1, 0, 0
#'   ),
#'   ncol = 6,
#'   byrow = TRUE
#' )
#'
#' # Check graph is connected:
#' is_graph_connected(adjacency_matrix = x)
#'
#' @export is_graph_connected
is_graph_connected <- function(adjacency_matrix) {
  
  # Record number of vertices for use in later conditionals:
  n_vertices <- dim(x = adjacency_matrix)[1]
  
  # If only one vertex then graph is considered connected (return TRUE)
  if (n_vertices == 1) return(TRUE)
  
  # If any individual vertex has degree zero then graph is not connected (return FALSE):
  if (any(x = apply(X = adjacency_matrix, MARGIN = 1, FUN = sum) == 0)) return(FALSE)
  
  # If only two vertices - and above test passed - then graph must be considered connected (return TRUE)
  if (n_vertices == 2) return(TRUE)
  
  # Set exponent for matrix power (i.e., maximum path length):
  k <- n_vertices - 1
  
  # Make diagonal one:
  diag(adjacency_matrix) <- 1
  
  # Store unadulterated matrix:
  start_matrix <- adjacency_matrix
  
  # Iteratively multiply matrix to return kth power:
  for(i in 1:(k - 1)) adjacency_matrix <- adjacency_matrix %*% start_matrix
  
  # Any zeroes mean matrix is not connected, is all positive then everything is connected:
  return(!any(x = adjacency_matrix == 0))
}
