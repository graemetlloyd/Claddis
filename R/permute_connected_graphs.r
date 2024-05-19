#' Permute all connected graphs
#'
#' @description
#'
#' Given a vertex count, permutes all connected graphs.
#'
#' @param n_vertices The number of vertices to connect.
#'
#' @details
#'
#' For the two vertex case there is only a single connected graph:
#'
#' \preformatted{A---B}
#'
#' (The labels A and B here simply indicate the two vertices and are not a true labelling.)
#'
#' If we add a third vertex, there are two connected graphs:
#'
#' \preformatted{A---B
#'  \ /
#'   C}
#'
#' And:
#'
#' \preformatted{A---B---C}
#'
#' This function permutes all such connected graphs for a given vertex count.
#'
#' Note that the output is in the form of a matrix of edges. For the three vertex case above these would be:
#'
#' \preformatted{     [,1] [,2]
#' [1,] "A"  "B"
#' [2,] "A"  "C"
#' [3,] "B"  "C"}
#'
#' And:
#'
#' \preformatted{     [,1] [,2]
#' [1,] "A"  "B"
#' [2,] "B"  "C"}
#'
#' Again, it is important to note that the labels A, B, and C here are purely "dummy" labels and should not be considered a graph labelling. To use the second graph as an example there are multiple labellings of this graph:
#'
#' \preformatted{A---B---C}
#'
#' And:
#'
#' \preformatted{B---A---C}
#'
#' And:
#'
#' \preformatted{A---C---B}
#'
#' However, these are all isomorphisms of the same unlabelled graph. Only the unique graphs themselves are returned here.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' A list of graphs (matrices of dummy labelled edges).
#'
#' @examples
#'
#' # Generate all connected graphs of four vertices:
#' permute_connected_graphs(n_vertices = 4)
#'
#' @export permute_connected_graphs
permute_connected_graphs <- function(n_vertices) {
  
  # Performs some checks:
  if (!is.numeric(x = n_vertices)) stop("n_vertices must be a number.")
  if (length(x = n_vertices) != 1) stop("n_vertices must be a single value.")
  if (n_vertices < 2) stop("n_vertices must be at least 2.")
  
  # Special case of two vertices (single possible graph) - return single fully conected graph):
  if (n_vertices == 2) return(value = list(data.frame(from = c("A", "B"), to = c("B", "A"), weight = c(1, 1))))
  
  # Generate vertex labels:
  vertices <- make_labels(N = n_vertices)

  # Create empty list to store graphs:
  graphs <- list()
  
  # Start with the fully connected graph:
  fully_connected_graph <- t(x = combn(x = vertices, m = 2))
  
  # Append fully connected graph to list:
  graphs[[1]] <- list(fully_connected_graph)
  
  # For each vertex removed down to the minimum for a connected graph:
  for(i in 1:length(x = (nrow(x = fully_connected_graph) - 1):(n_vertices - 1))) {
    
    # First make new graphs by knocking out one edge from previous graphs:
    new_graphs <- lapply(
      X = graphs[[length(x = graphs)]],
      FUN = function(graph) {
        lapply(
          X = as.list(x = 1:nrow(x = graph)),
          FUN = function(edge_to_remove) {
            graph[-edge_to_remove, , drop = FALSE]
          }
        )
      }
    )
    
    # Make new_graphs into single list:
    new_graphs <- do.call(what = c, args = new_graphs)
    
    # Check whether each new graph is connected:
    graph_is_connected <- unlist(x =
      lapply(
        X = new_graphs,
        FUN = function(new_graph) {
          
          # Make empty adjacency matrix for graph:
          adjacency_matrix <- matrix(
            data = 0,
            ncol = n_vertices,
            nrow = n_vertices,
            dimnames = list(vertices, vertices)
          )
              
          # Populate adjacency matrix
          for(i in 1:nrow(x = new_graph)) adjacency_matrix[new_graph[i, 1], new_graph[i, 2]] <- adjacency_matrix[new_graph[i, 2], new_graph[i, 1]] <- 1
              
          # Check graph is connected using adjacency matrix:
          is_graph_connected(adjacency_matrix = adjacency_matrix)
        }
      )
    )

    # Prune out any unconnected graphs:
    new_graphs <- new_graphs[graph_is_connected]
        
    # Get sorted vertex degrees to use to determine isomorphisms:
    sorted_vertex_degrees <- unlist(
      x = lapply(
        X = new_graphs,
        FUN = function(new_graph) {
          sorted_vector <- sort(x = as.vector(x = new_graph))
          rle_lengths <- rle(x = sorted_vector)$lengths
          paste(sort(x = rle_lengths, decreasing = TRUE), collapse = "-")
        }
      )
    )
        
    # Return just unique new graphs (dump isomorphisms):
    new_graphs <- new_graphs[!duplicated(x = sorted_vertex_degrees)]

    # Append new graphs to graphs:
    graphs[[(length(x = graphs) + 1)]] <- new_graphs
  }

  # Recompile into a single list of graphs:
  graphs <- do.call(what = c, args = graphs)
  
  # Construct as proper graphs (with return edges and dummy weights of one):
  graphs <- lapply(
    X = graphs,
    FUN = function(graph) {
      data.frame(
        from = c(graph[, 1], graph[, 2]),
        to = c(graph[, 2], graph[, 1]),
        weight = rep(x = 1, times = nrow(x = graph) * 2)
      )
    }
  )
  
  # Return (dummy labelled) graphs as a list:
  return(value = graphs)
}
