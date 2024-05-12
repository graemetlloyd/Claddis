#' Finds a minimum spanning tree of a stategraph
#'
#' @description
#'
#' Given a stategraph, returns a shortest tree connecting every state.
#'
#' @param stategraph An object of class \code{stateGraph}.
#'
#' @details
#'
#' The minimum parsimony length a phylogenetic hypothesis can have depends on both the stategraph of transition costs and the states actually sampled. If the stategraph vertices already represent sampled states (the assumption here) then this minimum length is reduced to a graph theory problem - the minimum spanning tree or, if a directed graph, the minimum weight spanning arboresence. (NB: if there are unsampled states then the \code{find_steiner_tree_of_stategraph} function should be used instead.) This function returns one such shortest tree (although others may exist). The sum of the weights of the edges or arcs returned is the minimum cost.
#'
#' As the algorithms used are graph theory based the function operates by simply calling \link{convert_costmatrix_to_stategraph} and \link{find_stategraph_minimum_span}. In practice, if the costmatrix represents a graph (transition costs are all symmetric) then Kruskal's algorithm is applied (Kruskal 1956). If costs are asymmetric, however, then the graph representation is a directed graph (or digraph) and so a version of Edmonds' algorithm is applied (Edmonds 1967).
#'
#' Note that Dollo characters represent a special case solution as although a penalty weight is applied to the edges intended to only ever be traversed once this weight should not be used when calculating tree lengths. The function catches this and returns the edges with the weight that would actually be counted for a minimum weight spanning tree.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Edmonds, J., 1967. Optimum branchings. \emph{Journal of Research of the National Bureau of Standards Section B}, \bold{71B}, 233-240.
#'
#' Kruskal, J. B., 1956. On the shortest spanning subtree of a graph and the traveling salesman problem. \emph{Proceedings of the American Mathematical Society}, \bold{7}, 48-50.
#'
#' @return
#'
#' A \code{data.frame} object describing a minimum spanning tree or minimum weight arboresence as a series of edges or arcs.
#'
#' @seealso
#'
#' \link{find_shortest_costmatrix_path}
#'
#' @examples
#'
#' # Make a four-state ordered character stategraph:
#' ordered_costmatrix <- make_costmatrix(
#'   min_state = "0",
#'   max_state = "3",
#'   character_type = "ordered"
#' )
#' ordered_stategraph <- convert_costmatrix_to_stategraph(costmatrix = ordered_costmatrix)
#'
#' # Find length of shortest spanning tree of stategraph:
#' find_stategraph_minimum_span(stategraph = ordered_stategraph)
#'
#' # Make a four-state unordered character stategraph:
#' unordered_costmatrix <- make_costmatrix(
#'   min_state = "0",
#'   max_state = "3",
#'   character_type = "unordered"
#' )
#' unordered_stategraph <- convert_costmatrix_to_stategraph(costmatrix = unordered_costmatrix)
#'
#' # Find length of shortest spanning tree of stategraph:
#' find_stategraph_minimum_span(stategraph = unordered_stategraph)
#'
#' # Make a four-state irreversible character stategraph:
#' irreversible_costmatrix <- make_costmatrix(
#'   min_state = "0",
#'   max_state = "3",
#'   character_type = "irreversible"
#' )
#' irreversible_stategraph <- convert_costmatrix_to_stategraph(costmatrix = irreversible_costmatrix)
#'
#' # Find length of shortest spanning tree of stategraph:
#' find_stategraph_minimum_span(stategraph = irreversible_stategraph)
#'
#' # Make a four-state Dollo character stategraph:
#' dollo_costmatrix <- make_costmatrix(
#'   min_state = "0",
#'   max_state = "3",
#'   character_type = "dollo"
#' )
#' dollo_stategraph <- convert_costmatrix_to_stategraph(costmatrix = dollo_costmatrix)
#'
#' # Find length of shortest spanning tree of stategraph:
#' find_stategraph_minimum_span(stategraph = dollo_stategraph)
#'
#' @export find_costmatrix_minimum_span
find_stategraph_minimum_span <- function(stategraph) {
  
  ### CHECK INPUT!!!
  
  # Special case of no edges (return empty matrix of edges):
  if (stategraph$n_arcs == 0) return(data.frame(from = "0", to = "1", weight = 1)[-1, ])
  
  # If the special case of a Dollo character:
  if (stategraph$type == "dollo") {
      
    # Make Dollo arcs from those that should be traversed just once:
    dollo_arcs <- stategraph$arcs[stategraph$arcs$weight == stategraph$dollo_penalty, ]
      
    # Set the weights of these arcs to one as this is what should actually be counted for a Dollo character:
    dollo_arcs$weight <- 1
      
    # Return Dollo arcs:
    return(value = dollo_arcs)
  }
  
  # If undirected case (apply Kruskal's algorithm):
  if (!stategraph$directed) {
      
    # Treat paired arcs as single weighted edges:
    weighted_edges <- stategraph$arcs[as.numeric(x = stategraph$arcs$from) < as.numeric(x = stategraph$arcs$to), ]
    
    # Sort from lowest to highest weight:
    weighted_edges <- weighted_edges[order(x = weighted_edges$weight), ]
      
    # Create starting (empty) mst edge list:
    mst_edges <- weighted_edges[-(1:nrow(x = weighted_edges)), ]
      
    # Create starting forest (single vertices for now):
    forest <- as.list(x = unique(x = c(weighted_edges$from, weighted_edges$to)))
      
    # As long as there are unconnected vertices:
    while(length(x = forest) > 1) {
        
      # Consider first (shortest) edge for consideration:
      edge_to_consider <- weighted_edges[1, ]
        
      # Find which tree(s) in forest are connected by edge:
      trees_connected_by_edge <- which(
        x = unlist(
          x = lapply(
            X = forest,
            FUN = function(tree) {
              any(
                c(
                  tree == edge_to_consider$from,
                  tree == edge_to_consider$to
                )
              )
            }
          )
        )
      )
        
      # Case if a new connection is made:
      if (length(x = trees_connected_by_edge) == 2) {
          
        # Add connected edge to forest:
        forest <- c(forest, list(unlist(x = forest[trees_connected_by_edge])))
          
        # Now remove unconnected versions from forest:
        forest <- forest[-trees_connected_by_edge]
          
        # Add edge to MST:
        mst_edges <- rbind(mst_edges, edge_to_consider)
      }
        
      # Now drop edge from list:
      weighted_edges <- weighted_edges[-1, ]
    }
      
    # Return MST edges:
    return(value = mst_edges)
  }
    
  # If directed case (Edmonds algorithm):
  if (stategraph$directed) {
    
    # First establish root options (must be able to reach all other vertices):
    root_options <- names(
      x = which(
        x = apply(
          X = convert_stategraph_to_costmatrix(stategraph = stategraph)$costmatrix < Inf,
          MARGIN = 1,
          FUN = all
        )
      )
    )
  
    # Little error check that should in theory never be hit if a valid stategraph was supplied:
    if (length(x = root_options) == 0) stop("No minimum-weight spanning arboresence is possible as no vertex exists that can reach all other vertices.")
    
    # Next reweight arcs (minimum arriving cost at vertex is zero):
    reweighted_arcs <- do.call(
      what = rbind,
      args = lapply(
        X = as.list(x = unique(x = stategraph$arcs$to)),
        FUN = function(unique_vertex) {
          arcs_arriving_at_unique_vertex <- stategraph$arcs[stategraph$arcs$to == unique_vertex, ]
          minimum_weight <- min(x = arcs_arriving_at_unique_vertex[, "weight"])
          arcs_arriving_at_unique_vertex$weight <- arcs_arriving_at_unique_vertex$weight - minimum_weight
          arcs_arriving_at_unique_vertex
        }
      )
    )
    
    # Generate a single shortest arboresence for each root option:
    shortest_arboresences <- lapply(
      X = as.list(x = root_options),
      FUN = function(root) {
          
        # First remove any arcs arriving at the root (as cannot be part of arboresence and hence are useless):
        if (any(x = reweighted_arcs$to == root)) reweighted_arcs <- reweighted_arcs[-which(x = reweighted_arcs$to == root), ]
          
        # Next establish what vertices need connecting:
        non_root_vertices <- setdiff(x = stategraph$vertices$label, y = root)
          
        # Permute all possible arboresences by considering every choice for an incoming arc of weight zero:
        permuted_arboresences <- expand.grid(
          lapply(
            X = as.list(x = non_root_vertices),
            FUN = function(non_root_vertex) {
              arriving_arcs <- reweighted_arcs[reweighted_arcs$to == non_root_vertex, ]
              zero_weight_arriving_arcs <- arriving_arcs[arriving_arcs$weight == 0, ]
              apply(X = zero_weight_arriving_arcs, MARGIN = 1, FUN = paste, collapse = "%")
            }
          )
        )
          
        # Reformat as proper from, to, weight data frames:
        permuted_arboresences <- apply(
          X = permuted_arboresences,
          MARGIN = 1,
          FUN = function(x) {
            possible_arboresence <- unname(obj = do.call(what = rbind, args = strsplit(x = x, split = "%")))
            data.frame(
              from = possible_arboresence[, 1],
              to = possible_arboresence[, 2],
              weight = possible_arboresence[, 3]
            )
          },
          simplify = FALSE
        )
          
        # Check for connectedness which by definition removes cycles meaning results are arboresences:
        actual_arboresences <- permuted_arboresences[which(
          x = unlist(
            x = lapply(
              X = permuted_arboresences,
              FUN = function(x) {
                length(
                  x = intersect(
                    x = unique(x = c(x$from, x$to)),
                    y = c(root, non_root_vertices)
                  )
                ) == length(
                  x = c(root, non_root_vertices)
                )
              }
            )
          )
        )]
        
        # Add in original weights:
        original_weight_arboresences <- lapply(
          X = actual_arboresences,
          FUN = function(x) {
            x$weight <- unlist(
              x = lapply(
                X = as.list(x = 1:nrow(x = x)),
                FUN = function(y) {
                  stategraph$arcs$weight[intersect(
                    x = which(x = stategraph$arcs$from == x$from[y]),
                    y = which(x = stategraph$arcs$to == x$to[y])
                  )]
                }
              )
            )
            x
          }
        )
        
        # As long as there are arboresences:
        if (length(x = original_weight_arboresences) > 0) {
          
          # Get total weight of each arboresence:
          arboresence_lengths <- unlist(x = lapply(X = original_weight_arboresences, FUN = function(x) sum(x = x$weight)))
          
          # Return first lowest weight arboresence:
          return(value = original_weight_arboresences[[which(x = arboresence_lengths == min(x = arboresence_lengths))[1]]])
        
        # If there re no arboresences:
        } else {
          
          # Create dummy infinite length arboresence:
          return(value = data.frame(from = root, to = non_root_vertices, weight = Inf))
        }
      }
    )
    
    # Little error check (shouldn't ever happen in practice!):
    if (all(unlist(x = lapply(X = shortest_arboresences, FUN = function(x) sum(x = x$weight))) == Inf)) stop("No valid shortest arboresences.")
    
    # Get total lengths of each root option:
    shortest_arboresence_lengths <- unlist(x = lapply(X = shortest_arboresences, FUN = function(x) sum(x = x$weight)))
    
    # Return first shortest arboresence:
    return(value = shortest_arboresences[[which(x = shortest_arboresence_lengths == min(x = shortest_arboresence_lengths))[1]]])
  }
}
