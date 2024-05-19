#' Permute costmatrices
#'
#' @description
#'
#' Given vectors of states and costs, permutes all possible costmatrices.
#'
#' @param states A vector of character states, e.g., "0", "1", "2".
#' @param costs A vector of numeric costs, e.g., 1, 2, Inf.
#' @param symmetry Must be one of \code{"symmetric"}, \code{"asymmetric"} or \code{"both"}.
#'
#' @details
#'
#' Costmatrices define the cost of each state-to-state transition, but they are restricted in what these costs can be (see \link{check_costMatrix}). Nevertheless, strictly speaking there are infinite possible costmatrices - even where costs are restricted to integer values (as TNT does; Goloboff et al. 2008; Goloboff and Catalano 2016), i.e., "stepmatrices" (Swofford and Maddison 1992). Thus this function operates on a finite system by requiring the user to specify a restricted set of states and individual cost values, with the function permuting every possible combination of finite costs. Note that not \emph{every} permutation will be returned as not all of these will be valid costmatrices (see \link{check_costMatrix} and \link{fix_costmatrix}). Others will not be returned because their cost \emph{ratio} can be considered redundant. For example, for a binary character (states "0", and "1") the following two costmatrices would be mutually redundant as the ratio of their costs is identical:
#'
#' \preformatted{  A B
#' A 0 1
#' B 2 0
#'
#'   A B
#' A 0 2
#' B 4 0}
#'
#' (If the user does want to consider these kinds of alternatives then a better solution is to simply weight the first matrix by two, or any other value, in any downstream analys(es).)
#'
#' For the function to work costs must be unique positive values. This includes infinity (\code{Inf} in R). Infinite costs can be used to denote a particular transition is impossible and allows defining (e.g.) irreversible characters, or those that force a particular root value.
#'
#' @return A list of unique costmatrices containing every possible combination of costs.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics/ \emph{Cladistics}, \bold{32}. 221-238
#'
#' Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. \emph{Cladistics}, \bold{24}, 774-786.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. \emph{In} R. L. Mayden (ed.) Systematics, Historical Ecology, and North American Freshwater Fishes. Stanford University Press, Stanford, p187-223.
#'
#' @examples
#'
#' # Permute all the ways to assign the costs 1 and 2 for a three state
#' # character:
#' permute_costmatrices(
#'   states = c("0", "1", "2"),
#'   costs = c(1, 2),
#'   symmetry = "both"
#' )
#'
#' @export permute_costmatrices
permute_costmatrices <- function(
  states = c("0", "1"),
  costs = c(1:3),
  symmetry = "both"
) {
  
  ### ADD CHECKS FOR SYMMETRY VARIABLE!!!!
  
  # Check states has positive length and stop and warn user if not:
  if (length(x = states) < 1) stop("states must have at least one value.")
  
  # Check states are in the form of characters and stop and warn user if not:
  if (!is.character(states)) stop("states must be a character vector. If they are 0, 1 etc. then place in quotes or use as.character().")
  
  # Check there are no duplicate states and stop and warn user if not:
  if (any(duplicated(x = states))) stop("states cannot include duplicate values.")
  
  # Check costs has positive length and stop and warn user if not:
  if (length(x = costs) < 1) stop("costs must have at least one value.")
  
  # Check costs are in the form of numbers and stop and warn user if not:
  if (!is.numeric(x = costs)) stop("costs must be a numeric vector.")
  
  # Check costs are all positive values and stop and war user if not:
  if (any(costs <= 0)) stop("costs must all be positive values.")
  
  # Check there are no duplicate costs and stop and warn user if not:
  if (any(duplicated(x = costs))) stop("costs cannot include duplicate values.")
  
  # Check costs is not just infinity and stop and warn user if not:
  if (all(costs == Inf)) stop("costs must include at least one non-infinite value.")
  
  # If only one state just return a single zero costmatrix:
  if (length(x = states) == 1) {
    return(
      list(
        make_costmatrix(
          min_state = states,
          max_state = states,
          character_type = "ordered"
        )
      )
    )
  }

  # Get number of vertices:
  n_vertices <- length(x = states)
  
  # If symmetric costmatries (undirected graphs) are requested:
  if (symmetry == "both" || symmetry == "symmetric") {
    
    # Cannot have infinite costs for symmetric graph as would not be connected:
    undirected_costs <- setdiff(x = costs, y = Inf)
    
    # Set directed costs length:
    n_undirected_costs <- length(x = undirected_costs)
    
    # First generate all possible connected graphs of n vertices:
    state_graphs <- permute_connected_graphs(n_vertices = n_vertices)
  
    # Now generate all possible edge weightings of these graphs using costs:
    undirected_graphs <- lapply(
      X = state_graphs,
      FUN = function(state_graph) {
      
        # Prune returning edges and weights from state graph:
        state_graph <- state_graph[1:(nrow(x = state_graph) / 2), c("from", "to")]
      
        # Create a matrix of costs (will be converted to a list for permutation below):
        costs_matrix <- matrix(
          data = rep(x = undirected_costs, times = nrow(x = state_graph)),
          ncol = n_undirected_costs,
          byrow = TRUE
        )
      
        # Permute all combos of edge weights (many will be redundant - this will be checked below):
        edge_weights <- t(
          expand.grid(
            apply(
              X = costs_matrix,
              MARGIN = 1,
              FUN = function(i) i,
              simplify = FALSE
            )
          )
        )
      
        # Get rid of annoying names:
        edge_weights <- unname(obj = edge_weights)
      
        # Make new state graphs using edge_weights:
        new_state_graphs <- apply(
          X = edge_weights,
          MARGIN = 2,
          function(weights) {
            data.frame(
              from = state_graph[, 1],
              to = state_graph[, 2],
              weight = weights
            )
          }
        )
      
        # Get graph encodings (vertex degrees plus weights for each edge, sorted):
        new_graph_encodings <- unlist(
          x = lapply(
            X = new_state_graphs,
            FUN = function(new_state_graph) {
              vertex_rle <- rle(x = sort(x = c(new_state_graph[, "from"], new_state_graph[, "to"])))
              vertex_degrees <- vertex_rle$lengths
              names(x = vertex_degrees) <- vertex_rle$values
              degrees_and_weights_matrix <- cbind(
                from = vertex_degrees[new_state_graph[, "from"]],
                weight = new_state_graph[, "weight"],
                to = vertex_degrees[new_state_graph[, "to"]]
              )
              edges_vector <- apply(
                X = degrees_and_weights_matrix,
                MARGIN = 1,
                FUN = paste,
                collapse = "-"
              )
              edges_vector <- unname(obj = edges_vector)
          
              # Return graph encoding:
              paste(sort(x = edges_vector), collapse = "%")
            }
          )
        )
      
        # Prune out duplicates to return just unique new state graphs:
        new_state_graphs <- new_state_graphs[!duplicated(x = new_graph_encodings)]
      }
    )
    
    # Remove any identical ratio graphs:
    undirected_graphs <- lapply(
      X = undirected_graphs,
      FUN = function(graph) {
        weight_sets <- do.call(
          what = rbind,
          args = lapply(
            X = graph,
            FUN = function(weight_sets) weight_sets[, "weight"] / min(x = weight_sets[, "weight"])
          )
        )
        ratio_costs <- apply(
          X = weight_sets,
          MARGIN = 1,
          FUN = paste,
          collapse = "-"
        )
        graph[!duplicated(x = ratio_costs)]
      }
    )
    
    # Recompile as a single list:
    undirected_graphs <- do.call(what = c, args = undirected_graphs)
    
    # Identify current graph labels (to replace with states):
    names(x = states) <- unique(
      x = sort(
        x = c(
          undirected_graphs[[1]][, "from"],
          undirected_graphs[[1]][, "to"]
        )
      )
    )
    
    # Use states to relabel graphs:
    undirected_graphs <- lapply(
      X = undirected_graphs,
      FUN = function(graph) {
        data.frame(
          from = unname(obj = states[graph[, "from"]]),
          to = unname(obj = states[graph[, "to"]]),
          weight = graph[, "weight"]
        )
      }
    )
    
    # Add complementary edges to make symmetric graph
    undirected_graphs <- lapply(
      X = undirected_graphs,
      FUN = function(graph) {
        data.frame(
          from = c(graph[, "from"], graph[, "to"]),
          to = c(graph[, "to"], graph[, "from"]),
          weight = c(graph[, "weight"], graph[, "weight"])
        )
      }
    )
    
    # Convert simple stategraphs into proper stategraph class:
    undirected_graphs <- lapply(
      X = undirected_graphs,
      FUN = function(stategraph) {
        
        # Make vertices:
        vertices <- data.frame(
          label = states,
          in_degree = 0,
          out_degree = 0,
          eccentricity = 0,
          periphery = 0,
          centre = 0
        )
        for(i in states) {
          vertices[vertices[, "label"] == i, "in_degree"] <- sum(x = stategraph[, "to"] == i)
          vertices[vertices[, "label"] == i, "out_degree"] <- sum(x = stategraph[, "from"] == i)
        }
        
        # Make adjacency matrix for graph:
        adjacency_matrix <- matrix(
          data = 0,
          nrow = n_vertices,
          ncol = n_vertices,
          dimnames = list(states, states)
        )
        for(i in 1:nrow(x = stategraph)) adjacency_matrix[stategraph[i, "from"], stategraph[i, "to"]] <- 1
        
        # Compile stategraph:
        stategraph <- list(
          n_vertices = n_vertices,
          n_arcs = nrow(x = stategraph),
          n_states = n_vertices,
          single_states = states,
          type = "custom",
          arcs = stategraph,
          vertices = vertices,
          radius = min(x = vertices[, "eccentricity"]),
          diameter = max(x = vertices[, "eccentricity"]),
          adjacency_matrix = adjacency_matrix,
          directed = FALSE,
          includes_polymorphisms = FALSE,
          polymorphism_costs = "additive",
          polymorphism_geometry = "simplex",
          polymorphism_distance = "euclidean",
          includes_uncertainties = FALSE,
          pruned = FALSE,
          dollo_penalty = 999,
          base_age = 1,
          weight = 1
        )
        
        # Make class stategraph:
        class(x = stategraph) <- "stateGraph"

        # Return now formatted stategraph:
        stategraph
      }
    )
    
    # Convert stategraph(s) to costmatri(ces):
    symmetric_costmatrices <- lapply(
      X = undirected_graphs,
      FUN = convert_stategraph_to_costmatrix
    )
    
    # If only want symmetric than can just return those now:
    if (symmetry == "symmetric") return(value = symmetric_costmatrices)
  }
  
  # If asymmetric costmatries (directed graphs) are requested:
  if (symmetry == "both" || symmetry == "asymmetric") {
    
    # Check for infinite costs which cannot currently be dealt with (they create a complex combinatoric problem that remians unsolved):
    if (any(costs == Inf)) stop("Function cannot currently handle infinite costs (which are treated as \"missing\" arcs).")
    
    # Set directed costs length:
    n_directed_costs <- length(x = costs)
    
    # First generate all possible connected graphs of n vertices:
    state_graphs <- permute_connected_graphs(n_vertices = n_vertices)
    
    # Permute asymmetric costmatrices from state graphs and costs:
    asymmetric_costmatrices <- lapply(
      X = state_graphs,
      FUN = function(state_graph) {
        
        # Set size of permutation (number of arcs):
        permutation_size <- nrow(x = state_graph)
        
        # Permute initial set of costs:
        permutations <- as.matrix(x = expand.grid(rep(x = list(costs), times = permutation_size)))
        
        # Set as ratios by dividing through by minimum cost:
        permutation_ratios <- apply(X = permutations, MARGIN = 1, FUN = function(i) i / min(x = i))
        
        # Find any duplicated cost ratios (can be eliminated):
        duplicate_ratios <- which(x = duplicated(x = apply(X = permutation_ratios, MARGIN = 2, FUN = paste, collapse = "%")))
        
        # If there are duplicatec cost ratios then remove these:
        if (length(x = duplicate_ratios) > 0) permutations <- permutations[-duplicate_ratios, , drop = FALSE]
        
        # Create permuted state graphs using permuted costs:
        permuted_state_graphs <- lapply(
          X = as.list(x = 1:nrow(x = permutations)),
          FUN = function(i) {
            state_graph[, "weight"] <- permutations[i, ]
            state_graph
          }
        )
        
        # Convert into proper state graphs:
        permuted_state_graphs <- lapply(
          X = permuted_state_graphs,
          FUN = function(permuted_state_graph) {
            
            # Make labels table:
            labels <- matrix(
              data = c(unique(x = c(permuted_state_graph[, "from"], permuted_state_graph[, "to"])), states),
              nrow = n_vertices
            )
            
            # Re-label graph with states:
            permuted_state_graph[, "from"] <- labels[match(x = permuted_state_graph[, "from"], table = labels[, 1]), 2]
            permuted_state_graph[, "to"] <- labels[match(x = permuted_state_graph[, "to"], table = labels[, 1]), 2]

            # Make vertices:
            vertices <- data.frame(
              label = states,
              in_degree = 0,
              out_degree = 0,
              eccentricity = 0,
              periphery = 0,
              centre = 0
            )
            for(i in states) {
              vertices[vertices[, "label"] == i, "in_degree"] <- sum(x = permuted_state_graph[, "to"] == i)
              vertices[vertices[, "label"] == i, "out_degree"] <- sum(x = permuted_state_graph[, "from"] == i)
            }
            
            # Make adjacency matrix for graph:
            adjacency_matrix <- matrix(
              data = 0,
              nrow = n_vertices,
              ncol = n_vertices,
              dimnames = list(states, states)
            )
            for(i in 1:nrow(x = permuted_state_graph)) adjacency_matrix[permuted_state_graph[i, "from"], permuted_state_graph[i, "to"]] <- 1
            
            # Compile stategraph:
            permuted_state_graph <- list(
              n_vertices = n_vertices,
              n_arcs = nrow(x = permuted_state_graph),
              n_states = n_vertices,
              single_states = states,
              type = "custom",
              arcs = permuted_state_graph,
              vertices = vertices,
              radius = min(x = vertices[, "eccentricity"]),
              diameter = max(x = vertices[, "eccentricity"]),
              adjacency_matrix = adjacency_matrix,
              directed = TRUE,
              includes_polymorphisms = FALSE,
              polymorphism_costs = "additive",
              polymorphism_geometry = "simplex",
              polymorphism_distance = "euclidean",
              includes_uncertainties = FALSE,
              pruned = FALSE,
              dollo_penalty = 999,
              base_age = 1,
              weight = 1
            )
            
            # Make class stategraph:
            class(x = permuted_state_graph) <- "stateGraph"

            # Return now formatted stategraph:
            permuted_state_graph
          }
        )
        
        # Convert permuted state graphs to costmatrices:
        asymmetric_costmatrices <- lapply(X = permuted_state_graphs, FUN = convert_stategraph_to_costmatrix)
        
        # Look for an symmetric costmatrices (need pruning):
        is_costmatrix_symmetric <- unlist(
          x = lapply(
            X = asymmetric_costmatrices,
            FUN = function(i) isSymmetric(object = i$costmatrix)
          )
        )

        # Prune any symmetric costmatrices:
        if (any(is_costmatrix_symmetric)) asymmetric_costmatrices <- asymmetric_costmatrices[-which(x = is_costmatrix_symmetric)]
        
        # Return asymmetric costmatrices:
        asymmetric_costmatrices
      }
    )
    
    # restructure as single list of costmatrices:
    asymmetric_costmatrices <- do.call( what = c, args = asymmetric_costmatrices)
    
    # Fix costmatrices:
    asymmetric_costmatrices <- lapply(X = asymmetric_costmatrices, FUN = fix_costmatrix)
    
    # Find any duplicate costmatrices:
    duplicated_costmatrices <- duplicated(x = unlist(x = lapply(X = asymmetric_costmatrices, FUN = function(i) paste(i$costmatrix, collapse = "%"))))
    
    # Remove any duplicated costmatrices:
    if (any(duplicated_costmatrices)) asymmetric_costmatrices <- asymmetric_costmatrices[-which(x = duplicated_costmatrices)]
    
    # If only want asymmetric than can just return those now:
    if (symmetry == "asymmetric") return(value = asymmetric_costmatrices)
  }
  
  # If still here must want both so return these:
  return(value = c(symmetric_costmatrices, asymmetric_costmatrices))
}
