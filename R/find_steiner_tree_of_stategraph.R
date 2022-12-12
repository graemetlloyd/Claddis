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
#' # XXXX:
#'
#' @export find_steiner_tree_of_stategraph
find_steiner_tree_of_stategraph <- function(stategraph, sampled_vertices) {

  ### EXAMPLE DATASET
  unordered_costmatrix <- make_costmatrix(
    min_state = 0,
    max_state= 5,
    character_type = "irreversible"
  )
  stategraph <- convert_costmatrix_to_stategraph(costmatrix = unordered_costmatrix)
  sampled_vertices <- c("1", "3", "4")
  
  ### CHECK SAMPLED VERTICES ARE UNIQUE AND PRESENT!
  ### CHECK STATEGRAPH IS A STATEGRAPH!
  ### ADD OPTION FOR SIMPLE SOLUTIONS FOR TWO VERTEX GRAPH (SHORTEST PATH) OR COMPLETE GRAPH/DIGRAPH (MST/MWSA) - WILL MAKE MUCH FASTER!!!!
  ### MAYBE ADD OTHER FAST OPTIONS FOR SIMPLE CASES (E.G. ORDERED BUT NO POLYMORPHISMS OR UNCERTAINTIES)
  ### WHAT TO RETURN? (COULD BE MULTILE SHORTEST STEINER TREES FOR EXAMPLE)
  ### ASYMMETRIC CHARACTERS CAN GET STUCK IF ROT IS NOT VALID (RUN OUT OF ARCS TO ADD BUT NOT BE ONGER THAN SHORTEST TREE
  
  # Case if only one sampled vertex then return tree as single vertex:
  if (length(x = sampled_vertices) == 1) return(sampled_vertices)
  
  # Case if at least two vertices:
  if (length(x = sampled_vertices) > 1) {
    
    # Set root options as every vertex that has at least one "from" arc:
    root_options <- as.character(x = unique(x = stategraph$arcs$from))
    
    # Generate all possible optimal Steiner trees:
    all_steiner_trees <- lapply(
      X = as.list(x = root_options),
      FUN = function(root_option) {
        
        # Initialise length threshold as Inf (used to avoid continuing searching for trees that are already too long):
        length_threshold <- Inf

        # Build starting trees:
        trees <- list(
          list(
            root = root_option,
            vertices_connected = root_option,
            used_arcs = stategraph$arcs[-(1:nrow(x = stategraph$arcs)), ],
            available_arcs = stategraph$arcs,
            tree_length = 0,
            sampled_vertices_connected = FALSE,
            active = TRUE
          )
        )
        
        # As long as there are active (unconnected, below threshold) trees:
        while(any(x = unlist(x = lapply(X = trees, FUN = function(tree) tree$active)))) {
          
          # Identify just active trees (i.e., those that still need to be connected and may be of minimum length):
          active_trees <- which(x = unlist(x = lapply(X = trees, FUN = function(tree) tree$active)))
          
          # Only keep working with active_trees
          for(active_tree in rev(x = active_trees)) {
            
            # Isolate active tree:
            tree <- trees[[active_tree]]
            
            # Find row numbers of new arcs (to add individually to create new trees):
            new_arcs <- unlist(
              x = lapply(
                X = as.list(x = tree$vertices_connected),
                FUN = function(vertex) which(x = tree$available_arcs$from == vertex)
              )
            )
            
            ### COULD BE NO NEW ARCS SO NEED CONDITIONALS TO DEAL WITH THIS (SET TO INACTIVE AND SET LENGTH TO INF?)
            
            # Duplicate tree to make new trees:
            new_trees <- rep(x = trees[active_tree], length(x = new_arcs))
            
            # Update each new tree by adding new arc and updating contents:
            new_trees <- mapply(
              FUN = function(new_tree, new_arc) {
                
                # Add new arc to used arcs:
                new_tree$used_arcs <- rbind(new_tree$used_arcs, new_tree$available_arcs[new_arc, ])
                
                # Update vertices connected:
                new_tree$vertices_connected <- unique(x = as.character(x = c(new_tree$used_arcs$from, new_tree$used_arcs$to)))
                
                # Update length of Steiner tree by adding weight of new arc:
                new_tree$tree_length <- new_tree$tree_length + new_tree$available_arcs[new_arc, "weight"]
                
                # Update available arcs by removing just added one:
                new_tree$available_arcs <- new_tree$available_arcs[-new_arc, , drop = FALSE]
                
                # Check if sampled vertices are now connected:
                if (length(x = setdiff(x = sampled_vertices, y = new_tree$vertices_connected)) == 0) {
                  new_tree$sampled_vertices_connected <- TRUE
                  new_tree$active <- FALSE
                }
                
                # If tree is too long to ever be the Steiner graph then make inactive:
                if (new_tree$tree_length >= length_threshold) new_tree$active <- FALSE
                
                # Return updated new tree:
                new_tree
              },
              new_tree = new_trees,
              new_arc = as.list(x = new_arcs),
              SIMPLIFY = FALSE
            )
            
            # Add new trees to trees:
            trees <- c(trees, new_trees)
          }
          
          # Remove the now superseded active trees:
          trees <- trees[-active_trees]

          # Update length_threshold to lowest connected value (or Inf if no trees are connected yet):
          length_threshold <- min(
            x = unlist(
              x = lapply(
                X = trees,
                FUN = function(tree) ifelse(
                  test = tree$sampled_vertices_connected,
                  yes = tree$tree_length,
                  no = Inf
                )
              )
            )
          )
          
          # Identify any dead trees (cannot be an optimal Steiner tree):
          dead_trees <- which(
            x = unlist(
              x = lapply(
                X = trees,
                FUN = function(tree) all(x = c(!tree$active, !tree$sampled_vertices_connected))
              )
            )
          )
          
          # Remove any dead trees found:
          if (length(x = dead_trees) > 0) trees <- trees[-dead_trees]
        }
        
        # Return trees:
        trees
      }
    )
    
    # Reconsitute into a single list (no nesting by root):
    all_steiner_trees <- do.call(
      what = c,
      args = all_steiner_trees
    )
    
    # Get just the lengths of the trees:
    steiner_tree_lengths <- unlist(x = lapply(X = all_steiner_trees, FUN = function(steiner_tree) steiner_tree$tree_length))
    
    # Collapse all_steiner_trees to just the shortest one(s):
    all_steiner_trees <- all_steiner_trees[steiner_tree_lengths == min(x = steiner_tree_lengths)]
    
    ### IF NOT DIRECTED CHECK FOR DUPLICATES (SAME BASIC ARCS JUST VARIED ORDER(S)) AND REMOVE (ROOT NO LONGER MATTERS?)
  }
}
