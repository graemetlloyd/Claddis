#' Edge matching function
#'
#' @description
#'
#' Given two trees where one is a pruned version of the other gives matching edges and nodes of pruned tree to original tree.
#'
#' @param original_tree A tree in phylo format.
#' @param pruned_tree A tree in phylo format that represents a pruned version of \code{original_tree}.
#'
#' @details
#'
#' Finds matching edge(s) and node(s) for a pruned tree in the original tree from which it was created. This is intended as an internal function, but may be of use to someone.
#'
#' @return
#'
#' \item{matching_edges}{A list of the matching edges.}
#' \item{matching_nodes}{A matrix of matching node numbers.}
#' \item{removed_edges}{A vector of the removed edges.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a random 10-taxon tree:
#' original_tree <- ape::rtree(n = 10)
#'
#' # Remove three leaves:
#' pruned_tree <- ape::drop.tip(phy = original_tree, tip = c("t1", "t3", "t8"))
#'
#' # Find matching edges:
#' X <- match_tree_edges(original_tree, pruned_tree)
#'
#' # Show matching edges:
#' X$matching_edges
#'
#' # Show matching nodes:
#' X$matching_nodes
#'
#' # Show removed edges:
#' X$removed_edges
#' @export match_tree_edges
match_tree_edges <- function(original_tree, pruned_tree) {

  # Get number of pruned tips:
  n_pruned_tips <- ape::Ntip(phy = pruned_tree)

  # Get number of pruned nodes:
  n_pruned_nodes <- ape::Nnode(phy = pruned_tree)

  # Conditional if pruned tree too small:
  if (n_pruned_tips < 3) stop("pruned_tree includes too few (<3) taxa to be used.")

  # Conditional in case where pruned tree taxa are not a subset of the original tree taxa:
  if (length(x = setdiff(x = pruned_tree$tip.label, y = original_tree$tip.label)) > 0) stop("pruned_tree cannot include taxa not present in original_tree.")

  # First find removed taxa (if any):
  removed_taxa <- setdiff(x = original_tree$tip.label, y = pruned_tree$tip.label)

  # If no taxa are removed:
  if (length(x = removed_taxa) == 0) {

    # Record removed edges as an empty vector:
    removed_edges <- numeric(0)

    # Return matching edges as list:
    matching_edges <- as.list(x = 1:nrow(original_tree$edge))

    # Create lists of nodes (which will be identical):
    clades <- corresponding_nodes <- c((n_pruned_tips + 1):(n_pruned_tips + n_pruned_nodes), 1:n_pruned_tips)

    # If taxa are removed:
  } else {

    # Create vector to store clades found in pruned tree:
    clades <- vector(mode = "character")

    # For each internal node of pruned tree:
    for (i in (n_pruned_tips + 1):(n_pruned_tips + n_pruned_nodes)) {

      # Get descendants of node (members of clade):
      clades <- c(clades, paste(pruned_tree$tip.label[strap::FindDescendants(n = i, tree = pruned_tree)], collapse = "%%SpLiTtEr%%"))

      # Update with node number:
      names(clades)[length(x = clades)] <- i
    }

    # Create vector to store corresponding node numbers in original tree:
    corresponding_nodes <- vector(mode = "numeric")

    # For each clade in pruned tree find and store node number in original tree:
    for (i in clades) corresponding_nodes <- c(corresponding_nodes, find_mrca(descendant_names = strsplit(i, "%%SpLiTtEr%%")[[1]], tree = original_tree))

    # Add tips to node numbers for original_tree:
    corresponding_nodes <- c(corresponding_nodes, match(pruned_tree$tip.label, original_tree$tip.label))

    # Add tips to node numbers for pruned_tree:
    clades <- c(as.numeric(names(clades)), 1:n_pruned_tips)

    # Make edge matrix for pruned tree using corresponding nodes in original tree:
    pruned_edges <- cbind(corresponding_nodes[match(pruned_tree$edge[, 1], clades)], corresponding_nodes[match(pruned_tree$edge[, 2], clades)])

    # Find edges that match EXACTLY with the original tree:
    matching_edges <- match(apply(pruned_edges, 1, paste, collapse = "%%"), apply(original_tree$edge, 1, paste, collapse = "%%"))

    # List non-matching edges for further searching:
    nonmatching_edges <- pruned_edges[is.na(matching_edges), ]

    # Only continue if there are non-matching edges (will be the case if only "outgroup(s)" are removed:
    if (length(x = nonmatching_edges) > 0) {

      # Correct stupid matrix to vector problem:
      if (!is.matrix(nonmatching_edges)) nonmatching_edges <- matrix(nonmatching_edges, ncol = 2)

      # For each non-matching edge:
      for (i in 1:nrow(nonmatching_edges)) {

        # Get start (ancestral) node:
        start_node <- nonmatching_edges[i, 1]

        # Get end (terminal) node:
        end_node <- nonmatching_edges[i, 2]

        # Create edges vector to store multiple edges that correspond to edge on pruned tree:
        edges <- match(end_node, original_tree$edge[, 2])

        # Keep going until start and end node are joined by contiguous edges:
        while (length(x = sort(x = match(original_tree$edge[edges, 1], start_node))) == 0) {

          # Update end node:
          end_node <- original_tree$edge[match(end_node, original_tree$edge[, 2]), 1]

          # Update edges:
          edges <- c(edges, match(end_node, original_tree$edge[, 2]))
        }

        # Update matching edges with multiple edges separated by a double-percent:
        matching_edges[which(x = is.na(matching_edges))[1]] <- paste(rev(edges), collapse = "%%")
      }

      # Get matching edges as list:
      matching_edges <- lapply(X = strsplit(matching_edges, "%%"), as.numeric)

      # If there are no non-matching edges:
    } else {

      # Get matching edges as list:
      matching_edges <- as.list(x = matching_edges)
    }

    # Get removed edges (branches from original tree missing in pruned tree:
    removed_edges <- setdiff(x = 1:nrow(original_tree$edge), y = as.numeric(unlist(x = matching_edges)))
  }

  # Add names to matching edges:
  names(matching_edges) <- 1:nrow(pruned_tree$edge)

  # Make list of matching nodes:
  node_matches <- cbind(pruned_node = clades, original_node = corresponding_nodes)

  # Return compiled output:
  list(matching_edges = matching_edges, matching_nodes = node_matches, removed_edges = removed_edges)
}
