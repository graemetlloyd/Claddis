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
#' \item{matching.nodes}{A matrix of matching node numbers.}
#' \item{removed_edges}{A vector of the removed edges.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a random 10-taxon tree:
#' original_tree <- ape::rtree(10)
#'
#' # Remove three leaves:
#' pruned_tree <- ape::drop.tip(original_tree, c("t1", "t3", "t8"))
#'
#' # Find matching edges:
#' X <- match_tree_edges(original_tree, pruned_tree)
#'
#' # Show matching edges:
#' X$matching_edges
#'
#' # Show removed edges:
#' X$removed_edges
#' @export match_tree_edges
match_tree_edges <- function(original_tree, pruned_tree) {

  # Conditional if pruned tree too small:
  if (ape::Ntip(pruned_tree) < 3) stop("pruned_tree includes too few (<3) taxa to be used.")

  # Conditional in case where pruned tree taxa are not a subset of the original tree taxa:
  if (length(setdiff(pruned_tree$tip.label, original_tree$tip.label)) > 0) stop("pruned_tree cannot include taxa not present in original_tree.")

  # First find removed taxa (if any):
  removed_taxa <- setdiff(original_tree$tip.label, pruned_tree$tip.label)

  # If no taxa are removed:
  if (length(removed_taxa) == 0) {

    # Record removed edges as an empty vector:
    removed_edges <- numeric(0)

    # Return matching edges as list:
    matching_edges <- as.list(1:nrow(original_tree$edge))

    # Create lists of nodes (which will be identical):
    clades <- corresponding_nodes <- c((ape::Ntip(pruned_tree) + 1):(ape::Ntip(pruned_tree) + ape::Nnode(pruned_tree)), 1:ape::Ntip(pruned_tree))

    # If taxa are removed:
  } else {

    # Create vector to store clades found in pruned tree:
    clades <- vector(mode = "character")

    # For each internal node of pruned tree:
    for (i in (ape::Ntip(pruned_tree) + 1):(ape::Ntip(pruned_tree) + ape::Nnode(pruned_tree))) {

      # Get descendants of node (members of clade):
      clades <- c(clades, paste(pruned_tree$tip.label[strap::FindDescendants(i, pruned_tree)], collapse = "%%SpLiTtEr%%"))

      # Update with node number:
      names(clades)[length(clades)] <- i
    }

    # Create vector to store corresponding node numbers in original tree:
    corresponding_nodes <- vector(mode = "numeric")

    # For each clade in pruned tree find and store node number in original tree:
    for (i in clades) corresponding_nodes <- c(corresponding_nodes, find_mrca(descendant_names = strsplit(i, "%%SpLiTtEr%%")[[1]], tree = original_tree))

    # Add tips to node numbers for original_tree:
    corresponding_nodes <- c(corresponding_nodes, match(pruned_tree$tip.label, original_tree$tip.label))

    # Add tips to node numbers for pruned_tree:
    clades <- c(as.numeric(names(clades)), 1:ape::Ntip(pruned_tree))

    # Make edge matrix for pruned tree using corresponding nodes in original tree:
    pruned.edges <- cbind(corresponding_nodes[match(pruned_tree$edge[, 1], clades)], corresponding_nodes[match(pruned_tree$edge[, 2], clades)])

    # Find edges that match EXACTLY with the original tree:
    matching_edges <- match(apply(pruned.edges, 1, paste, collapse = "%%"), apply(original_tree$edge, 1, paste, collapse = "%%"))

    # List non-matching edges for further searching:
    nonmatching_edges <- pruned.edges[is.na(matching_edges), ]

    # Only continue if there are non-matching edges (will be the case if only "outgroup(s)" are removed:
    if (length(nonmatching_edges) > 0) {

      # Correct stupid matrix to vector problem:
      if (!is.matrix(nonmatching_edges)) nonmatching_edges <- matrix(nonmatching_edges, ncol = 2)

      # For each non-matching edge:
      for (i in 1:nrow(nonmatching_edges)) {

        # Get start (ancestral) node:
        start.node <- nonmatching_edges[i, 1]

        # Get end (terminal) node:
        end.node <- nonmatching_edges[i, 2]

        # Create edges vector to store multiple edges that correspond to edge on pruned tree:
        edges <- match(end.node, original_tree$edge[, 2])

        # Keep going until start and end node are joined by contiguous edges:
        while (length(sort(x = match(original_tree$edge[edges, 1], start.node))) == 0) {

          # Update end node:
          end.node <- original_tree$edge[match(end.node, original_tree$edge[, 2]), 1]

          # Update edges:
          edges <- c(edges, match(end.node, original_tree$edge[, 2]))
        }

        # Update matching edges with multiple edges separated by a double-percent:
        matching_edges[which(is.na(matching_edges))[1]] <- paste(rev(edges), collapse = "%%")
      }

      # Get matching edges as list:
      matching_edges <- lapply(strsplit(matching_edges, "%%"), as.numeric)

      # If there are no non-matching edges:
    } else {

      # Get matching edges as list:
      matching_edges <- as.list(matching_edges)
    }

    # Get removed edges (branches from original tree missing in pruned tree:
    removed_edges <- setdiff(1:nrow(original_tree$edge), as.numeric(unlist(matching_edges)))
  }

  # Add names to matching edges:
  names(matching_edges) <- 1:nrow(pruned_tree$edge)

  # Make list of matching nodes:
  node.matches <- cbind(Pruned_node = clades, Original_node = corresponding_nodes)

  # Compile output:
  output <- list(matching_edges = matching_edges, matching.nodes = node.matches, removed_edges = removed_edges)

  # Return output:
  return(output)
}
