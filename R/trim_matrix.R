#' Trims a morphological distance matrix
#'
#' @description
#'
#' Trims a morphological distance matrix by removing objects that cause empty cells.
#'
#' @param distance_matrix A distance matrix in the format created by \link{calculate_morphological_distances}.
#' @param tree If the distance matrix includes ancestors this should be the tree (phylo object) used to estimate their states.
#'
#' @details
#'
#' Trims a morphological distance matrix by removing nodes (terminal or internal) that cause empty cells allowing it to be passed to an ordination function such as \link{cmdscale}.
#'
#' Some distances are not calculable from cladistic matrices if there are taxa that have no coded characters in common. This algorithm iteratively removes the taxa responsible for the most empty cells until the matrix is complete (no empty cells).
#'
#' If the matrix includes estimated ancestral states the user should also provide the tree used (as the \code{tree} argument). The function will then also remove the tips from the tree and where reconstructed ancestors also cause empty cells will prune the minimum number of descendants of that node. The function will then renumber the nodes in the distance matrix so they match the pruned tree.
#'
#' @return
#'
#' \item{distance_matrix}{A complete distance matrix with all cells filled. If there were no empty cells will return original.}
#' \item{tree}{A tree (if supplied) with the removed taxa (see below) pruned. If no taxa are dropped will return the same tree as inputted. If no tree is supplied this is set to NULL.}
#' \item{removed_taxa}{A character vector listing the taxa removed. If none are removed this will be set to NULL.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{calculate_morphological_distances}
#'
#' @examples
#'
#' # Get morphological distances for Michaux (1989) data set:
#' distances <- calculate_morphological_distances(cladistic_matrix = michaux_1989)
#'
#' # Attempt to trim max.distance_matrix:
#' trim_matrix(distance_matrix = distances$distance_matrix)
#' @export trim_matrix
trim_matrix <- function(distance_matrix, tree = NULL) {

  # If input is of class "dist" first convert to a regular matrix:
  if (inherits(distance_matrix, what = "dist")) distance_matrix <- as.matrix(distance_matrix)

  # Check the input is a distance matrix:
  if (!is.matrix(distance_matrix)) stop("Input must be a distance matrix (i.e., either an object of class \"dist\" or a square matrix).")

  # Case if there is no tree:
  if (is.null(tree)) {

    # Case if distance matrix is already complete:
    if (length(x = which(x = is.na(distance_matrix))) == 0) {

      # Warn user:
      print("There are no gaps in the distance matrix")

      # There are no taxa to be removed:
      removed_taxa <- NULL

      # Return compile data in single variable:
      return(list(distance_matrix = distance_matrix, tree = tree, removed_taxa = removed_taxa))

      # Case if distance matrix has gaps:
    } else {

      # Vector to store taxa that need removing:
      removes <- vector(mode = "character")

      # Whilst there are still gaps in the matrix:
      while (length(x = which(x = is.na(distance_matrix))) > 0) {

        # Vector to store number of NAs for each row of the matrix:
        na_lengths <- vector(mode = "numeric")

        # For each row of the matrix get the number of NAs:
        for (i in 1:length(x = distance_matrix[, 1])) na_lengths[i] <- length(x = which(x = is.na(distance_matrix[i, ])))

        # Take row with most NAs and make it the taxon to delete:
        taxon_to_delete <- rownames(x = distance_matrix)[which(x = na_lengths == max(na_lengths))[1]]

        # Find the taxon's row:
        row_to_delete <- which(x = rownames(x = distance_matrix) == taxon_to_delete)

        # Remove it from the distance matrix:
        distance_matrix <- distance_matrix[-row_to_delete, -row_to_delete]

        # Add taxon to removes list:
        removes <- c(removes, taxon_to_delete)
      }

      # Return compiled data:
      return(list(distance_matrix = distance_matrix, tree = tree, removed_taxa = removes))
    }

    # Case if there is a tree:
  } else {

    # Get tip count:
    n_tips <- ape::Ntip(phy = tree)

    # Get node count:
    n_nodes <- ape::Nnode(phy = tree)

    # Case if distance matrix is already complete:
    if (length(x = which(x = is.na(distance_matrix))) == 0) {

      # Warn user:
      print("There are no gaps in the distance matrix")

      # Create NULL vector for removed taxa:
      removed_taxa <- NULL

      # Return unmodified matrix and tree:
      return(list(distance_matrix = distance_matrix, tree = tree, removed_taxa = removed_taxa))

      # Case if distance matrix has gaps:
    } else {

      # Get list of node numbers as text:
      node_numbers <- as.character((n_tips + 1):(n_tips + n_nodes))

      # Rename distance_matrix rownames by descendant taxa that define them:
      for (i in match(node_numbers, rownames(x = distance_matrix))) colnames(x = distance_matrix)[i] <- rownames(x = distance_matrix)[i] <- paste(sort(x = tree$tip.label[strap::FindDescendants(n = rownames(x = distance_matrix)[i], tree = tree)]), collapse = "%%")

      # Vector to store taxa and nodes that need removing:
      removes <- vector(mode = "character")

      # Temporarily store complete distance matrix:
      temporary_distance_matrix <- distance_matrix

      # Whilst there are still gaps in the matrix:
      while (length(x = which(x = is.na(distance_matrix))) > 0) {

        # Vector to store number of NAs for each row of the matrix:
        na_lengths <- vector(mode = "numeric")

        # For each row of the matrix get the number of NAs:
        for (i in 1:length(x = distance_matrix[, 1])) na_lengths[i] <- length(x = which(x = is.na(distance_matrix[i, ])))

        # Take row with most NAs and make it the taxon to delete:
        taxon_to_delete <- rownames(x = distance_matrix)[which(x = na_lengths == max(na_lengths))[1]]

        # Find the taxons row:
        row_to_delete <- which(x = rownames(x = distance_matrix) == taxon_to_delete)

        # Remove it from the distance matrix:
        distance_matrix <- distance_matrix[-row_to_delete, -row_to_delete]

        # Add taxon to removes list:
        removes <- c(removes, taxon_to_delete)
      }

      # Restore complete distance matrix:
      distance_matrix <- temporary_distance_matrix

      # Find tips to remove:
      tips_to_remove <- removes[sort(x = match(tree$tip.label, removes))]

      # Find nodes to remove:
      nodes_to_remove <- removes[grep("%%", removes)]

      # Vector to store nodes to delete (i.e. from where tips to delete originate):
      tipname_nodes_to_remove <- vector(mode = "numeric")

      # For each tip to delete:
      for (i in tips_to_remove) {

        # Get originating node for tip:
        originating_node <- tree$edge[match(which(x = tree$tip.label == i), tree$edge[, 2]), 1]

        # Get name of originating node:
        originating_node_name <- paste(sort(x = tree$tip.label[strap::FindDescendants(n = originating_node, tree = tree)]), collapse = "%%")

        # Add to nodes to delete vector:
        tipname_nodes_to_remove <- unique(x = c(tipname_nodes_to_remove, originating_node_name))
      }

      # Check if nodes need to be deleted that are not related to tips to be deleted:
      if (length(x = setdiff(x = nodes_to_remove, y = tipname_nodes_to_remove)) > 0) {

        # Identify nodes left to remove:
        nodes_still_to_remove <- setdiff(x = nodes_to_remove, y = tipname_nodes_to_remove)

        # Go through each node:
        for (i in nodes_still_to_remove) {

          # Find node number:
          originating_node <- find_mrca(descendant_names = strsplit(i, "%%")[[1]], tree = tree)

          # Find descendants of that node:
          descendant_nodes <- tree$edge[which(x = tree$edge[, 1] == originating_node), 2]

          # Case if one of the descendants is a tip:
          if (length(x = which(x = descendant_nodes <= n_tips)) > 0) {

            # Get taxon to exclude:
            taxon_to_exclude <- tree$tip.label[min(descendant_nodes)]

            # Case if all descendants are internal nodes:
          } else {

            # If first descendant node has more taxa:
            if (length(x = strap::FindDescendants(n = descendant_nodes[1], tree = tree)) >= length(x = strap::FindDescendants(n = descendant_nodes[2], tree = tree))) {

              # Get taxon to exclude:
              taxon_to_exclude <- tree$tip.label[strap::FindDescendants(n = descendant_nodes[2], tree = tree)]

              # If second descendant node has more taxa:
            } else {

              # Get taxon to exclude:
              taxon_to_exclude <- tree$tip.label[strap::FindDescendants(n = descendant_nodes[1], tree = tree)]
            }
          }

          # Add to removes list:
          removes <- c(removes, taxon_to_exclude)

          # Add to tips to remove list:
          tips_to_remove <- c(tips_to_remove, taxon_to_exclude)
        }
      }

      # Prune matrix until complete:
      distance_matrix <- distance_matrix[-match(removes, rownames(x = distance_matrix)), -match(removes, rownames(x = distance_matrix))]

      # Loop to prune removed taxa from node names:
      for (i in 1:length(x = distance_matrix[, 1])) {

        # Split name into constituent taxa:
        taxa <- strsplit(rownames(x = distance_matrix)[i], "%%")[[1]]

        # Get only still present taxa and rename rows and columns:
        rownames(x = distance_matrix)[i] <- colnames(x = distance_matrix)[i] <- paste(sort(x = taxa[which(x = is.na(match(taxa, tips_to_remove)))]), collapse = "%%")
      }

      # Establish if there are redundant rows (nodes defined by a single tip or by no tips at all):
      redundant_rows <- c(which(x = duplicated(rownames(x = distance_matrix))), which(x = rownames(x = distance_matrix) == ""))

      # If there are any redundant rows remove them from distance matrix:
      if (length(x = redundant_rows) > 0) distance_matrix <- distance_matrix[-redundant_rows, -redundant_rows]

      # Store tree prior to any pruning:
      full_tree <- tree

      # Remove pruned taxa from tree:
      tree <- ape::drop.tip(phy = tree, tip = tips_to_remove)

      # Correct root time (if necessary):
      tree <- fix_root_time(original_tree = full_tree, pruned_tree = tree)

      # Find node names:
      node_names <- rownames(x = distance_matrix)[grep("%%", rownames(x = distance_matrix))]

      # Find node numbers:
      for (j in 1:length(x = node_names)) names(node_names)[j] <- find_mrca(descendant_names = strsplit(node_names[j], "%%")[[1]], tree = tree)

      # Replace node names with numbers:
      colnames(x = distance_matrix)[match(node_names, colnames(x = distance_matrix))] <- rownames(x = distance_matrix)[match(node_names, rownames(x = distance_matrix))] <- names(node_names)

      # Return data in single variable:
      return(list(distance_matrix = distance_matrix, tree = tree, removed_taxa = removes))
    }
  }
}
