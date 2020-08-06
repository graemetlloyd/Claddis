#' Plots character changes on branches
#'
#' @description
#'
#' Plots character changes in boxes on branches.
#'
#' @param character_changes A matrix of character changes.
#' @param time_tree Tree on which character changes occur.
#'
#' @details
#'
#' Takes the \code{character_changes} output from \link{test_rates} and plots it on the tree used to generate it.
#'
#' @return A plot of character changes on a tree.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Set random seed:
#' set.seed(17)
#'
#' # Generate a random tree for the Michaux data set:
#' time_tree <- ape::rtree(n = nrow(michaux_1989$matrix_1$matrix))
#'
#' # Update taxon names to match those in the data matrix:
#' time_tree$tip.label <- rownames(x = michaux_1989$matrix_1$matrix)
#'
#' # Set root time by making youngest taxon extant:
#' time_tree$root.time <- max(diag(x = ape::vcv(phy = time_tree)))
#'
#' # Get discrete character rates (includes changes):
#' out <- test_rates(time_tree, michaux_1989,
#'   seq(time_tree$root.time, 0, length.out = 3),
#'   branch_partitions = list(list(1)), alpha = 0.01
#' )
#'
#' # Plot character changes on the tree:
#' plot_changes_on_tree(
#'   out$inferred_character_changes,
#'   time_tree
#' )
#' @export plot_changes_on_tree
plot_changes_on_tree <- function(character_changes, time_tree) {

  # Update tree edge lengths to number of character changes:
  time_tree$edge.length <- rle(sort(x = c(character_changes[, "Edge"], 1:nrow(time_tree$edge))))$lengths - 1

  # Create empty edge labels vector:
  edge_labels <- rep(NA, nrow(time_tree$edge))

  # For each edge:
  for (i in 1:nrow(time_tree$edge)) {

    # Get rows for where changes occur:
    change_rows <- which(x = character_changes[, "Edge"] == i)

    # If there are changes on edge:
    if (length(x = change_rows) > 0) {

      # Compile all changes into edge label:
      edge_labels[i] <- paste(paste(character_changes[change_rows, "Character"], ": ", character_changes[change_rows, "From"], " -> ", character_changes[change_rows, "To"], sep = ""), collapse = "\n")
    }
  }

  # ADD DOT DOT DOT.....

  # Plot tree:
  plot(time_tree, direction = "upwards")

  # Add edge labels for changes:
  edgelabels(text = edge_labels, bg = "white", cex = 0.3)

  # NEED TO LADDERISE LEFT IF WRITING ON RIGHT OF BRANCHES...
}
