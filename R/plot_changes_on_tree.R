#' Plots character changes on branches
#'
#' @description
#'
#' Plots character changes in boxes on branches.
#'
#' @param character_changes A matrix of character changes.
#' @param time_tree Tree on which character changes occur.
#' @param label_size The size of the text for the barnch labels. Default is 0.5.
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
#' # Get first MPT for the Michaux data set:
#' time_tree <- ape::read.tree(text = paste0("(Ancilla:31.6,(Turrancilla:102.7,",
#'   "(Ancillista:1,Amalda:63.5):1):1);"))
#'
#' # Set root time for tree:
#' time_tree$root.time <- 103.7
#'
#' # Generate two equal length time bins:
#' time_bins <- matrix(data = c(seq(time_tree$root.time, 0, length.out = 3)[1:2],
#'   seq(time_tree$root.time, 0, length.out = 3)[2:3]), ncol = 2, dimnames = list(LETTERS[1:2],
#'   c("fad", "lad")))
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Get discrete character rates (includes changes):
#' out <- test_rates(
#'   time_tree = time_tree,
#'   cladistic_matrix = michaux_1989,
#'   time_bins = time_bins,
#'   branch_partitions = list(list(1)),
#'   alpha = 0.01
#' )
#'
#' # Plot character changes on the tree:
#' plot_changes_on_tree(
#'   character_changes = out$inferred_character_changes,
#'   time_tree = time_tree
#' )
#' @export plot_changes_on_tree
plot_changes_on_tree <- function(character_changes, time_tree, label_size = 0.5) {

  # Update tree edge lengths to number of character changes:
  time_tree$edge.length <- rle(sort(x = c(character_changes[, "edge"], 1:nrow(time_tree$edge))))$lengths - 0.5

  # Create empty edge labels vector:
  edge_labels <- rep(NA, nrow(time_tree$edge))

  # For each edge:
  for (i in 1:nrow(time_tree$edge)) {

    # Get rows for where changes occur:
    change_rows <- which(x = character_changes[, "edge"] == i)

    # If there are changes on edge:
    if (length(x = change_rows) > 0) {

      # Compile all changes into edge label:
      edge_labels[i] <- paste(paste(character_changes[change_rows, "character"], ": ", character_changes[change_rows, "from"], " -> ", character_changes[change_rows, "to"], sep = ""), collapse = "\n")
    }
  }

  # ADD DOT DOT DOT.....

  # Plot tree:
  plot(time_tree, direction = "upwards")

  # Add edge labels for changes:
  edgelabels(text = edge_labels, bg = "white", cex = label_size)

  # NEED TO LADDERISE LEFT IF WRITING ON RIGHT OF BRANCHES...
}
