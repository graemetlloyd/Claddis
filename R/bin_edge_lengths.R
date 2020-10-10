#' Edge-lengths present in time-bins
#'
#' @description
#'
#' Given a time-scaled tree and set of time bin boundaries will sum the edge-lengths present in each bin.
#'
#' @param time_tree A time-scaled tree in phylo format with a \code{$root.time} value.
#' @param time_bins An object of class \code{timeBins}.
#' @param pruned_tree A time-scaled tree in phylo format with a \code{$root.time} value that is a subset of \code{time_tree}.
#'
#' @details
#'
#' Calculates the total edge duration of a time-scaled tree present in a series of time bins. This is intended as an internal function for rate calculations, but may be of use to someone.
#'
#' The option of using a \code{pruned_tree} allows the user to correctly classify internal and terminal branches in a subtree of the larger tree. So for example, if taxa A and B are sisters then after pruning B the subtree branch leading to A is composed of an internal and a terminal branch on the complete tree.
#'
#' @return
#'
#' \item{binned_edge_lengths}{A vector giving the summed values in millions of years for each time bin. Names indicate the maximum and minimum values for each time bin.}
#' \item{binned_terminal_edge_lengths}{As above, but counting terminal edges only.}
#' \item{binned_internal_edge_lengths}{As above, but counting internal edges only.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a random 10-taxon tree:
#' time_tree <- ape::rtree(n = 10)
#'
#' # Add root age:
#' time_tree$root.time <- max(diag(ape::vcv(time_tree)))
#'
#' # Create time bins:
#' time_bins <- matrix(data = c(seq(from = time_tree$root.time, to = 0,
#'   length.out = 11)[1:10], seq(from = time_tree$root.time, to = 0,
#'   length.out = 11)[2:11]), ncol = 2, dimnames = list(LETTERS[1:10],
#'   c("fad", "lad")))
#'
#' # Set class:
#' class(time_bins) <- "timeBins"
#'
#' # Get edge lengths for each bin:
#' bin_edge_lengths(time_tree = time_tree, time_bins = time_bins)
#' @export bin_edge_lengths
bin_edge_lengths <- function(time_tree, time_bins, pruned_tree = NULL) {

  # Check time_bins is in a valid format and stop and warn user if not:
  if (!is.timeBins(x = time_bins)) stop(check_timeBins(time_bins = time_bins))

  # Find number of tips in tree:
  n_tips <- ape::Ntip(phy = time_tree)

  # Fine number of nodes in tree:
  n_nodes <- ape::Nnode(phy = time_tree)

  # Tree must have $root.time:
  if (is.null(time_tree$root.time)) stop("time_tree must have $root.time or function can not work.")

  # If tree is pruned from a larger version (where the terminal-internal dichotomy really applies):
  if (!is.null(pruned_tree)) {

    # Check pruned tree is subset of tree:
    if (!all(ape::drop.tip(phy = time_tree, tip = setdiff(x = time_tree$tip.label, y = pruned_tree$tip.label))$edge == pruned_tree$edge)) stop("ERROR: pruned_tree must be subtree of time_tree.")

    # Get dropped taxa:
    dropped_tips <- setdiff(x = time_tree$tip.label, y = pruned_tree$tip.label)

    # Collapse terminal branch lengths of dropped tips to zero:
    time_tree$edge.length[match(match(dropped_tips, time_tree$tip.label), time_tree$edge[, 2])] <- 0

    # Only continue if there are more branch lengths to collapse:
    if (sum(pruned_tree$edge.length) < sum(time_tree$edge.length)) {

      # Find descendant tips of each node:
      descendants <- lapply(X = as.list(x = (n_tips + 1):(n_tips + n_nodes)), strap::FindDescendants, tree = time_tree)

      # Add node numbers as names:
      names(descendants) <- (n_tips + 1):(n_tips + n_nodes)

      # Find edges to collapse (those with dropped descendants only):
      edges_to_collapse <- match(as.numeric(names(which(x = unlist(x = lapply(X = lapply(X = lapply(X = descendants, match, table = match(dropped_tips, time_tree$tip.label)), is.na), sum)) == 0))), time_tree$edge[, 2])

      # Collapse these branches to zero:
      time_tree$edge.length[edges_to_collapse] <- 0
    }
  }

  # Get terminal edge numbers:
  terminal_edges <- match(1:n_tips, time_tree$edge[, 2])

  # Get internal edge numbers:
  internal_edges <- setdiff(x = 1:nrow(time_tree$edge), y = terminal_edges)

  # Create all-zero vector to store ouput in:
  binned_internal_edge_lengths <- binned_terminal_edge_lengths <- binned_edge_lengths <- rep(0, nrow(x = time_bins))

  # Date nodes in tree:
  node_ages <- date_nodes(time_tree = time_tree)

  # Get maximum age for each edge:
  edge_maxima <- node_ages[time_tree$edge[, 1]]

  # Get minimum age for each edge:
  edge_minima <- node_ages[time_tree$edge[, 2]]

  # For each time bin:
  for (i in 1:nrow(x = time_bins)) {

    # Find out which edges (if any) are present in the bin:
    binned_edges <- intersect(which(x = edge_maxima > time_bins[i, "lad"]), which(x = edge_minima < time_bins[i, "fad"]))

    # If there is at least one edge in bin:
    if (length(x = binned_edges) > 0) {

      # Get maximum age for each edge in bin:
      binned_edge_maxima <- edge_maxima[binned_edges]

      # Get minimum age for each edge in bin:
      binned_edge_minima <- edge_minima[binned_edges]

      # Remove any part of edge that is found before the bin:
      if (sum(binned_edge_maxima > time_bins[i, "fad"]) > 0) binned_edge_maxima[binned_edge_maxima > time_bins[i, "fad"]] <- time_bins[i, "fad"]

      # Remove any part of edge that is found after the bin:
      if (sum(binned_edge_minima < time_bins[i, "lad"]) > 0) binned_edge_minima[binned_edge_minima < time_bins[i, "lad"]] <- time_bins[i, "lad"]

      # Get sum of edge lengths found in the bin:
      binned_edge_lengths[i] <- sum(binned_edge_maxima - binned_edge_minima)

      # Get sum of terminal edge lengths found in the bin:
      binned_terminal_edge_lengths[i] <- sum(binned_edge_maxima[match(setdiff(x = binned_edges, y = internal_edges), binned_edges)] - binned_edge_minima[match(setdiff(x = binned_edges, y = internal_edges), binned_edges)])

      # Get sum of internal edge lengths found in the bin:
      binned_internal_edge_lengths[i] <- sum(binned_edge_maxima[match(setdiff(x = binned_edges, y = terminal_edges), binned_edges)] - binned_edge_minima[match(setdiff(x = binned_edges, y = terminal_edges), binned_edges)])
    }
  }

  # Add time bin names to output:
  names(binned_terminal_edge_lengths) <- names(binned_internal_edge_lengths) <- names(binned_edge_lengths) <- rownames(time_bins)

  # Return compiled output:
  list(binned_edge_lengths = binned_edge_lengths, binned_terminal_edge_lengths = binned_terminal_edge_lengths, binned_internal_edge_lengths = binned_internal_edge_lengths)
}
