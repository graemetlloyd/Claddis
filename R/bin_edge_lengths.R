#' Edge-lengths present in time-bins
#'
#' @description
#'
#' Given a time-scaled tree and set of time bin boundaries will sum the edge-lengths present in each bin.
#'
#' @param time_tree A time-scaled tree in phylo format with a \code{$root.time} value.
#' @param time_bins A vector of ages in millions of years of time bin boundaries in old-to-young order.
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
#' \item{edge.length.in.bin}{A vector giving the summed values in millions of years for each time bin. Names indicate the maximum and minimum values for each time bin.}
#' \item{terminal.edge.length.in.bin}{As above, but counting terminal edges only.}
#' \item{internal.edge.length.in.bin}{As above, but counting internal edges only.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a random 10-taxon tree:
#' time_tree <- ape::rtree(10)
#'
#' # Add root age:
#' time_tree$root.time <- 100
#'
#' # Create time bins:
#' time_bins <- seq(100, 0, length.out = 11)
#'
#' # Get edge lengths for each bin:
#' bin_edge_lengths(time_tree, time_bins)
#' @export bin_edge_lengths
bin_edge_lengths <- function(time_tree, time_bins, pruned_tree = NULL) {

  # Tree must have $root.time:
  if (is.null(time_tree$root.time)) stop("time_tree must have $root.time or function can not work.")

  # If tree is pruned from a larger version (where the terminal-internal dichotomy really applies):
  if (!is.null(pruned_tree)) {

    # Check pruned tree is subset of tree:
    if (!all(ape::drop.tip(time_tree, setdiff(time_tree$tip.label, pruned_tree$tip.label))$edge == pruned_tree$edge)) stop("ERROR: pruned_tree must be subtree of time_tree.")

    # Get dropped taxa:
    dropped.tips <- setdiff(time_tree$tip.label, pruned_tree$tip.label)

    # Collapse terminal branch lengths of dropped tips to zero:
    time_tree$edge.length[match(match(dropped.tips, time_tree$tip.label), time_tree$edge[, 2])] <- 0

    # Only continue if there are more branch lengths to collapse:
    if (sum(pruned_tree$edge.length) < sum(time_tree$edge.length)) {

      # Find descendant tips of each node:
      descendant.tips <- lapply(as.list((ape::Ntip(time_tree) + 1):(ape::Ntip(time_tree) + ape::Nnode(time_tree))), strap::FindDescendants, tree = time_tree)

      # Add node numbers as names:
      names(descendant.tips) <- (ape::Ntip(time_tree) + 1):(ape::Ntip(time_tree) + ape::Nnode(time_tree))

      # Find edges to collapse (those with dropped descendants only):
      edges.to.collapse <- match(as.numeric(names(which(unlist(lapply(lapply(lapply(descendant.tips, match, table = match(dropped.tips, time_tree$tip.label)), is.na), sum)) == 0))), time_tree$edge[, 2])

      # Collapse these branches to zero:
      time_tree$edge.length[edges.to.collapse] <- 0
    }
  }

  # Get terminal edge numbers:
  terminal.edges <- match(1:ape::Ntip(time_tree), time_tree$edge[, 2])

  # Get internal edge numbers:
  internal.edges <- setdiff(1:nrow(time_tree$edge), terminal.edges)

  # Enforce old-to-young order of time bins:
  time_bins <- sort(x = time_bins, decreasing = TRUE)

  # Create all-zero vector to store ouput in:
  internal.edge.length.in.bin <- terminal.edge.length.in.bin <- edge.length.in.bin <- rep(0, length(time_bins) - 1)

  # Date nodes in tree:
  node.ages <- date_nodes(time_tree = time_tree)

  # Get maximum age for each edge:
  tree.edge.maxs <- node.ages[time_tree$edge[, 1]]

  # Get minimum age for each edge:
  tree.edge.mins <- node.ages[time_tree$edge[, 2]]

  # For each time bin:
  for (i in 2:length(time_bins)) {

    # Find out which edges (if any) are present in the bin:
    edges.in.bin <- intersect(which(tree.edge.maxs > time_bins[i]), which(tree.edge.mins < time_bins[(i - 1)]))

    # If there is at least one edge in bin:
    if (length(edges.in.bin) > 0) {

      # Get maximum age for each edge in bin:
      in.bin.edge.maxs <- tree.edge.maxs[edges.in.bin]

      # Get minimum age for each edge in bin:
      in.bin.edge.mins <- tree.edge.mins[edges.in.bin]

      # Remove any part of edge that is found before the bin:
      if (sum(in.bin.edge.maxs > time_bins[(i - 1)]) > 0) in.bin.edge.maxs[in.bin.edge.maxs > time_bins[(i - 1)]] <- time_bins[(i - 1)]

      # Remove any part of edge that is found after the bin:
      if (sum(in.bin.edge.mins < time_bins[i]) > 0) in.bin.edge.mins[in.bin.edge.mins < time_bins[i]] <- time_bins[i]

      # Get sum of edge lengths found in the bin:
      edge.length.in.bin[(i - 1)] <- sum(in.bin.edge.maxs - in.bin.edge.mins)

      # Get sum of terminal edge lengths found in the bin:
      terminal.edge.length.in.bin[(i - 1)] <- sum(in.bin.edge.maxs[match(setdiff(edges.in.bin, internal.edges), edges.in.bin)] - in.bin.edge.mins[match(setdiff(edges.in.bin, internal.edges), edges.in.bin)])

      # Get sum of internal edge lengths found in the bin:
      internal.edge.length.in.bin[(i - 1)] <- sum(in.bin.edge.maxs[match(setdiff(edges.in.bin, terminal.edges), edges.in.bin)] - in.bin.edge.mins[match(setdiff(edges.in.bin, terminal.edges), edges.in.bin)])
    }
  }

  # Add time bin max-mins as names:
  names(terminal.edge.length.in.bin) <- names(internal.edge.length.in.bin) <- names(edge.length.in.bin) <- apply(cbind(time_bins[1:(length(time_bins) - 1)], time_bins[2:length(time_bins)]), 1, paste, collapse = "-")

  # Compile output:
  output <- list(edge.length.in.bin = edge.length.in.bin, terminal.edge.length.in.bin = terminal.edge.length.in.bin, internal.edge.length.in.bin = internal.edge.length.in.bin)

  # Return output:
  return(output)
}
