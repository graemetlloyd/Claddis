#' Edge-lengths present in time-bins
#'
#' @description
#'
#' Given a time-scaled tree and set of time bin boundaries will sum the edge-lengths present in each bin.
#' 
#' @param time.tree A time-scaled tree in phylo format with a \code{$root.time} value.
#' @param time.bins A vector of ages in millions of years of time bin boundaries in old-to-young order.
#' @param pruned.tree A time-scaled tree in phylo format with a \code{$root.time} value that is a subset of \code{time.tree}.
#'
#' @details
#'
#' Calculates the total edge duration of a time-scaled tree present in a series of time bins. This is intended as an internal function for rate calculations, but may be of use to someone.
#'
#' The option of using a \code{pruned.tree} allows the user to correctly classify internal and terminal branches in a subtree of the larger tree. So for example, if taxa A and B are sisters then after pruning B the subtree branch leading to A is composed of an internal and a terminal branch on the complete tree.
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
#' time.tree <- rtree(10)
#' 
#' # Add root age:
#' time.tree$root.time <- 100
#' 
#' # Create time bins:
#' time.bins <- seq(100, 0, length.out = 11)
#' 
#' # Get edge lengths for each bin:
#' bin_edge_lengths(time.tree, time.bins)
#' 
#' @export bin_edge_lengths
bin_edge_lengths <- function(time.tree, time.bins, pruned.tree = NULL) {
	
	# Tree must have $root.time:
	if (is.null(time.tree$root.time)) stop("ERROR: time.tree must have $root.time or function can not work.")
	
	# If tree is pruned from a larger version (where the terminal-internal dichotomy really applies):
	if (!is.null(pruned.tree)) {
		
		# Check pruned tree is subset of tree:
		if (!all(drop.tip(time.tree, setdiff(time.tree$tip.label, pruned.tree$tip.label))$edge == pruned.tree$edge)) stop("ERROR: pruned.tree must be subtree of time.tree.")
		
		# Get dropped taxa:
		dropped.tips <- setdiff(time.tree$tip.label, pruned.tree$tip.label)
		
		# Collapse terminal branch lengths of dropped tips to zero:
		time.tree$edge.length[match(match(dropped.tips, time.tree$tip.label), time.tree$edge[, 2])] <- 0
		
		# Only continue if there are more branch lengths to collapse:
		if (sum(pruned.tree$edge.length) < sum(time.tree$edge.length)) {
			
			# Find descendant tips of each node:
			descendant.tips <- lapply(as.list((ape::Ntip(time.tree) + 1):(ape::Ntip(time.tree) + ape::Nnode(time.tree))), FindDescendants, tree = time.tree)
			
			# Add node numbers as names:
			names(descendant.tips) <- (ape::Ntip(time.tree) + 1):(ape::Ntip(time.tree) + ape::Nnode(time.tree))
			
			# Find edges to collapse (those with dropped descendants only):
			edges.to.collapse <- match(as.numeric(names(which(unlist(lapply(lapply(lapply(descendant.tips, match, table = match(dropped.tips, time.tree$tip.label)), is.na), sum)) == 0))), time.tree$edge[, 2])
			
			# Collapse these branches to zero:
			time.tree$edge.length[edges.to.collapse] <- 0
			
		}
		
	}
	
	# Get terminal edge numbers:
	terminal.edges <- match(1:ape::Ntip(time.tree), time.tree$edge[, 2])
	
	# Get internal edge numbers:
	internal.edges <- setdiff(1:nrow(time.tree$edge), terminal.edges)
	
	# Enforce old-to-young order of time bins:
	time.bins <- sort(time.bins, decreasing = TRUE)
	
	# Create all-zero vector to store ouput in:
	internal.edge.length.in.bin <- terminal.edge.length.in.bin <- edge.length.in.bin <- rep(0, length(time.bins) - 1)
	
	# Date nodes in tree:
	node.ages <- date_nodes(time.tree)
	
	# Get maximum age for each edge:
	tree.edge.maxs <- node.ages[time.tree$edge[, 1]]
	
	# Get minimum age for each edge:
	tree.edge.mins <- node.ages[time.tree$edge[, 2]]
	
	# For each time bin:
	for(i in 2:length(time.bins)) {
		
		# Find out which edges (if any) are present in the bin:
		edges.in.bin <- intersect(which(tree.edge.maxs > time.bins[i]), which(tree.edge.mins < time.bins[(i - 1)]))
		
		# If there is at least one edge in bin:
		if (length(edges.in.bin) > 0) {
			
			# Get maximum age for each edge in bin:
			in.bin.edge.maxs <- tree.edge.maxs[edges.in.bin]
			
			# Get minimum age for each edge in bin:
			in.bin.edge.mins <- tree.edge.mins[edges.in.bin]
			
			# Remove any part of edge that is found before the bin:
			if (sum(in.bin.edge.maxs > time.bins[(i - 1)]) > 0) in.bin.edge.maxs[in.bin.edge.maxs > time.bins[(i - 1)]] <- time.bins[(i - 1)]
			
			# Remove any part of edge that is found after the bin:
			if (sum(in.bin.edge.mins < time.bins[i]) > 0) in.bin.edge.mins[in.bin.edge.mins < time.bins[i]] <- time.bins[i]
			
			# Get sum of edge lengths found in the bin:
			edge.length.in.bin[(i - 1)] <- sum(in.bin.edge.maxs - in.bin.edge.mins)
			
			# Get sum of terminal edge lengths found in the bin:
			terminal.edge.length.in.bin[(i - 1)] <- sum(in.bin.edge.maxs[match(setdiff(edges.in.bin, internal.edges), edges.in.bin)] - in.bin.edge.mins[match(setdiff(edges.in.bin, internal.edges), edges.in.bin)])
			
			# Get sum of internal edge lengths found in the bin:
			internal.edge.length.in.bin[(i - 1)] <- sum(in.bin.edge.maxs[match(setdiff(edges.in.bin, terminal.edges), edges.in.bin)] - in.bin.edge.mins[match(setdiff(edges.in.bin, terminal.edges), edges.in.bin)])

		}
		
	}
	
	# Add time bin max-mins as names:
	names(terminal.edge.length.in.bin) <- names(internal.edge.length.in.bin) <- names(edge.length.in.bin) <- apply(cbind(time.bins[1:(length(time.bins) - 1)], time.bins[2:length(time.bins)]), 1, paste, collapse = "-")
	
	# Compile output:
	output <- list(edge.length.in.bin = edge.length.in.bin, terminal.edge.length.in.bin = terminal.edge.length.in.bin, internal.edge.length.in.bin = internal.edge.length.in.bin)

	# Return output:
	return(output)
	
}
