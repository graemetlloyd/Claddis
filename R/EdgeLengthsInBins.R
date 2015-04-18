EdgeLengthsInBins <- function(tree, time.bins, pruned.tree=NULL) {
	
	# Tree must have $root.time:
	if(is.null(tree$root.time)) stop("ERROR: Tree must have $root.time or function can not work.")
	
	# If tree is pruned from a larger version (where the terminal-internal dichotomy really applies):
	if(!is.null(pruned.tree)) {
		
		# Check pruned tree is subset of tree:
		if(!all(drop.tip(tree, setdiff(tree$tip.label, pruned.tree$tip.label))$edge == pruned.tree$edge)) stop("ERROR: Pruned tree must be subtree of tree.")
		
		# Get dropped taxa:
		dropped.tips <- setdiff(tree$tip.label, pruned.tree$tip.label)
		
		# Collapse terminal branch lengths of dropped tips to zero:
		tree$edge.length[match(match(dropped.tips, tree$tip.label), tree$edge[, 2])] <- 0
		
		# Only continue if there are more branch lengths to collapse:
		if(sum(pruned.tree$edge.length) < sum(tree$edge.length)) {
			
			# Find descendant tips of each node:
			descendant.tips <- lapply(as.list((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))), FindDescendants, tree=tree)
			
			# Add node numbers as names:
			names(descendant.tips) <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
			
			# Find edges to collapse (those with dropped descendants only):
			edges.to.collapse <- match(as.numeric(names(which(unlist(lapply(lapply(lapply(descendant.tips, match, table=match(dropped.tips, tree$tip.label)), is.na), sum)) == 0))), tree$edge[, 2])
			
			# Collapse these branches to zero:
			tree$edge.length[edges.to.collapse] <- 0
			
		}
		
	}
	
	# Get terminal edge numbers:
	terminal.edges <- match(1:Ntip(tree), tree$edge[, 2])
	
	# Get internal edge numbers:
	internal.edges <- setdiff(1:nrow(tree$edge), terminal.edges)
	
	# Enforce old-to-young order of time bins:
	time.bins <- sort(time.bins, decreasing=TRUE)
	
	# Create all-zero vector to store ouput in:
	internal.edge.length.in.bin <- terminal.edge.length.in.bin <- edge.length.in.bin <- rep(0, length(time.bins) - 1)
	
	# Date nodes in tree:
	node.ages <- GetNodeAges(tree)
	
	# Get maximum age for each edge:
	tree.edge.maxs <- node.ages[tree$edge[, 1]]
	
	# Get minimum age for each edge:
	tree.edge.mins <- node.ages[tree$edge[, 2]]
	
	# For each time bin:
	for(i in 2:length(time.bins)) {
		
		# Find out which edges (if any) are present in the bin:
		edges.in.bin <- intersect(which(tree.edge.maxs > time.bins[i]), which(tree.edge.mins < time.bins[(i - 1)]))
		
		# If there is at least one edge in bin:
		if(length(edges.in.bin) > 0) {
			
			# Get maximum age for each edge in bin:
			in.bin.edge.maxs <- tree.edge.maxs[edges.in.bin]
			
			# Get minimum age for each edge in bin:
			in.bin.edge.mins <- tree.edge.mins[edges.in.bin]
			
			# Remove any part of edge that is found before the bin:
			if(sum(in.bin.edge.maxs > time.bins[(i - 1)]) > 0) in.bin.edge.maxs[in.bin.edge.maxs > time.bins[(i - 1)]] <- time.bins[(i - 1)]
			
			# Remove any part of edge that is found after the bin:
			if(sum(in.bin.edge.mins < time.bins[i]) > 0) in.bin.edge.mins[in.bin.edge.mins < time.bins[i]] <- time.bins[i]
			
			# Get sum of edge lengths found in the bin:
			edge.length.in.bin[(i - 1)] <- sum(in.bin.edge.maxs - in.bin.edge.mins)
			
			# Get sum of terminal edge lengths found in the bin:
			terminal.edge.length.in.bin[(i - 1)] <- sum(in.bin.edge.maxs[match(setdiff(edges.in.bin, internal.edges), edges.in.bin)] - in.bin.edge.mins[match(setdiff(edges.in.bin, internal.edges), edges.in.bin)])
			
			# Get sum of internal edge lengths found in the bin:
			internal.edge.length.in.bin[(i - 1)] <- sum(in.bin.edge.maxs[match(setdiff(edges.in.bin, terminal.edges), edges.in.bin)] - in.bin.edge.mins[match(setdiff(edges.in.bin, terminal.edges), edges.in.bin)])

		}
		
	}
	
	# Add time bin max-mins as names:
	names(terminal.edge.length.in.bin) <- names(internal.edge.length.in.bin) <- names(edge.length.in.bin) <- apply(cbind(time.bins[1:(length(time.bins) - 1)], time.bins[2:length(time.bins)]), 1, paste, collapse="-")
	
	# Compile output:
	output <- list(edge.length.in.bin, terminal.edge.length.in.bin, internal.edge.length.in.bin)

	# Add names to output:
	names(output) <- c("edge.length.in.bin", "terminal.edge.length.in.bin", "internal.edge.length.in.bin")
	
	# Return output:
	return(output)
	
}
