EdgeLengthsInBins <- function(tree, time.bins) {
	
	# Tree must have $root.time:
	if(is.null(tree$root.time)) stop("ERROR: Tree must have $root.time or function can not work.")
	
	# Enforce old-to-young order of time bins:
	time.bins <- sort(time.bins, decreasing=T)
	
	# Create all-zero vector to store ouput in:
	edge.length.in.bin <- rep(0, length(time.bins) - 1)
	
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
			
		}
		
	}
	
	# Add time bin max-mins as names:
	names(edge.length.in.bin) <- apply(cbind(time.bins[1:(length(time.bins) - 1)], time.bins[2:length(time.bins)]), 1, paste, collapse="-")
	
	# Return edge lengths in bins:
	return(edge.length.in.bin)
	
}
