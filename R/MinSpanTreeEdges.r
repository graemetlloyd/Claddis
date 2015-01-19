MinSpanTreeEdges <- function(dist.matrix) {
	
	# Convert to matrix and set up links matrix:
	dist.matrix <- as.matrix(dist.matrix)

	# Get links matrix for minimum spanning tree:
	links.matrix <- mst(dist.matrix)
	
	# Update matrix class to NULL:
	class(links.matrix) <- NULL
	
	# Create empty matrix to store edges for minimum spanning tree:
	min.span.tree.edges <- matrix(nrow=0, ncol=2, dimnames=list(c(), c("From", "To")))
	
	# For each row:
	for(i in 1:(nrow(links.matrix) - 1)) {
		
		# For each column:
		for(j in (i + 1):ncol(links.matrix)) {
			
			# If there is a link then record it:
			if(links.matrix[i, j] == 1) min.span.tree.edges <- rbind(min.span.tree.edges, c(rownames(links.matrix)[i], colnames(links.matrix)[j]))
			
		}
		
	}
	
	# Get distances:
	distances <- diag(dist.matrix[min.span.tree.edges[, "From"], min.span.tree.edges[, "To"]])
	
	# Add names to distances:
	names(distances) <- apply(min.span.tree.edges, 1, paste, collapse="->")
	
	# Output distances for minimum soanning tree:
	return(distances)
		
}
