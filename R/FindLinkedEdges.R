FindLinkedEdges <- function(tree) {
	
	# Create matrix of linked edges (all starting as unlinked, 0):
	edge.link.matrix <- matrix(0, nrow=nrow(tree$edge), ncol=nrow(tree$edge), dimnames=list(1:nrow(tree$edge), 1:nrow(tree$edge)))
	
	# For each edge:
	for(i in 1:nrow(tree$edge)) {
		
		# Find linked edges:
		links <- setdiff(union(which(apply(tree$edge == tree$edge[i, 1], 1, sum) == 1), which(apply(tree$edge == tree$edge[i, 2], 1, sum) == 1)), i)
		
		# Code 1 (linked) for linked edges in matrix:
		edge.link.matrix[i, links] <- edge.link.matrix[links, i] <- 1
		
	}
	
	# Return matrix:
	return(edge.link.matrix)
	
}
