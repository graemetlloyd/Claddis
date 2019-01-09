#' Edge matching function
#'
#' @description
#'
#' Given two trees where one is a pruned version of the other gives matching edges and nodes of pruned tree to original tree.
#'
#' @param original.tree A tree in phylo format.
#' @param pruned.tree A tree in phylo format that represents a pruned version of \code{original.tree}.
#'
#' @details
#'
#' Finds matching edge(s) and node(s) for a pruned tree in the original tree from which it was created. This is intended as an internal function, but may be of use to someone.
#'
#' @return
#'
#' \item{matching.edges}{A list of the matching edges.}
#' \item{matching.nodes}{A matrix of matching node numbers.}
#' \item{removed.edges}{A vector of the removed edges.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#' 
#' # Create a random 10-taxon tree:
#' original.tree <- rtree(10)
#' 
#' # Remove three leaves:
#' pruned.tree <- drop.tip(original.tree, c("t1", "t3", "t8"))
#' 
#' # Find matching edges:
#' X <- EdgeMatch(original.tree, pruned.tree)
#' 
#' # Show matching edges:
#' X$matching.edges
#' 
#' # Show removed edges:
#' X$removed.edges
#' 
#' @export EdgeMatch
EdgeMatch <- function(original.tree, pruned.tree) {
	
	# Conditional if pruned tree too small:
	if(ape::Ntip(pruned.tree) < 3) stop("ERROR: pruned.tree includes too few (<3) taxa to be used.")
	
	# Conditional in case where pruned tree taxa are not a subset of the original tree taxa:
	if(length(setdiff(pruned.tree$tip.label, original.tree$tip.label)) > 0) stop("ERROR: pruned.tree cannot include taxa not present in original.tree.")
	
	# First find removed taxa (if any):
	removed.taxa <- setdiff(original.tree$tip.label, pruned.tree$tip.label)
	
	# If no taxa are removed:
	if(length(removed.taxa) == 0) {
		
		# Record removed edges as an empty vector:
		removed.edges <- numeric(0)
		
		# Return matching edges as list:
		matching.edges <- as.list(1:nrow(original.tree$edge))
		
		# Create lists of nodes (which will be identical):
		clades <- corresponding.nodes <- c((ape::Ntip(pruned.tree) + 1):(ape::Ntip(pruned.tree) + ape::Nnode(pruned.tree)), 1:ape::Ntip(pruned.tree))
		
	# If taxa are removed:
	} else {
		
		# Create vector to store clades found in pruned tree:
		clades <- vector(mode="character")
		
		# For each internal node of pruned tree:
		for(i in (ape::Ntip(pruned.tree) + 1):(ape::Ntip(pruned.tree) + ape::Nnode(pruned.tree))) {
			
			# Get descendants of node (members of clade):
			clades <- c(clades, paste(pruned.tree$tip.label[FindDescendants(i, pruned.tree)], collapse = "%%SpLiTtEr%%"))
			
			# Update with node number:
			names(clades)[length(clades)] <- i
			
		}
		
		# Create vector to store corresponding node numbers in original tree:
		corresponding.nodes <- vector(mode = "numeric")
		
		# For each clade in pruned tree find and store node number in original tree:
		for(i in clades) corresponding.nodes <- c(corresponding.nodes, FindAncestor(strsplit(i, "%%SpLiTtEr%%")[[1]], original.tree))
		
		# Add tips to node numbers for original.tree:
		corresponding.nodes <- c(corresponding.nodes, match(pruned.tree$tip.label, original.tree$tip.label))
		
		# Add tips to node numbers for pruned.tree:
		clades <- c(as.numeric(names(clades)), 1:ape::Ntip(pruned.tree))
		
		# Make edge matrix for pruned tree using corresponding nodes in original tree:
		pruned.edges <- cbind(corresponding.nodes[match(pruned.tree$edge[, 1], clades)], corresponding.nodes[match(pruned.tree$edge[, 2], clades)])
		
		# Find edges that match EXACTLY with the original tree:
		matching.edges <- match(apply(pruned.edges, 1, paste, collapse = "%%"), apply(original.tree$edge, 1, paste, collapse = "%%"))
		
		# List non-matching edges for further searching:
		nonmatching.edges <- pruned.edges[is.na(matching.edges), ]
		
		# Only continue if there are non-matching edges (will be the case if only "outgroup(s)" are removed:
		if(length(nonmatching.edges) > 0) {
		
			# Correct stupid matrix to vector problem:
			if(!is.matrix(nonmatching.edges)) nonmatching.edges <- matrix(nonmatching.edges, ncol = 2)
			
			# For each non-matching edge:
			for(i in 1:nrow(nonmatching.edges)) {
				
				# Get start (ancestral) node:
				start.node <- nonmatching.edges[i, 1]
				
				# Get end (terminal) node:
				end.node <- nonmatching.edges[i, 2]
				
				# Create edges vector to store multiple edges that correspond to edge on pruned tree:
				edges <- match(end.node, original.tree$edge[, 2])
				
				# Keep going until start and end node are joined by contiguous edges:
				while(length(sort(match(original.tree$edge[edges, 1], start.node))) == 0) {
					
					# Update end node:
					end.node <- original.tree$edge[match(end.node, original.tree$edge[, 2]), 1]
					
					# Update edges:
					edges <- c(edges, match(end.node, original.tree$edge[, 2]))
					
				}
				
				# Update matching edges with multiple edges separated by a double-percent:
				matching.edges[which(is.na(matching.edges))[1]] <- paste(rev(edges), collapse = "%%")
				
			}

			# Get matching edges as list:
			matching.edges <- lapply(strsplit(matching.edges, "%%"), as.numeric)

		# If there are no non-matching edges:
		} else {
			
			# Get matching edges as list:
			matching.edges <- as.list(matching.edges)
			
		}
		
		# Get removed edges (branches from original tree missing in pruned tree:
		removed.edges <- setdiff(1:nrow(original.tree$edge), as.numeric(unlist(matching.edges)))
		
	}
	
	# Add names to matching edges:
	names(matching.edges) <- 1:nrow(pruned.tree$edge)
	
	# Make list of matching nodes:
	node.matches <- cbind(clades, corresponding.nodes)
	
	# Add column names:
	colnames(node.matches) <- c("Pruned_node", "Original_node")
	
	# Compile output:
	output <- list(matching.edges, node.matches, removed.edges)
	
	# Add output names:
	names(output) <- c("matching.edges", "matching.nodes", "removed.edges")
	
	# Return output:
	return(output)
	
}
