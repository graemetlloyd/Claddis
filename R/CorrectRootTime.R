CorrectRootTime <- function(original.tree, pruned.tree) {
	
	# Conditional if pruned tree too small:
	if(Ntip(pruned.tree) < 3) stop("ERROR: pruned.tree includes too few (<3) taxa to be used.")
	
	# Conditional in case where pruned tree taxa are not a subset of the original tree taxa:
	if(length(setdiff(pruned.tree$tip.label, original.tree$tip.label)) > 0) stop("ERROR: pruned.tree cannot include taxa not present in original.tree.")
	
	# Update $root.time for pruned.tree:
	pruned.tree$root.time <- original.tree$root.time - mean(diag(vcv(original.tree))[names(diag(vcv(pruned.tree)))] - diag(vcv(pruned.tree)))
	
	# Return updated pruned tree:
	return(pruned.tree)
	
}
