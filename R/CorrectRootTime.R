#' Corrects root.time after taxa have been pruned from a tree
#'
#' @description
#'
#' Corrects root.time after taxa have been pruned from a tree using drop.tip
#' 
#' @param original.tree A tree in phylo format.
#' @param pruned.tree A tree in phylo format that represents a pruned version of \code{original.tree}.
#'
#' @details
#'
#' (NB: This function is designed to only cope with trees containing at least three tips.)
#'
#' When removing taxa from a time-scaled tree using \link{drop.tip} in \link{ape} \code{$root.time} is left unchanged. This can cause downstream problems if not corrected and that is what this function does.
#'
#' Note that \code{fixRootTime} in the \code{paleotree} package performs the same function, but is not called here to reduce the number of libraries on which \code{Claddis} is dependent. Interested users should also refer to the \code{dropPaleoTip} function in \code{paleotree}.
#'
#' @return Returns a tree (phylo object) with a corrected \code{$root.time}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#' 
#' # Create a simple four-taxon tree with branch lengths:
#' tree <- read.tree(text = "(A:1,(B:1,(C:1,D:1):1):1);")
#' 
#' # Set root age as 20 Ma:
#' tree$root.time <- 20
#' 
#' # Now prune taxon A:
#' pruned.tree <- drop.tip(tree, "A")
#' 
#' # Show that drop.tip has not updated the tree's root time:
#' pruned.tree$root.time
#' 
#' # Use the function to correct the root time:
#' pruned.tree <- CorrectRootTime(tree, pruned.tree)
#' 
#' # Show that the root time is now correct (19 Ma):
#' pruned.tree$root.time
#' 
#' @export CorrectRootTime
CorrectRootTime <- function(original.tree, pruned.tree) {
	
	# Conditional if pruned tree too small:
	if(ape::Ntip(pruned.tree) < 3) stop("ERROR: pruned.tree includes too few (<3) taxa to be used.")
	
	# Conditional in case where pruned tree taxa are not a subset of the original tree taxa:
	if(length(setdiff(pruned.tree$tip.label, original.tree$tip.label)) > 0) stop("ERROR: pruned.tree cannot include taxa not present in original.tree.")
	
	# Update $root.time for pruned.tree:
	pruned.tree$root.time <- original.tree$root.time - mean(diag(vcv(original.tree))[names(diag(vcv(pruned.tree)))] - diag(vcv(pruned.tree)))
	
	# Return updated pruned tree:
	return(pruned.tree)
	
}
