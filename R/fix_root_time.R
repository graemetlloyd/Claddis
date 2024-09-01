#' Fixes root.time after taxa have been pruned from a tree
#'
#' @description
#'
#' Fixes root.time after taxa have been pruned from a tree using ape::drop.tip
#'
#' @param original_tree A tree in phylo format.
#' @param pruned_tree A tree in phylo format that represents a pruned version of \code{original_tree}.
#'
#' @details
#'
#' (NB: This function is designed to only cope with trees containing at least three tips.)
#'
#' When removing taxa from a time-scaled tree using \link[ape]{drop.tip} in \link[ape]{ape} \code{$root.time} is left unchanged. This can cause downstream problems if not fixed and that is what this function does.
#'
#' Note that \code{fix_root_time} in the \code{paleotree} package performs the same function, but is not called here to reduce the number of libraries on which \code{Claddis} is dependent. Interested users should also refer to the \code{dropPaleoTip} function in \code{paleotree}.
#'
#' @return Returns a tree (phylo object) with a fixed \code{$root.time}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{drop_time_tip}
#'
#' @examples
#'
#' # Create a simple four-taxon tree with branch lengths:
#' tree <- ape::read.tree(text = "(A:1,(B:1,(C:1,D:1):1):1);")
#'
#' # Set root age as 20 Ma:
#' tree$root.time <- 20
#'
#' # Now prune taxon A:
#' pruned_tree <- ape::drop.tip(phy = tree, tip = "A")
#'
#' # Show that drop.tip has not updated the tree's root time:
#' pruned_tree$root.time
#'
#' # Use the function to fix the root time:
#' pruned_tree <- fix_root_time(original_tree = tree, pruned_tree = pruned_tree)
#'
#' # Show that the root time is now fixed (19 Ma):
#' pruned_tree$root.time
#' @export fix_root_time
fix_root_time <- function(original_tree, pruned_tree) {

  # Conditional if pruned tree too small:
  if (ape::Ntip(phy = pruned_tree) < 3) stop("pruned_tree includes too few (<3) taxa to be used.")

  # Conditional in case where pruned tree taxa are not a subset of the original tree taxa:
  if (length(x = setdiff(x = pruned_tree$tip.label, y = original_tree$tip.label)) > 0) stop("pruned_tree cannot include taxa not present in original_tree.")

  # Update $root.time for pruned_tree:
  pruned_tree$root.time <- original_tree$root.time - mean(diag(x = ape::vcv(phy = original_tree))[names(diag(x = ape::vcv(phy = pruned_tree)))] - diag(x = ape::vcv(pruned_tree)))

  # Return updated pruned tree:
  pruned_tree
}
