#' Drop tips from a time-scaled tree
#'
#' @description
#'
#' Drop tips from a time-scaled tree and update root.time accordingly
#'
#' @param time_tree A time-scaled tree in phylo format where branch lengths are durations and where a \code{$root.time} value indicates the root age.
#' @param tip_names A vector of tip names to be pruned from the tree.
#' @param ... Additional options to be passed to \code{ape::drop.tip}.
#'
#' @details
#'
#' (NB: This function is designed to only cope with trees containing at least three tips.)
#'
#' Usually ape formatted trees are pruned with the \link[ape]{drop.tip} function in \link[ape]{ape}. However, trees time-scaled using either the \code{paleotree} or \code{strap} packages have an additional important component, the root age (\code{$root.time}) that may need updating when tips are removed. (See \link{fix_root_time}.) Thus this function is a modified version of \link[ape]{drop.tip} that also performs the \link{fix_root_time} step.
#'
#' Note that \code{dropPaleoTip} in the \code{paleotree} package performs the exact same function, but is not called here to reduce the number of dependencies for \code{Claddis}.
#'
#' @return Returns a tree (phylo object) with pruned tips and corrected \code{$root.time}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{fix_root_time}
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
#' pruned_tree <- drop_time_tip(time_tree = tree, tip_names = "A")
#'
#' # Show that the root time is now fixed (19 Ma):
#' pruned_tree$root.time
#'
#' @export drop_time_tip
drop_time_tip <- function(time_tree, tip_names, ...) {
  
  # Add some top-level conditional checks in future (including checking that root.time is even there).
  
  # First generate pruned time tree:
  pruned_time_tree <- ape::drop.tip(phy = time_tree, tip = tip_names, ...)
  
  # Return tree with tips oruend and root rescaled:
  fix_root_time(original_tree = time_tree, pruned_tree = pruned_time_tree)
}
