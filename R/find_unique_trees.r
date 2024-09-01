#' Finds only the unique topologies amongst a set
#'
#' @description
#'
#' Given a set of trees with the same tip labels, returns just the unique topologies present.
#'
#' @param trees An object of class \code{multiPhylo}.
#'
#' @details
#'
#' Where labelled topologies are generated randomly or modified by (e.g.) removing a tip, it may be useful to isolate just those that are truly unique. The \code{ape} package already has a function for this (\link[ape]{unique.multiPhylo}), but it can be slow when the number of trees is large. This function is thus intended as a faster version.
#'
#' The function works by breaking down a tree into its' component bipartitions and treating the combination of these as the definition of the tree. It thus escapes problems due to the principle of free rotation. Specifically, these two trees are actually identical:
#'
#' \preformatted{A  B  C   D  E
#'  \/    \   \/
#'   \     \  /
#'    \     \/
#'     \    /
#'      \  /
#'       \/
#'
#' B  A  D  E   C
#'  \/    \/   /
#'   \     \  /
#'    \     \/
#'     \    /
#'      \  /
#'       \/}
#'
#' This becomes clearer if we decompose them into their bipartitions:
#'
#' AB, DE, CDE, ABCDE
#'
#' These correspond to the descendants of each internal node (branching point) and the last one is actually ignored (the root node) as it will be present in any tree.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' An object of class \code{"multiPhylo"}.
#'
#' @seealso
#'
#' \link[ape]{unique.multiPhylo}
#'
#' @examples
#'
#' # Make a set of three identical trees (differing only in "rotation" of nodes):
#' trees <- ape::read.tree(text = c(
#'   "((A,B),(C,(D,E)));",
#'   "((C,(D,E)),(A,B));",
#'   "((B,A),(C,(E,D)));")
#' )
#'
#' # Show that there is only one unique tree:
#' find_unique_trees(trees = trees)
#'
#' @export find_unique_trees
find_unique_trees <- function(trees) {
    
  # Checks to add:
  # - trees should be multiPhylo
  # - trees should have same tip labels
    
  # Get number of tips (assumes all trees have same tips):
  n_tips <- ape::Ntip(trees[[1]])
    
  # Return just unique topologies:
  trees[!duplicated(
    x = unlist(
      x = lapply(
        X = trees,
        FUN = function(tree) {
          paste(sort(
            x = unlist(x = lapply(
              X = as.list(x = n_tips + 2:tree$Nnode),
              FUN = function(node) {
                paste(sort(x = tree$tip.label[strap::FindDescendants(
                  n = node,
                  tree = tree
                )]), collapse = "%%")
              }
            ))
          ), collapse = "&&")
        }
      )
    )
  )]
}
