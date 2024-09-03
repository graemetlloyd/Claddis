#' Estimate ancestral states for a continuous character under squared-change parsimony
#'
#' @description
#'
#' Given a phylogeny and set of continuous tip values, returns a single estimate of ancestral states that minimise the total squared-change cost on the tree.
#'
#' @param tree A tree (object of class \code{phylo}) with or without branch lengths and having tip labels that match the elements of \code{tip_values}.
#' @param tip_values A labelled vector of continuous tip values, where the labels match the tip labels in \code{tree}. Tip values can be single numbers or a minimum-maximum range, see details.
#'
#' @details
#'
#' Multiple algorithms exist for estimating ancestral states for a univariate continuous character. This function specifically provides an ancestral state estimation that minimises the squared-change cost between each internal node and every other node that node is directly connected to. Note: in practice this approach can be identical to maximum likelihood ancestral state estimation under Brownian motion (e.g., Maddison 1991), leading Goloboff (2022) to argue that it isn't strictly a parsimony approach. However, this function extends squared-change parsimony to allow for ranges in tip values (see below) as well as implementing squared-change parsimony where branch lengths are variable and (hard) polytomies are permitted. As such it may still be of use to some users as distinct from other ancestral state estimation approaches for continuous characters.
#'
#' \bold{Algorithm}
#'
#' Although squared-change parsimony ancestral state estimates can be directly calculated using the approach of Maddison (1991), the algorithm used here is an optimisation approach as this allows tip ranges to be accommodated (see below). In practice this leads to minimal speed reduction in most cases as an initial estimate is made using the \code{fastAnc} function from the phytools package (Revell 2024) that will frequently be identical to the final solution (except where ranges in tip values are used).
#'
#' The optimisation used here is the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm which was independently developed by Broyden (1970), Fletcher (1970), Goldfarb (1970), and Shanno (1970) and is implemented as the \code{"L-BFGS-B"} method in the R \code{stats} function \code{optim()}. The choice of this method is what allows ranges in tip values to be accommodated by treating these as an additional parameter to be estimated but with box constraints (i.e., minimum and maximum possible values - the ranges for the tip value) applied.
#'
#' \bold{Ranges in tip values}
#'
#' Unlike other implementations of squared-change parsimony (to the best of my knowledge at least) this function removes the constraint that each tip must be represented by a single numeric value. Instead, any number of tip values can instead be represented by a range of values (i.e., a minimum and maximum value). These must be specified in the \code{tip_values} variable by separating the minimum and maximum values with an underscore character (_). E.g., \code{"0.33_0.53"} is the range 0.33 to 0.53. Note: an underscore is used rather than a dash as negative values are permitted and the dash symbol is reserved to indicate these.
#'
#' In practice ranges are treated as an additional paramater to be estimated (alongside the internal node estimates) with the same criterion of minimising the total length, or cost (cum of squares), for the tree as a whole. As such this is truly a parsimony algorithm and negates the criticism of Goloboff (2022) that squared-change parsimony is not a true parsimony algorithm and also differs from (say) a maximum likelihood estimate that assumes a central tendency (i.e., midpoint) to the range of values at a tip.
#'
#' The estimated tip values can also be examined a posteriori using the \code{tip_values} value from the function output.
#'
#' \bold{Polytomies in the input tree}
#'
#' The function can also handle polytomies in the input tree, albeit it only does so by treating these as "hard" polytomies. I.e., it does not permute the minimum sum of squares for every possible resolution of a polytomy.
#'
#' \bold{"Weighted" squared-change parsimony and branch lengths in the input tree}
#'
#' By default the function will apply so-called "unweighted" squared-change parsimony, meaning the user does not need to provide an input tree with branch lengths (only the topology). However, in practice "unweighted" parsimony really means treating every branch as having unit length (which is still \emph{a} weighting). However, a user may wish to estimate ancestral states where variable branch lengths are accounted for (e.g., such that a cherry with variable terminal branch lengths means the subtending node is more strongly influenced by the tip value at the end of the shorter terminal).
#'
#' In practice the function will implement weighted squared-change parsimony whenever the input tree already has edge-lengths applied. As such if the user does wish to implement "unweighted" squared-change parsimony they should be careful to supply an input tree where there are either no edge lengths or every edge is set to length one.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Broyden, C. G., 1970. The convergence of a class of double-rank minimization algorithms. \emph{Journal of the Institute of Mathematics and Its Applications}, \bold{6}, 76-90.
#'
#' Fletcher, R., 1970. A new approach to variable metric algorithms. \emph{Computer Journal}, \bold{13}, 317-322.
#'
#' Goloboff, P. A., 2022. \emph{Phylogenetic Analysis of Morphological Data, Volume 2: Refining Phylogenetic Analyses}. CRC Press, Boca Raton. 291 pp.
#'
#' Goldfarb, D., 1970. A family of variable metric updates derived by variational means. \emph{Mathematics of Computation}, \bold{24}, 23-26.
#'
#' Maddison, W. P., 1991. Squared-change parsimony reconstructions of ancestral states for continuous-valued characters on a phylogenetic tree. \emph{Systematic Zoology}, \bold{40}, 304-314.
#'
#' Revell, L. J., 2024. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
#'
#' Shanno, D. F., 1970. Conditioning of quasi-Newton methods for function minimization. \emph{Mathematics of Computation}, \bold{24}, 647-656.
#'
#' @return
#'
#' The function returns a list with three elements:
#'
#' \item{ancestral_state_estimates}{A single set of ancestral state estimates that represent a solution that minimises the total cost, or length, of the tree.}
#' \item{tip_values}{The (estimated) tip values that minimise the sum of squares. This will be identical to the \code{tip_values} input if no ranges in tip values were used, but otherwise will show a single tip value that minimised the total cost, or length, of the tree.}
#' \item{sum_of_squares}{A single numeric value indicating the total cost, or length, of the tree which is the sum of the squared length of each edge of the tree.}
#'
#' @examples
#'
#' # Estimate squared-change parsimony for four tips with single values:
#' scp <- estimate_squared_change_ancestors(
#'   tree = ape::read.tree(text = "((A,B),(C,D));"),
#'   tip_values = tip_values <- c("A" = 3.6, "B" = 4.1, "C" = 7.3, "D" = 8.8)
#' )
#'
#' # Plot results:
#' ape::plot.phylo(
#'   x = ape::read.tree(text = "((A,B),(C,D));"),
#'   main = paste("Sum of squares =", scp$sum_of_squares)
#' )
#' ape::tiplabels(text = scp$tip_values)
#' ape::nodelabels(text = round(x = scp$ancestral_state_estimates, digits = 2))
#'
#' # Estimate squared-change parsimony for four tips with a ranged value:
#' scp <- estimate_squared_change_ancestors(
#'   tree = ape::read.tree(text = "((A,B),(C,D));"),
#'   tip_values = tip_values <- c("A" = "3.6", "B" = "4.1", "C" = "7.3", "D" = "8.8_9.3")
#' )
#'
#' # Plot results:
#' ape::plot.phylo(
#'   x = ape::read.tree(text = "((A,B),(C,D));"),
#'   main = paste("Sum of squares =", scp$sum_of_squares)
#' )
#' ape::tiplabels(text = scp$tip_values)
#' ape::nodelabels(text = round(x = scp$ancestral_state_estimates, digits = 2))
#'
#' # Estimate weighted squared-change parsimony for four tips with single values:
#' scp <- estimate_squared_change_ancestors(
#'   tree = ape::read.tree(text = "((A:1,B:2):2,(C:2,D:3):4);"),
#'   tip_values = tip_values <- c("A" = 3.6, "B" = 4.1, "C" = 7.3, "D" = 8.8)
#' )
#'
#' # Plot results:
#' ape::plot.phylo(
#'   x = ape::read.tree(text = "((A:1,B:2):2,(C:2,D:3):4);"),
#'   main = paste("Sum of squares =", scp$sum_of_squares)
#' )
#' ape::tiplabels(text = scp$tip_values)
#' ape::nodelabels(text = round(x = scp$ancestral_state_estimates, digits = 2))
#'
#' # Estimate squared-change parsimony for four tips with single values and a single polytomy:
#' scp <- estimate_squared_change_ancestors(
#'   tree = ape::read.tree(text = "(A,(B,C,D));"),
#'   tip_values = tip_values <- c("A" = 3.6, "B" = 4.1, "C" = 7.3, "D" = 8.8)
#' )
#'
#' # Plot results:
#' ape::plot.phylo(
#'   x = ape::read.tree(text = "(A,(B,C,D));"),
#'   main = paste("Sum of squares =", scp$sum_of_squares)
#' )
#' ape::tiplabels(text = scp$tip_values)
#' ape::nodelabels(text = round(x = scp$ancestral_state_estimates, digits = 2))
#'
#' @export estimate_squared_change_ancestors
estimate_squared_change_ancestors <- function(tree, tip_values) {
  
  # CHECK NAMES MATCH BETWEEN X AND TREE AND ARE SAME LENGTH!
  
  # Reorder x by tree label:
  tip_values <- tip_values[tree$tip.label]
  
  # Set number of tips:
  n_tips <- ape::Ntip(tree)
  
  # Set number of edges:
  n_edges <- nrow(x = tree$edge)
  
  # If tree lacks edge lengths then set all edges to one ("unweighted" SCP):
  if (is.null(x = tree$edge.length)) tree$edge.length[1:nrow(x = tree$edge)] <- 1
  
  # Check edge lengths are all positive:
  if (any(tree$edge.length <= 0)) stop("Branch lengths on tree must all be positive!")
  
  # Make a new x out of the means of x:
  tip_value_means <- unlist(x = lapply(X = strsplit(x = sapply(X = tip_values, FUN = as.character), split = "_"), FUN = function(i) mean(x = as.numeric(x = i))))
  
  # If only one node (star tree special case):
  if (tree$Nnode == 1) {
    
    # Set root as mean of means:
    init <- c(mean(x = tip_value_means))
    
    # Add root number:
    names(x = init) <- n_tips + 1
  }
  
  # If at least two internal nodes use phytools fasAnc to get a quick initial value for each internal node:
  if (tree$Nnode > 1) init <- phytools::fastAnc(tree = tree, x = tip_value_means)

  # Initialise tips with ranges as FALSE:
  tips_with_ranges <- FALSE

  # Initialise range minima and maxima:
  range_minima <- range_maxima <- c()
  
  # If there are tips with ranges:
  if (length(x = grep(pattern = "_", x = tip_values)) > 0) {
    
    # Overwrite tips with ranges with TRUE:
    tips_with_ranges <- TRUE
    
    # Get names of tips with ranges:
    ranged_tip_names <- names(x = tip_values[grep(pattern = "_", x = tip_values)])
    
    # Get numbers of ranged tip values:
    ranged_tip_numbers <- unname(obj = sapply(X = ranged_tip_names, FUN = function(i) which(x = tree$tip.label == i)))
    
    # Get just the ranged values:
    ranged_values <- tip_values[grep(pattern = "_", x = tip_values)]
    
    # Get number of ranged tips:
    n_ranged_values <- length(x = ranged_tip_numbers)
    
    # Create initial tip values for ranged tips:
    initial_xs <- tip_value_means[ranged_tip_names]
    
    # Add tip numbers as names to tips:
    names(x = initial_xs) <- ranged_tip_numbers
    
    # Add ranged tips to intialised values:
    init <- c(init, initial_xs)
    
    # Get range minima:
    range_minima <- unname(obj = sapply(X = ranged_values, FUN = function(i) min(x = as.numeric(x = strsplit(x = i, split = "_")[[1]]))))
    
    # Get range maxima
    range_maxima <- unname(obj = sapply(X = ranged_values, FUN = function(i) max(x = as.numeric(x = strsplit(x = i, split = "_")[[1]]))))
  }
  
  # Subfunction to calculate sum of squares:
  sum_of_squares <- function(parameter_estimates, tip_values, tip_value_means, tree, n_tips, n_edges, tips_with_ranges, ranged_tip_numbers) {
    
    # Build intial edge matrix of node values:
    A <- matrix(
      data = c(
        tip_value_means[tree$tip.label],
        parameter_estimates[as.character(1:tree$Nnode + n_tips)]
      )[tree$edge],
      nrow = n_edges,
      ncol = 2
    )
    
    # If some tips are ranges:
    if (tips_with_ranges) {
      
      # Add tip values before calculating sum of sqaures:
      A[match(ranged_tip_numbers, tree$edge[, 2]), 2] <- parameter_estimates[as.character(x = ranged_tip_numbers)]
    }
    
    # Return sum of squares:
    sum((A[, 2] - A[, 1]) ^ 2 / tree$edge.length)
  }
  
  # Optimise fit to minimum sum of squares:
  fit <- stats::optim(
    par = init,
    fn = sum_of_squares,
    method = "L-BFGS-B",
    tip_values = tip_values,
    tip_value_means = tip_value_means,
    tree = tree,
    n_tips = n_tips,
    n_edges = n_edges,
    tips_with_ranges = tips_with_ranges,
    ranged_tip_numbers = ranged_tip_numbers,
    control = list(maxit = 5000),
    lower = c(rep(x = -Inf, times = tree$Nnode), range_minima),
    upper = c(rep(x = Inf, times = tree$Nnode), range_maxima)
  )
  
  # Update tip va;ues with optimal value if ranges used:
  if (tips_with_ranges) tip_values[ranged_tip_numbers] <- unclass(x = fit$par)[as.character(x = ranged_tip_numbers)]
  
  # Return list of values:
  list(
    ancestral_state_estimates = unclass(x = fit$par)[1:tree$Nnode],
    tip_values = sapply(X = tip_values, FUN = as.numeric),
    sum_of_squares = fit$value
  )
}
