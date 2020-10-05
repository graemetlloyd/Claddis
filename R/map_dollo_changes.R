#' Stochastic Character Map For Dollo Character
#'
#' @description
#'
#' Given a tree with binary tip states produces a stochastic Dollo character map.
#'
#' @param time_tree A tree in phylo format with positive branch lengths and a value for \code{$root.time}.
#' @param tip_states A named vector of tip states (must be 0 or 1), where the names match \code{tree$tip.label}.
#'
#' @details
#'
#' The non-ideal solution from Tarver et al. (2018) to the problem of generating a stochastic character map for a Dollo character (i.e., a single gain of the derived state, 1) with any number of losses (1 -> 0).
#'
#' The function operates as follows:
#'
#' 1) Establishes the least inclusive clade exhibiting the derived state (1).
#' 2) Assumes a single gain occurred with equal probability along the branch subtending this clade.
#' 3) Prunes the inclusive clade to generate a subtree with a strong root prior of the derived state (1).
#' 4) Calls \code{make.simmap} from the \code{phytools} package to generate a stochastic character map using a model where only losses are possible.
#' 5) Outputs both the stochastic character map (time spent in each state on each branch) and a matrix of state changes.
#'
#' NB: As the map is stochastic the answer will be different each time the function is run and multiple replicates are strongly advised in order to characterise this uncertainty.
#'
#' @return
#'
#' \item{changes}{A matrix of all changes (gains and losses).}
#' \item{stochastic_character_map}{The stochastic character map.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Tarver, J. E., Taylor, R. S., Puttick, M. N., Lloyd, G. T., Pett, W., Fromm, B., Schirrmeister, B. E., Pisani, D., Peterson, K. J. and Donoghue, P. C. J., 2018. Well-annotated microRNAomes do not evidence pervasive miRNA loss. \emph{Genome Biology and Evolution}, \bold{6}, 1457-1470.
#'
#' @examples
#'
#' # Build example ten-tip tree:
#' time_tree <- ape::read.tree(text = paste0("(A:1,(B:1,((C:1,(D:1,(E:1,F:1):1):1):1,",
#'   "((G:1,H:1):1,(I:1,J:1):1):1):1):1);"))
#'
#' # Arbitrarily add a root.time value of 100 Ma:
#' time_tree$root.time <- 100
#'
#' # Build example tip state values:
#' tip_states <- c(A = 0, B = 0, C = 1, D = 1, E = 0, F = 1, G = 1, H = 1, I = 0, J = 1)
#'
#' # Run map_dollo_changes on data and store output:
#' out <- map_dollo_changes(time_tree, tip_states)
#'
#' # View matrix of changes:
#' out$changes
#'
#' # View stochastic character map (time spent in each state on each branch):
#' out$stochastic_character_map
#' @export map_dollo_changes
map_dollo_changes <- function(time_tree, tip_states) {

  # Output of this should match other stochastic character map function and ultimately be a format that can be handed to test_rates or disparity functions

  # Find number of tips:
  n_tips <- ape::Ntip(phy = time_tree)

  # Find number of nodes:
  n_nodes <- ape::Nnode(phy = time_tree)

  # Check tree has branch lengths:
  if (is.null(time_tree$edge.length)) stop("time_tree must have branch lengths.")

  # Check tree has a root time:
  if (is.null(time_tree$root.time)) stop("time_tree must have a $root.time value.")

  # Record non-matching names:
  nonmatching_names <- c(setdiff(x = names(tip_states), y = time_tree$tip.label), setdiff(x = time_tree$tip.label, y = names(tip_states)))

  # If there are any non-matching names stop and notify user:
  if (length(x = nonmatching_names) > 0) stop(paste("Non-matching names between tree and tips states found:", paste(nonmatching_names, collapse = ", ")))

  # Check for states other than zero or one:
  nonbinary_states <- setdiff(x = unique(x = tip_states), y = c(0, 1))

  # If states other than zero or one found warn user:
  if (length(x = nonbinary_states) > 0) stop("States other than 0 or 1 found. All states must be 0 or 1.")

  # Make Dollo-like model where only losses (1 -> 0) are allowed (will later force root state to be one):
  dollo_model <- matrix(c(0, 1, 0, 0), nrow = 2, ncol = 2, dimnames = list(c("0", "1"), c("0", "1")))

  # Find new root (indicating least inclusive clade of
  new_root <- find_mrca(descendant_names = names(which(x = tip_states == 1)), tree = time_tree)

  # If new root is really the old root:
  if ((n_tips + 1) == new_root && sum(tip_states == 1) > 1) {

    # Set acqusition branch to zero as precedes root:
    acquisition_branch <- 0

    # Set acquisition time to root time:
    acquisition_time <- time_tree$root.time

    # If tip states vary:
    if (length(x = unique(x = tip_states)) > 1) {

      # Get stochastic character map for full tree using Dollo model and a strong root prior of one:
      stochastic_character_map <- phytools::make.simmap(tree = time_tree, x = tip_states[time_tree$tip.label], model = dollo_model, pi = c(0, 1), message = FALSE)$maps

      # If tip states are constant:
    } else {

      # Create SCM with no losses:
      stochastic_character_map <- as.list(x = time_tree$edge.length)

      # Label all states as derived:
      for (i in 1:length(x = stochastic_character_map)) names(stochastic_character_map[[i]]) <- as.character(unique(x = tip_states))
    }

    # If new root reflects a subtree:
  } else {

    # Case if character is an autapomorphy:
    if (sum(tip_states == 1) == 1) {

      # Create base stochastic character map (SCM:
      stochastic_character_map <- as.list(x = c(1:nrow(time_tree$edge)))

      # For each branch in the SCM:
      for (i in 1:length(x = stochastic_character_map)) {

        # Create a null value (vector equal to branch length):
        x <- time_tree$edge.length[i]

        # Label null vector with null (root) state of zero:
        names(x) <- "0"

        # Update base SCM with null value:
        stochastic_character_map[[i]] <- x
      }

      # Get branch along which (single) acqusition of derived state occurs:
      acquisition_branch <- match(which(x = tip_states == 1), time_tree$edge[, 2])

      # Get bounds for acquisition times:
      acquisition_bounds <- date_nodes(time_tree = time_tree)[time_tree$edge[acquisition_branch, ]]

      # Draw an acusition time from a uniform distribution between the bounds:
      acquisition_time <- stats::runif(n = 1, min = acquisition_bounds[2], max = acquisition_bounds[1])

      # Create an SCM branch for the acquisition:
      acquisition_branch_stochastic_character_map <- c(acquisition_bounds[1] - acquisition_time, acquisition_time - acquisition_bounds[2])

      # Add labels to the acquisition branch:
      names(acquisition_branch_stochastic_character_map) <- c("0", "1")

      # Add acquisition branch to the SCM:
      stochastic_character_map[[acquisition_branch]] <- acquisition_branch_stochastic_character_map

      # Case if at least two taxa exhibit the derived state:
    } else {

      # Find members of the least inclusive clade:
      clade_members <- time_tree$tip.label[strap::FindDescendants(n = new_root, tree = time_tree)]

      # Find non-members of the least inclusive clade:
      nonclade_members <- setdiff(x = time_tree$tip.label, y = clade_members)

      # Case if only two members in clade:
      if (length(x = clade_members) == 2) {

        # Get branch along which (single) acqusition of derived state occurs:
        acquisition_branch <- match(new_root, time_tree$edge[, 2])

        # Get bounds for acquisition times:
        acquisition_bounds <- date_nodes(time_tree = time_tree)[time_tree$edge[acquisition_branch, ]]

        # Draw an acusition time from a uniform distribution between the bounds:
        acquisition_time <- stats::runif(n = 1, min = acquisition_bounds[2], max = acquisition_bounds[1])

        # Create base stochastic character map (SCM):
        stochastic_character_map <- as.list(x = c(1:nrow(time_tree$edge)))

        # For each branch in the SCM:
        for (i in 1:length(x = stochastic_character_map)) {

          # Create a null value (vector equal to branch length):
          x <- time_tree$edge.length[i]

          # Label null vector with null (root) state of zero:
          names(x) <- "0"

          # Update base SCM with null value:
          stochastic_character_map[[i]] <- x
        }

        # Create an SCM branch for the acquisition:
        acquisition_branch_stochastic_character_map <- c(acquisition_bounds[1] - acquisition_time, acquisition_time - acquisition_bounds[2])

        # Add labels to the acqusition branch:
        names(acquisition_branch_stochastic_character_map) <- c("0", "1")

        # Add acquisition branch to the SCM:
        stochastic_character_map[[acquisition_branch]] <- acquisition_branch_stochastic_character_map

        # Get two descending edges:
        descendant_edges <- find_descendant_edges(new_root, time_tree)

        # Update state of descending edges to derived (1):
        for (i in descendant_edges) names(SCM[[i]]) <- "1"

        # Case if at least three members in clade:
      } else {
        
        # Prune taxa external to least inclusive clade to create a pruned time tree:
        new_tree <- drop_time_tip(time_tree = time_tree, tip_names = nonclade_members)
        
        # Get tip and node counts for new tree:
        n_new_tips <- ape::Ntip(phy = new_tree)
        n_new_nodes <- ape::Nnode(phy = new_tree)

        # Update tip states for new pruned tree:
        new_tips <- tip_states[clade_members]

        # Get branch along which (single) acqusition of derived state occurs:
        acquisition_branch <- match(new_root, time_tree$edge[, 2])

        # Get bounds for acquisition times:
        acquisition_bounds <- date_nodes(time_tree = time_tree)[time_tree$edge[acquisition_branch, ]]

        # Draw an acusition time from a uniform distribution between the bounds:
        acquisition_time <- stats::runif(n = 1, min = acquisition_bounds[2], max = acquisition_bounds[1])

        # Create base stochastic character map (SCM):
        stochastic_character_map <- as.list(x = c(1:nrow(time_tree$edge)))

        # For each branch in the SCM:
        for (i in 1:length(x = stochastic_character_map)) {

          # Create a null value (vector equal to branch length):
          x <- time_tree$edge.length[i]

          # Label null vector with null (root) state of zero:
          names(x) <- "0"

          # Update base SCM with null value:
          stochastic_character_map[[i]] <- x
        }

        # Create an SCM branch for the acquisition:
        acquisition_branch_stochastic_character_map <- c(acquisition_bounds[1] - acquisition_time, acquisition_time - acquisition_bounds[2])

        # Add labels to the acqusition branch:
        names(acquisition_branch_stochastic_character_map) <- c("0", "1")

        # Add acquisition branch to the SCM:
        stochastic_character_map[[acquisition_branch]] <- acquisition_branch_stochastic_character_map

        # If tip states vary:
        if (length(x = unique(x = new_tips)) > 1) {

          # Now do real SCM on pruned tree using Dollo model and a strong root prior of one:
          real_stochastic_character_map <- phytools::make.simmap(tree = new_tree, x = new_tips[new_tree$tip.label], model = dollo_model, pi = c(0, 1), message = FALSE)$maps

          # If tip state is constant:
        } else {

          # Create SCM with no losses:
          real_stochastic_character_map <- as.list(x = new_tree$edge.length)

          # Label all states as derived:
          for (i in 1:length(x = real_stochastic_character_map)) names(real_stochastic_character_map[[i]]) <- "1"
        }

        # Create edge matrix for original tree:
        original_edges <- time_tree$edge

        # Create edge matrix for pruned tree:
        new_edges <- new_tree$edge

        # Update tip names for original tree edge matrix:
        for (i in 1:n_tips) original_edges[which(x = original_edges[, 2] == i), 2] <- time_tree$tip.label[i]

        # Update tip names for pruned tree edge matrix:
        for (i in 1:n_new_tips) new_edges[which(x = new_edges[, 2] == i), 2] <- new_tree$tip.label[i]

        # Update node names for original tree edge matrix:
        for (i in (n_tips + 1):(n_tips + n_nodes)) original_edges[which(x = original_edges == i)] <- paste(sort(x = time_tree$tip.label[strap::FindDescendants(n = i, tree = time_tree)]), collapse = "")

        # Update node names for pruned tree edge matrix:
        for (i in (n_new_tips + 1):(n_new_tips + n_new_nodes)) new_edges[which(x = new_edges == i)] <- paste(sort(x = new_tree$tip.label[strap::FindDescendants(n = i, tree = new_tree)]), collapse = "")

        # Collapse original edge matrix to from-to straings for matching:
        original_edges <- apply(original_edges, 1, paste, collapse = "%%TO%%")

        # Collapse pruned edge matrix to from-to straings for matching:
        new_edges <- apply(new_edges, 1, paste, collapse = "%%TO%%")

        # Match edges between pruned and original trees:
        edge_matches <- match(new_edges, original_edges)

        # Store real SCM in SCM for full tree:
        stochastic_character_map[edge_matches] <- real_stochastic_character_map
      }
    }
  }

  # Get node ages for tree in advance of establishing character change times:
  all_node_ages <- date_nodes(time_tree = time_tree)

  # Create changes matrix:
  changes_matrix <- matrix(nrow = 0, ncol = 4, dimnames = list(c(), c("branch", "from", "to", "time")))

  # Get branches that record changes:
  change_branches <- which(x = lapply(X = stochastic_character_map, length) > 1)

  # Special case of acquisition occurring prior to root (add extra line to changes matrix):
  if (acquisition_branch == 0) changes_matrix <- rbind(changes_matrix, c(acquisition_branch, 0, 1, acquisition_time))

  # As long as there are changes to record (i.e., it is not a constant character):
  if (length(x = change_branches) > 0) {

    # For each change branch:
    for (i in change_branches) {

      # Record branch number:
      change_branch <- i

      # Record from and to states:
      change_states <- as.numeric(names(stochastic_character_map[[i]]))

      # Record change time:
      change_time <- all_node_ages[time_tree$edge[i, 1]] - stochastic_character_map[[i]][1]

      # Make changes line ready for inserting into matrix:
      change_line <- c(change_branch, change_states, change_time)

      # Add changes to matrix:
      changes_matrix <- rbind(changes_matrix, change_line)
    }
  }

  # Clean up row names:
  rownames(x = changes_matrix) <- NULL

  # Return output as list:
  list(changes = changes_matrix, stochastic_character_map = stochastic_character_map)
}
