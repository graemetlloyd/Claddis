#' Ancestral Character State Estimation
#'
#' @description
#'
#' Given a tree and a cladistic matrix uses likelihood to estimate the ancestral states for every character.
#'
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @param time_tree A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{cladistic_matrix}.
#' @param estimate_all_nodes Logical that allows the user to make estimates for all ancestral values. The default (\code{FALSE}) will only make estimates for nodes that link coded terminals (recommended).
#' @param estimate_tip_values Logical that allows the user to make estimates for tip values. The default (\code{FALSE}) will only makes estimates for internal nodes (recommended).
#' @param inapplicables_as_missing Logical that decides whether or not to treat inapplicables as missing (TRUE) or not (FALSE, the default and recommended option).
#' @param polymorphism_behaviour One of either "equalp" or "treatasmissing".
#' @param uncertainty_behaviour One of either "equalp" or "treatasmissing".
#' @param threshold The threshold value to use when collapsing marginal likelihoods to discrete state(s).
#' @param all_missing_allowed Logical to allow all missing character values (generally not recommended, hence default is FALSE).
#'
#' @details
#'
#' At its' core the function uses either the \link{rerootingMethod} (Yang et al. 1995) as implemented in the \link{phytools} package (for discrete characters) or the \link{ace} function in the \link{ape} package (for continuous characters) to make ancestral state estimates. For discrete characters these are collapsed to the most likely state (or states, given equal likelihoods or likelihood within a defined \code{threshold} value). In the latter case the resulting states are represented as an uncertainty (i.e., states separated by a slash, e.g., 0/1). This is the method developed for Brusatte et al. (2014).
#'
#' The function can deal with ordered or unordered characters and does so by allowing only indirect transitions (from 0 to 2 must pass through 1) or direct transitions (from 0 straight to 2), respectively. However, more complex step matrix transitions are not currently supported.
#'
#' Ancestral state estimation is complicated where polymorphic or uncertain tip values exist. These are not currently well handled here, although see the \code{fitpolyMk} function in \link{phytools} for a way these could be dealt with in future. The only available options right now are to either treat multiple states as being equally probable of the "true" tip state (i.e., a uniform prior) or to avoid dealing with them completely by treating them as missing (NA) values.
#'
#' It is also possible to try to use phylogenetic information to infer missing states, both for internal nodes (e.g., those leading to missing tip states) and for tips. This is captured by the \code{estimate_all_nodes} and \code{estimate_tip_values} options. These have been partially explored by Lloyd (2018), who cuationed against their use.
#'
#' @return
#'
#' The function will return the same \code{cladistic_matrix}, but with two key additions: 1. Internal nodes (numbered by \link{ape} formatting) will appear after taxa in each matrix block with estimated states coded for them, and 2. The time-scaled tree used will be added to \code{cladistic_matrix} as \code{cladistic_matrix$topper$tree}. Note that if using the \code{estimate_tip_values = TRUE} option then tip values may also be changed from those provided as input.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Thomas Guillerme \email{guillert@@tcd.ie}
#'
#' @references
#'
#' Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. \emph{Current Biology}, 24, 2386-2392.
#'
#' Lloyd, G. T., 2018. Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. \emph{Palaeontology}, \bold{61}, 637-645.
#'
#' Yang, Z., Kumar, S. and Nei, M., 1995. A new method of inference of ancestral nucleotide and amino acid sequences. \emph{Genetics}, \bold{141}, 1641-1650.
#'
#' @examples
#'
#' # Set random seed:
#' set.seed(4)
#'
#' # Generate a random tree for the Day data set:
#' time_tree <- ape::rtree(n = nrow(day_2016$matrix_1$matrix))
#'
#' # Update taxon names to match those in the data matrix:
#' time_tree$tip.label <- rownames(x = day_2016$matrix_1$matrix)
#'
#' # Set root time by making youngest taxon extant:
#' time_tree$root.time <- max(diag(x = ape::vcv(phy = time_tree)))
#'
#' # Use Day matrix as cladistic matrix:
#' cladistic_matrix <- day_2016
#'
#' # Prune most characters out to make example run fast:
#' cladistic_matrix <- prune_cladistic_matrix(cladistic_matrix,
#'   characters2prune = c(2:3, 5:37)
#' )
#'
#' # Estimate ancestral states:
#' estimate_ancestral_states(
#'   cladistic_matrix = cladistic_matrix,
#'   time_tree = time_tree
#' )
#' @export estimate_ancestral_states
estimate_ancestral_states <- function(cladistic_matrix, time_tree, estimate_all_nodes = FALSE, estimate_tip_values = FALSE, inapplicables_as_missing = FALSE, polymorphism_behaviour = "equalp", uncertainty_behaviour = "equalp", threshold = 0.01, all_missing_allowed = FALSE) {

  # How to get predicted tip states for a continuous character? (Phytools answer: http://blog.phytools.org/2013/11/reconstructed-ancestral-tip-states-for.html)
  #   - So basically under ML just inherit state from ancestral node (really this is mean of distribution where sd would grow with duration of branch so to allow the possibility of variance this could also be sampled stochastically
  # How to deal with step matrices?
  # How to deal with models where intermediate tip states are not even in sample
  # Change help file to explain interactions between all options, e.g., if doing all chars then polymorphisms used for discrete, midpoint for continuous etc.
  # Handle all missing/inapplicable case properly
  # Handle only two tips case properly
  # Add Liam Revell polymorphism options

  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")

  # Get number of tips in tree:
  n_tips <- ape::Ntip(phy = time_tree)

  # Get number of nodes in tree:
  n_nodes <- ape::Nnode(phy = time_tree)

  # Catch problem with trees with no branch lengths:
  if (is.null(time_tree$edge.length)) stop("time_tree must have branch lengths.")

  # Catch problem with polytomies:
  if (time_tree$Nnode < (n_tips - 1)) stop("time_tree must be fully bifurcating.")

  # Catch problem with zero-length branches:
  if (any(time_tree$edge.length == 0)) stop("time_tree must not have zero-length branches.")

  # Check for step matrices and stop and warn if found:
  if (length(x = cladistic_matrix$topper$step_matrices) > 0) stop("Function can not currently deal with step matrices.")

  # Check estimate_all_nodes is a logical:
  if (!is.logical(estimate_all_nodes)) stop("estimate_all_nodes must be a logical (TRUE or FALSE).")

  # Check estimate_tip_values is a logical:
  if (!is.logical(estimate_tip_values)) stop("estimate_tip_values must be a logical (TRUE or FALSE).")

  # Check inapplicables_as_missing is a logical:
  if (!is.logical(inapplicables_as_missing)) stop("inapplicables_as_missing must be a logical (TRUE or FALSE).")

  # Check polymorphism_behaviour is a single allowable value:
  if (length(x = polymorphism_behaviour) != 1 || !any(c("equalp", "treatasmissing") == polymorphism_behaviour)) stop("polymorphism_behaviour must be a single value of either, \"equalp\" or \"treatasmissing\".")

  # Check uncertainty_behaviour is a single allowable value:
  if (length(x = uncertainty_behaviour) != 1 || !any(c("equalp", "treatasmissing") == uncertainty_behaviour)) stop("uncertainty_behaviour must be a single value of either, \"equalp\" or \"treatasmissing\".")

  # Check threshold is a numeric value between the limits of zero and one:
  if (!is.numeric(threshold) || threshold > 0.5 || threshold < 0) stop("threshold must be a numeric value between 0 and 0.5.")

  # Collapse matrix to vectors for each character (state and ordering combination):
  collapsed_matrix <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], function(x) apply(rbind(x$matrix, x$ordering), 2, paste, collapse = ""))))

  # Isolate ordering elements:
  ordering <- unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering"))

  # Isolate minimum values:
  minimum_values <- unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values"))

  # Isolate maximum values:
  maximum_values <- unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values"))

  # Store raw original matrix:
  raw_cladistic_matrix <- cladistic_matrix

  # Combine matrix blocks into a single matrix:
  cladistic_matrix <- original_matrix <- do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"))

  # Find any failed name matches:
  failed_name_matches <- c(setdiff(x = rownames(x = cladistic_matrix), y = time_tree$tip.label), setdiff(x = time_tree$tip.label, y = rownames(x = cladistic_matrix)))

  # Check there are no failed name matches and stop and report if found:
  if (length(x = failed_name_matches) > 0) stop(paste("The following names do not match between the tree and matrix: ", paste(sort(x = failed_name_matches), collapse = ", "), ". Check spelling and try again.", sep = ""))

  # If treating inapplicables as missing (and there is at least one inapplicable) replace with NA:
  if (inapplicables_as_missing && length(x = which(x = cladistic_matrix == "")) > 0) cladistic_matrix[which(x = cladistic_matrix == "")] <- NA

  # If treating polymorphisms as missing:
  if (polymorphism_behaviour == "treatasmissing" && length(x = grep("&", cladistic_matrix)) > 0) cladistic_matrix[grep("&", cladistic_matrix)] <- NA

  # If treating uncertainties as missing:
  if (uncertainty_behaviour == "treatasmissing" && length(x = grep("/", cladistic_matrix)) > 0) cladistic_matrix[grep("/", cladistic_matrix)] <- NA

  # Get vector of character numbers where all values are NA:
  dataless_characters <- which(x = apply(cladistic_matrix, 2, function(x) all(is.na(x))))

  # Look for all missing characters and stop and wanr user if found:
  if (!all_missing_allowed && length(x = dataless_characters) > 0) stop(paste0("The following characters are coded as missing across all tips: ", paste0(dataless_characters, collapse = ", "), ". This can arise either because of the input data (in which case it is recommended that the user prune these characters using prune_cladistic_matrix) or because of the chosen options for inapplicables_as_missing, polymorphism_behaviour, and/or uncertainty_behaviour (in which case the user may wish to chose different values for these)."))

  # Convert tip states into a list:
  data_list <- apply(cladistic_matrix, 2, function(x) list(tip_states = x))

  # For each character:
  for (i in 1:length(x = data_list)) {

    # Add minimum value to list:
    data_list[[i]]$minimum_values <- unname(minimum_values[i])

    # Add maximum value to list:
    data_list[[i]]$maximum_values <- unname(maximum_values[i])

    # Add ordering to list:
    data_list[[i]]$ordering <- unname(ordering[i])

    # Add tree to list:
    data_list[[i]]$tree <- time_tree
  }

  # If estimating values for all characters (need to set dummy tip states for missing values):
  if (estimate_all_nodes) {

    # Subfunction to fill missing values (and inapplicables if desired):
    fill_missing <- function(tip_states) {

      # Find which rows correspond to missing states:
      missingRows <- which(x = is.na(tip_states$tip_states))

      # If missing states found:
      if (length(x = missingRows) > 0) {

        # Build missing state by either forming a polymorphism of all possible tip states, or if continuous the midpoint value:
        fill_states <- ifelse(tip_states$ordering == "cont", (tip_states$minimum_values + tip_states$maximum_values) / 2, paste(tip_states$minimum_values:tip_states$maximum_values, collapse = "/"))

        # Insert missing values:
        tip_states$tip_states[missingRows] <- fill_states
      }

      # Return tip states with missing values replaced:
      tip_states
    }

    # Apply fill missing function across all characters:
    data_list <- lapply(X = data_list, fill_missing)
  }

  # Subfunction to prune tips with missing or inapplicable values:
  prune_tips <- function(x) {

    # Find all missing or inapplicable value tip names:
    missing <- names(sort(x = c(which(x = x$tip_states == ""), which(x = is.na(x$tip_states)))))

    # Work out how many tips will be left after pruning:
    n_tips_remaining <- length(x = setdiff(x = names(x$tip_states), y = missing))

    # If there is at least one missing value:
    if (length(x = missing) > 0) {

      # If less than two tips will remain then set tree as NULL:
      if (n_tips_remaining < 2) x$tree <- NULL

      # If at least two tips will remain prune missing values from tree:
      if (n_tips_remaining > 1) x$tree <- ape::drop.tip(phy = x$tree, tip = missing)

      # Collapse tip states:
      x$tip_states <- x$tip_states[setdiff(x = names(x$tip_states), y = missing)]
    }

    # Return pruned output:
    x
  }

  # Prune out missing and inapplicable tips:
  data_list <- lapply(X = data_list, prune_tips)

  # Subfunction to build tip state matrices:
  convert_tip_states_to_matrix <- function(x) {

    # As long asthere is at least one tip state:
    if (length(x = x$tip_states) > 0) {

      # If the character is not continuous (i.e., it is some form of discrete character):
      if (x$ordering != "cont") {

        # Temporarily store tip states so matrix format can overwrite the stored version below:
        tip_states <- x$tip_states

        # Create matrix of tip state probabilities:
        x$tip_states <- matrix(0, nrow = length(x = x$tip_states), ncol = x$maximum_values - x$minimum_values + 1, dimnames = list(names(x$tip_states), x$minimum_values:x$maximum_values))

        # For each character state if a single state is coded store probability as 1:
        for (i in colnames(x = x$tip_states)) x$tip_states[tip_states == i, i] <- 1

        # If there are polymorphisms and/or uncertainties:
        if (length(x = grep("&|/", tip_states)) > 0) {

          # Get polymorphism locations:
          polymorphisms <- grep("&", tip_states)

          # Get uncertainty locations:
          uncertainties <- grep("/", tip_states)

          # If there are polymorphisms and using the "equalp" (equal probability of each state) option:
          if (length(x = polymorphisms) > 0 && polymorphism_behaviour == "equalp") {

            # For each polymorphisms set each state as equally probable:
            for (i in polymorphisms) x$tip_states[i, strsplit(tip_states[i], split = "&")[[1]]] <- 1 / length(x = strsplit(tip_states[i], split = "&")[[1]])
          }

          # If there are uncertainties and using the "equalp" (equal probability of each state) option:
          if (length(x = uncertainties) > 0 && uncertainty_behaviour == "equalp") {

            # For each uncertainty set each state as equally probable:
            for (i in uncertainties) x$tip_states[i, strsplit(tip_states[i], split = "/")[[1]]] <- 1 / length(x = strsplit(tip_states[i], split = "/")[[1]])
          }
        }

        # If a continuous character:
      } else {

        # Simply make tip states the numeric values (should never be a polymorphism) as a vector:
        x$tip_states <- as.numeric(x$tip_states)
      }

      # If tip state has no length (all values are missing):
    } else {

      # Create row-less tip states matrix:
      x$tip_states <- matrix(nrow = 0, ncol = 1, dimnames = list(c(), "0"))
    }

    # Return the revised input in the same list format:
    x
  }

  # Reformat tip states ready for ancestral estimation:
  data_list <- lapply(X = data_list, convert_tip_states_to_matrix)

  # Subfunction to build character model:
  build_character_model <- function(x) {

    # Set default model to equal rates (works for all binary or unordered characters):
    x$model <- "ER"

    # If a character is both ordered and has at least three states:
    if ((x$maximum_values - x$minimum_values) > 1 && x$ordering == "ord") {

      # Get number of states:
      n_states <- (x$maximum_values - x$minimum_values) + 1

      # Build all zero matrix to begin with:
      x$model <- matrix(0, nrow = n_states, ncol = n_states, dimnames = list(x$minimum_values:x$maximum_values, x$minimum_values:x$maximum_values))

      # for each (just) off-diagonal value store 1 (i.e., N steps to move between adjacent states):
      for (i in 2:n_states) x$model[(i - 1), i] <- x$model[i, (i - 1)] <- 1
    }

    # Return full output:
    x
  }

  # Add ancestral state model for each character:
  data_list <- lapply(X = data_list, build_character_model)

  # Subfunction to get ancestral states:
  estimate_ancestral_state <- function(x, estimate_tip_values, threshold) {

    # As long as there is a tree:
    if (!is.null(x$tree)) {

      # If character is continuous:
      if (x$ordering == "cont") {

        # Get ancestral states using ace:
        x$ancestral_states <- ace(x = x$tip_states, phy = x$tree)$ace

        # If character is discrete:
      } else {

        # If invariant character:
        if (ncol(x$tip_states) == 1) {

          # Get number of tips:
          n_tips <- ape::Ntip(phy = x$tree)

          # Set ancestral states as all the same:
          x$ancestral_states <- matrix(rep(x = 1, times = (n_tips + x$tree$Nnode)), ncol = 1, dimnames = list(c(x$tree$tip.label, (n_tips + 1):(n_tips + x$tree$Nnode)), colnames(x = x$tip_states)))
        }

        # If variant character then get ancestral states using rerooting method:
        if (ncol(x$tip_states) > 1) x$ancestral_states <- phytools::rerootingMethod(tree = x$tree, x = x$tip_states, model = x$model)$marginal.anc

        # Reformat to most likely state
        x$ancestral_states <- unlist(x = lapply(X = lapply(X = apply(x$ancestral_states, 1, list), unlist), function(x) {
          paste(names(x[x > (max(x) - threshold)]), collapse = "/")
        }))

        # If not estimating tip values then prune these:
        if (!estimate_tip_values) x$ancestral_states <- x$ancestral_states[-match(x$tree$tip.label, names(x$ancestral_states))]
      }

      # If no tree:
    } else {

      # Set ancestral states as NULL:
      x$ancestral_states <- vector(mode = "character")
    }

    # Return full output of x:
    x
  }

  # Get ancestral states for each character:
  data_list <- lapply(X = data_list, estimate_ancestral_state, estimate_tip_values, threshold)

  # Get Newick strings of all sampled subtrees (to use to avoid redundancy in tree node mapping):
  newick_strings <- unlist(x = lapply(X = data_list, function(x) ifelse(is.null(x$tree), NA, ape::write.tree(x$tree))))

  # Get just unique strings (i.e., the minimum set needded to map everything to the full tree):
  unique_newick_strings <- unique(x = newick_strings[!is.na(newick_strings)])

  # Convert unique Newick strings to unique trees:
  unique_trees <- ape::read.tree(text = unique_newick_strings)

  # If only a single tree reformat as a list:
  if (inherits(unique_trees, what = "phylo")) unique_trees <- list(unique_trees)

  # Subfunction to map nodes from pruned tree to full tree:
  map_to_full_tree <- function(pruned_tree, full_tree) {

    # Get number of tips of pruned tree:
    n_tips <- ape::Ntip(phy = pruned_tree)

    # Get number of nodes of pruned tree:
    n_nodes <- ape::Nnode(phy = pruned_tree)

    # Get all internal node numbers for pruned tree:
    node_numbers <- (n_tips + 1):(n_tips + n_nodes)

    # If the pruned tree is different to the full tree:
    if (write.tree(pruned_tree) != write.tree(full_tree)) {

      # Get descendants of each node in pruned tree:
      node_descendants <- lapply(X = as.list(x = node_numbers), function(x) pruned_tree$tip.label[strap::FindDescendants(n = x, tree = pruned_tree)])

      # Get corresponding ancestral node in full tree:
      ancestral_nodes <- unlist(x = lapply(X = node_descendants, function(x) find_mrca(descendant_names = x, tree = full_tree)))

      # If pruned tree is identical to full tree (not pruned at all):
    } else {

      # Set ancestors as node numbers:
      ancestral_nodes <- node_numbers
    }

    # Output matrix matching node numbers of pruned tree to full tree:
    matrix(c(node_numbers, ancestral_nodes), ncol = 2, dimnames = list(c(), c("pruned_node", "full_node")))
  }

  # Get pruned node to full node for each unique tree:
  node_maps <- lapply(X = unique_trees, map_to_full_tree, full_tree = time_tree)

  # Build out for all trees (adds in any duplicated trees):
  node_maps <- node_maps[match(newick_strings, unique_newick_strings)]

  # Add node maps to data list:
  for (i in 1:length(x = data_list)) data_list[[i]]$node_maps <- node_maps[[i]]

  # Get all node names and numbers:
  nodes <- c(rownames(x = original_matrix), (n_tips + 1):(n_tips + n_nodes))

  # Renumber nodes of ancestral states:
  data_list <- lapply(X = data_list, function(x) {

    # Renumber nodes:
    names(x$ancestral_states)[match(as.character(x$node_maps[, "pruned_node"]), names(x$ancestral_states))] <- as.character(x$node_maps[, "full_node"])

    # Return renumbered nodes:
    x
  })

  # Collapse down to an ancestral state matrix ready for output:
  ancestral_state_matrix <- do.call(what = cbind, args = lapply(X = data_list, function(x) {

    # Get ancestral states for nodes:
    x$ancestral_states <- x$ancestral_states[nodes]

    # Add node names:
    names(x$ancestral_states) <- nodes

    # Return output:
    x$ancestral_states
  }))

  # Isolate estimated tip values:
  tip_matrix <- ancestral_state_matrix[rownames(x = original_matrix), ]

  # If there are any missing values:
  if (any(is.na(tip_matrix))) {

    # Isolate missing values:
    missing_tip_states <- which(x = is.na(tip_matrix))

    # Replace missing values with original (unmodified) input values:
    tip_matrix[missing_tip_states] <- original_matrix[missing_tip_states]

    # Add tip values back into full output:
    ancestral_state_matrix[rownames(x = original_matrix), ] <- tip_matrix
  }

  # Get column (character) count for each matrix block:
  matrix_columns <- unlist(x = lapply(X = lapply(X = raw_cladistic_matrix[2:length(x = raw_cladistic_matrix)], "[[", "matrix"), ncol))

  # For each matrix block:
  for (i in 1:length(x = matrix_columns)) {

    # Insert portion of ancestral state estimate into block:
    raw_cladistic_matrix[[(i + 1)]]$matrix <- ancestral_state_matrix[, 1:matrix_columns[i], drop = FALSE]

    # Remove that portion from the block:
    ancestral_state_matrix <- ancestral_state_matrix[, -(1:matrix_columns[i]), drop = FALSE]
  }

  # Overwrite ancestral state output with updated raw input:
  ancestral_state_matrix <- raw_cladistic_matrix

  # Add tree to output:
  ancestral_state_matrix$topper$tree <- time_tree

  # Return ancestral state matrix:
  ancestral_state_matrix
}
