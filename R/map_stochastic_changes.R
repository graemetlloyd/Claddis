#' Finds all state changes on a tree using stochastic character mapping
#'
#' @description
#'
#' Takes a cladistic matrix and time-scaled tree and makes point estimates for every character change using stochastic character mapping.
#'
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @param time_tree A time-scaled tree (phylo object) that represents the relationships of the taxa in \code{cladistic_matrix}.
#' @param time_bins A vector of ages representing the boundaries of a series of time bins.
#' @param n_simulations The number of simulations to perform (passed to \link{make.simmap}.
#' @param polymorphism_behaviour What to do with polymorphic (&) characters. One of "equalp", "missing", or "random". See details.
#' @param uncertainty_behaviour What to do with uncertain (/) characters. One of "equalp", "missing", or "random". See details.
#' @param inapplicable_behaviour What to do with inapplicable characters. Only one option currently ("missing"). See details.
#'
#' @details
#'
#' Important: this function is not yet complete and should not be used.
#'
#' A wrapper function for \link{make.simmap} in the \link{phytools} package.
#'
#' This function is intended to enumerate all possible changes on a tree (including to and from missing or inapplicable states) under the assumptions of stochastic character mapping as an alternative means of establishing branch-lengths (for rate analyses) or recording the state occupied at a particular point in time for disparity analyses.
#'
#' @return
#'
#' \item{all_state_changes}{A matrix of rows for each change with columns corresponding to the character, the simulation number, the edge number, the time the change occurred, and the start and end states.}
#' \item{character_times}{A vector of the sampled tree-length (in Ma) for each character.}
#' \item{binned_edge_lengths}{A matrix of time bins (columns) and characters (rows) indicating the sampled tree-length (in Ma).}
#' \item{binned_terminal_edge_lengths}{As above, but for terminal edges only.}
#' \item{binned_internal_edge_lengths}{As above, but for internal edges only.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Set random seed:
#' set.seed(2)
#'
#' # Use Day 2016 as source matrix:
#' cladistic_matrix <- day_2016
#'
#' # Prune out continuous characters:
#' cladistic_matrix <- prune_cladistic_matrix(
#'   cladistic_matrix =
#'     cladistic_matrix, blocks2prune = 1
#' )
#'
#' # Prune out majority of characters so
#' # example runs quickly:
#' cladistic_matrix <- prune_cladistic_matrix(
#'   cladistic_matrix =
#'     cladistic_matrix, characters2prune = 1:32
#' )
#'
#' # Generete random tree for matrix taxa:
#' time_tree <- ape::rtree(n = nrow(day_2016$matrix_1$matrix))
#'
#' # Add taxon names to tree:
#' time_tree$tip.label <- rownames(x = day_2016$matrix_1$matrix)
#'
#' # Add root age to tree:
#' time_tree$root.time <- max(diag(x = ape::vcv(phy = time_tree)))
#'
#' # Get all state changes for two simulations:
#' state_changes <-
#'   map_stochastic_changes(
#'     cladistic_matrix = cladistic_matrix,
#'     time_tree = time_tree, time_bins = seq(time_tree$root.time, 0,
#'       length.out = 3
#'     ), n_simulations = 2
#'   )
#'
#' # View matrix of all stochstic character changes:
#' state_changes$all_state_changes
#'
#' # View vector of sampled time for each character:
#' state_changes$character_times
#'
#' # View matrix of edge lengths in each time bin:
#' state_changes$binned_edge_lengths
#'
#' # View matrix of terminal edge lengths in each time bin:
#' state_changes$binned_terminal_edge_lengths
#'
#' # View matrix of internal edge lengths in each time bin:
#' state_changes$binned_internal_edge_lengths
#' @export map_stochastic_changes
map_stochastic_changes <- function(cladistic_matrix, time_tree, time_bins, n_simulations = 10, polymorphism_behaviour = "equalp", uncertainty_behaviour = "equalp", inapplicable_behaviour = "missing") {

  # IMPROVE CUSTOMISATION OF MAKE.SIMMAP WITH OPTIONS FOR PI, Q ETC. (FLAT PRIOR ON ROOT MAY BE PARTICULARLY BAD? ALLOW MAYBE SKEWING TOWARDS OUTGROUP STATE AS SOME KIND OF SLIDING VALUE?).
  # AND ONLY PERFORM SCM ON UNIQUE STATE DISTRIBUTON-CHARACTER TYPE COMBOS
  # MAJOR ISSUE IS NO EASY WAY TO MAKE MODEL FOR ORDERED MULTISTATE CHARACTER WHEN NOT ALL STATES ARE FOUND AT TIPS (E.G., 0 and 2 sampled, but not 1)
  # CHECK FOR ALL NA CHARACTERS AS THESE WILL NEED TO BE REMOVED.
  # MOVE MISSING AS POLYMORPHISM/UNCERTIANT BEHAVIOUR TO TOP AS MORE EFFICIENT.
  # ADD INAPPLICABLE OPTION THAT TIES SWITCH TO "" TO DEPENDENT CHARACTER A LA calculate_morphological_distances APPROACH.
  # MODEL AND TIP STATE DIMENSIONS MAY VARY IF SAY ONLY VARIANCE IS A POLYMORPHIC CHARACTER BUT "RANDOM" IS USED???

  # ANY REMAINING POLYMORPHISMS ARE FOR EQUAL P
  # IF USING EQUALP OR RANDOM AT END THEN NEED TO RECORD WEIRD CHANGE OF, SAY, 0 TO 0&1
  # IF USING MISSING NEED TO RECORD NA TO 0&1 CHANGE
  # GONNA HAVE TO DEAL WITH N SIMULATIONS SPREAD OVER MULTIPLE TIP STATES (SOME WILL BE UNIQUE, OTHERS WILL NEED DUPLICATION)
  # SMEARING BACK AND SMEARING FORWARD SOMETHIG TO NITHING AND NOTHING TO SOMETHING CHANGES MAY BE DIFFERENT. PERHAPS HAVE OPTION TO TREAT TIMESTAMP FOR THESE TO VARY.
  # IF ADDING TREES TO OUTPUT CONVERT THEM TO CLASS MULTIPHYLO, E.G. STOCHASTIC CHARACTER MAPS WITH NAS - NOTE THIS IN MANUAL TOO AS AN EXTENSION OF WHAT PHTTOOLS DOES.

  # Check for continuous and step matrices and stop and warn user if found:
  if (length(x = setdiff(x = unique(x = unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], function(x) x$ordering))), y = c("unord", "ord"))) > 0) stop("cladistic_matrix can only contain characters of type \"ord\" or \"unord\" (i.e., no step matrices or continuous characters).")

  # Check tree has branch lengths:
  if (is.null(time_tree$edge.length)) stop("time_tree does not have branch lengths (durations). Try timescaling the tree, e.g., with DatePhylo.")

  # Check branches all have positive length:
  if (any(time_tree$edge.length == 0)) stop("All branch lengths must be positive (no zero-length branches).")

  # Check time_tree has root age:
  if (is.null(time_tree$root.time)) stop("time_tree is missing $root.time. Try setting this before continuing, e.g., time_tree$root.time <- 104.2.")

  # Check polymorphism_behaviour is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = polymorphism_behaviour, y = c("equalp", "missing", "random"))) > 0) stop("polymorphism_behaviour must be one of must be one of \"equalp\", \"missing\" or \"random\".")

  # Check uncertainty_behaviour is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = uncertainty_behaviour, y = c("equalp", "missing", "random"))) > 0) stop("uncertainty_behaviour must be one of \"equalp\", \"missing\" or \"random\".")

  # Check inapplicable_behaviour is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = inapplicable_behaviour, y = c("missing"))) > 0) stop("inapplicable_behaviour must be \"missing\".")

  # Ensure time bins are in correct order:
  time_bins <- sort(x = unique(x = time_bins), decreasing = TRUE)

  # Get tree node ages:
  node_ages <- date_nodes(time_tree = time_tree)

  # Build all data into single matrix:
  matrix_block <- do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"))

  # If inapplicable_behaviour is missing replace inaplicables with NAs:
  if (inapplicable_behaviour == "missing" && any(matrix_block[!is.na(matrix_block)] == "")) matrix_block[which(x = matrix_block == "")] <- NA

  # Assemble all ordering values into a single vector:
  ordering <- unname(do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering")))

  # Assemble all minimum values into a single vector:
  minimum_values <- unname(do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values")))

  # Assemble all maximum values into a single vector:
  maximum_values <- unname(do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values")))

  # Assemble all maximum values into a single vector:
  character_weights <- unname(do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "character_weights")))

  # Build each character into list values starting with tip state lists (of N Simulations in length):
  character_list <- lapply(X = lapply(X = apply(matrix_block, 2, list), unlist), function(x) {
    y <- list()
    y$tip_states <- lapply(X = apply(matrix(rep(x, times = n_simulations), ncol = n_simulations, dimnames = list(names(x), c())), 2, list), unlist)
    y
  })

  # Add weight to list for each character:
  for (i in 1:length(x = character_list)) character_list[[i]]$weight <- character_weights[i]

  # Add ordering to list for each character:
  for (i in 1:length(x = character_list)) character_list[[i]]$ordering <- ordering[i]

  # Add full tree to list for each character:
  for (i in 1:length(x = character_list)) character_list[[i]]$full_tree <- time_tree

  # Subfunction to perform polymorphism and uncertainty changes:
  edit_polymorphisms <- function(x, polymorphism_behaviour, uncertainty_behaviour) {

    # Isolate tip states list:
    tip_state_list <- x$tip_states

    # If polymorphism_behaviour is missing replace all polymorphisms with NAs:
    if (polymorphism_behaviour == "missing") {
      tip_state_list <- lapply(X = tip_state_list, function(x) {
        cells_to_alter <- grep("&", x)
        if (length(x = cells_to_alter) > 0) x[cells_to_alter] <- NA
        return(x)
      })
    }

    # If uncertainty_behaviour is missing replace all uncertainties with NAs:
    if (uncertainty_behaviour == "missing") {
      tip_state_list <- lapply(X = tip_state_list, function(x) {
        cells_to_alter <- grep("/", x)
        if (length(x = cells_to_alter) > 0) x[cells_to_alter] <- NA
        return(x)
      })
    }

    # If polymorphism_behaviour is random replace all polymorphisms with one included state at random:
    if (polymorphism_behaviour == "random") {
      tip_state_list <- lapply(X = tip_state_list, function(x) {
        cells_to_alter <- grep("&", x)
        if (length(x = cells_to_alter) > 0) for (i in 1:length(x = cells_to_alter)) x[cells_to_alter[i]] <- sample(strsplit(x[cells_to_alter[i]], split = "&")[[1]], size = 1)
        return(x)
      })
    }

    # If uncertainty_behaviour is random replace all uncertainties with one included state at random:
    if (uncertainty_behaviour == "random") {
      tip_state_list <- lapply(X = tip_state_list, function(x) {
        cells_to_alter <- grep("/", x)
        if (length(x = cells_to_alter) > 0) for (i in 1:length(x = cells_to_alter)) x[cells_to_alter[i]] <- sample(strsplit(x[cells_to_alter[i]], split = "/")[[1]], size = 1)
        return(x)
      })
    }

    # If uncertainty_behaviour is equalp replace all uncertainties withpolymorphism equivalent (makes grep simpler later):
    if (uncertainty_behaviour == "equalp") tip_state_list <- lapply(X = tip_state_list, function(x) gsub(pattern = "/", replacement = "&", x = x))

    # Reinsert modified tip states:
    x$tip_states <- tip_state_list

    # Return list with modified tip states added:
    x
  }

  # Alter polymorphisms and uncertainties according to polymorphism_behaviour and uncertainty_behaviour settings:
  character_list <- lapply(X = character_list, edit_polymorphisms, polymorphism_behaviour = polymorphism_behaviour, uncertainty_behaviour = uncertainty_behaviour)

  # Create pruned tipstates by removing all missing values:
  character_list <- lapply(X = character_list, function(x) {
    x$pruned_tip_states <- lapply(X = x$tip_states, function(y) y[!is.na(y)])
    x
  })

  # Get minimum and maximum values for each character (has to be done post polymorphism and inapplicable steps or will not work correctly):
  character_list <- lapply(X = character_list, function(x) {
    ranges <- range(as.numeric(unlist(x = strsplit(unlist(x = x$pruned_tip_states), split = "&"))))
    x$minimum_values <- ranges[1]
    x$maximum_values <- ranges[2]
    x
  })

  # Subfunction to form tip states matrices ready for simmap function:
  build_tip_state_matrices <- function(x) {

    # Isolate pruned tip states list:
    pruned_tip_state_list <- x$pruned_tip_states

    # Subfunction to build tip state matrix:
    build_tip_state_matrix <- function(y, minimum_values, maximum_values) {

      # Create empty tip state matrix:
      tip_state_matrix <- matrix(0, nrow = length(x = y), ncol = length(x = minimum_values:maximum_values), dimnames = list(names(y), as.character(minimum_values:maximum_values)))

      # Get row and column data as list:
      rows_and_columns <- lapply(X = lapply(X = as.list(x = y), strsplit, split = "&"), unlist)

      # Isolate rows (taxon names):
      rows <- unname(unlist(x = mapply(rep, names(rows_and_columns), unlist(x = lapply(X = rows_and_columns, length)))))

      # Isolate columns:
      columns <- unname(unlist(x = rows_and_columns))

      # Fill tip state matrix:
      for (i in 1:length(x = rows)) tip_state_matrix[rows[i], columns[i]] <- 1

      # Convert each row to a probability (unnecessary if only one column):
      if (ncol(tip_state_matrix) > 1) tip_state_matrix <- t(apply(tip_state_matrix, 1, function(x) x / sum(x)))

      # Return tip state matrix:
      return(tip_state_matrix)
    }

    # Build tip state matrices:
    pruned_tip_state_list <- lapply(X = pruned_tip_state_list, build_tip_state_matrix, minimum_values = x$minimum_values, maximum_values = x$maximum_values)

    # Overwrite pruned tip states with tip state matrices:
    x$pruned_tip_states <- pruned_tip_state_list

    # Return full output:
    x
  }

  # Build tip state matrices:
  character_list <- lapply(X = character_list, build_tip_state_matrices)

  # Subfunction to build character model:
  build_character_model <- function(x) {

    # If character is ordered:
    if (x$ordering == "ord") {

      # Create character model with all transitions impossible:
      character_model <- matrix(0, ncol = length(x = x$minimum_values:x$maximum_values), nrow = length(x = x$minimum_values:x$maximum_values), dimnames = list(x$minimum_values:x$maximum_values, x$minimum_values:x$maximum_values))

      # Create character matrix with only off-diagonal (adjacent) states possible:
      if (ncol(character_model) > 1) for (i in 2:ncol(character_model)) character_model[i, (i - 1)] <- character_model[(i - 1), i] <- 1
    }

    # If character is unordered:
    if (x$ordering == "unord") {

      # Create character model with all transitions possible:
      character_model <- matrix(1, ncol = length(x = x$minimum_values:x$maximum_values), nrow = length(x = x$minimum_values:x$maximum_values), dimnames = list(x$minimum_values:x$maximum_values, x$minimum_values:x$maximum_values))

      # Set diagonals to zero:
      diag(x = character_model) <- 0
    }

    # Store character model:
    x$character_model <- character_model

    # Return full output:
    x
  }

  # Build character models for each character:
  character_list <- lapply(X = character_list, build_character_model)

  # Get total tips in tree:
  n_tips <- ape::Ntip(phy = time_tree)

  # Subfunction to prune tree to just tips with data:
  prune_tree <- function(x, n_tips) {

    # Find any tips to drop:
    tips_to_drop <- setdiff(x = x$full_tree$tip.label, y = rownames(x = x$pruned_tip_states[[1]]))

    # If there are tips to drop:
    if (length(x = tips_to_drop) > 0) {

      # If exactly one tip will remain:
      if ((n_tips - length(x = tips_to_drop)) == 1) {

        # Set pruned tree initially as full tree:
        x$pruned_tree <- x$full_tree

        # Set tip number (of single coded taxon):
        tip_number <- which(x = time_tree$tip.label == setdiff(x = time_tree$tip.label, y = tips_to_drop))

        # Find single tip edge:
        tip_edge <- which(x = time_tree$edge[, 2] == tip_number)

        # Set all other edge lengths to zero:
        x$pruned_tree$edge.length[setdiff(x = 1:length(x = time_tree$edge), y = tip_edge)] <- 0

        # Reset root time to beginning of tip edge:
        x$pruned_tree$root.time <- unname(node_ages[time_tree$edge[tip_edge, 1]])
      }

      # If exactly two tips will remain:
      if ((n_tips - length(x = tips_to_drop)) == 2) {

        # Prune tree and store:
        x$pruned_tree <- ape::drop.tip(phy = x$full_tree, tip = tips_to_drop)

        # Correct root time manually:
        x$pruned_tree$root.time <- unname(node_ages[find_mrca(descendant_names = setdiff(x = time_tree$tip.label, y = tips_to_drop), tree = time_tree)])
      }

      # If at least three tips will remain:
      if ((n_tips - length(x = tips_to_drop)) > 2) {

        # Prune tree and store:
        x$pruned_tree <- ape::drop.tip(phy = x$full_tree, tip = tips_to_drop)

        # Ensure pruned trees $root.time value is correct:
        x$pruned_tree <- fix_root_time(time_tree, x$pruned_tree)
      }

      # If no tips need to be dropped:
    } else {

      # Store full tree as pruned version:
      x$pruned_tree <- x$full_tree
    }

    # Return full output:
    return(x)
  }

  # Get pruned trees with only tips from pruned tip states returned:
  character_list <- lapply(X = character_list, FUN = prune_tree, n_tips = n_tips)

  # Subfunction to get node ages for pruned tree:
  date_pruned_nodes <- function(x) {

    # Get number of (scorable) tips for pruned tree:
    n_tips <- lapply(X = x$pruned_tip_states, nrow)[[1]]

    # If exactly one tip:
    if (n_tips == 1) {

      # Set pruned node ages as beginning and end of single branch:
      x$pruned_date_nodes <- c(x$pruned_tree$root.time - sum(x$pruned_tree$edge.length), x$pruned_tree$root.time)

      # Set node numbers as 1 for the tip and 2 for the node subtending the sole branch:
      names(x$pruned_date_nodes) <- as.character(1:2)
    }

    # If exactly two tips:
    if (n_tips == 2) {

      # Set node ages by subtracting two branch lengths from root time and adding root time at end:
      x$pruned_date_nodes <- c(x$pruned_tree$root.time - x$pruned_tree$edge.length[order(x$pruned_tree$edge[, 2])], x$pruned_tree$root.time)

      # Set node ages as 1 and 2 (for tips) and 3 for root:
      names(x$pruned_date_nodes) <- as.character(1:3)
    }

    # If more than two tips just apply get node ages function:
    if (n_tips > 2) x$pruned_date_nodes <- date_nodes(time_tree = x$pruned_tree)

    # Return full output:
    return(x)
  }

  # Get node ages for pruned tree:
  character_list <- lapply(X = character_list, date_pruned_nodes)

  # Subfunction to perform edge matches between pruned and full tree:
  match_edges <- function(x, tree) {

    # Get number of (scorable) tips for pruned tree:
    n_tips <- lapply(X = x$pruned_tip_states, nrow)[[1]]

    # If exactly one tip:
    if (n_tips == 1) {

      # Compile single branch output (i.e., single branch of tree with positive length):
      x$edge_matches <- list(which(x = x$pruned_tree$edge.length > 0))

      # Add name (has to be one because only one edge):
      names(x$edge_matches) <- "1"
    }

    # If exactly two tips:
    if (n_tips == 2) {

      # Get two tip names (in order):
      tip_names <- x$pruned_tree$tip.label

      # Get shared ancestor node on full tree:
      ancestor_node <- find_mrca(descendant_names = tip_names, tree = tree)

      # Get two tip numbers for two tips on full tree:
      tip_numbers <- unlist(x = lapply(X = lapply(X = as.list(x = tip_names), "==", tree$tip.label), which))

      # Create empty list ready to store edge matches:
      x$edge_matches <- list()

      # For each tip number:
      for (i in tip_numbers) {

        # Get first edge found (terminal branch):
        edges_found <- which(x = tree$edge[, 2] == i)

        # Set current node to ith node:
        current_node <- i

        # While the ancestor has not been hit:
        while (tree$edge[edges_found[1], 1] != ancestor_node) {

          # Reset current node to start of current branch:
          current_node <- tree$edge[edges_found[1], 1]

          # Add new edges found to vector:
          edges_found <- c(which(x = tree$edge[, 2] == current_node), edges_found)
        }

        # Store edges found in root-to-tip order:
        x$edge_matches[[which(x = tip_numbers == i)]] <- edges_found
      }

      # Add edge names from pruned tree:
      names(x$edge_matches) <- as.character(1:2)
    }

    # If more than two tips simply use function normally:
    if (n_tips > 2) x$edge_matches <- match_tree_edges(tree, x$pruned_tree)$matching_edges

    # Return full output:
    return(x)
  }

  # Get edge matches between pruned and full trees (for recording true edge changes later):
  character_list <- lapply(X = character_list, match_edges, tree = time_tree)

  # Subfunction to get edge lengths in bins (usable later for rate calculations):
  get_binned_edge_lengths <- function(x) {

    # Set temporary tree as full tree (as modifying branch lengths but will rant to retain these later:
    temporary_tree <- x$full_tree

    # Find any branch lengths on full tree to set as zero (effectively excluding them from the edge lengths as they correspond to missing values):
    branches_to_set_as_zeroes <- setdiff(x = 1:nrow(temporary_tree$edge), y = unname(unlist(x = x$edge_matches)))

    # If there are nranches to set lengtsh to zero then do so and store:
    if (length(x = branches_to_set_as_zeroes) > 0) temporary_tree$edge.length[branches_to_set_as_zeroes] <- 0

    # Get edge lengths in bins:
    x$edge_lengths_in_bins <- bin_edge_lengths(time_tree = temporary_tree, time_bins = time_bins)

    # Return full output:
    return(x)
  }

  # Get edge lengths in bins:
  character_list <- lapply(X = character_list, get_binned_edge_lengths)

  # Subfunction to perform actual stochastic character maps:
  build_character_map_trees <- function(x) {

    # Get number of (scorable) tips for pruned tree:
    n_tips <- lapply(X = x$pruned_tip_states, nrow)[[1]]

    # If exactly one tip:
    if (n_tips == 1) {

      # Subfunction to generate stochastic character map like output but for the single taxon case:
      build_single_taxon_map <- function(y, tree) {

        # Set output as the tree initially:
        output <- tree

        # Add maps to output of just the single branch's edge length:
        output$maps <- list(sum(tree$edge.length))

        # Add state name to maps:
        output$maps <- lapply(X = output$maps, function(z) {
          names(z) <- colnames(x = y)
          z
        })

        # Return output:
        output
      }

      # Perform stochastic character mapping and store output:
      x$stochastic_character_map_trees <- lapply(X = x$pruned_tip_states, build_single_taxon_map, tree = x$pruned_tree)
    }

    # If exactly two tips:
    if (n_tips == 2) {

      # Subfunction to generate stochastic character map like output but for the two taxon case:
      build_two_taxon_map <- function(y, tree) {

        # Set output as the tree initially:
        output <- tree

        # Add maps to output of just the single branch's edge length:
        Output$maps <- as.list(x = unname(tree$edge.length))

        # If only one state (character is constant) add state name to maps:
        if (ncol(y) == 1) {
          output$maps <- lapply(X = output$maps, function(z) {
            names(z) <- colnames(x = y)
            z
          })
        }

        # If character is variant:
        if (ncol(y) == 2) {

          # Get raw data for calculating root state probabiity (incorporates tip sattes and reciprocal of branch lengths):
          root_state_raw_data <- apply(sweep(y[tree$tip.label, ], MARGIN = 1, matrix(1 / (tree$edge.length / sum(tree$edge.length))), "*"), 2, sum)

          # Divide through by sum to get true probability:
          root_state_probability <- root_state_raw_data / sum(root_state_raw_data)

          # Sample root state with probabilities based on tip state(s) and reciprocal of branch lengths:
          root_state <- sample(names(root_state_probability), size = 1, prob = root_state_probability)

          # Initially set map names to root state (one may change later):
          output$maps <- lapply(X = output$maps, function(z) {
            names(z) <- root_state
            z
          })

          # Get tip states with root column pruned (helps identify branch with changes):
          pruned_root_tip_states <- y[tree$tip.label, -which(x = colnames(x = y) == root_state), drop = FALSE]

          # Get tips with terminal changes (could concievably be empty if there is a polymorphism):
          tips_with_terminal_changes <- names(which(x = apply(pruned_root_tip_states, 1, "==", 1)))

          # If there is a tip with a change to record:
          if (length(x = tips_with_terminal_changes) > 0) {

            # Get the tip number:
            tip_number <- which(x = tree$tip.label == tips_with_terminal_changes)

            # Get the edge the tip corresponds to:
            edge_number <- which(x = tree$edge[, 2] == tip_number)

            # Get tip state:
            tip_state <- colnames(x = pruned_root_tip_states[, which(x = pruned_root_tip_states[tips_with_terminal_changes, ] == 1), drop = FALSE])

            # Get edge length:
            edge_length <- unname(output$maps[[edge_number]])

            # Pick a split point (proportion along branch where state change occurs):
            split_point <- stats::runif(n = 1, min = 0, max = 1)

            # Cnvert proportions into absolute branch length values:
            output$maps[[edge_number]] <- c(split_point, 1 - split_point) * edge_length

            # Add state names to output:
            names(output$maps[[edge_number]]) <- c(root_state, tip_state)
          }
        }

        # Return output:
        output
      }

      # Perform stochastic character mapping and store output:
      x$stochastic_character_map_trees <- lapply(X = x$pruned_tip_states, build_two_taxon_map, tree = x$pruned_tree)
    }

    # If more than two tips:
    if (n_tips > 2) {

      # Subfunction to perform stochastic character mapping:
      build_map_tree <- function(y, tree, model) {

        # If character is constant (invariant):
        if (ncol(y) == 1) {

          # Set initial output as tree:
          output <- tree

          # Set maps as edge lengths:
          output$maps <- as.list(x = tree$edge.length)

          # Add output names to maps as invariant character:
          output$maps <- lapply(X = output$maps, function(z) {
            names(z) <- colnames(x = y)
            z
          })
        }

        # If character is variable, perform regular stochastic character mapping:
        if (ncol(y) > 1) output <- phytools::make.simmap(tree = tree, x = y, nsim = 1, model = model, pi = "estimated", message = FALSE)

        # Return output:
        output
      }

      # Perform stochastic character mapping and store output:
      x$stochastic_character_map_trees <- lapply(X = x$pruned_tip_states, build_map_tree, tree = x$pruned_tree, model = x$character_model)
    }

    # Return all output:
    return(x)
  }

  # Generate initial stochastic character map trees:
  character_list <- lapply(X = character_list, build_character_map_trees)

  # Subfunction to map stochastic chracter maps of pruned tree to full trees:
  map_characters_to_full_tree <- function(x) {

    # Only proceed if pruned tree is actually smaller than full tree (otherwise data are fine as is):
    if (lapply(X = x$pruned_tip_states, nrow)[[1]] < ape::Ntip(phy = x$full_tree)) {

      # Create empty list to store stochastic character maps for full trees:
      y <- list()

      # Start by filling out list with full tree:
      for (i in 1:n_simulations) y[[i]] <- x$full_tree

      # Generate null stochatsic character map from edge lengths:
      null_map <- as.list(x = x$full_tree$edge.length)

      # Add NA as default state (i.e., missing data) - will want to correct this for inapplicables later:
      null_map <- lapply(X = null_map, function(z) {
        names(z) <- NA
        z
      })

      # Add null stochastic character map to list:
      for (i in 1:n_simulations) y[[i]]$maps <- null_map

      # Find any single edge matches (where pruned edge matches a single edge in the full tree):
      single_edge_matches <- unname(which(x = unlist(x = lapply(X = x$edge_matches, length)) == 1))

      # Find any multiple edge matches (where a pruned edge matches more than one edge in the full tree):
      multiple_edge_matches <- unname(which(x = unlist(x = lapply(X = x$edge_matches, length)) > 1))

      # If there is at least one single edge match then map pruned edges to full tree for them:
      if (length(x = single_edge_matches) > 0) for (i in 1:length(x = y)) y[[i]]$maps[unname(unlist(x = x$edge_matches[single_edge_matches]))] <- x$stochastic_character_map_trees[[i]]$maps[single_edge_matches]

      # If there is at least one multiple edge match:
      if (length(x = multiple_edge_matches) > 0) {

        # For each simulation:
        for (i in 1:length(x = y)) {

          # Find multiple edge matches where pruned edge has no changes (single state persists):
          multiple_edge_single_state_matches <- multiple_edge_matches[which(x = lapply(X = x$stochastic_character_map_trees[[i]]$maps[multiple_edge_matches], length) == 1)]

          # Find multiple edge matches where pruned edge has at least one change (multiple states sampled):
          multiple_edge_multiple_state_matches <- multiple_edge_matches[which(x = lapply(X = x$stochastic_character_map_trees[[i]]$maps[multiple_edge_matches], length) > 1)]

          # If there are multiple edge but single state matches:
          if (length(x = multiple_edge_single_state_matches) > 0) {

            # For each match simply replace NA with the single character state:
            for (j in multiple_edge_single_state_matches) {
              y[[i]]$maps[unname(unlist(x = x$edge_matches[j]))] <- lapply(X = y[[i]]$maps[unname(unlist(x = x$edge_matches[j]))], function(z) {
                names(z) <- names(x$stochastic_character_map_trees[[i]]$maps[j][[1]])
                z
              })
            }
          }

          # If there are multiple edge multiple state matches:
          if (length(x = multiple_edge_multiple_state_matches) > 0) {

            # For each multiple state multiple edge match:
            for (j in multiple_edge_multiple_state_matches) {

              # Get matching edges on full tree:
              matching_edges <- unname(unlist(x = x$edge_matches[j]))

              # Get edge lengths of matching full tree edges:
              matching_edge_lengths <- unname(unlist(x = y[[i]]$maps[unname(unlist(x = x$edge_matches[j]))]))

              # Get pruned stochastic maps for current pruned edge:
              pruned_stochastic_maps <- x$stochastic_character_map_trees[[i]]$maps[j][[1]]

              # Get starting age of current pruned edge (will use to get change times that can be binned by full tree edge later):
              start_age_of_pruned_edge <- unname(x$pruned_date_nodes[x$pruned_tree$edge[j, 1]])

              # Get times at which character changes occur:
              change_times <- unname(start_age_of_pruned_edge - cumsum(pruned_stochastic_maps[1:(length(x = pruned_stochastic_maps) - 1)]))

              # Build matrix of from-to changes:
              from_to_changes <- cbind(names(pruned_stochastic_maps[1:(length(x = pruned_stochastic_maps) - 1)]), names(pruned_stochastic_maps[2:length(x = pruned_stochastic_maps)]))

              # Set matching edges on full tree as time bins:
              matching_edges_as_time_bins <- c(start_age_of_pruned_edge, start_age_of_pruned_edge - cumsum(matching_edge_lengths))

              # Get edge on whicb change occurs:
              change_edges <- unlist(x = lapply(X = as.list(x = change_times), function(z) min(which(x = z > matching_edges_as_time_bins)) - 1))

              # Create full tree stochastic character map of correct size:
              full_stochastic_maps <- lapply(X = as.list(x = rle(sort(x = c(change_edges, 1:length(x = matching_edges))))$lengths), function(z) rep(0, z))

              # Set current time as beginning of pruned edge:
              current_time <- start_age_of_pruned_edge

              # Set current edge as 1 (will increment through while loop below):
              current_edge <- 1

              # Set current state as first from value in from-to matrix:
              current_state <- from_to_changes[1, 1]

              # Make vector of edge switch times:
              edge_switch_times <- matching_edges_as_time_bins[2:length(x = matching_edges_as_time_bins)]

              # Whilst there are still changes or edge switches left to deal with:
              while (length(x = c(edge_switch_times, change_times)) > 0) {

                # Set next event time:
                next_event <- max(c(edge_switch_times, change_times))

                # If next event is to switch edges:
                if (edge_switch_times[1] == next_event) {

                  # Find current position on stochastic map:
                  current_position <- which(x = full_stochastic_maps[[current_edge]] == 0)[1]

                  # Store edge length:
                  full_stochastic_maps[[current_edge]][current_position] <- current_time - next_event

                  # Store current state:
                  names(full_stochastic_maps[[current_edge]])[current_position] <- current_state

                  # Update current time:
                  current_time <- next_event

                  # Update current edge:
                  current_edge <- current_edge + 1

                  # Update edge_switch_times by removing last change:
                  edge_switch_times <- edge_switch_times[-1]

                  # If next event is a change of character state:
                } else {

                  # Find current position on stochastic map:
                  current_position <- which(x = full_stochastic_maps[[current_edge]] == 0)[1]

                  # Store edge length:
                  full_stochastic_maps[[current_edge]][current_position] <- current_time - next_event

                  # Store current state:
                  names(full_stochastic_maps[[current_edge]])[current_position] <- current_state

                  # As long as the from-to matrix still exists:
                  if (nrow(from_to_changes) > 0) {

                    # Update current state:
                    current_state <- from_to_changes[1, 2]

                    # Prune change from from-to matrix:
                    from_to_changes <- from_to_changes[-1, , drop = FALSE]
                  }

                  # Update current time:
                  current_time <- next_event

                  # Prune change time from vector:
                  change_times <- change_times[-1]
                }
              }

              # Store full stochstic character map in y:
              y[[i]]$maps[unname(unlist(x = x$edge_matches[j]))] <- full_stochastic_maps
            }
          }
        }
      }

      # Overwrite pruned tree stochastic character maps with full tree stochastic character maps:
      x$stochastic_character_map_trees <- y
    }

    # Return full output:
    return(x)
  }

  # Map stochastic characters to full tree in preparation for recording changes:
  character_list <- lapply(X = character_list, map_characters_to_full_tree)

  # Subfunction to extract character changes matrix from trees:
  reformat_changes <- function(x) {

    # Find the root edge for the tree:
    root_edge <- match(ape::Ntip(phy = x$full_tree) + 1, x$full_tree$edge[, 1])

    # Subfunction to get character changes and root state:
    extract_changes <- function(y) {

      # Get root state:
      root_state <- as.numeric(names(y$maps[[root_edge]][1]))

      # Find any edges with changes on them:
      edges_with_changes <- which(x = unlist(x = lapply(X = y$maps, length)) > 1)

      # As long as there is at least one change:
      if (length(x = edges_with_changes) > 0) {

        # Get ages at start of edges with changes (subtracting from this will give change times later:
        age_at_start_of_edges_with_changes <- node_ages[x$full_tree$edge[edges_with_changes, 1]]

        # Get from and to states for each change:
        froms_and_tos <- matrix(as.numeric(unlist(x = strsplit(unlist(x = lapply(X = y$maps[edges_with_changes], function(z) paste(names(z[1:(length(x = z) - 1)]), names(z[2:length(x = z)]), sep = "%%"))), split = "%%"))), ncol = 2, byrow = TRUE)

        # Get character change times:
        character_change_times <- unname(unlist(x = mapply("-", age_at_start_of_edges_with_changes, lapply(X = y$maps[edges_with_changes], function(z) cumsum(z[1:(length(x = z) - 1)])))))

        # Build changes matrix:
        changes_matrix <- cbind(froms_and_tos, unlist(x = mapply(rep, edges_with_changes, unlist(x = lapply(X = y$maps[edges_with_changes], length)) - 1)), character_change_times)

        # Add column names to matrix:
        colnames(x = changes_matrix) <- c("from", "to", "edge", "time")

        # If no changes occur:
      } else {

        # Create changes matrix with no changes:
        changes_matrix <- matrix(nrow = 0, ncol = 4, dimnames = list(c(), c("from", "to", "edge", "time")))
      }

      # Build matrix of from to states for each edge (i.e., state at start of edge and state at end of edge):
      edge_from_to <- cbind(as.numeric(unlist(x = lapply(X = y$maps, function(y) names(y[1])))), as.numeric(unlist(x = lapply(X = y$maps, function(y) names(y[length(x = y)])))))

      # Get preceding states for each edge (will help identify edge-to-edge changes:
      preceding_state <- edge_from_to[match(time_tree$edge[, 1], time_tree$edge[, 2]), 2]

      # Get following states for each edge (will help identify edge-to-edge changes:
      following_state <- edge_from_to[match(time_tree$edge[, 2], time_tree$edge[, 1]), 1]

      # Correct following edge to eliminate terminal changes whic do not need to be recorded:
      following_state[match(1:n_tips, time_tree$edge[, 2])] <- edge_from_to[match(1:n_tips, time_tree$edge[, 2]), 2]

      # Get preceding changes (will want to add to changes matrix):
      preceding_changes <- which(x = apply(cbind(is.na(preceding_state), !is.na(edge_from_to[, 1])), 1, all))

      # Get following changes (will want to add to changes matrix):
      following_changes <- which(x = apply(cbind(!is.na(edge_from_to[, 2]), is.na(following_state)), 1, all))

      # For each unique preceding change:
      for (i in preceding_changes[!duplicated(time_tree$edge[preceding_changes, 1])]) {

        # If simply the root state:
        if (n_tips + 1 == time_tree$edge[i, 1]) {

          # Add root change to matrix as edge zero:
          changes_matrix <- rbind(changes_matrix, c(NA, edge_from_to[i, 1], 0, time_tree$root.time))

          # If some other state:
        } else {

          # Add NA to something state change to changes matrix:
          changes_matrix <- rbind(changes_matrix, c(NA, edge_from_to[i, 1], which(x = time_tree$edge[, 2] == time_tree$edge[i, 1]), unname(node_ages[time_tree$edge[i, 1]])))
        }
      }

      # If there are any following changes (i.e., is there at least one missing value at the tips):
      if (length(x = following_changes) > 0) {

        # Add following changes to changes matrix:
        changes_matrix <- rbind(changes_matrix, cbind(edge_from_to[following_changes, 2], following_state[following_changes], unlist(x = lapply(X = as.list(x = time_tree$edge[following_changes, 2]), function(z) {
          following_edges <- which(x = time_tree$edge[, 1] == z)
          following_edges[is.na(edge_from_to[following_edges, 1])]
        })), unname(node_ages[time_tree$edge[following_changes, 2]])))
      }

      # Return output:
      changes_matrix
    }

    # Get root and changes for each stochastic character map and store:
    x$stochastic_character_changes <- lapply(X = x$stochastic_character_map_trees, extract_changes)

    # Return output:
    return(x)
  }

  # Get root state and character changes for each stochastic character map:
  character_list <- lapply(X = character_list, reformat_changes)

  # Subfunction to collapse changes across simulations into a single matrix:
  compile_changes_as_matrix <- function(x) {

    # Get simulation numbers (will be new column in matrix):
    simulation_numbers <- unlist(x = mapply(rep, 1:n_simulations, unlist(x = lapply(X = x$stochastic_character_changes, function(y) nrow(y)))))

    # Collapse changes into single matrix:
    changes_matrix <- do.call(what = rbind, args = lapply(X = x$stochastic_character_changes, function(y) y))

    # Add simulation numbers:
    changes_matrix <- cbind(matrix(simulation_numbers), changes_matrix)

    # Add column name:
    colnames(x = changes_matrix)[1] <- "simulation_number"

    # Overwrite stochastic character matrices with new collapsed format:
    x$stochastic_character_changes <- changes_matrix

    # Return full output:
    return(x)
  }

  # Collapse stochastic character matrices to single matrix for each character:
  character_list <- lapply(X = character_list, compile_changes_as_matrix)

  # Compile all state changes into a single matrix:
  all_state_changes <- cbind(matrix(unlist(x = mapply(rep, 1:ncol(matrix_block), lapply(X = character_list, function(x) nrow(x$stochastic_character_changes))))), do.call(what = rbind, args = lapply(X = character_list, function(x) x$stochastic_character_changes)))

  # Add column name for character:
  colnames(x = all_state_changes)[1] <- "character"

  # Get rid of pesky state rownames:
  rownames(x = all_state_changes) <- NULL

  # Sort all state changes by edge number:
  all_state_changes <- all_state_changes[order(all_state_changes[, "edge"]), ]

  # Sort all state changes by character number:
  all_state_changes <- all_state_changes[order(all_state_changes[, "character"]), ]

  # Sort all state changes by simulation number:
  all_state_changes <- all_state_changes[order(all_state_changes[, "simulation_number"]), ]

  # Get character times (length of subtrees for each character):
  character_times <- unlist(x = lapply(X = character_list, function(x) sum(x$edge_lengths_in_bins$binned_edge_lengths)))

  # Get edge length for each bin by character:
  binned_edge_lengths <- do.call(what = rbind, args = lapply(X = character_list, function(x) x$edge_lengths_in_bins$binned_edge_lengths))

  # Get edge length for each bin by character (terminal branches only):
  binned_terminal_edge_lengths <- do.call(what = rbind, args = lapply(X = character_list, function(x) x$edge_lengths_in_bins$binned_terminal_edge_lengths))

  # Get edge length for each bin by character (internal branches only):
  binned_internal_edge_lengths <- do.call(what = rbind, args = lapply(X = character_list, function(x) x$edge_lengths_in_bins$binned_internal_edge_lengths))

  # Compile output as list:
  output <- list(all_state_changes = all_state_changes, character_times = character_times, binned_edge_lengths = binned_edge_lengths, binned_terminal_edge_lengths = binned_terminal_edge_lengths, binned_internal_edge_lengths = binned_internal_edge_lengths)

  # Return output:
  return(invisible(output))
}
