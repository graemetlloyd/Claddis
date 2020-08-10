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
#' \item{root_states}{A matrix of the root states for each character (column) and simulation (rows).}
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
#' StateChanges <-
#'   map_stochastic_changes(
#'     cladistic_matrix = cladistic_matrix,
#'     time_tree = time_tree, time_bins = seq(time_tree$root.time, 0,
#'       length.out = 3
#'     ), n_simulations = 2
#'   )
#'
#' # View matrix of all stochstic character changes:
#' StateChanges$all_state_changes
#'
#' # View vector of sampled time for each
#' # character:
#' StateChanges$character_times
#'
#' # View matrix of edge lengths in each time bin:
#' StateChanges$binned_edge_lengths
#'
#' # View matrix of termnial edge lengths in each time bin:
#' StateChanges$binned_terminal_edge_lengths
#'
#' # View matrix of internal edge lengths in each time bin:
#' StateChanges$binned_internal_edge_lengths
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
  MatrixBlock <- do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"))

  # If inapplicable_behaviour is missing replace inaplicables with NAs:
  if (inapplicable_behaviour == "missing" && any(MatrixBlock[!is.na(MatrixBlock)] == "")) MatrixBlock[which(x = MatrixBlock == "")] <- NA

  # Assemble all ordering values into a single vector:
  ordering <- unname(do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering")))

  # Assemble all minimum values into a single vector:
  minimum_values <- unname(do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values")))

  # Assemble all maximum values into a single vector:
  maximum_values <- unname(do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values")))

  # Assemble all maximum values into a single vector:
  character_weights <- unname(do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "character_weights")))

  # Build each character into list values starting with tip state lists (of N Simulations in length):
  CharacterList <- lapply(X = lapply(X = apply(MatrixBlock, 2, list), unlist), function(x) {
    y <- list()
    y$TipStates <- lapply(X = apply(matrix(rep(x, times = n_simulations), ncol = n_simulations, dimnames = list(names(x), c())), 2, list), unlist)
    return(y)
  })

  # Add weight to list for each character:
  for (i in 1:length(x = CharacterList)) CharacterList[[i]]$Weight <- character_weights[i]

  # Add ordering to list for each character:
  for (i in 1:length(x = CharacterList)) CharacterList[[i]]$ordering <- ordering[i]

  # Add full tree to list for each character:
  for (i in 1:length(x = CharacterList)) CharacterList[[i]]$full_tree <- time_tree

  # Subfunction to perform polymorphism and uncertainty changes:
  edit_polymorphisms <- function(x, polymorphism_behaviour, uncertainty_behaviour) {

    # Isolate tip states list:
    TipStateList <- x$TipStates

    # If polymorphism_behaviour is missing replace all polymorphisms with NAs:
    if (polymorphism_behaviour == "missing") {
      TipStateList <- lapply(X = TipStateList, function(x) {
        CellsToAlter <- grep("&", x)
        if (length(x = CellsToAlter) > 0) x[CellsToAlter] <- NA
        return(x)
      })
    }

    # If uncertainty_behaviour is missing replace all uncertainties with NAs:
    if (uncertainty_behaviour == "missing") {
      TipStateList <- lapply(X = TipStateList, function(x) {
        CellsToAlter <- grep("/", x)
        if (length(x = CellsToAlter) > 0) x[CellsToAlter] <- NA
        return(x)
      })
    }

    # If polymorphism_behaviour is random replace all polymorphisms with one included state at random:
    if (polymorphism_behaviour == "random") {
      TipStateList <- lapply(X = TipStateList, function(x) {
        CellsToAlter <- grep("&", x)
        if (length(x = CellsToAlter) > 0) for (i in 1:length(x = CellsToAlter)) x[CellsToAlter[i]] <- sample(strsplit(x[CellsToAlter[i]], split = "&")[[1]], size = 1)
        return(x)
      })
    }

    # If uncertainty_behaviour is random replace all uncertainties with one included state at random:
    if (uncertainty_behaviour == "random") {
      TipStateList <- lapply(X = TipStateList, function(x) {
        CellsToAlter <- grep("/", x)
        if (length(x = CellsToAlter) > 0) for (i in 1:length(x = CellsToAlter)) x[CellsToAlter[i]] <- sample(strsplit(x[CellsToAlter[i]], split = "/")[[1]], size = 1)
        return(x)
      })
    }

    # If uncertainty_behaviour is equalp replace all uncertainties withpolymorphism equivalent (makes grep simpler later):
    if (uncertainty_behaviour == "equalp") TipStateList <- lapply(X = TipStateList, function(x) gsub(pattern = "/", replacement = "&", x = x))

    # Reinsert modified tip states:
    x$TipStates <- TipStateList

    # Return list with modified tip states added:
    return(x)
  }

  # Alter polymorphisms and uncertainties according to polymorphism_behaviour and uncertainty_behaviour settings:
  CharacterList <- lapply(X = CharacterList, edit_polymorphisms, polymorphism_behaviour = polymorphism_behaviour, uncertainty_behaviour = uncertainty_behaviour)

  # Create pruned tipstates by removing all missing values:
  CharacterList <- lapply(X = CharacterList, function(x) {
    x$pruned_tip_states <- lapply(X = x$TipStates, function(y) y[!is.na(y)])
    x
  })

  # Get minimum and maximum values for each character (has to be done post polymorphism and inapplicable steps or will not work correctly):
  CharacterList <- lapply(X = CharacterList, function(x) {
    Ranges <- range(as.numeric(unlist(x = strsplit(unlist(x = x$pruned_tip_states), split = "&"))))
    x$minimum_values <- Ranges[1]
    x$maximum_values <- Ranges[2]
    x
  })

  # Subfunction to form tip states matrices ready for simmap function:
  build_tip_state_matrices <- function(x) {

    # Isolate pruned tip states list:
    PrunedTipStateList <- x$pruned_tip_states

    # Subfunction to build tip state matrix:
    build_tip_state_matrix <- function(y, minimum_values, maximum_values) {

      # Create empty tip state matrix:
      TipStateMatrix <- matrix(0, nrow = length(x = y), ncol = length(x = minimum_values:maximum_values), dimnames = list(names(y), as.character(minimum_values:maximum_values)))

      # Get row and column data as list:
      RowsAndColumns <- lapply(X = lapply(X = as.list(x = y), strsplit, split = "&"), unlist)

      # Isolate rows (taxon names):
      Rows <- unname(unlist(x = mapply(rep, names(RowsAndColumns), unlist(x = lapply(X = RowsAndColumns, length)))))

      # Isolate columns:
      Columns <- unname(unlist(x = RowsAndColumns))

      # Fill tip state matrix:
      for (i in 1:length(x = Rows)) TipStateMatrix[Rows[i], Columns[i]] <- 1

      # Convert each row to a probability (unnecessary if only one column):
      if (ncol(TipStateMatrix) > 1) TipStateMatrix <- t(apply(TipStateMatrix, 1, function(x) x / sum(x)))

      # Return tip state matrix:
      return(TipStateMatrix)
    }

    # Build tip state matrices:
    PrunedTipStateList <- lapply(X = PrunedTipStateList, build_tip_state_matrix, minimum_values = x$minimum_values, maximum_values = x$maximum_values)

    # Overwrite pruned tip states with tip state matrices:
    x$pruned_tip_states <- PrunedTipStateList

    # Return full output:
    return(x)
  }

  # Build tip state matrices:
  CharacterList <- lapply(X = CharacterList, build_tip_state_matrices)

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
    return(x)
  }

  # Build character models for each character:
  CharacterList <- lapply(X = CharacterList, build_character_model)

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
        TipEdge <- which(x = time_tree$edge[, 2] == tip_number)

        # Set all other edge lengths to zero:
        x$pruned_tree$edge.length[setdiff(x = 1:length(x = time_tree$edge), y = TipEdge)] <- 0

        # Reset root time to beginning of tip edge:
        x$pruned_tree$root.time <- unname(node_ages[time_tree$edge[TipEdge, 1]])
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
  CharacterList <- lapply(X = CharacterList, FUN = prune_tree, n_tips = n_tips)

  # Subfunction to get node ages for pruned tree:
  date_pruned_nodes <- function(x) {

    # Get number of (scorable) tips for pruned tree:
    n_tips <- lapply(X = x$pruned_tip_states, nrow)[[1]]

    # If exactly one tip:
    if (n_tips == 1) {

      # Set pruned node ages as beginning and end of single branch:
      x$Pruneddate_nodes <- c(x$pruned_tree$root.time - sum(x$pruned_tree$edge.length), x$pruned_tree$root.time)

      # Set node numbers as 1 for the tip and 2 for the node subtending the sole branch:
      names(x$Pruneddate_nodes) <- as.character(1:2)
    }

    # If exactly two tips:
    if (n_tips == 2) {

      # Set node ages by subtracting two branch lengths from root time and adding root time at end:
      x$Pruneddate_nodes <- c(x$pruned_tree$root.time - x$pruned_tree$edge.length[order(x$pruned_tree$edge[, 2])], x$pruned_tree$root.time)

      # Set node ages as 1 and 2 (for tips) and 3 for root:
      names(x$Pruneddate_nodes) <- as.character(1:3)
    }

    # If more than two tips just apply get node ages function:
    if (n_tips > 2) x$Pruneddate_nodes <- date_nodes(time_tree = x$pruned_tree)

    # Return full output:
    return(x)
  }

  # Get node ages for pruned tree:
  CharacterList <- lapply(X = CharacterList, date_pruned_nodes)

  # Subfunction to perform edge matches between pruned and full tree:
  match_edges <- function(x, tree) {

    # Get number of (scorable) tips for pruned tree:
    n_tips <- lapply(X = x$pruned_tip_states, nrow)[[1]]

    # If exactly one tip:
    if (n_tips == 1) {

      # Compile single branch output (i.e., single branch of tree with positive length):
      x$PrunedTofull_treeEdgeMatches <- list(which(x = x$pruned_tree$edge.length > 0))

      # Add name (has to be one because only one edge):
      names(x$PrunedTofull_treeEdgeMatches) <- "1"
    }

    # If exactly two tips:
    if (n_tips == 2) {

      # Get two tip names (in order):
      TipNames <- x$pruned_tree$tip.label

      # Get shared ancestor node on full tree:
      AncestorNode <- find_mrca(descendant_names = TipNames, tree = tree)

      # Get two tip numbers for two tips on full tree:
      tip_numbers <- unlist(x = lapply(X = lapply(X = as.list(x = TipNames), "==", tree$tip.label), which))

      # Create empty list ready to store edge matches:
      x$PrunedTofull_treeEdgeMatches <- list()

      # For each tip number:
      for (i in tip_numbers) {

        # Get first edge found (terminal branch):
        EdgesFound <- which(x = tree$edge[, 2] == i)

        # Set current node to ith node:
        CurrentNode <- i

        # While the ancestor has not been hit:
        while (tree$edge[EdgesFound[1], 1] != AncestorNode) {

          # Reset current node to start of current branch:
          CurrentNode <- tree$edge[EdgesFound[1], 1]

          # Add new edges found to vector:
          EdgesFound <- c(which(x = tree$edge[, 2] == CurrentNode), EdgesFound)
        }

        # Store edges found in root-to-tip order:
        x$PrunedTofull_treeEdgeMatches[[which(x = tip_numbers == i)]] <- EdgesFound
      }

      # Add edge names from pruned tree:
      names(x$PrunedTofull_treeEdgeMatches) <- as.character(1:2)
    }

    # If more than two tips simply use function normally:
    if (n_tips > 2) x$PrunedTofull_treeEdgeMatches <- match_tree_edges(tree, x$pruned_tree)$matching_edges

    # Return full output:
    return(x)
  }

  # Get edge matches between pruned and full trees (for recording true edge changes later):
  CharacterList <- lapply(X = CharacterList, match_edges, tree = time_tree)

  # Subfunction to get edge lengths in bins (usable later for rate calculations):
  get_binned_edge_lengths <- function(x) {

    # Set temporary tree as full tree (as modifying branch lengths but will rant to retain these later:
    TemporaryTree <- x$full_tree

    # Find any branch lengths on full tree to set as zero (effectively excluding them from the edge lengths as they correspond to missing values):
    BranchesToSetAsZeroes <- setdiff(x = 1:nrow(TemporaryTree$edge), y = unname(unlist(x = x$PrunedTofull_treeEdgeMatches)))

    # If there are nranches to set lengtsh to zero then do so and store:
    if (length(x = BranchesToSetAsZeroes) > 0) TemporaryTree$edge.length[BranchesToSetAsZeroes] <- 0

    # Get edge lengths in bins:
    x$EdgeLengthsInBins <- bin_edge_lengths(time_tree = TemporaryTree, time_bins = time_bins)

    # Return full output:
    return(x)
  }

  # Get edge lengths in bins:
  CharacterList <- lapply(X = CharacterList, get_binned_edge_lengths)

  # Subfunction to perform actual stochastic character maps:
  build_character_map_trees <- function(x) {

    # Get number of (scorable) tips for pruned tree:
    n_tips <- lapply(X = x$pruned_tip_states, nrow)[[1]]

    # If exactly one tip:
    if (n_tips == 1) {

      # Subfunction to generate stochastic character map like output but for the single taxon case:
      build_single_taxon_map <- function(y, tree) {

        # Set output as the tree initially:
        Output <- tree

        # Add maps to output of just the single branch's edge length:
        Output$maps <- list(sum(tree$edge.length))

        # Add state name to maps:
        Output$maps <- lapply(X = Output$maps, function(z) {
          names(z) <- colnames(x = y)
          z
        })

        # Return output:
        Output
      }

      # Perform stochastic character mapping and store output:
      x$StochasticCharacterMapTrees <- lapply(X = x$pruned_tip_states, build_single_taxon_map, tree = x$pruned_tree)
    }

    # If exactly two tips:
    if (n_tips == 2) {

      # Subfunction to generate stochastic character map like output but for the two taxon case:
      build_two_taxon_map <- function(y, tree) {

        # Set output as the tree initially:
        Output <- tree

        # Add maps to output of just the single branch's edge length:
        Output$maps <- as.list(x = unname(tree$edge.length))

        # If only one state (character is constant) add state name to maps:
        if (ncol(y) == 1) {
          Output$maps <- lapply(X = Output$maps, function(z) {
            names(z) <- colnames(x = y)
            z
          })
        }

        # If character is variant:
        if (ncol(y) == 2) {

          # Get raw data for calculating root state probabiity (incorporates tip sattes and reciprocal of branch lengths):
          RootStateRawData <- apply(sweep(y[tree$tip.label, ], MARGIN = 1, matrix(1 / (tree$edge.length / sum(tree$edge.length))), "*"), 2, sum)

          # Divide through by sum to get true probability:
          RootStateProbability <- RootStateRawData / sum(RootStateRawData)

          # Sample root state with probabilities based on tip state(s) and reciprocal of branch lengths:
          RootState <- sample(names(RootStateProbability), size = 1, prob = RootStateProbability)

          # Initially set map names to root state (one may change later):
          Output$maps <- lapply(X = Output$maps, function(z) {
            names(z) <- RootState
            z
          })

          # Get tip states with root column pruned (helps identify branch with changes):
          PrunedRootTipStates <- y[tree$tip.label, -which(x = colnames(x = y) == RootState), drop = FALSE]

          # Get tips with terminal changes (could concievably be empty if there is a polymorphism):
          TipsWithTerminalChanges <- names(which(x = apply(PrunedRootTipStates, 1, "==", 1)))

          # If there is a tip with a change to record:
          if (length(x = TipsWithTerminalChanges) > 0) {

            # Get the tip number:
            tip_number <- which(x = tree$tip.label == TipsWithTerminalChanges)

            # Get the edge the tip corresponds to:
            EdgeNumber <- which(x = tree$edge[, 2] == tip_number)

            # Get tip state:
            TipState <- colnames(x = PrunedRootTipStates[, which(x = PrunedRootTipStates[TipsWithTerminalChanges, ] == 1), drop = FALSE])

            # Get edge length:
            EdgeLength <- unname(Output$maps[[EdgeNumber]])

            # Pick a split point (proportion along branch where state change occurs):
            SplitPoint <- stats::runif(n = 1, min = 0, max = 1)

            # Cnvert proportions into absolute branch length values:
            Output$maps[[EdgeNumber]] <- c(SplitPoint, 1 - SplitPoint) * EdgeLength

            # Add state names to output:
            names(Output$maps[[EdgeNumber]]) <- c(RootState, TipState)
          }
        }

        # Return output:
        return(Output)
      }

      # Perform stochastic character mapping and store output:
      x$StochasticCharacterMapTrees <- lapply(X = x$pruned_tip_states, build_two_taxon_map, tree = x$pruned_tree)
    }

    # If more than two tips:
    if (n_tips > 2) {

      # Subfunction to perform stochastic character mapping:
      build_map_tree <- function(y, tree, model) {

        # If character is constant (invariant):
        if (ncol(y) == 1) {

          # Set initial output as tree:
          Output <- tree

          # Set maps as edge lengths:
          Output$maps <- as.list(x = tree$edge.length)

          # Add output names to maps as invariant character:
          Output$maps <- lapply(X = Output$maps, function(z) {
            names(z) <- colnames(x = y)
            z
          })
        }

        # If character is variable, perform regular stochastic character mapping:
        if (ncol(y) > 1) Output <- phytools::make.simmap(tree = tree, x = y, nsim = 1, model = model, pi = "estimated", message = FALSE)

        # Return output:
        return(Output)
      }

      # Perform stochastic character mapping and store output:
      x$StochasticCharacterMapTrees <- lapply(X = x$pruned_tip_states, build_map_tree, tree = x$pruned_tree, model = x$character_model)
    }

    # Return all output:
    return(x)
  }

  # Generate initial stochastic character map trees:
  CharacterList <- lapply(X = CharacterList, build_character_map_trees)

  # Subfunction to map stochastic chracter maps of pruned tree to full trees:
  map_characters_to_full_tree <- function(x) {

    # Only proceed if pruned tree is actually smaller than full tree (otherwise data are fine as is):
    if (lapply(X = x$pruned_tip_states, nrow)[[1]] < ape::Ntip(phy = x$full_tree)) {

      # Create empty list to store stochastic character maps for full trees:
      y <- list()

      # Start by filling out list with full tree:
      for (i in 1:n_simulations) y[[i]] <- x$full_tree

      # Generate null stochatsic character map from edge lengths:
      NullMap <- as.list(x = x$full_tree$edge.length)

      # Add NA as default state (i.e., missing data) - will want to correct this for inapplicables later:
      NullMap <- lapply(X = NullMap, function(z) {
        names(z) <- NA
        z
      })

      # Add null stochastic character map to list:
      for (i in 1:n_simulations) y[[i]]$maps <- NullMap

      # Find any single edge matches (where pruned edge matches a single edge in the full tree):
      SingleEdgeMatches <- unname(which(x = unlist(x = lapply(X = x$PrunedTofull_treeEdgeMatches, length)) == 1))

      # Find any multiple edge matches (where a pruned edge matches more than one edge in the full tree):
      MultipleEdgeMatches <- unname(which(x = unlist(x = lapply(X = x$PrunedTofull_treeEdgeMatches, length)) > 1))

      # If there is at least one single edge match then map pruned edges to full tree for them:
      if (length(x = SingleEdgeMatches) > 0) for (i in 1:length(x = y)) y[[i]]$maps[unname(unlist(x = x$PrunedTofull_treeEdgeMatches[SingleEdgeMatches]))] <- x$StochasticCharacterMapTrees[[i]]$maps[SingleEdgeMatches]

      # If there is at least one multiple edge match:
      if (length(x = MultipleEdgeMatches) > 0) {

        # For each simulation:
        for (i in 1:length(x = y)) {

          # Find multiple edge matches where pruned edge has no changes (single state persists):
          MultipleEdgeSingleStateMatches <- MultipleEdgeMatches[which(x = lapply(X = x$StochasticCharacterMapTrees[[i]]$maps[MultipleEdgeMatches], length) == 1)]

          # Find multiple edge matches where pruned edge has at least one change (multiple states sampled):
          MultipleEdgeMultipleStateMatches <- MultipleEdgeMatches[which(x = lapply(X = x$StochasticCharacterMapTrees[[i]]$maps[MultipleEdgeMatches], length) > 1)]

          # If there are multiple edge but single state matches:
          if (length(x = MultipleEdgeSingleStateMatches) > 0) {

            # For each match simply replace NA with the single character state:
            for (j in MultipleEdgeSingleStateMatches) {
              y[[i]]$maps[unname(unlist(x = x$PrunedTofull_treeEdgeMatches[j]))] <- lapply(X = y[[i]]$maps[unname(unlist(x = x$PrunedTofull_treeEdgeMatches[j]))], function(z) {
                names(z) <- names(x$StochasticCharacterMapTrees[[i]]$maps[j][[1]])
                z
              })
            }
          }

          # If there are multiple edge multiple state matches:
          if (length(x = MultipleEdgeMultipleStateMatches) > 0) {

            # For each multiple state multiple edge match:
            for (j in MultipleEdgeMultipleStateMatches) {

              # Get matching edges on full tree:
              MatchingEdges <- unname(unlist(x = x$PrunedTofull_treeEdgeMatches[j]))

              # Get edge lengths of matching full tree edges:
              MatchingEdgeLengths <- unname(unlist(x = y[[i]]$maps[unname(unlist(x = x$PrunedTofull_treeEdgeMatches[j]))]))

              # Get pruned stochastic maps for current pruned edge:
              PrunedStochasticMaps <- x$StochasticCharacterMapTrees[[i]]$maps[j][[1]]

              # Get starting age of current pruned edge (will use to get change times that can be binned by full tree edge later):
              StartAgeOfPrunedEdge <- unname(x$Pruneddate_nodes[x$pruned_tree$edge[j, 1]])

              # Get times at which character changes occur:
              change_times <- unname(StartAgeOfPrunedEdge - cumsum(PrunedStochasticMaps[1:(length(x = PrunedStochasticMaps) - 1)]))

              # Build matrix of from-to changes:
              FromToChanges <- cbind(names(PrunedStochasticMaps[1:(length(x = PrunedStochasticMaps) - 1)]), names(PrunedStochasticMaps[2:length(x = PrunedStochasticMaps)]))

              # Set matching edges on full tree as time bins:
              MatchingEdgesAstime_bins <- c(StartAgeOfPrunedEdge, StartAgeOfPrunedEdge - cumsum(MatchingEdgeLengths))

              # Get edge on whicb change occurs:
              ChangeEdges <- unlist(x = lapply(X = as.list(x = change_times), function(z) min(which(x = z > MatchingEdgesAstime_bins)) - 1))

              # Create full tree stochastic character map of correct size:
              FullStochasticMaps <- lapply(X = as.list(x = rle(sort(x = c(ChangeEdges, 1:length(x = MatchingEdges))))$lengths), function(z) rep(0, z))

              # Set current time as beginning of pruned edge:
              CurrentTime <- StartAgeOfPrunedEdge

              # Set current edge as 1 (will increment through while loop below):
              CurrentEdge <- 1

              # Set current state as first from value in from-to matrix:
              CurrentState <- FromToChanges[1, 1]

              # Make vector of edge switch times:
              EdgeSwitchTimes <- MatchingEdgesAstime_bins[2:length(x = MatchingEdgesAstime_bins)]

              # Whilst there are still changes or edge switches left to deal with:
              while (length(x = c(EdgeSwitchTimes, change_times)) > 0) {

                # Set next event time:
                NextEvent <- max(c(EdgeSwitchTimes, change_times))

                # If next event is to switch edges:
                if (EdgeSwitchTimes[1] == NextEvent) {

                  # Find current position on stochastic map:
                  CurrentPosition <- which(x = FullStochasticMaps[[CurrentEdge]] == 0)[1]

                  # Store edge length:
                  FullStochasticMaps[[CurrentEdge]][CurrentPosition] <- CurrentTime - NextEvent

                  # Store current state:
                  names(FullStochasticMaps[[CurrentEdge]])[CurrentPosition] <- CurrentState

                  # Update current time:
                  CurrentTime <- NextEvent

                  # Update current edge:
                  CurrentEdge <- CurrentEdge + 1

                  # Update EdgeSwitchTimes by removing last change:
                  EdgeSwitchTimes <- EdgeSwitchTimes[-1]

                  # If next event is a change of character state:
                } else {

                  # Find current position on stochastic map:
                  CurrentPosition <- which(x = FullStochasticMaps[[CurrentEdge]] == 0)[1]

                  # Store edge length:
                  FullStochasticMaps[[CurrentEdge]][CurrentPosition] <- CurrentTime - NextEvent

                  # Store current state:
                  names(FullStochasticMaps[[CurrentEdge]])[CurrentPosition] <- CurrentState

                  # As long as the from-to matrix still exists:
                  if (nrow(FromToChanges) > 0) {

                    # Update current state:
                    CurrentState <- FromToChanges[1, 2]

                    # Prune change from from-to matrix:
                    FromToChanges <- FromToChanges[-1, , drop = FALSE]
                  }

                  # Update current time:
                  CurrentTime <- NextEvent

                  # Prune change time from vector:
                  change_times <- change_times[-1]
                }
              }

              # Store full stochstic character map in y:
              y[[i]]$maps[unname(unlist(x = x$PrunedTofull_treeEdgeMatches[j]))] <- FullStochasticMaps
            }
          }
        }
      }

      # Overwrite pruned tree stochastic character maps with full tree stochastic character maps:
      x$StochasticCharacterMapTrees <- y
    }

    # Return full output:
    return(x)
  }

  # Map stochastic characters to full tree in preparation for recording changes:
  CharacterList <- lapply(X = CharacterList, map_characters_to_full_tree)

  # Subfunction to extract character changes matrix from trees:
  reformat_changes <- function(x) {

    # Find the root edge for the tree:
    RootEdge <- match(ape::Ntip(phy = x$full_tree) + 1, x$full_tree$edge[, 1])

    # Subfunction to get character changes and root state:
    extract_changes <- function(y) {

      # Get root state:
      RootState <- as.numeric(names(y$maps[[RootEdge]][1]))

      # Find any edges with changes on them:
      EdgesWithChanges <- which(x = unlist(x = lapply(X = y$maps, length)) > 1)

      # As long as there is at least one change:
      if (length(x = EdgesWithChanges) > 0) {

        # Get ages at start of edges with changes (subtracting from this will give change times later:
        AgeAtStartOfEdgesWithChanges <- node_ages[x$full_tree$edge[EdgesWithChanges, 1]]

        # Get from and to states for each change:
        FromsAndTos <- matrix(as.numeric(unlist(x = strsplit(unlist(x = lapply(X = y$maps[EdgesWithChanges], function(z) paste(names(z[1:(length(x = z) - 1)]), names(z[2:length(x = z)]), sep = "%%"))), split = "%%"))), ncol = 2, byrow = TRUE)

        # Get character change times:
        Characterchange_times <- unname(unlist(x = mapply("-", AgeAtStartOfEdgesWithChanges, lapply(X = y$maps[EdgesWithChanges], function(z) cumsum(z[1:(length(x = z) - 1)])))))

        # Build changes matrix:
        ChangesMatrix <- cbind(FromsAndTos, unlist(x = mapply(rep, EdgesWithChanges, unlist(x = lapply(X = y$maps[EdgesWithChanges], length)) - 1)), Characterchange_times)

        # Add column names to matrix:
        colnames(x = ChangesMatrix) <- c("From", "To", "Edge", "Time")

        # If no changes occur:
      } else {

        # Create changes matrix with no changes:
        ChangesMatrix <- matrix(nrow = 0, ncol = 4, dimnames = list(c(), c("From", "To", "Edge", "Time")))
      }

      # Build matrix of from to states for each edge (i.e., state at start of edge and state at end of edge):
      EdgeFromTo <- cbind(as.numeric(unlist(x = lapply(X = y$maps, function(y) names(y[1])))), as.numeric(unlist(x = lapply(X = y$maps, function(y) names(y[length(x = y)])))))

      # Get preceding states for each edge (will help identify edge-to-edge changes:
      PrecedingState <- EdgeFromTo[match(time_tree$edge[, 1], time_tree$edge[, 2]), 2]

      # Get following states for each edge (will help identify edge-to-edge changes:
      FollowingState <- EdgeFromTo[match(time_tree$edge[, 2], time_tree$edge[, 1]), 1]

      # Correct following edge to eliminate terminal changes whic do not need to be recorded:
      FollowingState[match(1:n_tips, time_tree$edge[, 2])] <- EdgeFromTo[match(1:n_tips, time_tree$edge[, 2]), 2]

      # Get preceding changes (will want to add to changes matrix):
      PrecedingChanges <- which(x = apply(cbind(is.na(PrecedingState), !is.na(EdgeFromTo[, 1])), 1, all))

      # Get following changes (will want to add to changes matrix):
      FollowingChanges <- which(x = apply(cbind(!is.na(EdgeFromTo[, 2]), is.na(FollowingState)), 1, all))

      # For each unique preceding change:
      for (i in PrecedingChanges[!duplicated(time_tree$edge[PrecedingChanges, 1])]) {

        # If simply the root state:
        if (n_tips + 1 == time_tree$edge[i, 1]) {

          # Add root change to matrix as edge zero:
          ChangesMatrix <- rbind(ChangesMatrix, c(NA, EdgeFromTo[i, 1], 0, time_tree$root.time))

          # If some other state:
        } else {

          # Add NA to something state change to changes matrix:
          ChangesMatrix <- rbind(ChangesMatrix, c(NA, EdgeFromTo[i, 1], which(x = time_tree$edge[, 2] == time_tree$edge[i, 1]), unname(node_ages[time_tree$edge[i, 1]])))
        }
      }

      # If there are any following changes (i.e., is there at least one missing value at the tips):
      if (length(x = FollowingChanges) > 0) {

        # Add following changes to changes matrix:
        ChangesMatrix <- rbind(ChangesMatrix, cbind(EdgeFromTo[FollowingChanges, 2], FollowingState[FollowingChanges], unlist(x = lapply(X = as.list(x = time_tree$edge[FollowingChanges, 2]), function(z) {
          FollowingEdges <- which(x = time_tree$edge[, 1] == z)
          FollowingEdges[is.na(EdgeFromTo[FollowingEdges, 1])]
        })), unname(node_ages[time_tree$edge[FollowingChanges, 2]])))
      }

      # Return output:
      return(ChangesMatrix)
    }

    # Get root and changes for each stochastic character map and store:
    x$Stochasticcharacter_changes <- lapply(X = x$StochasticCharacterMapTrees, extract_changes)

    # Return output:
    return(x)
  }

  # Get root state and character changes for each stochastic character map:
  CharacterList <- lapply(X = CharacterList, reformat_changes)

  # Subfunction to collapse changes across simulations into a single matrix:
  compile_changes_as_matrix <- function(x) {

    # Get simulation numbers (will be new column in matrix):
    SimulationNumbers <- unlist(x = mapply(rep, 1:n_simulations, unlist(x = lapply(X = x$Stochasticcharacter_changes, function(y) nrow(y)))))

    # Collapse changes into single matrix:
    ChangesMatrix <- do.call(what = rbind, args = lapply(X = x$Stochasticcharacter_changes, function(y) y))

    # Add simulation numbers:
    ChangesMatrix <- cbind(matrix(SimulationNumbers), ChangesMatrix)

    # Add column name:
    colnames(x = ChangesMatrix)[1] <- "SimulationNumber"

    # Overwrite stochastic character matrices with new collapsed format:
    x$Stochasticcharacter_changes <- ChangesMatrix

    # Return full output:
    return(x)
  }

  # Collapse stochastic character matrices to single matrix for each character:
  CharacterList <- lapply(X = CharacterList, compile_changes_as_matrix)

  # Compile all state changes into a single matrix:
  all_state_changes <- cbind(matrix(unlist(x = mapply(rep, 1:ncol(MatrixBlock), lapply(X = CharacterList, function(x) nrow(x$Stochasticcharacter_changes))))), do.call(what = rbind, args = lapply(X = CharacterList, function(x) x$Stochasticcharacter_changes)))

  # Add column name for character:
  colnames(x = all_state_changes)[1] <- "Character"

  # Get rid of pesky state rownames:
  rownames(x = all_state_changes) <- NULL

  # Sort all state changes by edge number:
  all_state_changes <- all_state_changes[order(all_state_changes[, "Edge"]), ]

  # Sort all state changes by character number:
  all_state_changes <- all_state_changes[order(all_state_changes[, "Character"]), ]

  # Sort all state changes by simulation number:
  all_state_changes <- all_state_changes[order(all_state_changes[, "SimulationNumber"]), ]

  # Get character times (length of subtrees for each character):
  character_times <- unlist(x = lapply(X = CharacterList, function(x) sum(x$EdgeLengthsInBins$binned_edge_lengths)))

  # Get edge length for each bin by character:
  binned_edge_lengths <- do.call(what = rbind, args = lapply(X = CharacterList, function(x) x$EdgeLengthsInBins$binned_edge_lengths))

  # Get edge length for each bin by character (terminal branches only):
  binned_terminal_edge_lengths <- do.call(what = rbind, args = lapply(X = CharacterList, function(x) x$EdgeLengthsInBins$binned_terminal_edge_lengths))

  # Get edge length for each bin by character (internal branches only):
  binned_internal_edge_lengths <- do.call(what = rbind, args = lapply(X = CharacterList, function(x) x$EdgeLengthsInBins$binned_internal_edge_lengths))

  # Compile output as list:
  output <- list(all_state_changes = all_state_changes, character_times = character_times, binned_edge_lengths = binned_edge_lengths, binned_terminal_edge_lengths = binned_terminal_edge_lengths, binned_internal_edge_lengths = binned_internal_edge_lengths)

  # Return output:
  return(invisible(output))
}
