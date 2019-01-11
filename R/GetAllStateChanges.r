#' Finds all state changes on a tree using stochastic character mapping
#'
#' @description
#'
#' Takes a cladistic matrix and time-scaled tree and makes point estimates for every character change using stochastic character mapping.
#'
#' @param CladisticMatrix A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @param Tree A time-scaled tree (phylo object) that represents the relationships of the taxa in \code{CladisticMatrix}.
#' @param TimeBins A vector of ages representing the boundaries of a series of time bins.
#' @param NSimulations The number of simulations to perform (passed to \code{make.simmap}.
#' @param PolymorphismBehaviour What to do with polymorphic (&) characters. One of "equalp", "missing", or "random". See details.
#' @param UncertaintyBehaviour What to do with uncertain (/) characters. One of "equalp", "missing", or "random". See details.
#' @param InapplicableBehaviour What to do with inapplicable characters. Only one option currently ("missing"). See details.
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
#' \item{RootStates}{A matrix of the root states for each character (column) and simulation (rows).}
#' \item{AllStateChanges}{A matrix of rows for each change with columns corresponding to the character, the simulation number, the edge number, the time the change occurred, and the start and end states.}
#' \item{CharacterTimes}{A vector of the sampled tree-length (in Ma) for each character.}
#' \item{EdgeLengthsPerBin}{A matrix of time bins (columns) and characters (rows) indicating the sampled tree-length (in Ma).}
#' \item{TerminalEdgeLengthsPerBin}{As above, but for terminal edges only.}
#' \item{InternalEdgeLengthsPerBin}{As above, but for internal edges only.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Set random seed:
#' set.seed(17)
#'
#' # Generate a random tree for the Michaux data set:
#' Tree <- rtree(nrow(Michaux1989$Matrix_1$Matrix))
#'
#' # Update taxon names to match those in the data matrix:
#' Tree$tip.label <- rownames(Michaux1989$Matrix_1$Matrix)
#'
#' # Set root time by making youngest taxon extant:
#' Tree$root.time <- max(diag(vcv(Tree)))
#'
#' # Get all state changes:
#' GetAllStateChanges(Michaux1989, Tree,
#'   seq(Tree$root.time, 0, length.out = 3), NSimulations = 2)
#'
#' @export GetAllStateChanges
GetAllStateChanges <- function(CladisticMatrix, Tree, TimeBins, NSimulations = 10, PolymorphismBehaviour = "equalp", UncertaintyBehaviour = "equalp", InapplicableBehaviour = "missing") {
  
  # IMPROVE CUSTOMISATION OF MAKE.SIMMAP WITH OPTIONS FOR PI, Q ETC.
  # AND ONLY PERFORM SCM ON UNIQUE STATE DISTRIBUTON-CHARACTER TYPE COMBOS
  # MAJOR ISSUE IS NO EASY WAY TO MAKE MODEL FOR ORDERED MULTISTATE CHARACTER WHEN NOT ALL STATES ARE FOUND AT TIPS (E.G., 0 and 2 sampled, but not 1)
  # CHECK FOR ALL NA CHARACTERS AS THESE WILL NEED TO BE REMOVED.
  
  
  
  #CladisticMatrix <- Day2016
  #CladisticMatrix <- MatrixPruner(CladisticMatrix = CladisticMatrix, blocks2prune = 1)
  #Tree <- rtree(nrow(Day2016$Matrix_1$Matrix))
  #Tree$tip.label <- rownames(Day2016$Matrix_1$Matrix)
  #Tree$root.time <- max(diag(vcv(Tree)))
  #TimeBins <- seq(Tree$root.time, 0, length.out = 3)
  #NSimulations <- 2
  #PolymorphismBehaviour = "equalp"
  #UncertaintyBehaviour = "equalp"
  #InapplicableBehaviour = "missing"
  
  
  
  # Check for continuous and step matrices and stop and warn user if found:
  if(length(setdiff(unique(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], function(x) x$Ordering))), c("unord", "ord"))) > 0) stop("CladisticMatrix can only contain characters of type \"ord\" or \"unord\" (i.e., no step matrices or continuous characters).")
  
  # Check tree has branch lengths:
  if(is.null(Tree$edge.length)) stop("Tree does not have branch lengths (durations). Try timescaling the Tree, e.g., with DatePhylo.")
  
  # Check Tree has root age:
  if(is.null(Tree$root.time)) stop("Tree is missing $root.time. Try setting this before continuing, e.g., Tree$root.time <- 104.2.")
  
  # Check PolymorphismBehaviour is correctly formatted or stop and warn user:
  if(length(setdiff(PolymorphismBehaviour, c("equalp", "missing", "random"))) > 0) stop("PolymorphismBehaviour must be one of must be one of \"equalp\", \"missing\" or \"random\".")
  
  # Check UncertaintyBehaviour is correctly formatted or stop and warn user:
  if(length(setdiff(UncertaintyBehaviour, c("equalp", "missing", "random"))) > 0) stop("UncertaintyBehaviour must be one of \"equalp\", \"missing\" or \"random\".")
  
  # Check InapplicableBehaviour is correctly formatted or stop and warn user:
  if(length(setdiff(InapplicableBehaviour, c("missing"))) > 0) stop("InapplicableBehaviour must be \"missing\".")
  
  # Ensure time bins are in correct order:
  TimeBins <- sort(unique(TimeBins), decreasing = TRUE)

  # Get tree node ages:
  TreeNodeAges <- GetNodeAges(Tree)
  
  # Build all data into single matrix:
  MatrixBlock <- do.call(cbind, lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"))
  
  # Assemble all ordering values into a single vector:
  Ordering <- unname(do.call(c, lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Ordering")))
  
  # Assemble all minimum values into a single vector:
  MinVals <- unname(do.call(c, lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MinVals")))
  
  # Assemble all maximum values into a single vector:
  MaxVals <- unname(do.call(c, lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MaxVals")))
  
  # Assemble all maximum values into a single vector:
  Weights <- unname(do.call(c, lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Weights")))
  
  # Build each character into list values starting with tip state lists (of N Simulations in length):
  CharacterList <- lapply(lapply(apply(MatrixBlock, 2, list), unlist), function(x) {y <- list(); y$TipStates <- lapply(apply(matrix(rep(x, times = NSimulations), ncol = NSimulations, dimnames = list(names(x), c())), 2, list), unlist); return(y)})
  
  # Add weight to list for each character:
  for(i in 1:length(CharacterList)) CharacterList[[i]]$Weight <- Weights[i]
  
  # Add ordering to list for each character:
  for(i in 1:length(CharacterList)) CharacterList[[i]]$Ordering <- Ordering[i]
  
  # Add full tree to list for each character:
  for(i in 1:length(CharacterList)) CharacterList[[i]]$FullTree <- Tree
  
  # Subfunction to perform polymorphism and uncertainty changes:
  EditPolymorphisms <- function(x, PolymorphismBehaviour, UncertaintyBehaviour) {
    
    # Isolate tip states list:
    TipStateList <- x$TipStates
    
    # If PolymorphismBehaviour is missing replace all polymorphisms with NAs:
    if(PolymorphismBehaviour == "missing") TipStateList <- lapply(TipStateList, function(x) {CellsToAlter <- grep("&", x); if(length(CellsToAlter) > 0) x[CellsToAlter] <- NA; return(x)})
    
    # If UncertaintyBehaviour is missing replace all uncertainties with NAs:
    if(UncertaintyBehaviour == "missing") TipStateList <- lapply(TipStateList, function(x) {CellsToAlter <- grep("/", x); if(length(CellsToAlter) > 0) x[CellsToAlter] <- NA; return(x)})
    
    # If PolymorphismBehaviour is random replace all polymorphisms with one included state at random:
    if(PolymorphismBehaviour == "random") TipStateList <- lapply(TipStateList, function(x) {CellsToAlter <- grep("&", x); if(length(CellsToAlter) > 0) for(i in 1:length(CellsToAlter)) x[CellsToAlter[i]] <- sample(strsplit(x[CellsToAlter[i]], split = "&")[[1]], size = 1); return(x)})
    
    # If UncertaintyBehaviour is random replace all uncertainties with one included state at random:
    if(UncertaintyBehaviour == "random") TipStateList <- lapply(TipStateList, function(x) {CellsToAlter <- grep("/", x); if(length(CellsToAlter) > 0) for(i in 1:length(CellsToAlter)) x[CellsToAlter[i]] <- sample(strsplit(x[CellsToAlter[i]], split = "/")[[1]], size = 1); return(x)})
    
    # If UncertaintyBehaviour is equalp replace all uncertainties withpolymorphism equivalent (makes grep simpler later):
    if(UncertaintyBehaviour == "equalp") TipStateList <- lapply(TipStateList, function(x) gsub("/", "&", x))

    # Reinsert modified tip states:
    x$TipStates <- TipStateList
    
    # Return list with modified tip states added:
    return(x)
    
  }
  
  # Alter polymorphisms and uncertainties according to PolymorphismBehaviour and UncertaintyBehaviour settings:
  CharacterList <- lapply(CharacterList, EditPolymorphisms, PolymorphismBehaviour = PolymorphismBehaviour, UncertaintyBehaviour = UncertaintyBehaviour)
  
  # If InapplicableBehaviour is missing replace inaplicables with NAs:
  if(InapplicableBehaviour == "missing") CharacterList <- lapply(CharacterList, function(x) {x$TipStates <- lapply(x$TipStates, function(y) {Inapplicables <- which(y == ""); if(length(Inapplicables) > 0) y[Inapplicables] <- NA; y}); x})
  
  # Create pruned tipstates by removing all missing values:
  CharacterList <- lapply(CharacterList, function(x) {x$PrunedTipStates <- lapply(x$TipStates, function(y) y[!is.na(y)]); x})
  
  # Get minimum and maximum values for each character (has to be done post polymorphism and inapplicable steps or will not work correctly):
  CharacterList <- lapply(CharacterList, function(x) {Ranges <- range(as.numeric(unlist(strsplit(unlist(x$PrunedTipStates), split = "&")))); x$MinVal <- Ranges[1]; x$MaxVal <- Ranges[2]; x})
  
  # Subfunction to form tip states matrices ready for simmap function:
  BuildTipStateMatrices <- function(x) {
    
    # Isolate pruned tip states list:
    PrunedTipStateList <- x$PrunedTipStates
    
    # Subfunction to build tip state matrix:
    BuildTipStateMatrix <- function(y, MinVal, MaxVal) {
      
      # Create empty tip state matrix:
      TipStateMatrix <- matrix(0, nrow = length(y), ncol = length(MinVal:MaxVal), dimnames = list(names(y), as.character(MinVal:MaxVal)))
      
      # Get row and column data as list:
      RowsAndColumns <- lapply(lapply(as.list(y), strsplit, split = "&"), unlist)
      
      # Isolate rows (taxon names):
      Rows <- unname(unlist(mapply(rep, names(RowsAndColumns), unlist(lapply(RowsAndColumns, length)))))
      
      # Isolate columns:
      Columns <- unname(unlist(RowsAndColumns))
      
      # Fill tip state matrix:
      for(i in 1:length(Rows)) TipStateMatrix[Rows[i], Columns[i]] <- 1
      
      # Convert each row to a probability
      TipStateMatrix <- t(apply(TipStateMatrix, 1, function(x) x / sum(x)))
      
      # Return tip state matrix:
      return(TipStateMatrix)

    }
    
    # Build tip state matrices:
    PrunedTipStateList <- lapply(PrunedTipStateList, BuildTipStateMatrix, MinVal = x$MinVal, MaxVal = x$MaxVal)
    
    # Overwrite pruned tip states with tip state matrices:
    x$PrunedTipStates <- PrunedTipStateList
    
    # Return full output:
    return(x)
    
  }
  
  # Build tip state matrices:
  CharacterList <- lapply(CharacterList, BuildTipStateMatrices)
  
  # Subfunction to build character model:
  BuildCharacterModel <- function(x) {
    
    # If character is ordered:
    if(x$Ordering == "ord") {
      
      # Create character model with all transitions impossible:
      CharacterModel <- matrix(0, ncol = length(x$MinVal:x$MaxVal), nrow = length(x$MinVal:x$MaxVal), dimnames = list(x$MinVal:x$MaxVal, x$MinVal:x$MaxVal))
      
      # Create character matrix with only off-diagonal (adjacent) states possible:
      if(ncol(CharacterModel) > 1) for(i in 2:ncol(CharacterModel)) CharacterModel[i, (i - 1)] <- CharacterModel[(i - 1), i] <- 1
      
    }
    
    # If character is unordered:
    if(x$Ordering == "unord") {
      
      # Create character model with all transitions possible:
      CharacterModel <- matrix(1, ncol = length(x$MinVal:x$MaxVal), nrow = length(x$MinVal:x$MaxVal), dimnames = list(x$MinVal:x$MaxVal, x$MinVal:x$MaxVal))
      
      # Set diagonals to zero:
      diag(CharacterModel) <- 0
      
    }
    
    # Store character model:
    x$CharacterModel <- CharacterModel
    
    # Return full output:
    return(x)
  
  }
  
  # Build character models for each character:
  CharacterList <- lapply(CharacterList, BuildCharacterModel)
  
  # Subfunction to prune tree to just tips with data:
  PruneTree <- function(x) {
    
    # Find any tips to drop:
    TipsToDrop <- setdiff(x$FullTree$tip.label, rownames(x$PrunedTipStates[[1]]))
    
    # MODIFY BELOW TO DEAL WITH CASES OF PRUNING ALL BUT ONE TIP OR ALL BUT TWO TIPS:
    
    # If there are tips to drop:
    if(length(TipsToDrop) > 0) {
      
      # Prune tree and store:
      x$PrunedTree <- drop.tip(x$FullTree, TipsToDrop)
      
      # Ensure pruned trees $root.time value is correct:
      x$PrunedTree <- CorrectRootTime(Tree, x$PrunedTree)

    # If no tips need to be dropped:
    } else {
      
      # Store full tree as pruned version:
      x$PrunedTree <- x$FullTree
      
    }
    
    # Return full output:
    return(x)
  
  }
  
  # Get pruned trees with only tips from pruned tip states returned:
  CharacterList <- lapply(CharacterList, PruneTree)
  
  # Get node ages for pruned tree:
  CharacterList <- lapply(CharacterList, function(x) {x$PrunedNodeAges <- GetNodeAges(x$PrunedTree); x})
  
  # Subfunction to perform actual stochastic character maps:
  BuildStochasticCharacterMapTrees <- function(x) {
    
    # CONDITIONAL IF TOO FEW TIPS? NEED TO FILL THESE OUT LATER:
    if(ape::Ntip(x$PrunedTree) == 1) stop("Only one tip scored for a character.") # WILL HAVE TO GENERATE MAPS LIKE OUTPUT FOR THI LATER PLUS MODIFY BRANCH LENGTHS IN PRUNED TREES TO GET EDGE LENGTSH RIGHT?
    if(ape::Ntip(x$PrunedTree) == 2) stop("Only two tips scored for a character.") # WILL HAVE TO GENERATE MAPS LIKE OUTPUT FOR THI LATER PLUS MODIFY BRANCH LENGTHS IN PRUNED TREES TO GET EDGE LENGTSH RIGHT?
    
    # Apply make.simmap function to pruned data:
    if(ape::Ntip(x$PrunedTree) > 2) x$StochasticCharacterMapTrees <- lapply(x$PrunedTipStates, function(y) make.simmap(y, tree = x$PrunedTree, nsim = 1, model = x$CharacterModel, pi = "estimated", message = FALSE))
    
    # Return all output:
    return(x)
    
  }
  
  # Generate initial stochastic character map trees:
  CharacterList <- lapply(CharacterList, BuildStochasticCharacterMapTrees)
  
  # Subfunction to extract character changes matrix from trees:
  ReformatCharacterChanges <- function(x) {
    
    # Find the root edge for the pruned tree:
    RootEdge <- match(ape::Ntip(x$PrunedTree) + 1, x$PrunedTree$edge[, 1])
    
    # Subfunction to get character changes and root state:
    ExtractCharacterChanges <- function(y) {
      
      # Get root state:
      RootState <- as.numeric(names(y$maps[[RootEdge]][1]))
      
      # Find any edges with changes on them:
      EdgesWithChanges <- which(unlist(lapply(y$maps, length)) > 1)
      
      # As long as there is at least one change:
      if(length(EdgesWithChanges) > 0) {
        
        # Get ages at start of edges with changes (subtracting from this will give change times later:
        AgeAtStartOfEdgesWithChanges <- x$PrunedNodeAges[x$PrunedTree$edge[EdgesWithChanges, 1]]
        
        # Get from and to states for each change:
        FromsAndTos <- matrix(as.numeric(unlist(strsplit(unlist(lapply(y$maps[EdgesWithChanges], function(z) paste(names(z[1:(length(z) - 1)]), names(z[2:length(z)]), sep = "%%"))), split = "%%"))), ncol = 2, byrow = TRUE)
        
        # Get character change times:
        CharacterChangeTimes <- unname(unlist(mapply('-', AgeAtStartOfEdgesWithChanges, lapply(y$maps[EdgesWithChanges], function(z) cumsum(z[1:(length(z) - 1)])))))
        
        # Build changes matrix:
        ChangesMatrix <- cbind(FromsAndTos, unlist(mapply(rep, EdgesWithChanges, unlist(lapply(y$maps[EdgesWithChanges], length)) - 1)), CharacterChangeTimes)
        
        # Add column names to matrix:
        colnames(ChangesMatrix) <- c("From", "To", "Edge", "Time")
        
      # If no changes occur:
      } else {
        
        # Create changes matrix with no changes:
        ChangesMatrix <- matrix(nrow = 0, ncol = 4, dimnames = list(c(), c("From", "To", "Edge", "Time")))
        
      }
      
      # Compile output (root and changes):
      StochasticCharacterMapOutput <- list(RootState, ChangesMatrix)
      
      # Name output (root and changes):
      names(StochasticCharacterMapOutput) <- c("RootState", "ChangesMatrix")
      
      # Return output:
      return(StochasticCharacterMapOutput)
      
    }
    
    # Get root and changes for each stochastic character map and store:
    x$StochasticCharacterChanges <- lapply(x$StochasticCharacterMapTrees, ExtractCharacterChanges)
    
    # Return output:
    return(x)
    
  }
  
  # Get root state and character changes for each stochastic character map:
  CharacterList <- lapply(CharacterList, ReformatCharacterChanges)
  
  # Get edge matches between pruned and full trees (for recording true edge changes later):
  CharacterList <- lapply(CharacterList, function(x) {x$PrunedToFullTreeEdgeMatches <- EdgeMatch(Tree, x$PrunedTree)$matching.edges; x})
  
  # Get edge lengths in bins:
  CharacterList <- lapply(CharacterList, function(x) {x$EdgeLengthsInBins <- EdgeLengthsInBins(Tree, time.bins = TimeBins, pruned.tree = x$PrunedTree); x})
  
  # Subfunction to collapse changes across simulations into a single matrix (plus vector of root states):
  CompileChangesIntoSingleMatrix <- function(x) {
    
    # Get simulation numbers (will be new column in matrix):
    SimulationNumbers <- unlist(mapply(rep, 1:NSimulations, unlist(lapply(x$StochasticCharacterChanges, function(y) nrow(y$ChangesMatrix)))))
    
    # Collapse changes into single matrix:
    ChangesMatrix <- do.call(rbind, lapply(x$StochasticCharacterChanges, function(y) y$ChangesMatrix))
    
    # Add simulation numbers:
    ChangesMatrix <- cbind(matrix(SimulationNumbers), ChangesMatrix)
    
    # Add column name:
    colnames(ChangesMatrix)[1] <- "SimulationNumber"
    
    # Get vector of root states:
    RootStates <- unlist(lapply(x$StochasticCharacterChanges, function(y) y$RootState))
    
    # Compile output:
    Output <- list(RootStates, ChangesMatrix)
    
    # Add names to output:
    names(Output) <- c("RootStates", "ChangesMatrix")
    
    # Overwrite stochastic character matrices with new collapsed format:
    x$StochasticCharacterChanges <- Output
    
    # Return full output:
    return(x)
    
  }
  
  # Collapse stochastic character matrices to single matrix for each character:
  CharacterList <- lapply(CharacterList, CompileChangesIntoSingleMatrix)
  
  
  
  
  
  
  
  
  # NEED NA TO CODED TRANISTIONS AND NA TO SOMETHING TRANSITIONS
  # IF UNORDERED AND NOT ALL STATES SAMPLED CAN JUST REMOVE THOSE STATES
  # IF ORDERED AND UNSAMPLED STATES THEN NEED TO STOP AND WARN USER
  # ANY REMAINING POLYMORPHISMS ARE FOR EQUAL P
  # IF USING EQUALP OR RANDOM AT END THEN NEED TO RECORD WEIRD CHANGE OF, SAY, 0 TO 0&1
  # IF USING MISSING NEED TO RECORD NA TO 0&1 CHANGE
  # GONNA HAVE TO DEAL WITH N SIMULATIONS SPREAD OVER MULTIPLE TIP STATES (SOME WILL BE UNIQUE, OTHERS WILL NEED DUPLICATION)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Get rows that correspond to current characters only:
  #char.rows <- sort(unlist(lapply(lapply(as.list(MissingStringCharacters[[i]]), '==', allstatechanges[, "Character"]), which)))
  
  # Find edges in character tree that correspond to multiple edges in original tree:
  #edgeswithmultiplematches <- as.numeric(names(which(unlist(lapply(edgematches$matching.edges[unique(allstatechanges[char.rows, "Edge"])], length)) > 1)))
  
  # Find edges in character tree that correspond to single edges in original tree:
  #edgeswithsinglematches <- as.numeric(names(which(unlist(lapply(edgematches$matching.edges[unique(allstatechanges[char.rows, "Edge"])], length)) == 1)))
  
  # Get single edge match replacemnt edge numbers:
  #singlematchreplacements <- unlist(edgematches$matching.edges[edgeswithsinglematches])
  
  # Create vector of length one to store edge replacements:
  #edgereplacements <- NA
  
  # If there is more than one edge to replace create longer vector of NAs to store edge replacements:
  #if(length(char.rows) > 1) edgereplacements <- rep(NA, nrow(allstatechanges[char.rows, ]))
  
  # As long as there are single edge replacements find single match replacements and update them:
  #if(length(edgeswithsinglematches) > 0) for(j in edgeswithsinglematches) edgereplacements[which(allstatechanges[char.rows, "Edge"] == j)] <- singlematchreplacements[match(j, edgeswithsinglematches)]
  
  # As long as there are multiple replacement edges:
  #if(length(edgeswithmultiplematches) > 0) {
  
  # For each edge with multiple replacements:
  #for(j in edgeswithmultiplematches) {
  
  # Get matching edges as edge numbers:
  #matching.edges <- edgematches$matching.edges[[as.character(j)]]
  
  # Get matching edges as to-from node numbers:
  #matching.edges.matrix <- Tree$edge[matching.edges, ]
  
  # Get node ages for matching edges:
  #matching.edges.ages <- cbind(tree.nodeages[matching.edges.matrix[, 1]], tree.nodeages[matching.edges.matrix[, 2]])
  
  # Get rows from allstatechanges which will need to be changed:
  #rowstochange <- which(allstatechanges[char.rows, "Edge"] == j)
  
  # For each matching edge:
  #for(k in 1:nrow(matching.edges.matrix)) {
  
  # Get maximum edge age:
  #max.edge.age <- matching.edges.ages[k, 1]
  
  # Get minimum edge age:
  #min.edge.age <- matching.edges.ages[k, 2]
  
  # Get changes on current edge:
  #changes.on.current.edge <- intersect(which(allstatechanges[char.rows[rowstochange], "Age"] > min.edge.age), which(allstatechanges[char.rows[rowstochange], "Age"] <= max.edge.age))
  
  # Update edge replacements:
  #edgereplacements[which(allstatechanges[char.rows, "Edge"] == j)[changes.on.current.edge]] <- matching.edges[k]
  
  
  # Update state changes with edges from complete tree:
  #allstatechanges[char.rows, "Edge"] <- edgereplacements
  
  # Anything with a row sum of zero in this list can be ignored:
  #edgelinks <- matrix(FindLinkedEdges(Tree)[edgematches$removed.edges, -edgematches$removed.edges], nrow = length(edgematches$removed.edges), dimnames = list(edgematches$removed.edges, colnames(FindLinkedEdges(Tree))[-edgematches$removed.edges]))
  
  # Get list of missing edges that will have a node on a sampled edge:
  #missingedges <- as.numeric(rownames(edgelinks)[which(apply(edgelinks, 1, sum) > 0)])
  
  # Case if missing edge terminates at root of pruned tree:
  #if(length(sort(match(Tree$edge[missingedges, 2], setdiff(Tree$edge[sort(unlist(edgematches$matching.edges)), 1], Tree$edge[sort(unlist(edgematches$matching.edges)), 2])))) > 0) {
  
  # Find node number of root in original tree:
  #pruned.root.node <- setdiff(Tree$edge[sort(unlist(edgematches$matching.edges)), 1], Tree$edge[sort(unlist(edgematches$matching.edges)), 2])
  
  # Find root terminating edge:
  #root.terminating.edge <- missingedges[match(pruned.root.node, Tree$edge[missingedges, 2])]
  
  # Add root edge to allstatechanges:
  #allstatechanges <- rbind(allstatechanges, cbind(sort(rep(MissingStringCharacters[[i]], NSimulations)), rep(1:NSimulations, ncol(root.states)), rep(root.terminating.edge, length(root.states)), rep(pruned.tree$root.time, length(root.states)), rep(NA, length(root.states)), as.vector(root.states)))
  
  # Update missing edges (remove edge prior to root or pruned tree):
  #missingedges <- setdiff(missingedges, root.terminating.edge)
  
  # If root is sampled in pruned tree:
  #} else {
  
  # Ensure root state (on edge "0") is recorded even if no other changes occur:
  #allstatechanges <- rbind(allstatechanges, cbind(sort(rep(MissingStringCharacters[[i]], NSimulations)), rep(1:NSimulations, ncol(root.states)), rep(0, length(root.states)), rep(Tree$root.time, length(root.states)), rep(NA, length(root.states)), as.vector(root.states)))
  
  #}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Compile all root states into a single matrix:
  RootStates <- do.call(cbind, lapply(CharacterList, function(x) x$StochasticCharacterChanges$RootStates))
  
  # Compile all state changes into a singel matrix:
  AllStateChanges <- cbind(matrix(unlist(mapply(rep, 1:ncol(MatrixBlock), lapply(CharacterList, function(x) nrow(x$StochasticCharacterChanges$ChangesMatrix))))), do.call(rbind, lapply(CharacterList, function(x) x$StochasticCharacterChanges$ChangesMatrix)))
  
  # Add column name for character:
  colnames(AllStateChanges)[1] <- "Character"
  
  # Get rid of pesky state rownames:
  rownames(AllStateChanges) <- NULL
  
  # Sort all state changes by edge number:
  AllStateChanges <- AllStateChanges[order(AllStateChanges[, "Edge"]), ]
  
  # Sort all state changes by character number:
  AllStateChanges <- AllStateChanges[order(AllStateChanges[, "Character"]), ]
  
  # Sort all state changes by simulation number:
  AllStateChanges <- AllStateChanges[order(AllStateChanges[, "SimulationNumber"]), ]
  
  # Get character times (length of subtrees for each character):
  CharacterTimes <- unlist(lapply(CharacterList, function(x) sum(x$EdgeLengthsInBins$edge.length.in.bin)))
  
  # Get edge length for each bin by character:
  EdgeLengthsPerBin <- do.call(rbind, lapply(CharacterList, function(x) x$EdgeLengthsInBins$edge.length.in.bin))
  
  # Get edge length for each bin by character (terminal branches only):
  TerminalEdgeLengthsPerBin <- do.call(rbind, lapply(CharacterList, function(x) x$EdgeLengthsInBins$terminal.edge.length.in.bin))
  
  # Get edge length for each bin by character (internal branches only):
  InternalEdgeLengthsPerBin <- do.call(rbind, lapply(CharacterList, function(x) x$EdgeLengthsInBins$internal.edge.length.in.bin))
  
  # Compile output as list:
  output <- list(RootStates, AllStateChanges, CharacterTimes, EdgeLengthsPerBin, TerminalEdgeLengthsPerBin, InternalEdgeLengthsPerBin)
  
  # Add names to output:
  names(output) <- c("RootStates", "AllStateChanges", "CharacterTimes", "EdgeLengthsPerBin", "TerminalEdgeLengthsPerBin", "InternalEdgeLengthsPerBin")
  
  # Return output:
  return(output)
  
}
