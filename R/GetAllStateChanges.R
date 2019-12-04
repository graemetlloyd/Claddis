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
#' set.seed(2)
#'
#' # Use Day 2016 as source matrix:
#' CladisticMatrix <- Day2016
#'
#' # Prune out continuous characters:
#' CladisticMatrix <- MatrixPruner(CladisticMatrix =
#'   CladisticMatrix, blocks2prune = 1)
#'
#' # Prune out majority of characters so
#' # example runs quickly:
#' CladisticMatrix <- MatrixPruner(CladisticMatrix =
#'   CladisticMatrix, characters2prune = 1:32)
#'
#' # Generete random tree for matrix taxa:
#' Tree <- rtree(nrow(Day2016$Matrix_1$Matrix))
#'
#' # Add taxon names to tree:
#' Tree$tip.label <- rownames(Day2016$Matrix_1$Matrix)
#'
#' # Add root age to tree:
#' Tree$root.time <- max(diag(vcv(Tree)))
#'
#' # Get all state changes for two simulations:
#' StateChanges <-
#'   GetAllStateChanges(CladisticMatrix = CladisticMatrix,
#'   Tree = Tree, TimeBins = seq(Tree$root.time, 0,
#'   length.out = 3), NSimulations = 2)
#'
#' # View matrix of all stochstic character changes:
#' StateChanges$AllStateChanges
#'
#' # View vector of sampled time for each
#' # character:
#' StateChanges$CharacterTimes
#'
#' # View matrix of edge lengths in each time bin:
#' StateChanges$EdgeLengthsPerBin
#'
#' # View matrix of termnial edge lengths in each time bin:
#' StateChanges$TerminalEdgeLengthsPerBin
#'
#' # View matrix of internal edge lengths in each time bin:
#' StateChanges$InternalEdgeLengthsPerBin
#'
#' @export GetAllStateChanges
GetAllStateChanges <- function(CladisticMatrix, Tree, TimeBins, NSimulations = 10, PolymorphismBehaviour = "equalp", UncertaintyBehaviour = "equalp", InapplicableBehaviour = "missing") {
  
  # IMPROVE CUSTOMISATION OF MAKE.SIMMAP WITH OPTIONS FOR PI, Q ETC. (FLAT PRIOR ON ROOT MAY BE PARTICULARLY BAD? ALLOW MAYBE SKEWING TOWARDS OUTGROUP STATE AS SOME KIND OF SLIDING VALUE?).
  # AND ONLY PERFORM SCM ON UNIQUE STATE DISTRIBUTON-CHARACTER TYPE COMBOS
  # MAJOR ISSUE IS NO EASY WAY TO MAKE MODEL FOR ORDERED MULTISTATE CHARACTER WHEN NOT ALL STATES ARE FOUND AT TIPS (E.G., 0 and 2 sampled, but not 1)
  # CHECK FOR ALL NA CHARACTERS AS THESE WILL NEED TO BE REMOVED.
  # MOVE MISSING AS POLYMORPHISM/UNCERTIANT BEHAVIOUR TO TOP AS MORE EFFICIENT.
  # ADD INAPPLICABLE OPTION THAT TIES SWITCH TO "" TO DEPENDENT CHARACTER A LA MORPHDISTMATRIX APPROACH.
  # MODEL AND TIP STATE DIMENSIONS MAY VARY IF SAY ONLY VARIANCE IS A POLYMORPHIC CHARACTER BUT "RANDOM" IS USED???
  
  # ANY REMAINING POLYMORPHISMS ARE FOR EQUAL P
  # IF USING EQUALP OR RANDOM AT END THEN NEED TO RECORD WEIRD CHANGE OF, SAY, 0 TO 0&1
  # IF USING MISSING NEED TO RECORD NA TO 0&1 CHANGE
  # GONNA HAVE TO DEAL WITH N SIMULATIONS SPREAD OVER MULTIPLE TIP STATES (SOME WILL BE UNIQUE, OTHERS WILL NEED DUPLICATION)
  # SMEARING BACK AND SMEARING FORWARD SOMETHIG TO NITHING AND NOTHING TO SOMETHING CHANGES MAY BE DIFFERENT. PERHAPS HAVE OPTION TO TREAT TIMESTAMP FOR THESE TO VARY.
  # IF ADDING TREES TO OUTPUT CONVERT THEM TO CLASS MULTIPHYLO, E.G. STOCHASTIC CHARACTER MAPS WITH NAS - NOTE THIS IN MANUAL TOO AS AN EXTENSION OF WHAT PHTTOOLS DOES.
  
  # Check for continuous and step matrices and stop and warn user if found:
  if(length(setdiff(unique(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], function(x) x$Ordering))), c("unord", "ord"))) > 0) stop("CladisticMatrix can only contain characters of type \"ord\" or \"unord\" (i.e., no step matrices or continuous characters).")
  
  # Check tree has branch lengths:
  if(is.null(Tree$edge.length)) stop("Tree does not have branch lengths (durations). Try timescaling the Tree, e.g., with DatePhylo.")
  
  # Check branches all have positive length:
  if(any(Tree$edge.length == 0)) stop("All branch lengths must be positive (no zero-length branches).")
  
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
  
  # If InapplicableBehaviour is missing replace inaplicables with NAs:
  if(InapplicableBehaviour == "missing" && any(MatrixBlock[!is.na(MatrixBlock)] == "")) MatrixBlock[which(MatrixBlock == "")] <- NA
  
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
      
      # Convert each row to a probability (unnecessary if only one column):
      if(ncol(TipStateMatrix) > 1) TipStateMatrix <- t(apply(TipStateMatrix, 1, function(x) x / sum(x)))
      
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
  
  # Get total tips in tree:
  TipsInTree <- ape::Ntip(Tree)

  # Subfunction to prune tree to just tips with data:
  PruneTree <- function(x) {
    
    # Find any tips to drop:
    TipsToDrop <- setdiff(x$FullTree$tip.label, rownames(x$PrunedTipStates[[1]]))
    
    # If there are tips to drop:
    if(length(TipsToDrop) > 0) {
      
      # If exactly one tip will remain:
      if((TipsInTree - length(TipsToDrop)) == 1) {
        
        # Set pruned tree initially as full tree:
        x$PrunedTree <- x$FullTree
        
        # Set tip number (of single coded taxon):
        TipNumber <- which(Tree$tip.label == setdiff(Tree$tip.label, TipsToDrop))
        
        # Find single tip edge:
        TipEdge <- which(Tree$edge[, 2] == TipNumber)
        
        # Set all other edge lengths to zero:
        x$PrunedTree$edge.length[setdiff(1:length(Tree$edge), TipEdge)] <- 0
        
        # Reset root time to beginning of tip edge:
        x$PrunedTree$root.time <- unname(TreeNodeAges[Tree$edge[TipEdge, 1]])
        
      }
      
      # If exactly two tips will remain:
      if((TipsInTree - length(TipsToDrop)) == 2) {
        
        # Prune tree and store:
        x$PrunedTree <- drop.tip(x$FullTree, TipsToDrop)
        
        # Correct root time manually:
        x$PrunedTree$root.time <- unname(TreeNodeAges[FindAncestor(setdiff(Tree$tip.label, TipsToDrop), Tree)])
        
      }
      
      # If at least three tips will remain:
      if((TipsInTree - length(TipsToDrop)) > 2) {
        
        # Prune tree and store:
        x$PrunedTree <- drop.tip(x$FullTree, TipsToDrop)
        
        # Ensure pruned trees $root.time value is correct:
        x$PrunedTree <- CorrectRootTime(Tree, x$PrunedTree)
        
      }
      
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
  
  # Subfunction to get node ages for pruned tree:
  GetPrunedNodeAges <- function(x) {
    
    # Get number of (scorable) tips for pruned tree:
    NumberOfTips <- lapply(x$PrunedTipStates, nrow)[[1]]
    
    # If exactly one tip:
    if(NumberOfTips == 1) {
      
      # Set pruned node ages as beginning and end of single branch:
      x$PrunedNodeAges <- c(x$PrunedTree$root.time - sum(x$PrunedTree$edge.length), x$PrunedTree$root.time)
      
      # Set node numbers as 1 for the tip and 2 for the node subtending the sole branch:
      names(x$PrunedNodeAges) <- as.character(1:2)
      
    }
    
    # If exactly two tips:
    if(NumberOfTips == 2) {
      
      # Set node ages by subtracting two branch lengths from root time and adding root time at end:
      x$PrunedNodeAges <- c(x$PrunedTree$root.time - x$PrunedTree$edge.length[order(x$PrunedTree$edge[, 2])], x$PrunedTree$root.time)
      
      # Set node ages as 1 and 2 (for tips) and 3 for root:
      names(x$PrunedNodeAges) <- as.character(1:3)
      
    }
    
    # If more than two tips just apply get node ages function:
    if(NumberOfTips > 2) x$PrunedNodeAges <- GetNodeAges(x$PrunedTree)
    
    # Return full output:
    return(x)
    
  }
  
  # Get node ages for pruned tree:
  CharacterList <- lapply(CharacterList, GetPrunedNodeAges)
  
  # Subfunction to perform edge matches between pruned and full tree:
  PerformEdgeMatches <- function(x, tree) {
    
    # Get number of (scorable) tips for pruned tree:
    NumberOfTips <- lapply(x$PrunedTipStates, nrow)[[1]]
    
    # If exactly one tip:
    if(NumberOfTips == 1) {
      
      # Compile single branch output (i.e., single branch of tree with positive length):
      x$PrunedToFullTreeEdgeMatches <- list(which(x$PrunedTree$edge.length > 0))
      
      # Add name (has to be one because only one edge):
      names(x$PrunedToFullTreeEdgeMatches) <- "1"
      
    }
    
    # If exactly two tips:
    if(NumberOfTips == 2) {
      
      # Get two tip names (in order):
      TipNames <- x$PrunedTree$tip.label
      
      # Get shared ancestor node on full tree:
      AncestorNode <- FindAncestor(descs = TipNames, tree = tree)
      
      # Get two tip numbers for two tips on full tree:
      TipNumbers <- unlist(lapply(lapply(as.list(TipNames), '==', tree$tip.label), which))
      
      # Create empty list ready to store edge matches:
      x$PrunedToFullTreeEdgeMatches <- list()
      
      # For each tip number:
      for(i in TipNumbers) {
        
        # Get first edge found (terminal branch):
        EdgesFound <- which(tree$edge[, 2] == i)
        
        # Set current node to ith node:
        CurrentNode <- i
        
        # While the ancestor has not been hit:
        while(tree$edge[EdgesFound[1], 1] != AncestorNode) {
          
          # Reset current node to start of current branch:
          CurrentNode <- tree$edge[EdgesFound[1], 1]
          
          # Add new edges found to vector:
          EdgesFound <- c(which(tree$edge[, 2] == CurrentNode), EdgesFound)
          
        }
        
        # Store edges found in root-to-tip order:
        x$PrunedToFullTreeEdgeMatches[[which(TipNumbers == i)]] <- EdgesFound
        
      }
      
      # Add edge names from pruned tree:
      names(x$PrunedToFullTreeEdgeMatches) <- as.character(1:2)
      
    }
    
    # If more than two tips simply use function normally:
    if(NumberOfTips > 2) x$PrunedToFullTreeEdgeMatches <- EdgeMatch(Tree, x$PrunedTree)$matching.edges
    
    # Return full output:
    return(x)
    
  }
  
  # Get edge matches between pruned and full trees (for recording true edge changes later):
  CharacterList <- lapply(CharacterList, PerformEdgeMatches, tree = Tree)
  
  # Subfunction to get edge lengths in bins (usable later for rate calculations):
  GetEdgeLengthsInBins <- function(x) {
    
    # Set temporary tree as full tree (as modifying branch lengths but will rant to retain these later:
    TemporaryTree <- x$FullTree
    
    # Find any branch lengtsh on full tree to set as zero (effectively excluding them from the edge lengths as they correspond to missing values):
    BranchesToSetAsZeroes <- setdiff(1:nrow(TemporaryTree$edge), unname(unlist(x$PrunedToFullTreeEdgeMatches)))
    
    # If there are nranches to set lengtsh to zero then do so and store:
    if(length(BranchesToSetAsZeroes) > 0) TemporaryTree$edge.length[BranchesToSetAsZeroes] <- 0
    
    # Get edge lengths in bins:
    x$EdgeLengthsInBins <- EdgeLengthsInBins(tree = TemporaryTree, time.bins = TimeBins)
    
    # Return full output:
    return(x)
    
  }
  
  # Get edge lengths in bins:
  CharacterList <- lapply(CharacterList, GetEdgeLengthsInBins)
  
  # Subfunction to perform actual stochastic character maps:
  BuildStochasticCharacterMapTrees <- function(x) {
    
    # Get number of (scorable) tips for pruned tree:
    NumberOfTips <- lapply(x$PrunedTipStates, nrow)[[1]]
    
    # If exactly one tip:
    if(NumberOfTips == 1) {
      
      # Subfunction to generate stochastic character map like output but for the single taon case:
      BuildOneTaxonStochasticCharacterMap <- function(y, tree) {
        
        # Set output as the tree initially:
        Output <- tree
        
        # Add maps to output of just the single branch's edge length:
        Output$maps <- list(sum(tree$edge.length))
        
        # Add state name to maps:
        Output$maps <- lapply(Output$maps, function(z) {names(z) <- colnames(y); z})
        
        # Return output:
        return(Output)
        
      }
      
      # Perform stochastic character mapping and store output:
      x$StochasticCharacterMapTrees <- lapply(x$PrunedTipStates, BuildOneTaxonStochasticCharacterMap, tree = x$PrunedTree)
      
    }
    
    # If exactly two tips:
    if(NumberOfTips == 2) {
      
      # Subfunction to generate stochastic character map like output but for the two taxon case:
      BuildTwoTaxonStochasticCharacterMap <- function(y, tree) {
        
        # Set output as the tree initially:
        Output <- tree
        
        # Add maps to output of just the single branch's edge length:
        Output$maps <- as.list(unname(tree$edge.length))
        
        # If only one state (character is constant) add state name to maps:
        if(ncol(y) == 1) Output$maps <- lapply(Output$maps, function(z) {names(z) <- colnames(y); z})
        
        # If character is variant:
        if(ncol(y) == 2) {
          
          # Get raw data for calculating root state probabiity (incorporates tip sattes and reciprocal of branch lengths):
          RootStateRawData <- apply(sweep(y[tree$tip.label, ], MARGIN = 1, matrix(1 / (tree$edge.length / sum(tree$edge.length))), '*'), 2, sum)
          
          # Divide through by sum to get true probability:
          RootStateProbability <- RootStateRawData / sum(RootStateRawData)
          
          # Sample root state with probabilities based on tip state(s) and reciprocal of branch lengths:
          RootState <- sample(names(RootStateProbability), size = 1, prob = RootStateProbability)
          
          # Initially set map names to root state (one may change later):
          Output$maps <- lapply(Output$maps, function(z) {names(z) <- RootState; z})
          
          # Get tip states with root column pruned (helps identify branch with changes):
          PrunedRootTipStates <- y[tree$tip.label, -which(colnames(y) == RootState), drop = FALSE]
          
          # Get tips with terminal changes (could concievably be empty if there is a polymorphism):
          TipsWithTerminalChanges <- names(which(apply(PrunedRootTipStates, 1, '==', 1)))
          
          # If there is a tip with a change to record:
          if(length(TipsWithTerminalChanges) > 0) {
            
            # Get the tip number:
            TipNumber <- which(tree$tip.label == TipsWithTerminalChanges)
            
            # Get the edge the tip corresponds to:
            EdgeNumber <- which(tree$edge[, 2] == TipNumber)
            
            # Get tip state:
            TipState <- colnames(PrunedRootTipStates[, which(PrunedRootTipStates[TipsWithTerminalChanges, ] == 1), drop = FALSE])
            
            # Get edge length:
            EdgeLength <- unname(Output$maps[[EdgeNumber]])
            
            # Pick a split point (proportion along branch where state change occurs):
            SplitPoint <- runif(n = 1, min = 0, max = 1)
            
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
      x$StochasticCharacterMapTrees <- lapply(x$PrunedTipStates, BuildTwoTaxonStochasticCharacterMap, tree = x$PrunedTree)
      
    }
    
    # If more than two tips:
    if(NumberOfTips > 2) {
      
      # Subfunction to perform stochastic character mapping:
      GetStochasticCharacterMapTree <- function(y, tree, model) {
        
        # If character is constant (invariant):
        if(ncol(y) == 1) {
          
          # Set initial output as tree:
          Output <- tree
          
          # Set maps as edge lengths:
          Output$maps <- as.list(tree$edge.length)
          
          # Add output names to maps as invariant character:
          Output$maps <- lapply(Output$maps, function(z) {names(z) <- colnames(y); z})
          
        }
        
        # If character is variable, perform regular stochastic character mapping:
        if(ncol(y) > 1) Output <- make.simmap(y, tree = tree, nsim = 1, model = model, pi = "estimated", message = FALSE)
        
        # Return output:
        return(Output)
        
      }
      
      # Perform stochastic character mapping and store output:
      x$StochasticCharacterMapTrees <- lapply(x$PrunedTipStates, GetStochasticCharacterMapTree, tree = x$PrunedTree, model = x$CharacterModel)
      
    }
    
    # Return all output:
    return(x)
    
  }
  
  # Generate initial stochastic character map trees:
  CharacterList <- lapply(CharacterList, BuildStochasticCharacterMapTrees)
  
  # Subfunction to map stochastic chracter maps of pruned tree to full trees:
  MapStochasticCharactersToFullTree <- function(x) {
    
    # Only proceed if pruned tree is actually smaller than full tree (otherwise data are fine as is):
    if(lapply(x$PrunedTipStates, nrow)[[1]] < ape::Ntip(x$FullTree)) {
      
      # Create empty list to store stochastic character maps for full trees:
      y <- list()
      
      # Start by filling out list with full tree:
      for(i in 1:NSimulations) y[[i]] <- x$FullTree
      
      # Generate null stochatsic character map from edge lengths:
      NullMap <- as.list(x$FullTree$edge.length)
      
      # Add NA as default state (i.e., missing data) - will want to correct this for inapplicables later:
      NullMap <- lapply(NullMap, function(z) {names(z) <- NA; z})
      
      # Add null stochastic character map to list:
      for(i in 1:NSimulations) y[[i]]$maps <- NullMap
      
      # Find any single edge matches (where pruned edge matches a single edge in the full tree):
      SingleEdgeMatches <- unname(which(unlist(lapply(x$PrunedToFullTreeEdgeMatches, length)) == 1))
      
      # Find any multiple edge matches (where a pruned edge matches more than one edge in the full tree):
      MultipleEdgeMatches <- unname(which(unlist(lapply(x$PrunedToFullTreeEdgeMatches, length)) > 1))
      
      # If there is at least one single edge match then map pruned edges to full tree for them:
      if(length(SingleEdgeMatches) > 0) for(i in 1:length(y)) y[[i]]$maps[unname(unlist(x$PrunedToFullTreeEdgeMatches[SingleEdgeMatches]))] <- x$StochasticCharacterMapTrees[[i]]$maps[SingleEdgeMatches]
      
      # If there is at least one multiple edge match:
      if(length(MultipleEdgeMatches) > 0) {
        
        # For each simulation:
        for(i in 1:length(y)) {
          
          # Find multiple edge matches where pruned edge has no changes (single state persists):
          MultipleEdgeSingleStateMatches <- MultipleEdgeMatches[which(lapply(x$StochasticCharacterMapTrees[[i]]$maps[MultipleEdgeMatches], length) == 1)]
          
          # Find multiple edge matches where pruned edge has at least one change (multiple states sampled):
          MultipleEdgeMultipleStateMatches <- MultipleEdgeMatches[which(lapply(x$StochasticCharacterMapTrees[[i]]$maps[MultipleEdgeMatches], length) > 1)]
          
          # If there are multiple edge but single state matches:
          if(length(MultipleEdgeSingleStateMatches) > 0) {
            
            # For each match simply replace NA with the single character state:
            for(j in MultipleEdgeSingleStateMatches) y[[i]]$maps[unname(unlist(x$PrunedToFullTreeEdgeMatches[j]))] <- lapply(y[[i]]$maps[unname(unlist(x$PrunedToFullTreeEdgeMatches[j]))], function(z) {names(z) <- names(x$StochasticCharacterMapTrees[[i]]$maps[j][[1]]); z})
            
          }
          
          # If there are multiple edge multiple state matches:
          if(length(MultipleEdgeMultipleStateMatches) > 0) {
            
            # For each multiple state multiple edge match:
            for(j in MultipleEdgeMultipleStateMatches) {
              
              # Get matching edges on full tree:
              MatchingEdges <- unname(unlist(x$PrunedToFullTreeEdgeMatches[j]))
              
              # Get edge lengths of matching full tree edges:
              MatchingEdgeLengths <- unname(unlist(y[[i]]$maps[unname(unlist(x$PrunedToFullTreeEdgeMatches[j]))]))
              
              # Get pruned stochastic maps for current pruned edge:
              PrunedStochasticMaps <- x$StochasticCharacterMapTrees[[i]]$maps[j][[1]]
              
              # Get starting age of current pruned edge (will use to get change times that can be binned by full tree edge later):
              StartAgeOfPrunedEdge <- unname(x$PrunedNodeAges[x$PrunedTree$edge[j, 1]])
              
              # Get times at which character changes occur:
              ChangeTimes <- unname(StartAgeOfPrunedEdge - cumsum(PrunedStochasticMaps[1:(length(PrunedStochasticMaps) - 1)]))
              
              # Build matrix of from-to changes:
              FromToChanges <- cbind(names(PrunedStochasticMaps[1:(length(PrunedStochasticMaps) - 1)]), names(PrunedStochasticMaps[2:length(PrunedStochasticMaps)]))
              
              # Set matching edges on full tree as time bins:
              MatchingEdgesAsTimeBins <- c(StartAgeOfPrunedEdge, StartAgeOfPrunedEdge - cumsum(MatchingEdgeLengths))
              
              # Get edge on whicb change occurs:
              ChangeEdges <- unlist(lapply(as.list(ChangeTimes), function(z) min(which(z > MatchingEdgesAsTimeBins)) - 1))
              
              # Create full tree stochastic character map of correct size:
              FullStochasticMaps <- lapply(as.list(rle(sort(c(ChangeEdges, 1:length(MatchingEdges))))$lengths), function(z) rep(0, z))
              
              # Set current time as beginning of pruned edge:
              CurrentTime <- StartAgeOfPrunedEdge
              
              # Set current edge as 1 (will increment through while loop below):
              CurrentEdge <- 1
              
              # Set current state as first from value in from-to matrix:
              CurrentState <- FromToChanges[1, 1]
              
              # Make vector of edge switch times:
              EdgeSwitchTimes <- MatchingEdgesAsTimeBins[2:length(MatchingEdgesAsTimeBins)]
              
              # Whilst there are still changes or edge switches left to deal with:
              while(length(c(EdgeSwitchTimes, ChangeTimes)) > 0) {
                
                # Set next event time:
                NextEvent <- max(c(EdgeSwitchTimes, ChangeTimes))
                
                # If next event is to switch edges:
                if(EdgeSwitchTimes[1] == NextEvent) {
                  
                  # Find current position on stochastic map:
                  CurrentPosition <- which(FullStochasticMaps[[CurrentEdge]] == 0)[1]
                  
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
                  CurrentPosition <- which(FullStochasticMaps[[CurrentEdge]] == 0)[1]
                  
                  # Store edge length:
                  FullStochasticMaps[[CurrentEdge]][CurrentPosition] <- CurrentTime - NextEvent

                  # Store current state:
                  names(FullStochasticMaps[[CurrentEdge]])[CurrentPosition] <- CurrentState
                  
                  # As long as the from-to matrix still exists:
                  if(nrow(FromToChanges) > 0) {
                    
                    # Update current state:
                    CurrentState <- FromToChanges[1, 2]
                    
                    # Prune change from from-to matrix:
                    FromToChanges <- FromToChanges[-1, , drop = FALSE]
                    
                  }
                  
                  # Update current time:
                  CurrentTime <- NextEvent
                  
                  # Prune change time from vector:
                  ChangeTimes <- ChangeTimes[-1]
                  
                }
                
              }
              
              # Store full stochstic character map in y:
              y[[i]]$maps[unname(unlist(x$PrunedToFullTreeEdgeMatches[j]))] <- FullStochasticMaps
              
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
  CharacterList <- lapply(CharacterList, MapStochasticCharactersToFullTree)

  # Subfunction to extract character changes matrix from trees:
  ReformatCharacterChanges <- function(x) {
    
    # Find the root edge for the tree:
    RootEdge <- match(ape::Ntip(x$FullTree) + 1, x$FullTree$edge[, 1])
    
    # Subfunction to get character changes and root state:
    ExtractCharacterChanges <- function(y) {
      
      # Get root state:
      RootState <- as.numeric(names(y$maps[[RootEdge]][1]))
      
      # Find any edges with changes on them:
      EdgesWithChanges <- which(unlist(lapply(y$maps, length)) > 1)
      
      # As long as there is at least one change:
      if(length(EdgesWithChanges) > 0) {
        
        # Get ages at start of edges with changes (subtracting from this will give change times later:
        AgeAtStartOfEdgesWithChanges <- TreeNodeAges[x$FullTree$edge[EdgesWithChanges, 1]]
        
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
      
      # Build matrix of from to states for each edge (i.e., state at start of edge and state at end of edge):
      EdgeFromTo <- cbind(as.numeric(unlist(lapply(y$maps, function(y) names(y[1])))), as.numeric(unlist(lapply(y$maps, function(y) names(y[length(y)])))))
      
      # Get preceding states for each edge (will help identify edge-to-edge changes:
      PrecedingState <- EdgeFromTo[match(Tree$edge[, 1], Tree$edge[, 2]), 2]
      
      # Get following states for each edge (will help identify edge-to-edge changes:
      FollowingState <- EdgeFromTo[match(Tree$edge[, 2], Tree$edge[, 1]), 1]
      
      # Correct following edge to eliminate terminal changes whic do not need to be recorded:
      FollowingState[match(1:ape::Ntip(Tree), Tree$edge[, 2])] <- EdgeFromTo[match(1:ape::Ntip(Tree), Tree$edge[, 2]), 2]
      
      # Get preceding changes (will want to add to changes matrix):
      PrecedingChanges <- which(apply(cbind(is.na(PrecedingState), !is.na(EdgeFromTo[, 1])), 1, all))
      
      # Get following changes (will want to add to changes matrix):
      FollowingChanges <- which(apply(cbind(!is.na(EdgeFromTo[, 2]), is.na(FollowingState)), 1, all))
      
      # For each unique preceding change:
      for(i in PrecedingChanges[!duplicated(Tree$edge[PrecedingChanges, 1])]) {
        
        # If simply the root state:
        if(ape::Ntip(Tree) + 1 == Tree$edge[i, 1]) {
          
          # Add root change to matrix as edge zero:
          ChangesMatrix <- rbind(ChangesMatrix, c(NA, EdgeFromTo[i, 1], 0, Tree$root.time))
          
        # If some other state:
        } else {
          
          # Add NA to something state change to changes matrix:
          ChangesMatrix <- rbind(ChangesMatrix, c(NA, EdgeFromTo[i, 1], which(Tree$edge[, 2] == Tree$edge[i, 1]), unname(TreeNodeAges[Tree$edge[i, 1]])))
          
        }
      
      }
      
      # If there are any following changes (i.e., is there at least one missing value at the tips):
      if(length(FollowingChanges) > 0) {
        
        # Add following changes to changes matrix:
        ChangesMatrix <- rbind(ChangesMatrix, cbind(EdgeFromTo[FollowingChanges, 2], FollowingState[FollowingChanges], unlist(lapply(as.list(Tree$edge[FollowingChanges, 2]), function(z) {FollowingEdges <- which(Tree$edge[, 1] == z); FollowingEdges[is.na(EdgeFromTo[FollowingEdges, 1])]})), unname(TreeNodeAges[Tree$edge[FollowingChanges, 2]])))
        
      }
      
      # Return output:
      return(ChangesMatrix)
      
    }
    
    # Get root and changes for each stochastic character map and store:
    x$StochasticCharacterChanges <- lapply(x$StochasticCharacterMapTrees, ExtractCharacterChanges)
    
    # Return output:
    return(x)
    
  }
  
  # Get root state and character changes for each stochastic character map:
  CharacterList <- lapply(CharacterList, ReformatCharacterChanges)
  
  # Subfunction to collapse changes across simulations into a single matrix:
  CompileChangesIntoSingleMatrix <- function(x) {
    
    # Get simulation numbers (will be new column in matrix):
    SimulationNumbers <- unlist(mapply(rep, 1:NSimulations, unlist(lapply(x$StochasticCharacterChanges, function(y) nrow(y)))))
    
    # Collapse changes into single matrix:
    ChangesMatrix <- do.call(rbind, lapply(x$StochasticCharacterChanges, function(y) y))
    
    # Add simulation numbers:
    ChangesMatrix <- cbind(matrix(SimulationNumbers), ChangesMatrix)
    
    # Add column name:
    colnames(ChangesMatrix)[1] <- "SimulationNumber"
    
    # Overwrite stochastic character matrices with new collapsed format:
    x$StochasticCharacterChanges <- ChangesMatrix
    
    # Return full output:
    return(x)
    
  }
  
  # Collapse stochastic character matrices to single matrix for each character:
  CharacterList <- lapply(CharacterList, CompileChangesIntoSingleMatrix)
  
  # Compile all state changes into a single matrix:
  AllStateChanges <- cbind(matrix(unlist(mapply(rep, 1:ncol(MatrixBlock), lapply(CharacterList, function(x) nrow(x$StochasticCharacterChanges))))), do.call(rbind, lapply(CharacterList, function(x) x$StochasticCharacterChanges)))
  
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
  output <- list(AllStateChanges, CharacterTimes, EdgeLengthsPerBin, TerminalEdgeLengthsPerBin, InternalEdgeLengthsPerBin)
  
  # Add names to output:
  names(output) <- c("AllStateChanges", "CharacterTimes", "EdgeLengthsPerBin", "TerminalEdgeLengthsPerBin", "InternalEdgeLengthsPerBin")
  
  # Return output:
  return(invisible(output))
  
}
