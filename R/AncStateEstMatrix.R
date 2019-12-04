#' Ancestral Character State Estimation
#'
#' @description
#'
#' Given a tree and a cladistic matrix uses likelihood to estimate the ancestral states for every character.
#'
#' @param CladisticMatrix A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @param Tree A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{CladisticMatrix}.
#' @param EstimateAllNodes Logical that allows the user to make estimates for all ancestral values. The default (\code{FALSE}) will only make estimates for nodes that link coded terminals (recommended).
#' @param EstimateTipValues Logical that allows the user to make estimates for tip values. The default (\code{FALSE}) will only makes estimates for internal nodes (recommended).
#' @param InapplicablesAsMissing Logical that decides whether or not to treat inapplicables as missing (TRUE) or not (FALSE, the default and recommended option).
#' @param PolymorphismBehaviour One of either "equalp" or "treatasmissing".
#' @param UncertaintyBehaviour One of either "equalp" or "treatasmissing".
#' @param Threshold The threshold value to use when collapsing marginal likelihoods to discrete state(s).
#'
#' @details
#'
#' Uses either the \link{rerootingMethod} (Yang et al. 1995) as implemented in the \link{phytools} package (discrete characters) or the \link{ace} function in the \link{ape} package (continuous characters) to make ancestral state estimates. For discrete characters these are collapsed to the most likely state (or states, given equal likelihoods or likelihood within a defined threshold value). In the latter case the resulting states are represented as an uncertainty (i.e., states separated by a slash, e.g., 0/1. This is the method used by Brusatte et al. (2014).
#'
#' @return \item{anc.lik.matrix}{A matrix of nodes (hypothetical ancestors; rows) against characters (columns) listing the reconstructed ancestral states.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Thomas Guillerme \email{guillert@@tcd.ie}
#'
#' @references
#'
#' Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. Current Biology, 24, 2386-2392.
#' 
#' Yang, Z., Kumar, S. and Nei, M., 1995. A new method of inference of ancestral nucleotide and amino acid sequences. Genetics, 141, 1641-1650.
#'
#' @examples
#' 
#' # Set random seed:
#' set.seed(4)
#' 
#' # Generate a random tree for the Day data set:
#' tree <- rtree(n = nrow(Day2016$Matrix_1$Matrix))
#' 
#' # Update taxon names to match those in the data matrix:
#' tree$tip.label <- rownames(Day2016$Matrix_1$Matrix)
#' 
#' # Set root time by making youngest taxon extant:
#' tree$root.time <- max(diag(vcv(tree)))
#'
#' # Use Day matrix as cladistic matrix:
#' CladisticMatrix <- Day2016
#'
#' # Prune most characters out to make example run fast:
#' CladisticMatrix <- MatrixPruner(CladisticMatrix,
#'   characters2prune = c(2:3, 5:37))
#'
#' # Estimate ancestral states:
#' AncStateEstMatrix(CladisticMatrix = CladisticMatrix,
#'   Tree = tree)
#' 
#' @export AncStateEstMatrix
AncStateEstMatrix <- function(CladisticMatrix, Tree, EstimateAllNodes = FALSE, EstimateTipValues = FALSE, InapplicablesAsMissing = FALSE, PolymorphismBehaviour = "equalp", UncertaintyBehaviour = "equalp", Threshold = 0.01) {
  
  # How to get tip states for a continuous character?
  # How to deal with step matrices?
  # Change help file to explain interactions between all options, e.g., if doing all chars then polymorphsims used for discrete, midpoint for continuous etc.
  
  # Catch problem with trees with no branch lengths:
  if(is.null(Tree$edge.length)) stop("Tree must have branch lengths.")
  
  # Catch problem with polytomies:
  if(Tree$Nnode < (ape::Ntip(Tree) - 1)) stop("Tree must be fully bifurcating.")
  
  # Catch problem with zero-length branches:
  if(any(Tree$edge.length == 0)) stop("Tree must not have zero-length branches.")
  
  # Check for step matrices and stop and warn if found:
  if(length(CladisticMatrix$Topper$StepMatrices) > 0) stop("Function can not currently deal with step matrices.")
  
  # Check EstimateAllNodes is a logical:
  if(!is.logical(EstimateAllNodes)) stop("EstimateAllNodes must be a logical (TRUE or FALSE).")
  
  # Check EstimateTipValues is a logical:
  if(!is.logical(EstimateTipValues)) stop("EstimateTipValues must be a logical (TRUE or FALSE).")
  
  # Check InapplicablesAsMissing is a logical:
  if(!is.logical(InapplicablesAsMissing)) stop("InapplicablesAsMissing must be a logical (TRUE or FALSE).")
  
  # Collapse matrix to vectors for each character (state and ordering combination):
  collapse.matrix <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], function(x) apply(rbind(x$Matrix, x$Ordering), 2, paste, collapse = ""))))
  
  # Isolate ordering elements:
  ordering <- unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Ordering"))
  
  # Isolate minimum values:
  min.vals <- unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MinVals"))
  
  # Isolate maximum values:
  max.vals <- unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MaxVals"))
  
  # Store raw original matrix:
  RawCladisticMatrix <- CladisticMatrix
  
  # Combine matrix blocks into a single matrix:
  CladisticMatrix <- OriginalMatrix <- do.call(cbind, lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"))
  
  # Find any failed name matches:
  FailedNameMatches <- c(setdiff(rownames(CladisticMatrix), Tree$tip.label), setdiff(Tree$tip.label, rownames(CladisticMatrix)))
  
  # Check there are no failed name matches and stop and report if found:
  if(length(FailedNameMatches) > 0) stop(paste("The following names do not match between the tree and matrix: ", paste(sort(FailedNameMatches), collapse = ", "), ". Check spelling and try again.", sep = ""))
  
  # Check PolymorphismBehaviour is a valid option and stop and report if not:
  if(length(setdiff(PolymorphismBehaviour, c("equalp", "treatasmissing"))) > 0) stop("PolymorphismBehaviour must be one of \"equalp\" or \"treatasmissing\".")
  
  # Check UncertaintyBehaviour is a valid option and stop and report if not:
  if(length(setdiff(UncertaintyBehaviour, c("equalp", "treatasmissing"))) > 0) stop("UncertaintyBehaviour must be one of \"equalp\" or \"treatasmissing\".")
  
  # Check threshold is in correct window of possible values:
  if(Threshold < 0 || Threshold > 0.5) stop("Threshold must be between 0 and 0.5.")
  
  # If treating inapplicables as missing (and there is at least one inapplicable) replace with NA:
  if(InapplicablesAsMissing && length(which(CladisticMatrix == "")) > 0) CladisticMatrix[which(CladisticMatrix == "")] <- NA
  
  # If treating polymorphisms as missing:
  if(PolymorphismBehaviour == "treatasmissing" && length(grep("&", CladisticMatrix)) > 0) CladisticMatrix[grep("&", CladisticMatrix)] <- NA
  
  # If treating uncertainties as missing:
  if(UncertaintyBehaviour == "treatasmissing" && length(grep("/", CladisticMatrix)) > 0) CladisticMatrix[grep("/", CladisticMatrix)] <- NA
  
  # Convert tip states into a list:
  DataAsList <- apply(CladisticMatrix, 2, list)
  
  # Add Tipsattes name to list:
  DataAsList <- lapply(DataAsList, function(x) {names(x) <- "TipStates"; return(x)})
  
  # For each character:
  for(i in 1:length(DataAsList)) {
    
    # Add minimum value to list:
    DataAsList[[i]]$MinVal <- unname(min.vals[i])
    
    # Add maximum value to list:
    DataAsList[[i]]$MaxVal <- unname(max.vals[i])
    
    # Add ordering to list:
    DataAsList[[i]]$Ordering <- unname(ordering[i])
    
    # Add tree to list:
    DataAsList[[i]]$Tree <- Tree
    
  }
  
  # If estimating values for all characters (need to set dummy tip states for missing values):
  if(EstimateAllNodes) {
    
    # Subfunction to fill missing values (and inapplicables if desired):
    FillMissing <- function(TipStates) {
      
      # Find which rows correspond to missing states:
      MissingRows <- which(is.na(TipStates$TipStates))
      
      # If missing states found:
      if(length(MissingRows) > 0) {
        
        # Build missing state by either forming a polymorphism of all possible tip states, or if continuous the midpoint value:
        FillStates <- ifelse(TipStates$Ordering == "cont", (TipStates$MinVal + TipStates$MaxVal) / 2, paste(TipStates$MinVal:TipStates$MaxVal, collapse = "/"))
        
        # Insert missing values:
        TipStates$TipStates[MissingRows] <- FillStates
        
      }
      
      # Return tip states with missing values replaced:
      return(TipStates)
      
    }
    
    # Apply fill missing fucntion across all characters:
    DataAsList <- lapply(DataAsList, FillMissing)
    
  }
  
  # Subfunction to prune tips with missing or inapplicable values:
  PruneTips <- function(x) {
    
    # Find all missing or inapplicable value tip names:
    Missing <- names(sort(c(which(x$TipStates == ""), which(is.na(x$TipStates)))))
    
    # If there is at least one:
    if(length(Missing) > 0) {
      
      # Remove tips from tree:
      x$Tree <- drop.tip(phy = x$Tree, tip = Missing)
      
      # Remove tips from tip states:
      x$TipStates <- x$TipStates[setdiff(names(x$TipStates), Missing)]
      
    }
    
    # Return pruend output:
    return(x)
    
  }
  
  # Prune out missing and inapplicable tips:
  DataAsList <- lapply(DataAsList, PruneTips)
  
  # Subfunction to build tip state matrices:
  TipStateVectorToMatrix <- function(x) {
    
    # If the character is not continuous (i.e., it is some form of discrete character):
    if(x$Ordering != "cont") {
      
      # Temporarily store tip states so matrix format can overwrite the stored version below:
      TipStates <- x$TipStates
      
      # Create marix of tip state probabilities:
      x$TipStates <- matrix(0, nrow = length(x$TipStates), ncol = x$MaxVal - x$MinVal + 1, dimnames = list(names(x$TipStates), x$MinVal:x$MaxVal))
      
      # For each character state if a single state is coded store probability as 1:
      for(i in colnames(x$TipStates)) x$TipStates[TipStates == i, i] <- 1
      
      # If there are polymorphisms and/or uncertainties:
      if(length(grep("&|/", TipStates)) > 0) {
        
        # Get polymorphism locations:
        Polymorphisms <- grep("&", TipStates)
        
        # Get uncertainty locations:
        Uncertainties <- grep("/", TipStates)
        
        # If there are polymorphisms and using the "equalp" (equal probability of each state) option:
        if(length(Polymorphisms) > 0 && PolymorphismBehaviour == "equalp") {
          
          # For each polymorphisms set each state as equally probable:
          for(i in Polymorphisms) x$TipStates[i, strsplit(TipStates[i], split = "&")[[1]]] <- 1 / length(strsplit(TipStates[i], split = "&")[[1]])
          
        }
        
        # If there are uncertainties and using the "equalp" (equal probability of each state) option:
        if(length(Uncertainties) > 0 && UncertaintyBehaviour == "equalp") {
          
          # For each uncertainty set each state as equally probable:
          for(i in Uncertainties) x$TipStates[i, strsplit(TipStates[i], split = "/")[[1]]] <- 1 / length(strsplit(TipStates[i], split = "/")[[1]])
          
        }
        
      }
      
    # If a continuous character:
    } else {
      
      # Simply make tip states the numeric values (should never be a polymorphism) as a vector:
      x$TipStates <- as.numeric(x$TipStates)
      
    }
    
    # Return the revised input in the same list format:
    return(x)
    
  }
  
  # Reformat tip states ready for ancestral estimation:
  DataAsList <- lapply(DataAsList, TipStateVectorToMatrix)
  
  # Subfunction to build tip state matrices:
  ModelBuilder <- function(x) {
    
    # Set default model to equal rates (works for all binary or unordered characters):
    x$Model <- "ER"
    
    # If a character is both ordered and has at least three states:
    if((x$MaxVal - x$MinVal) > 1 && x$Ordering == "ord") {
      
      # Get number of states:
      NStates <- (x$MaxVal - x$MinVal) + 1
      
      # Build all zero matrix to begin with:
      x$Model <- matrix(0, nrow = NStates, ncol = NStates, dimnames = list(x$MinVal:x$MaxVal, x$MinVal:x$MaxVal))
      
      # for each (just) off-diagonal value store 1 (i.e., N steps to move between adjacent states):
      for(i in 2:NStates) x$Model[(i - 1), i] <- x$Model[i, (i - 1)] <- 1
      
    }
    
    # Return full output:
    return(x)
  
  }
  
  # Add ancestral state model for each character:
  DataAsList <- lapply(DataAsList, ModelBuilder)

  # Sunfunctio to get ancestral states:
  GetAncStates <- function(x, EstimateTipValues, Threshold) {
    
    # If character is continuous:
    if(x$Ordering == "cont") {
      
      # Get ancestral states using ace:
      x$AncestralStates <- ace(x = x$TipStates, phy = x$Tree)$ace
      
    # If character is discrete:
    } else {
      
      # Get ancestral states using rerooting method:
      x$AncestralStates <- rerootingMethod(tree = x$Tree, x = x$TipStates, model = x$Model)$marginal.anc
      
      # Reformat to most likely state
      x$AncestralStates <- unlist(lapply(lapply(apply(x$AncestralStates, 1, list), unlist), function(x) {paste(names(x[x > (max(x) - Threshold)]), collapse = "/")}))
      
      # If not estimating tip values then prune these:
      if(!EstimateTipValues) x$AncestralStates <- x$AncestralStates[-match(x$Tree$tip.label, names(x$AncestralStates))]
      
    }
    
    # Return full output of x:
    return(x)
    
  }
  
  # Get ancestral states for each character:
  DataAsList <- lapply(DataAsList, GetAncStates, EstimateTipValues, Threshold)
  
  # Get Newick strings of all sampled subtrees (to use to avoid reudunancy in tree node mapping):
  NewickStrings <- unlist(lapply(lapply(DataAsList, '[[', "Tree"), ape::write.tree))
  
  # Get just unique strings (i.e., just those trees that need toa ctyal map to full tree):
  UniqueNewickStrings <- unique(NewickStrings)
  
  # Convert unique Newcik strings to unique trees:
  UniqueTrees <- read.tree(text = UniqueNewickStrings)
  
  # If only a single tree reformat as a list:
  if(inherits(UniqueTrees, what = "phylo")){
    UniqueTrees <- list(UniqueTrees)
    }
  
  # Subfunction map nodes from pruned tree to full tree:
  MapPrunedTreeNodesToFullTreeNodes <- function(tree, fulltree) {
    
    # Get number of tips of pruned tree:
    NTips <- ape::Ntip(tree)
    
    # Get number of nodes of pruend tree:
    NNodes <- ape::Nnode(tree)
    
    # Get all internal node numbers for peruend tree:
    NodeNumbers <- (NTips + 1):(NTips + NNodes)
    
    # If the pruned tree is different to the full tree:
    if(write.tree(tree) != write.tree(fulltree)) {
      
      # Get descendants of each node in pruned tree:
      Descendants <- lapply(as.list(NodeNumbers), function(x) tree$tip.label[strap::FindDescendants(x, tree = tree)])
      
      # Get corresponding ancestral node in full tree:
      Ancestors <- unlist(lapply(Descendants, function(x) Claddis::FindAncestor(descs = x, tree = fulltree)))
      
    # If pruned tree is identical to full tree (not pruned at all):
    } else {
      
      # Set ancestors as node numbers:
      Ancestors <- NodeNumbers
      
    }
    
    # Output matrix matching node numbers of pruned tree to full tree:
    return(matrix(c(NodeNumbers, Ancestors), ncol = 2, dimnames = list(c(), c("PrunedNode", "FullNode"))))
    
  }
  
  # Get pruned node to full node for each unique tree:
  NodeMapsList <- lapply(UniqueTrees, MapPrunedTreeNodesToFullTreeNodes, fulltree = Tree)
  
  # Build out for all trees (adds in any duplicated trees):
  NodeMapsList <- NodeMapsList[match(NewickStrings, UniqueNewickStrings)]
  
  # Add node maps to data list:
  for(i in 1:length(DataAsList)) DataAsList[[i]]$NodeMaps <- NodeMapsList[[i]]
  
  # Get number of tips in tree:
  NTips <- ape::Ntip(Tree)
  
  # Get number of ndoes in tree:
  NNodes <- ape::Nnode(Tree)
  
  # Get all node names and numbers:
  Nodes <- c(rownames(OriginalMatrix), (NTips + 1):(NTips + NNodes))
  
  # Renumber nodes of ancestral states:
  DataAsList <- lapply(DataAsList, function(x) { names(x$AncestralStates)[match(as.character(x$NodeMaps[, "PrunedNode"]), names(x$AncestralStates))] <- as.character(x$NodeMaps[, "FullNode"]); return(x) })
  
  # Collapse down to an ancestral state matrix ready for output:
  AncestralStateMatrix <- do.call(cbind, lapply(DataAsList, function(x) {x$AncestralStates <- x$AncestralStates[Nodes]; names(x$AncestralStates) <- Nodes; return(x$AncestralStates)}))
  
  # Isolate estimated tip values:
  TipMatrix <- AncestralStateMatrix[rownames(OriginalMatrix), ]
  
  # If there are any missing values:
  if(any(is.na(TipMatrix))) {
    
    # Isolate misisng values:
    MissingTipStates <- which(is.na(TipMatrix))
    
    # Replace missing values with original (unmodified) input values:
    TipMatrix[MissingTipStates] <- OriginalMatrix[MissingTipStates]
    
    # Add tip values back into full output:
    AncestralStateMatrix[rownames(OriginalMatrix), ] <- TipMatrix
    
  }
  
  # Get column (character) count for each matrix block:
  MatrixColumns <- unlist(lapply(lapply(RawCladisticMatrix[2:length(RawCladisticMatrix)], '[[', "Matrix"), ncol))
  
  # For each matrix block:
  for(i in 1:length(MatrixColumns)) {
    
    # Insert portion of ancestral state estimate into block:
    RawCladisticMatrix[[(i + 1)]]$Matrix <- AncestralStateMatrix[, 1:MatrixColumns[i], drop = FALSE]
    
    # Remove that portion from the block:
    AncestralStateMatrix <- AncestralStateMatrix[, -(1:MatrixColumns[i]), drop = FALSE]
    
  }
  
  # Overwrite ancestral state output with updated raw input:
  AncestralStateMatrix <- RawCladisticMatrix
  
  # Add tree to output:
  AncestralStateMatrix$Topper$Tree <- Tree
  
  # Return ancestral state matrix:
  return(AncestralStateMatrix)

}
