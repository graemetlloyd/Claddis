#' Ancestral Character State Estimation
#'
#' @description
#'
#' Given a tree and a cladistic matrix uses likelihood to estimate the ancestral states for every character.
#'
#' @param cladistic.matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @param time.tree A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{cladistic.matrix}.
#' @param estimate.all.nodes Logical that allows the user to make estimates for all ancestral values. The default (\code{FALSE}) will only make estimates for nodes that link coded terminals (recommended).
#' @param estimate.tip.values Logical that allows the user to make estimates for tip values. The default (\code{FALSE}) will only makes estimates for internal nodes (recommended).
#' @param inapplicables.as.missing Logical that decides whether or not to treat inapplicables as missing (TRUE) or not (FALSE, the default and recommended option).
#' @param polymorphism.behaviour One of either "equalp" or "treatasmissing".
#' @param uncertainty.behaviour One of either "equalp" or "treatasmissing".
#' @param threshold The threshold value to use when collapsing marginal likelihoods to discrete state(s).
#' @param allow.all.missing Logical to allow all missing character values (generally not recommended, hence default is FALSE).
#'
#' @details
#'
#' At its' core the function uses either the \link{rerootingMethod} (Yang et al. 1995) as implemented in the \link{phytools} package (for discrete characters) or the \link{ace} function in the \link{ape} package (for continuous characters) to make ancestral state estimates. For discrete characters these are collapsed to the most likely state (or states, given equal likelihoods or likelihood within a defined threshold value). In the latter case the resulting states are represented as an uncertainty (i.e., states separated by a slash, e.g., 0/1. This is the method used by Brusatte et al. (2014).
#'
#' The method can deal with ordered or unordered characters and does so by allowing only indirect transitions (from 0 to 2 must pass through 1) or direct transitions (from 0 straight to 2), respectively. However, more complex step matrix transitions are not currently supported.
#'
#' Ancestral state estimation is complicated where polymorphic or
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
#' time.tree <- rtree(n = nrow(Day2016$Matrix_1$Matrix))
#' 
#' # Update taxon names to match those in the data matrix:
#' time.tree$tip.label <- rownames(Day2016$Matrix_1$Matrix)
#' 
#' # Set root time by making youngest taxon extant:
#' time.tree$root.time <- max(diag(vcv(time.tree)))
#'
#' # Use Day matrix as cladistic matrix:
#' cladistic.matrix <- Day2016
#'
#' # Prune most characters out to make example run fast:
#' cladistic.matrix <- prune_cladistic_matrix(cladistic.matrix,
#'   characters2prune = c(2:3, 5:37))
#'
#' # Estimate ancestral states:
#' estimate_ancestral_states(cladistic.matrix = cladistic.matrix,
#'   time.tree = time.tree)
#' 
#' @export estimate_ancestral_states
estimate_ancestral_states <- function(cladistic.matrix, time.tree, estimate.all.nodes = FALSE, estimate.tip.values = FALSE, inapplicables.as.missing = FALSE, polymorphism.behaviour = "equalp", uncertainty.behaviour = "equalp", threshold = 0.01, allow.all.missing = FALSE) {
  
  # How to get tip states for a continuous character? (Phytools answer: http://blog.phytools.org/2013/11/reconstructed-ancestral-tip-states-for.html)
  #   - So basically under ML just inherit state from ancestral node (really this is mean of distribution where sd would grow with duration of branch so to allow the possibility of variance this could also be sampled stochastically
  # How to deal with step matrices?
  # How to deal with models where intermediate tip states are not even in sample
  # Change help file to explain interactions between all options, e.g., if doing all chars then polymorphisms used for discrete, midpoint for continuous etc.
  # Handle all missing/inapplicable case properly
  # Handle only two tips case properly
  # Add Liam Revell polymorphism options
  # Finish trailing off help file!
  # Add Lloyd citation
  
  # Catch problem with trees with no branch lengths:
  if(is.null(time.tree$edge.length)) stop("time.tree must have branch lengths.")
  
  # Catch problem with polytomies:
  if(time.tree$Nnode < (ape::Ntip(time.tree) - 1)) stop("time.tree must be fully bifurcating.")
  
  # Catch problem with zero-length branches:
  if(any(time.tree$edge.length == 0)) stop("time.tree must not have zero-length branches.")
  
  # Check for step matrices and stop and warn if found:
  if(length(cladistic.matrix$Topper$StepMatrices) > 0) stop("Function can not currently deal with step matrices.")
  
  # Check estimate.all.nodes is a logical:
  if(!is.logical(estimate.all.nodes)) stop("estimate.all.nodes must be a logical (TRUE or FALSE).")
  
  # Check estimate.tip.values is a logical:
  if(!is.logical(estimate.tip.values)) stop("estimate.tip.values must be a logical (TRUE or FALSE).")
  
  # Check inapplicables.as.missing is a logical:
  if(!is.logical(inapplicables.as.missing)) stop("inapplicables.as.missing must be a logical (TRUE or FALSE).")
  
  # Check polymorphism.behaviour is a single allowable value:
  if(length(polymorphism.behaviour) != 1 || !any(c("equalp", "treatasmissing") == polymorphism.behaviour)) stop("polymorphism.behaviour must be a single value of either, \"equalp\" or \"treatasmissing\".")
  
  # Check uncertainty.behaviour is a single allowable value:
  if(length(uncertainty.behaviour) != 1 || !any(c("equalp", "treatasmissing") == uncertainty.behaviour)) stop("uncertainty.behaviour must be a single value of either, \"equalp\" or \"treatasmissing\".")
  
  # Check threshold is a numeric value between the limits of zero and one:
  if(!is.numeric(threshold) || threshold > 0.5 || threshold < 0) stop("threshold must be a numeric value between 0 and 0.5.")
  
  # Collapse matrix to vectors for each character (state and ordering combination):
  collapse.matrix <- unname(unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], function(x) apply(rbind(x$Matrix, x$Ordering), 2, paste, collapse = ""))))
  
  # Isolate ordering elements:
  ordering <- unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "Ordering"))
  
  # Isolate minimum values:
  min.vals <- unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "MinVals"))
  
  # Isolate maximum values:
  max.vals <- unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "MaxVals"))
  
  # Store raw original matrix:
  Rawcladistic.matrix <- cladistic.matrix
  
  # Combine matrix blocks into a single matrix:
  cladistic.matrix <- OriginalMatrix <- do.call(cbind, lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "Matrix"))
  
  # Find any failed name matches:
  FailedNameMatches <- c(setdiff(rownames(cladistic.matrix), time.tree$tip.label), setdiff(time.tree$tip.label, rownames(cladistic.matrix)))
  
  # Check there are no failed name matches and stop and report if found:
  if(length(FailedNameMatches) > 0) stop(paste("The following names do not match between the tree and matrix: ", paste(sort(FailedNameMatches), collapse = ", "), ". Check spelling and try again.", sep = ""))
  
  # If treating inapplicables as missing (and there is at least one inapplicable) replace with NA:
  if(inapplicables.as.missing && length(which(cladistic.matrix == "")) > 0) cladistic.matrix[which(cladistic.matrix == "")] <- NA
  
  # If treating polymorphisms as missing:
  if(polymorphism.behaviour == "treatasmissing" && length(grep("&", cladistic.matrix)) > 0) cladistic.matrix[grep("&", cladistic.matrix)] <- NA
  
  # If treating uncertainties as missing:
  if(uncertainty.behaviour == "treatasmissing" && length(grep("/", cladistic.matrix)) > 0) cladistic.matrix[grep("/", cladistic.matrix)] <- NA
  
  # Get vector of character numbers where all values are NA:
  AllMissingCharacters <- which(apply(cladistic.matrix, 2, function(x) all(is.na(x))))
  
  # Look for all missing characters and stop and wanr user if found:
  if(!allow.all.missing && length(AllMissingCharacters) > 0) stop(paste0("The following characters are coded as missing across all tips: ", paste0(AllMissingCharacters, collapse = ", "), ". This can arise either because of the input data (in which case it is recommended that the user prune these characters using Claddis::prune_cladistic_matrix) or because of the chosen options for inapplicables.as.missing, polymorphism.behaviour, and/or uncertainty.behaviour (in which case the user may wish to chose different values for these)."))
  
  # Convert tip states into a list:
  DataAsList <- apply(cladistic.matrix, 2, function(x) list(TipStates = x))
  
  # For each character:
  for(i in 1:length(DataAsList)) {
    
    # Add minimum value to list:
    DataAsList[[i]]$MinVal <- unname(min.vals[i])
    
    # Add maximum value to list:
    DataAsList[[i]]$MaxVal <- unname(max.vals[i])
    
    # Add ordering to list:
    DataAsList[[i]]$Ordering <- unname(ordering[i])
    
    # Add tree to list:
    DataAsList[[i]]$Tree <- time.tree
    
  }
  
  # If estimating values for all characters (need to set dummy tip states for missing values):
  if(estimate.all.nodes) {
    
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
    
    # Apply fill missing function across all characters:
    DataAsList <- lapply(DataAsList, FillMissing)
    
  }
  
  # Subfunction to prune tips with missing or inapplicable values:
  PruneTips <- function(x) {
    
    # Find all missing or inapplicable value tip names:
    Missing <- names(sort(c(which(x$TipStates == ""), which(is.na(x$TipStates)))))
    
    # Work out how many tips will be left after pruning:
    NTipsRemaining <- length(setdiff(names(x$TipStates), Missing))
    
    # If there is at least one missing value:
    if(length(Missing) > 0) {
      
      # If less than two tips will remain then set tree as NULL:
      if(NTipsRemaining < 2) x$Tree <- NULL

      # If at least two tips will remain prune missing values from tree:
      if(NTipsRemaining > 1) x$Tree <- drop.tip(phy = x$Tree, tip = Missing)
      
      # Collapse tip states:
      x$TipStates <- x$TipStates[setdiff(names(x$TipStates), Missing)]

    }

    # Return pruned output:
    return(x)
    
  }
  
  # Prune out missing and inapplicable tips:
  DataAsList <- lapply(DataAsList, PruneTips)
  
  # Subfunction to build tip state matrices:
  TipStateVectorToMatrix <- function(x) {
    
    # As long asthere is at least one tip state:
    if(length(x$TipStates) > 0) {
      
      # If the character is not continuous (i.e., it is some form of discrete character):
      if(x$Ordering != "cont") {
        
        # Temporarily store tip states so matrix format can overwrite the stored version below:
        TipStates <- x$TipStates
        
        # Create matrix of tip state probabilities:
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
          if(length(Polymorphisms) > 0 && polymorphism.behaviour == "equalp") {
            
            # For each polymorphisms set each state as equally probable:
            for(i in Polymorphisms) x$TipStates[i, strsplit(TipStates[i], split = "&")[[1]]] <- 1 / length(strsplit(TipStates[i], split = "&")[[1]])
            
          }
          
          # If there are uncertainties and using the "equalp" (equal probability of each state) option:
          if(length(Uncertainties) > 0 && uncertainty.behaviour == "equalp") {
            
            # For each uncertainty set each state as equally probable:
            for(i in Uncertainties) x$TipStates[i, strsplit(TipStates[i], split = "/")[[1]]] <- 1 / length(strsplit(TipStates[i], split = "/")[[1]])
            
          }
          
        }
        
      # If a continuous character:
      } else {
        
        # Simply make tip states the numeric values (should never be a polymorphism) as a vector:
        x$TipStates <- as.numeric(x$TipStates)
        
      }
      
    # If tip state has no length (all values are missing):
    } else {
      
      # Create row-less tip states matrix:
      x$TipStates <- matrix(nrow = 0, ncol = 1, dimnames = list(c(), "0"))
      
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

  # Subfunction to get ancestral states:
  GetAncStates <- function(x, estimate.tip.values, threshold) {
    
    # As long as there is a tree:
    if(!is.null(x$Tree)) {
      
      # If character is continuous:
      if(x$Ordering == "cont") {
        
        # Get ancestral states using ace:
        x$AncestralStates <- ace(x = x$TipStates, phy = x$Tree)$ace
        
      # If character is discrete:
      } else {
        
        # If invariant character:
        if(ncol(x$TipStates) == 1) {
          
          # Get number of tips:
          NTreeTips <- ape::Ntip(x$Tree)
          
          # Set ancestral states as all the same:
          x$AncestralStates <- matrix(rep(x = 1, times = (NTreeTips + x$Tree$Nnode)), ncol = 1, dimnames = list(c(x$Tree$tip.label, (NTreeTips + 1):(NTreeTips + x$Tree$Nnode)), colnames(x$TipStates)))
          
        }
        
        # If variant character then get ancestral states using rerooting method:
        if(ncol(x$TipStates) > 1) x$AncestralStates <- rerootingMethod(tree = x$Tree, x = x$TipStates, model = x$Model)$marginal.anc
        
        # Reformat to most likely state
        x$AncestralStates <- unlist(lapply(lapply(apply(x$AncestralStates, 1, list), unlist), function(x) {paste(names(x[x > (max(x) - threshold)]), collapse = "/")}))
        
        # If not estimating tip values then prune these:
        if(!estimate.tip.values) x$AncestralStates <- x$AncestralStates[-match(x$Tree$tip.label, names(x$AncestralStates))]
        
      }
     
    # If no tree:
    } else {
      
      # Set ancestral states as NULL:
      x$AncestralStates <- vector(mode = "character")
      
    }
    
    # Return full output of x:
    return(x)
    
  }
  
  # Get ancestral states for each character:
  DataAsList <- lapply(DataAsList, GetAncStates, estimate.tip.values, threshold)
  
  # Get Newick strings of all sampled subtrees (to use to avoid redundancy in tree node mapping):
  NewickStrings <- unlist(lapply(DataAsList, function(x) ifelse(is.null(x$Tree), NA, ape::write.tree(x$Tree))))
  
  # Get just unique strings (i.e., the minimum set needded to map everything to the full tree):
  UniqueNewickStrings <- unique(NewickStrings[!is.na(NewickStrings)])
  
  # Convert unique Newick strings to unique trees:
  UniqueTrees <- read.tree(text = UniqueNewickStrings)
  
  # If only a single tree reformat as a list:
  if(inherits(UniqueTrees, what = "phylo")) UniqueTrees <- list(UniqueTrees)
  
  # Subfunction map nodes from pruned tree to full tree:
  MapPrunedTreeNodesToFullTreeNodes <- function(tree, fulltree) {
    
    # Get number of tips of pruned tree:
    NTips <- ape::Ntip(tree)
    
    # Get number of nodes of pruend tree:
    NNodes <- ape::Nnode(tree)
    
    # Get all internal node numbers for pruned tree:
    NodeNumbers <- (NTips + 1):(NTips + NNodes)
    
    # If the pruned tree is different to the full tree:
    if(write.tree(tree) != write.tree(fulltree)) {
      
      # Get descendants of each node in pruned tree:
      Descendants <- lapply(as.list(NodeNumbers), function(x) tree$tip.label[strap::FindDescendants(x, tree = tree)])
      
      # Get corresponding ancestral node in full tree:
      Ancestors <- unlist(lapply(Descendants, function(x) Claddis::find_mrca(descs = x, tree = fulltree)))
      
    # If pruned tree is identical to full tree (not pruned at all):
    } else {
      
      # Set ancestors as node numbers:
      Ancestors <- NodeNumbers
      
    }
    
    # Output matrix matching node numbers of pruned tree to full tree:
    return(matrix(c(NodeNumbers, Ancestors), ncol = 2, dimnames = list(c(), c("PrunedNode", "FullNode"))))
    
  }
  
  # Get pruned node to full node for each unique tree:
  NodeMapsList <- lapply(UniqueTrees, MapPrunedTreeNodesToFullTreeNodes, fulltree = time.tree)
  
  # Build out for all trees (adds in any duplicated trees):
  NodeMapsList <- NodeMapsList[match(NewickStrings, UniqueNewickStrings)]
  
  # Add node maps to data list:
  for(i in 1:length(DataAsList)) DataAsList[[i]]$NodeMaps <- NodeMapsList[[i]]
  
  # Get number of tips in tree:
  NTips <- ape::Ntip(time.tree)
  
  # Get number of nodes in tree:
  NNodes <- ape::Nnode(time.tree)
  
  # Get all node names and numbers:
  Nodes <- c(rownames(OriginalMatrix), (NTips + 1):(NTips + NNodes))
  
  # Renumber nodes of ancestral states:
  DataAsList <- lapply(DataAsList, function(x) {
    
    # Renumber nodes:
    names(x$AncestralStates)[match(as.character(x$NodeMaps[, "PrunedNode"]), names(x$AncestralStates))] <- as.character(x$NodeMaps[, "FullNode"])
    
    # Return renumbered nodes:
    return(x)
  
  })
  
  # Collapse down to an ancestral state matrix ready for output:
  AncestralStateMatrix <- do.call(cbind, lapply(DataAsList, function(x) {
    
    # Get acestral states for nodes:
    x$AncestralStates <- x$AncestralStates[Nodes]
    
    # Add node names:
    names(x$AncestralStates) <- Nodes
    
    # Return output:
    return(x$AncestralStates)
  
  }))
  
  # Isolate estimated tip values:
  TipMatrix <- AncestralStateMatrix[rownames(OriginalMatrix), ]
  
  # If there are any missing values:
  if(any(is.na(TipMatrix))) {
    
    # Isolate missing values:
    MissingTipStates <- which(is.na(TipMatrix))
    
    # Replace missing values with original (unmodified) input values:
    TipMatrix[MissingTipStates] <- OriginalMatrix[MissingTipStates]
    
    # Add tip values back into full output:
    AncestralStateMatrix[rownames(OriginalMatrix), ] <- TipMatrix
    
  }
  
  # Get column (character) count for each matrix block:
  MatrixColumns <- unlist(lapply(lapply(Rawcladistic.matrix[2:length(Rawcladistic.matrix)], '[[', "Matrix"), ncol))
  
  # For each matrix block:
  for(i in 1:length(MatrixColumns)) {
    
    # Insert portion of ancestral state estimate into block:
    Rawcladistic.matrix[[(i + 1)]]$Matrix <- AncestralStateMatrix[, 1:MatrixColumns[i], drop = FALSE]
    
    # Remove that portion from the block:
    AncestralStateMatrix <- AncestralStateMatrix[, -(1:MatrixColumns[i]), drop = FALSE]
    
  }
  
  # Overwrite ancestral state output with updated raw input:
  AncestralStateMatrix <- Rawcladistic.matrix
  
  # Add tree to output:
  AncestralStateMatrix$Topper$Tree <- time.tree
  
  # Return ancestral state matrix:
  return(AncestralStateMatrix)

}
