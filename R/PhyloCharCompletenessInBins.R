#' Phylogenetic character completeness in time-bins
#'
#' @description
#'
#' Given a cladistic matrix, time-scaled tree, and set of time bin boundaries will return the proportional character completeness in each bin.
#'
#' @param CladisticMatrix A cladistic matrix in the form imported by \link{ReadMorphNexus}.
#' @param TimeTree A time-scaled phylogenetic tree containing all the taxa in \code{CladisticMatrix}.
#' @param TimeBins A set of time bin boundaries (oldest to youngest) in millions of years.
#' @param plot An optional choice to plot the results (default is \code{FALSE}).
#' @param CI The confidence interval to be used as a proportion (0 to 1). Default is 0.95 (i.e., 95\%).
#'
#' @details
#'
#' Character completeness metrics have been used as an additional metric for comparing fossil record quality across time, space, and taxa. However, these only usually refer to point samples of fossils in bins, and not our ability to infer information along the branches of a phylogenetic tree.
#'
#' This function returns the proportional phylogenetic character completeness for a set of time bins.
#'
#' @return
#'
#' A list summarising the mean, upper and lower 95% confidence interval, and per character proportional character completeness in each time bin.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a random tree for the Day et al. 2016 data set:
#' Day2016tree <- rtree(nrow(Day2016$Matrix_1$Matrix))
#' Day2016tree$tip.label <- rownames(Day2016$Matrix_1$Matrix)
#' Day2016tree$root.time <- max(diag(vcv(Day2016tree)))
#'
#' # Get proportional phylogenetic character completeness in ten equal-length
#' # time bins:
#' PhyloCharCompletenessInBins(CladisticMatrix = Day2016,
#'   TimeTree = Day2016tree, TimeBins = seq(from =
#'   Day2016tree$root.time, to = Day2016tree$root.time -
#'   max(diag(vcv(Day2016tree))), length.out = 11))
#'
#' @export PhyloCharCompletenessInBins
PhyloCharCompletenessInBins <- function(CladisticMatrix, TimeTree, TimeBins, plot = FALSE, CI = 0.95) {
  
  # Subfunction for getting missing and inapplicable characters:
  MissingAndInapplicables <- function(x) {
    
    # Get all inapplicables:
    x_ina <- apply(apply(x, 2, '==', ""), 2, as.numeric)
    
    # Replace NAs with 0:
    x_ina[is.na(x_ina)] <- 0
    
    # Get missing:
    x_mis <- apply(apply(x, 2, is.na), 2, as.numeric)
    
    # Return numeric matrix (1 for missing or inapplicable):
    return(x_mis + x_ina)
    
  }
  
  # Create vector of time bin mid-points:
  TimeBinMidpoints <- TimeBins[1:(length(TimeBins) - 1)] + abs(diff(TimeBins)) / 2
  
  # Create time bin names:
  TimeBinNames <- paste(round(TimeBins[1:(length(TimeBins) - 1)], 1), round(TimeBins[2:length(TimeBins)], 1), sep = "-")
  
  # Get total number of characters:
  NCharacters <- sum(unlist(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), ncol)))
  
  # Get edge lengths in bins for complete tree (measure of a complete character):
  CompleteEdgesInBins <- EdgeLengthsInBins(TimeTree, TimeBins)$edge.length.in.bin
  
  # Set uo missing values vector (no missing values are set as empty characters (""):
  MissingValues <- rep("", NCharacters)
  
  # If there are missing or inapplicable values collapse row numbers for them with double percentage:
  if(any(unlist(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), is.na))) || any(unlist(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), '==', "")))) MissingValues <- unname(unlist(lapply(lapply(lapply(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), MissingAndInapplicables), apply, 2, '==', 1), apply, 2, which), lapply, paste, collapse = "%%")))
  
  # Set up matrix to store edge lengths in each character bin (columns) per character (rows):
  EdgeLengthsInBinsByCharacter <- matrix(0, ncol = length(TimeBins) - 1, nrow = NCharacters)
  
  # For each unique missing value combination (no point repeating those that are the same):
  for(i in unique(MissingValues)) {
    
    # If there is missing data:
    if(nchar(i) > 0) {
      
      # List taxa to prune:
      TaxaToPrune <- rownames(CladisticMatrix$Matrix_1$Matrix)[as.numeric(strsplit(i, "%%")[[1]])]
      
      # Check that there are still enough taxa left for a tree to exist:
      if(length(setdiff(TimeTree$tip.label, TaxaToPrune)) > 1) {
        
        # Remove tips with missing data from tree:
        PrunedTree <- drop.tip(TimeTree, TaxaToPrune)
        
        # Need to correct root time to make sure time binning makes sense:
        PrunedTree <- CorrectRootTime(TimeTree, PrunedTree)
        
      # If there is one or fewer taxa:
      } else {
        
        # Set pruned tree as NA:
        PrunedTree <- NA
        
      }
      
    # If there is not missing data:
    } else {
      
      # Set complete tree as pruned tree:
      PrunedTree <- TimeTree
      
    }
    
    # As long as the tree exists (i.e., it is not pruend down to one or zero taxa) store edge lengths in bin:
    if(!is.na(PrunedTree)[1]) EdgeLengthsInBinsByCharacter[which(MissingValues == i), ] <- matrix(rep(EdgeLengthsInBins(PrunedTree, TimeBins)$edge.length.in.bin, length(which(MissingValues == i))), ncol = ncol(EdgeLengthsInBinsByCharacter), byrow = TRUE)
    
  }
  
  # Calculate and store proportional character completeness:
  ProportionalCompletenessInBinsByCharacter <- EdgeLengthsInBinsByCharacter / matrix(rep(CompleteEdgesInBins, NCharacters), nrow = NCharacters, byrow = TRUE)
  
  # Calculate mean proportional character completeness:
  MeanProportionalCompletenessInBins <- apply(ProportionalCompletenessInBinsByCharacter, 2, mean)
  
  # Calculate upper 95% CI proportional character completeness:
  Upper95ProportionalCompletenessInBins <- apply(ProportionalCompletenessInBinsByCharacter, 2, sort)[ceiling((CI + ((1 - CI) / 2)) * NCharacters), ]
  
  # Calculate lower 95% CI proportional character completeness:
  Lower95ProportionalCompletenessInBins <- apply(ProportionalCompletenessInBinsByCharacter, 2, sort)[max(c(1, floor(((1 - CI) / 2) * NCharacters))), ]
  
  # If plotting output:
  if(plot) {
    
    # Set plot environment for two plots (one on top of the other):
    par(mfrow = c(2, 1))
    
    # Create empty plot first (y-axis limits set from 0 to 1):
    plot(x = TimeBinMidpoints, y = apply(ProportionalCompletenessInBinsByCharacter, 2, mean), type = "n", ylim = c(0, 1), xlim = c(max(TimeBinMidpoints), min(TimeBinMidpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness")
    
    # Add 95% CI as shaded polygon:
    polygon(x = c(TimeBinMidpoints, rev(TimeBinMidpoints)), y = c(Upper95ProportionalCompletenessInBins, rev(Lower95ProportionalCompletenessInBins)), col = "grey", border = 0)
    
    # Plot mean character completeness on top:
    points(x = TimeBinMidpoints, y = MeanProportionalCompletenessInBins, type = "l", ylim = c(0, 1), xlim = c(max(TimeBinMidpoints), min(TimeBinMidpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness", lwd = 2)
    
    # Add legend:
    graphics::legend(x = max(TimeBinMidpoints), y = 0.4, c(paste(CI, "% CI", sep = ""), "Mean"), col = c("grey", "black"), lwd = c(8, 2), merge = TRUE, bg = "white")
    
    # Create empty plot first (y-axis limits set from 0 to 1):
    plot(x = TimeBinMidpoints, y = apply(ProportionalCompletenessInBinsByCharacter, 2, mean), type = "n", ylim = c(0, 1), xlim = c(max(TimeBinMidpoints), min(TimeBinMidpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness")
    
    # Plot each individual characters proportional completeness:
    for(i in 1:NCharacters) points(x = TimeBinMidpoints, y = ProportionalCompletenessInBinsByCharacter[i, ], type = "l", col = "grey")
    
    # Plot mean character completeness on top:
    points(x = TimeBinMidpoints, y = MeanProportionalCompletenessInBins, type = "l", ylim = c(0, 1), xlim = c(max(TimeBinMidpoints), min(TimeBinMidpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness", lwd = 2)
    
    # Add legend:
    graphics::legend(x = max(TimeBinMidpoints), y = 0.4, c("Individual characters", "Mean"), col = c("grey", "black"), lwd = c(1, 2), merge = TRUE, bg = "white")
    
    # Reset plotting environment:
    par(mfrow = c(1, 1))
    
  }
  
  # Add time bin names to output:
  names(MeanProportionalCompletenessInBins) <- names(Upper95ProportionalCompletenessInBins) <- names(Lower95ProportionalCompletenessInBins) <- colnames(ProportionalCompletenessInBinsByCharacter) <- TimeBinNames
  
  # Compile output variables:
  output <- list(MeanProportionalCompletenessInBins, Upper95ProportionalCompletenessInBins, Lower95ProportionalCompletenessInBins, ProportionalCompletenessInBinsByCharacter)
  
  # Add names to output:
  names(output) <- c("MeanProportionalCharacterCompletenessPerTimeBin", "Upper95PercentCIProportionalCharacterCompletenessPerTimeBin", "Lower95PercentCIProportionalCharacterCompletenessPerTimeBin", "ProportionalCharacterCompletenessPerTimeBinByCharacter")
  
  # Return output:
  return(output)
  
}
