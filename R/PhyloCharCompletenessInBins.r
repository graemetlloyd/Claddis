#' Phylogenetic character completeness in time-bins
#'
#' Given a cladistic matrix, time-scaled tree, and set of time bin boundaries will return the proportional character completeness in each bin.
#'
#' Character completeness metrics have been used as an additional metric for comparing fossil record quality across time, space, and taxa. However, these only usually refer to point samples of fossils in bins, and not our ability to infer information along the branches of a phylogenetic tree.
#'
#' This function returns the proportional phylogenetic character completeness for a set of time bins.
#'
#' @param CladMatrix  A cladistic matrix in the form imported by \link{ReadMorphNexus}.
#' @param TimeTree  A time-scaled phylogenetic tree containing all the taxa in \code{CladMatrix}.
#' @param TimeBins  A set of time bin boundaries (oldest to youngest) in millions of years.
#' @param plot An optional choice to plot the results (default is \code{FALSE}).
#'
#' @return
#'
#' A list summarising the mean, upper and lower 95% conifidence interval, and per character proportional character completeness in each time bin.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a random tree for the Gauthier 1986 data set:
#' Gauthier1986tree <- rtree(nrow(Gauthier1986$matrix))
#' Gauthier1986tree$tip.label <- rownames(Gauthier1986$matrix)
#' Gauthier1986tree$root.time <- max(diag(vcv(Gauthier1986tree)))
#'
#' # Get proportional phylogenetic character completeness in ten equal-length
#' # time bins:
#' PhyloCharCompletenessInBins(CladMatrix = Gauthier1986,
#'   TimeTree = Gauthier1986tree, TimeBins = seq(from =
#'   Gauthier1986tree$root.time, to = Gauthier1986tree$root.time -
#'   max(diag(vcv(Gauthier1986tree))), length.out = 11))
#'
#' @export PhyloCharCompletenessInBins
PhyloCharCompletenessInBins <- function(CladMatrix, TimeTree, TimeBins, plot = FALSE) {
    
    # Create vector of time bin mid-points:
    TimeBinMidpoints <- TimeBins[1:(length(TimeBins) - 1)] + abs(diff(TimeBins)) / 2
    
    # Create time bin names:
    TimeBinNames <- paste(round(TimeBins[1:(length(TimeBins) - 1)], 1), round(TimeBins[2:length(TimeBins)], 1), sep = "-")
    
    # Get total number of characters:
    NCharacters <- ncol(CladMatrix$matrix)
    
    # Get edge lengths in bins for complete tree (measure of a complete character):
    CompleteEdgesInBins <- EdgeLengthsInBins(TimeTree, TimeBins)$edge.length.in.bin
    
    # Set uo missing values vector (no missing values are set as empty characters (""):
    MissingValues <- rep("", ncol(CladMatrix$matrix))
    
    # If there are missing values collapse row numbers for them with double percentage:
    if(any(is.na(CladMatrix$matrix))) MissingValues <- unlist(lapply(apply(apply(CladMatrix$matrix, 2, is.na), 2, which), paste, collapse = "%%"))
    
    # Set up matrix to store edge lengths in each character bin (columns) per character (rows):
    EdgeLengthsInBinsByCharacter <- matrix(0, ncol = length(TimeBins) - 1, nrow = NCharacters)
    
    # For each unique missing value combination (no point repeating those that are the same):
    for(i in unique(MissingValues)) {
        
        # If there is missing data:
        if(nchar(i) > 0) {
            
            # List taxa to prune:
            TaxaToPrune <- rownames(CladMatrix$matrix)[as.numeric(strsplit(i, "%%")[[1]])]
            
            # Check
            
            if(length(setdiff(TimeTree$tip.label, TaxaToPrune)) > 1) {
            
                PrunedTree <- drop.tip(TimeTree, TaxaToPrune)
            
                PrunedTree <- CorrectRootTime(TimeTree, PrunedTree)
            
            }
            
        # If there is not missing data:
        } else {
            
            # Set compete tree as pruned tree:
            PrunedTree <- TimeTree
            
        }
        
        # Store data for all characters with same configuration (set of taxa) of missing tip data:
        for(j in which(MissingValues == i)) EdgeLengthsInBinsByCharacter[j, ] <- EdgeLengthsInBins(PrunedTree, TimeBins)$edge.length.in.bin
        
    }
    
    # Calculate and store proportional character completeness:
    ProportionalCompletenessInBinsByCharacter <- EdgeLengthsInBinsByCharacter / matrix(rep(CompleteEdgesInBins, NCharacters), nrow = NCharacters, byrow = TRUE)
    
    # Calculate mean proportional character completeness:
    MeanProportionalCompletenessInBins <- apply(ProportionalCompletenessInBinsByCharacter, 2, mean)
    
    # Calculate upper 95% CI proportional character completeness:
    Upper95ProportionalCompletenessInBins <- apply(ProportionalCompletenessInBinsByCharacter, 2, sort)[ceiling(0.975 * NCharacters), ]
    
    # Calculate lower 95% CI proportional character completeness:
    Lower95ProportionalCompletenessInBins <- apply(ProportionalCompletenessInBinsByCharacter, 2, sort)[max(c(1, floor(0.025 * NCharacters))), ]

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
        legend(x = max(TimeBinMidpoints), y = 0.4, c("95% CI", "Mean"), col = c("grey", "black"), lwd = c(8, 2), merge = TRUE, bg = "white")

        # Create empty plot first (y-axis limits set from 0 to 1):
        plot(x = TimeBinMidpoints, y = apply(ProportionalCompletenessInBinsByCharacter, 2, mean), type = "n", ylim = c(0, 1), xlim = c(max(TimeBinMidpoints), min(TimeBinMidpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness")
        
        # Plot each individual characters proportional completeness:
        for(i in 1:NCharacters) points(x = TimeBinMidpoints, y = ProportionalCompletenessInBinsByCharacter[i, ], type = "l", col = "grey")
        
        # Plot mean character completeness on top:
        points(x = TimeBinMidpoints, y = MeanProportionalCompletenessInBins, type = "l", ylim = c(0, 1), xlim = c(max(TimeBinMidpoints), min(TimeBinMidpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness", lwd = 2)
        
        # Add legend:
        legend(x = max(TimeBinMidpoints), y = 0.4, c("Individual characters", "Mean"), col = c("grey", "black"), lwd = c(1, 2), merge = TRUE, bg = "white")
        
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
