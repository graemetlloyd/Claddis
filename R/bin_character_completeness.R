#' Phylogenetic character completeness in time-bins
#'
#' @description
#'
#' Given a cladistic matrix, time-scaled tree, and set of time bin boundaries will return the proportional character completeness in each bin.
#'
#' @param cladistic_matrix A cladistic matrix in the form imported by \link{read_nexus_matrix}.
#' @param time_tree A time-scaled phylogenetic tree containing all the taxa in \code{cladistic_matrix}.
#' @param time_bins A set of time bin boundaries (oldest to youngest) in millions of years.
#' @param plot An optional choice to plot the results (default is \code{FALSE}).
#' @param confidence.interval The confidence interval to be used as a proportion (0 to 1). Default is 0.95 (i.e., 95\%).
#'
#' @details
#'
#' Character completeness metrics have been used as an additional metric for comparing fossil record quality across time, space, and taxa. However, these only usually refer to point samples of fossils in bins, and not our ability to infer information along the branches of a phylogenetic tree.
#'
#' This function returns the proportional phylogenetic character completeness for a set of time bins.
#'
#' @return
#'
#' A list summarising the mean, upper and lower confidence interval, and per character proportional character completeness in each time bin.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a random tree for the Day et al. 2016 data set:
#' day_2016tree <- ape::rtree(n = nrow(day_2016$matrix_1$matrix))
#' day_2016tree$tip.label <- rownames(x = day_2016$matrix_1$matrix)
#' day_2016tree$root.time <- max(diag(x = ape::vcv(phy = day_2016tree)))
#'
#' # Get proportional phylogenetic character completeness in ten equal-length
#' # time bins:
#' bin_character_completeness(
#'   cladistic_matrix = day_2016,
#'   time_tree = day_2016tree, time_bins = seq(
#'     from =
#'       day_2016tree$root.time, to = day_2016tree$root.time -
#'       max(diag(x = ape::vcv(phy = day_2016tree))), length.out = 11
#'   )
#' )
#' @export bin_character_completeness
bin_character_completeness <- function(cladistic_matrix, time_tree, time_bins, plot = FALSE, confidence.interval = 0.95) {

  # TO DO:
  #
  # - Add an aphylogenetic option
  # - Allow binning in other ways. E.g., by character grouping.

  # Subfunction for getting missing and inapplicable characters:
  find_missing_and_inapplicable <- function(x) {

    # Get all inapplicables:
    x.inapplicables <- apply(apply(x, 2, "==", ""), 2, as.numeric)

    # Replace NAs with 0:
    x.inapplicables[is.na(x.inapplicables)] <- 0

    # Get missing:
    x.missing <- apply(apply(x, 2, is.na), 2, as.numeric)

    # Return numeric matrix (1 for missing or inapplicable):
    return(x.missing + x.inapplicables)
  }

  # Create vector of time bin mid-points:
  time.bin.midpoints <- time_bins[1:(length(x = time_bins) - 1)] + abs(diff(time_bins)) / 2

  # Create time bin names:
  time.bin.names <- paste(round(time_bins[1:(length(x = time_bins) - 1)], 1), round(time_bins[2:length(x = time_bins)], 1), sep = "-")

  # Get total number of characters:
  n.characters <- sum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol)))

  # Get edge lengths in bins for complete tree (measure of a complete character):
  complete.edges.in.bins <- bin_edge_lengths(time_tree = time_tree, time_bins = time_bins)$binned_edge_lengths

  # Set uo missing values vector (no missing values are set as empty characters (""):
  missing.values <- rep("", n.characters)

  # If there are missing or inapplicable values collapse row numbers for them with double percentage:
  if (any(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), is.na))) || any(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), "==", "")))) missing.values <- unname(unlist(x = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), find_missing_and_inapplicable), apply, 2, "==", 1), apply, 2, which), lapply, paste, collapse = "%%")))

  # Set up matrix to store edge lengths in each character bin (columns) per character (rows):
  edge.lengths.in.bins.by.character <- matrix(0, ncol = length(x = time_bins) - 1, nrow = n.characters)

  # For each unique missing value combination (no point repeating those that are the same):
  for (i in unique(x = missing.values)) {

    # If there is missing data:
    if (nchar(i) > 0) {

      # List taxa to prune:
      taxa.to.prune <- rownames(x = cladistic_matrix$matrix_1$matrix)[as.numeric(strsplit(i, "%%")[[1]])]

      # Check that there are still enough taxa left for a tree to exist:
      if (length(x = setdiff(x = time_tree$tip.label, y = taxa.to.prune)) > 1) {

        # Remove tips with missing data from tree:
        pruned_tree <- ape::drop.tip(phy = time_tree, tip = taxa.to.prune)

        # Need to correct root time to make sure time binning makes sense:
        pruned_tree <- fix_root_time(time_tree, pruned_tree)

        # If there is one or fewer taxa:
      } else {

        # Set pruned tree as NA:
        pruned_tree <- NA
      }

      # If there is not missing data:
    } else {

      # Set complete tree as pruned tree:
      pruned_tree <- time_tree
    }

    # As long as the tree exists (i.e., it is not pruend down to one or zero taxa) store edge lengths in bin:
    if (!is.na(pruned_tree)[1]) edge.lengths.in.bins.by.character[which(x = missing.values == i), ] <- matrix(rep(bin_edge_lengths(time_tree = pruned_tree, time_bins = time_bins)$binned_edge_lengths, length(x = which(x = missing.values == i))), ncol = ncol(edge.lengths.in.bins.by.character), byrow = TRUE)
  }

  # Calculate and store proportional character completeness:
  proportional.completeness.in.bins.by.character <- edge.lengths.in.bins.by.character / matrix(rep(complete.edges.in.bins, n.characters), nrow = n.characters, byrow = TRUE)

  # Calculate mean proportional character completeness:
  mean.proportional.completeness.in.bins <- apply(proportional.completeness.in.bins.by.character, 2, mean)

  # Calculate upper 95% confidence interval proportional character completeness:
  upper.proportional.completeness.in.bins <- apply(proportional.completeness.in.bins.by.character, 2, sort)[ceiling((confidence.interval + ((1 - confidence.interval) / 2)) * n.characters), ]

  # Calculate lower 95% confidence interval proportional character completeness:
  lower.proportional.completeness.in.bins <- apply(proportional.completeness.in.bins.by.character, 2, sort)[max(c(1, floor(((1 - confidence.interval) / 2) * n.characters))), ]

  # If plotting output:
  if (plot) {

    # Set plot environment for two plots (one on top of the other):
    par(mfrow = c(2, 1))

    # Create empty plot first (y-axis limits set from 0 to 1):
    plot(x = time.bin.midpoints, y = apply(proportional.completeness.in.bins.by.character, 2, mean), type = "n", ylim = c(0, 1), xlim = c(max(time.bin.midpoints), min(time.bin.midpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness")

    # Add 95% confidence interval as shaded polygon:
    polygon(x = c(time.bin.midpoints, rev(time.bin.midpoints)), y = c(upper.proportional.completeness.in.bins, rev(lower.proportional.completeness.in.bins)), col = "grey", border = 0)

    # Plot mean character completeness on top:
    points(x = time.bin.midpoints, y = mean.proportional.completeness.in.bins, type = "l", ylim = c(0, 1), xlim = c(max(time.bin.midpoints), min(time.bin.midpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness", lwd = 2)

    # Add legend:
    graphics::legend(x = max(time.bin.midpoints), y = 0.4, c(paste(confidence.interval, "% confidence interval", sep = ""), "Mean"), col = c("grey", "black"), lwd = c(8, 2), merge = TRUE, bg = "white")

    # Create empty plot first (y-axis limits set from 0 to 1):
    plot(x = time.bin.midpoints, y = apply(proportional.completeness.in.bins.by.character, 2, mean), type = "n", ylim = c(0, 1), xlim = c(max(time.bin.midpoints), min(time.bin.midpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness")

    # Plot each individual characters proportional completeness:
    for (i in 1:n.characters) points(x = time.bin.midpoints, y = proportional.completeness.in.bins.by.character[i, ], type = "l", col = "grey")

    # Plot mean character completeness on top:
    points(x = time.bin.midpoints, y = mean.proportional.completeness.in.bins, type = "l", ylim = c(0, 1), xlim = c(max(time.bin.midpoints), min(time.bin.midpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness", lwd = 2)

    # Add legend:
    graphics::legend(x = max(time.bin.midpoints), y = 0.4, c("Individual characters", "Mean"), col = c("grey", "black"), lwd = c(1, 2), merge = TRUE, bg = "white")

    # Reset plotting environment:
    par(mfrow = c(1, 1))
  }

  # Add time bin names to output:
  names(mean.proportional.completeness.in.bins) <- names(upper.proportional.completeness.in.bins) <- names(lower.proportional.completeness.in.bins) <- colnames(x = proportional.completeness.in.bins.by.character) <- time.bin.names

  # Compile output variables:
  output <- list(mean.proportional.character.completeness.per.time.bin = mean.proportional.completeness.in.bins, upper.confidence.interval.proportional.character.completeness.per.time.bin = upper.proportional.completeness.in.bins, lower.confidence.intervalproportional.character.completeness.per.time.bin = lower.proportional.completeness.in.bins, proportional.character.completeness.per.time.binByCharacter = proportional.completeness.in.bins.by.character)

  # Return output:
  return(output)
}
