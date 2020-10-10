#' Phylogenetic character completeness in time-bins
#'
#' @description
#'
#' Given a cladistic matrix, time-scaled tree, and set of time bin boundaries will return the proportional character completeness in each bin.
#'
#' @param cladistic_matrix A cladistic matrix in the form imported by \link{read_nexus_matrix}.
#' @param time_tree A time-scaled phylogenetic tree containing all the taxa in \code{cladistic_matrix}.
#' @param time_bins An object of class \code{timeBins}.
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
#' # Build ten equal-length time bins spanning the tree:
#' time_bins <- matrix(data = c(seq(from = day_2016tree$root.time,
#'   to = day_2016tree$root.time - max(diag(x = ape::vcv(phy = day_2016tree))),
#'   length.out = 11)[1:10], seq(from = day_2016tree$root.time,
#'   to = day_2016tree$root.time - max(diag(x = ape::vcv(phy = day_2016tree))),
#'   length.out = 11)[2:11]), ncol = 2, dimnames = list(LETTERS[1:10], c("fad", "lad")))
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Get proportional phylogenetic character completeness in ten equal-length
#' # time bins:
#' bin_character_completeness(
#'   cladistic_matrix = day_2016,
#'   time_tree = day_2016tree,
#'   time_bins = time_bins
#' )
#'
#' # Same, but with a plot:
#' bin_character_completeness(
#'   cladistic_matrix = day_2016,
#'   time_tree = day_2016tree,
#'   time_bins = time_bins,
#'   plot = TRUE
#' )
#' @export bin_character_completeness
bin_character_completeness <- function(cladistic_matrix, time_tree, time_bins, plot = FALSE, confidence.interval = 0.95) {

  # TO DO:
  #
  # - Add an aphylogenetic option
  # - Allow binning in other ways. E.g., by character grouping.

  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")
  
  # Check time_bins is in a valid format and stop and warn user if not:
  if (!is.timeBins(x = time_bins)) stop(check_timeBins(time_bins = time_bins))
  
  # Subfunction for getting missing and inapplicable characters:
  find_missing_and_inapplicable <- function(x) {

    # Get all inapplicables:
    inapplicables <- apply(apply(x, 2, "==", ""), 2, as.numeric)

    # Replace NAs with 0:
    inapplicables[is.na(inapplicables)] <- 0

    # Get missing:
    missings <- apply(apply(x, 2, is.na), 2, as.numeric)

    # Return numeric matrix (1 for missing or inapplicable):
    return(missings + inapplicables)
  }

  # Create vector of time bin mid-points:
  time_bin_midpoints <- find_time_bin_midpoints(time_bins = time_bins)

  # Set time bin names:
  time_bin_names <- rownames(time_bins)

  # Get total number of characters:
  n_characters <- sum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = "[[", "matrix"), FUN = ncol)))

  # Get edge lengths in bins for complete tree (measure of a complete character):
  binned_edge_lengths <- bin_edge_lengths(time_tree = time_tree, time_bins = time_bins)$binned_edge_lengths

  # Set uo missing values vector (no missing values are set as empty characters (""):
  missing_values <- rep("", n_characters)

  # If there are missing or inapplicable values collapse row numbers for them with double percentage:
  if (any(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), is.na))) || any(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), "==", "")))) missing_values <- unname(unlist(x = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), find_missing_and_inapplicable), apply, 2, "==", 1), apply, 2, which), lapply, paste, collapse = "%%")))

  # Set up matrix to store edge lengths in each character bin (columns) per character (rows):
  binned_edges_by_character <- matrix(data = 0, ncol = nrow(x = time_bins), nrow = n_characters)

  # For each unique missing value combination (no point repeating those that are the same):
  for (i in unique(x = missing_values)) {

    # If there is missing data:
    if (nchar(x = i) > 0) {

      # List taxa to prune:
      taxa_to_prune <- rownames(x = cladistic_matrix$matrix_1$matrix)[as.numeric(strsplit(i, "%%")[[1]])]

      # Check that there are still enough taxa left for a tree to exist:
      if (length(x = setdiff(x = time_tree$tip.label, y = taxa_to_prune)) > 1) {
        
        # Remove tips with missing data from time tree:
        pruned_tree <- drop_time_tip(time_tree = time_tree, tip_names = taxa_to_prune)

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

    # As long as the tree exists (i.e., it is not pruned down to one or zero taxa) store edge lengths in bin:
    if (!is.na(pruned_tree)[1]) binned_edges_by_character[which(x = missing_values == i), ] <- matrix(rep(bin_edge_lengths(time_tree = pruned_tree, time_bins = time_bins)$binned_edge_lengths, length(x = which(x = missing_values == i))), ncol = ncol(binned_edges_by_character), byrow = TRUE)
  }
  
  # Calculate and store proportional character completeness:
  binned_character_completeness <- binned_edges_by_character / matrix(rep(binned_edge_lengths, n_characters), nrow = n_characters, byrow = TRUE)

  # Calculate mean proportional character completeness:
  mean_binned_completeness <- apply(binned_character_completeness, 2, mean)

  # Calculate upper 95% confidence interval proportional character completeness:
  upper_binned_completeness <- apply(binned_character_completeness, 2, sort)[ceiling((confidence.interval + ((1 - confidence.interval) / 2)) * n_characters), ]

  # Calculate lower 95% confidence interval proportional character completeness:
  lower_binned_completeness <- apply(binned_character_completeness, 2, sort)[max(c(1, floor(((1 - confidence.interval) / 2) * n_characters))), ]

  # If plotting output:
  if (plot) {

    # Set plot environment for two plots (one on top of the other):
    par(mfrow = c(2, 1))

    # Create empty plot first (y-axis limits set from 0 to 1):
    graphics::plot(x = time_bin_midpoints, y = apply(binned_character_completeness, 2, mean), type = "n", ylim = c(0, 1), xlim = c(max(time_bin_midpoints), min(time_bin_midpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness")

    # Add 95% confidence interval as shaded polygon:
    graphics::polygon(x = c(time_bin_midpoints, rev(time_bin_midpoints)), y = c(upper_binned_completeness, rev(lower_binned_completeness)), col = "grey", border = 0)

    # Plot mean character completeness on top:
    graphics::points(x = time_bin_midpoints, y = mean_binned_completeness, type = "l", ylim = c(0, 1), xlim = c(max(time_bin_midpoints), min(time_bin_midpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness", lwd = 2)

    # Add legend:
    graphics::legend(x = max(time_bin_midpoints), y = 0.4, c(paste(confidence.interval, "% confidence interval", sep = ""), "Mean"), col = c("grey", "black"), lwd = c(8, 2), merge = TRUE, bg = "white")

    # Create empty plot first (y-axis limits set from 0 to 1):
    graphics::plot(x = time_bin_midpoints, y = apply(binned_character_completeness, 2, mean), type = "n", ylim = c(0, 1), xlim = c(max(time_bin_midpoints), min(time_bin_midpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness")

    # Plot each individual characters proportional completeness:
    for (i in 1:n_characters) graphics::points(x = time_bin_midpoints, y = binned_character_completeness[i, ], type = "l", col = "grey")

    # Plot mean character completeness on top:
    graphics::points(x = time_bin_midpoints, y = mean_binned_completeness, type = "l", ylim = c(0, 1), xlim = c(max(time_bin_midpoints), min(time_bin_midpoints)), xlab = "Time (Ma)", ylab = "Proportional Character Completeness", lwd = 2)

    # Add legend:
    graphics::legend(x = max(time_bin_midpoints), y = 0.4, legend = c("Individual characters", "Mean"), col = c("grey", "black"), lwd = c(1, 2), merge = TRUE, bg = "white")

    # Reset plotting environment:
    graphics::par(mfrow = c(1, 1))
  }

  # Add time bin names to output:
  names(mean_binned_completeness) <- names(upper_binned_completeness) <- names(lower_binned_completeness) <- colnames(x = binned_character_completeness) <- time_bin_names

  # Return compiled output:
  list(mean_completeness = mean_binned_completeness, upper_completeness = upper_binned_completeness, lower_completeness = lower_binned_completeness, completeness_by_character = binned_character_completeness)
}
