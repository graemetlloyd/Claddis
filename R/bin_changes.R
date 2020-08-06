#' Counts the changes in a series of time bins
#'
#' @description
#'
#' Given a vector of dates for a series of time bins and another for the times when a character change occurred will return the total number of changes in each bin.
#'
#' @param change_times A vector of ages in millions of years at which character changes are hypothesised to have occurred.
#' @param time_bins A vector of ages in millions of years of time bin boundaries in old-to-young order.
#'
#' @details
#'
#' Calculates the total number of evolutionary changes in a series of time bins. This is intended as an internal function for rate calculations, but could be used for other purposes (e.g., counting any point events in a series of time bins).
#'
#' @return
#'
#' A vector giving the number of changes for each time bin. Names indicate the maximum and minimum (bottom and top) values for each time bin.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a random dataset of 100 changes (between 100 and 0 Ma):
#' change_times <- stats::runif(n = 100, min = 0, max = 100)
#'
#' # Create 10 equal-length time bins:
#' time_bins <- seq(100, 0, length.out = 11)
#'
#' # Get N changes for each bin:
#' bin_changes(change_times, time_bins)
#' @export bin_changes
bin_changes <- function(change_times, time_bins) {

  # EXPLAIN HOW TIMES ON BOUDNARIES WORK AND CHECK TOTAL COUNTS MAKE SENSE AT THE END
  # MAYBE SWITCH TO LAPPLY INSTEAD OF FOR LOOP
  # ADD BOUNDARY TIME OPTION? I.E., WHICH BIN SHOULD THEY BE ASSIGNED TO?

  # Enforce old-to-young order of time bins:
  time_bins <- sort(x = time_bins, decreasing = TRUE)

  # Create all-zero vector to store ouput in:
  binned_changes <- rep(x = 0, times = length(x = time_bins) - 1)

  # For each time bin:
  for (i in 2:length(x = time_bins)) {

    # Find out which edges (if any) are present in the bin:
    binned_changes[(i - 1)] <- length(x = intersect(which(x = change_times > time_bins[i]), which(x = change_times <= time_bins[(i - 1)])))
  }

  # Add time bin max-mins as names:
  names(binned_changes) <- apply(cbind(time_bins[1:(length(x = time_bins) - 1)], time_bins[2:length(x = time_bins)]), 1, paste, collapse = "-")

  # Return edge lengths in bins:
  binned_changes
}
