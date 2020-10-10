#' Counts the changes in a series of time bins
#'
#' @description
#'
#' Given a vector of dates for a series of time bins and another for the times when a character change occurred will return the total number of changes in each bin.
#'
#' @param change_times A vector of ages in millions of years at which character changes are hypothesised to have occurred.
#' @param time_bins An object of class \code{timeBins}.
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
#' time_bins <- matrix(data = c(seq(from = 100, to = 10, length.out = 10),
#'   seq(from = 90, to = 0, length.out = 10)), ncol = 2,
#'   dimnames = list(LETTERS[1:10], c("fad", "lad")))
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Get N changes for each bin:
#' bin_changes(change_times = change_times, time_bins = time_bins)
#' @export bin_changes
bin_changes <- function(change_times, time_bins) {

  # EXPLAIN HOW TIMES ON BOUDNARIES WORK AND CHECK TOTAL COUNTS MAKE SENSE AT THE END
  # MAYBE SWITCH TO LAPPLY INSTEAD OF FOR LOOP
  # ADD BOUNDARY TIME OPTION? I.E., WHICH BIN SHOULD THEY BE ASSIGNED TO?

  # Check time_bins is in a valid format and stop and warn user if not:
  if (!is.timeBins(x = time_bins)) stop(check_timeBins(time_bins = time_bins))

  # Create all-zero vector to store ouput in:
  binned_changes <- rep(x = 0, times = nrow(x = time_bins))

  # For each time bin:
  for (i in 1:nrow(x = time_bins)) {

    # Find out which edges (if any) are present in the bin:
    binned_changes[i] <- length(x = intersect(which(x = change_times > time_bins[i, "lad"]), which(x = change_times <= time_bins[i, "fad"])))
  }

  # Add time bin names to binned changes:
  names(binned_changes) <- rownames(time_bins)

  # Return edge lengths in bins:
  binned_changes
}
