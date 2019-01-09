#' Counts the changes in a series of time bins
#'
#' @description
#'
#' Given a vector of dates for a series of time bins and another for the times when a character change occurred will return the total number of changes in each bin.
#' 
#' @param change.times A vector of ages in millions of years at which character changes are hypothesised to have occurred.
#' @param time.bins A vector of ages in millions of years of time bin boundaries in old-to-young order.
#'
#' @details
#'
#' Calculates the total number of evolutionary changes in a series of time bins. This is intended as an internal function for rate calculations, but could be used for other purposes (e.g., counting any point events in a series of time bins).
#'
#' @return A vector giving the number of changes for each time bin. Names indicate the maximum and minimum (bottom and top) values for each time bin.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#' 
#' # Create a random dataset of 100 changes:
#' change.times <- runif(100, 0, 100)
#' 
#' # Create time bins:
#' time.bins <- seq(100, 0, length.out=11)
#' 
#' # Get N changes for each bin:
#' ChangesInBins(change.times, time.bins)
#' 
#' @export ChangesInBins
ChangesInBins <- function(change.times, time.bins) {
	
	# Enforce old-to-young order of time bins:
	time.bins <- sort(time.bins, decreasing = TRUE)
	
	# Create all-zero vector to store ouput in:
	changes.in.bin <- rep(0, length(time.bins) - 1)
	
	# For each time bin:
	for(i in 2:length(time.bins)) {
		
		# Find out which edges (if any) are present in the bin:
		changes.in.bin[(i - 1)] <- length(intersect(which(change.times > time.bins[i]), which(change.times <= time.bins[(i - 1)])))
		
	}
	
	# Add time bin max-mins as names:
	names(changes.in.bin) <- apply(cbind(time.bins[1:(length(time.bins) - 1)], time.bins[2:length(time.bins)]), 1, paste, collapse = "-")
	
	# Return edge lengths in bins:
	return(changes.in.bin)
	
}
