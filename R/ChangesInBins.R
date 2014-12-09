ChangesInBins <- function(change.times, time.bins) {
	
	# Enforce old-to-young order of time bins:
	time.bins <- sort(time.bins, decreasing=T)
	
	# Create all-zero vector to store ouput in:
	changes.in.bin <- rep(0, length(time.bins) - 1)
	
	# For each time bin:
	for(i in 2:length(time.bins)) {
		
		# Find out which edges (if any) are present in the bin:
		changes.in.bin[(i - 1)] <- length(intersect(which(change.times > time.bins[i]), which(change.times <= time.bins[(i - 1)])))
		
	}
	
	# Add time bin max-mins as names:
	names(changes.in.bin) <- apply(cbind(time.bins[1:(length(time.bins) - 1)], time.bins[2:length(time.bins)]), 1, paste, collapse="-")
	
	# Return edge lengths in bins:
	return(changes.in.bin)
	
}
