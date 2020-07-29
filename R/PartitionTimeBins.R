#' Time bin partitioner
#'
#' @description
#'
#' Generates all possible contiguous partitions of N time bins.
#'
#' @param NTimeBins The number of time bins.
#' @param NPartitonsToInclude Either "All" or a vector of requested partition sizes.
#'
#' @details
#'
#' This function is designed for use with the \link{DiscreteCharacterRate} function and generates all possible contiguous partitions of N time bins. This allows use of an information criterion like AIC to pick a "best" partition, weighing fit and partition number simultaneously.
#'
#' You can also ask for only partitions of a specific number using the \code{NPartitonsToInclude} option. For example, \code{NPartitonsToInclude = c(0, 1, 2)} will only return partitions of 1, 2, or 3 sets of elements.
#'
#' @return
#'
#' Returns a list of lists of vectors ready for use in \link{DiscreteCharacterRate}.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get all partitions for four time bins:
#' PartitionTimeBins(4)
#'
#' @export PartitionTimeBins
PartitionTimeBins <- function(NTimeBins, NPartitonsToInclude = "All") {
  
  # Build time bins vector:
  TimeBins <- 1:NTimeBins
  
  # Check there are multiple time bins:
  if(length(TimeBins) < 2) stop("There must be at least two time bins.")
  
  # Format NPartitonsToInclude as a vector of all possible numbers if input is "All":
  if(any(NPartitonsToInclude == "All")) NPartitonsToInclude <- 0:(length(TimeBins) - 1)
  
  # Get number of possible "switches", i.e., positions where a partition can (1) or cannot (0) be:
  NSwitches <- length(TimeBins) - 1
  
  # Calculate expeted number of paritions (if large can stop adn warn user):
  NPartitions <- sum(unlist(lapply(as.list(NPartitonsToInclude), function(x) ncol(combn(NSwitches, x)))))
  
  # Small correction if zero is included (no partitions or rather a partition of one):
  if(any(NPartitonsToInclude == 0)) NPartitions<- NPartitions + 1
  
  # Generate starting splitsiwtches vector:
  SplitSwitches <- as.character(0:1)
  
  # Generate all possible combinations of switches iteratively:
  while(nchar(SplitSwitches[1]) < (length(TimeBins) - 1)) SplitSwitches <- apply(expand.grid(SplitSwitches, as.character(0:1)), 1, paste, collapse = "")
  
  # Work out partition sizes for subsetting below:
  PartitionSizes <- unlist(lapply(strsplit(SplitSwitches, split = ""), function(x) sum(as.numeric(x))))
  
  # Collpase switches vector to just those of the reequired input size:
  SplitSwitches <- SplitSwitches[!is.na(match(PartitionSizes, NPartitonsToInclude))]
  
  # Subfnction to generate parition positions:
  PartitionPositions <- function(SwitchSequence) {
    
    # How long should the output vector be?:
    VectorLength <- nchar(SwitchSequence) + 1
    
    # Turn siwtches into split after points (adding end in to complete sequence):
    SplitAfter <- c(which(as.numeric(strsplit(SwitchSequence, split = "")[[1]]) == 1), VectorLength)
    
    # Seed start position:
    StartPosition <- 1
    
    # Set empty output list:
    Output <- list()
    
    # For each split:
    for(i in SplitAfter) {
      
      # Add to output:
      Output[[(length(Output) + 1)]] <- StartPosition:i
      
      # Update start position:
      StartPosition <- i + 1
      
    }
    
    # Return list of vectors:
    return(Output)
    
  }
  
  # Convert split siwtch strings to lists of vectors of partition elements:
  ParitionPositionList <- lapply(as.list(SplitSwitches), function(x) PartitionPositions(x))
  
  # Return list of vectors of partition elements:
  return(ParitionPositionList)
  
}
