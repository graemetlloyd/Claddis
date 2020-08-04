#' Time bin partitioner
#'
#' @description
#'
#' Generates all possible contiguous partitions of N time bins.
#'
#' @param n_time_bins The number of time bins.
#' @param partitions_to_include Either "All" or a vector of requested partition sizes.
#'
#' @details
#'
#' This function is designed for use with the \link{test_rates} function and generates all possible contiguous partitions of N time bins. This allows use of an information criterion like AIC to pick a "best" partition, weighing fit and partition number simultaneously.
#'
#' You can also ask for only partitions of a specific number using the \code{partitions_to_include} option. For example, \code{partitions_to_include = c(0, 1, 2)} will only return partitions of 1, 2, or 3 sets of elements.
#'
#' @return
#'
#' Returns a list of lists of vectors ready for use in \link{test_rates}.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get all partitions for four time bins:
#' partition_time_bins(n_time_bins = 4)
#' @export partition_time_bins
partition_time_bins <- function(n_time_bins, partitions_to_include = "All") {

  # Build time bins vector:
  time_bins <- 1:n_time_bins

  # Check there are multiple time bins:
  if (length(time_bins) < 2) stop("There must be at least two time bins.")

  # Format partitions_to_include as a vector of all possible numbers if input is "All":
  if (any(partitions_to_include == "All")) partitions_to_include <- 0:(length(time_bins) - 1)

  # Get number of possible "switches", i.e., positions where a partition can (1) or cannot (0) be:
  n_switches <- length(time_bins) - 1

  # Calculate expeted number of paritions (if large can stop adn warn user):
  n_partitions <- sum(unlist(lapply(as.list(partitions_to_include), function(x) ncol(combn(n_switches, x)))))

  # Small correction if zero is included (no partitions or rather a partition of one):
  if (any(partitions_to_include == 0)) n_partitions <- n_partitions + 1

  # Generate starting splitsiwtches vector:
  split_switches <- as.character(0:1)

  # Generate all possible combinations of switches iteratively:
  while (nchar(split_switches[1]) < (length(time_bins) - 1)) split_switches <- apply(expand.grid(split_switches, as.character(0:1)), 1, paste, collapse = "")

  # Work out partition sizes for subsetting below:
  partition_sizes <- unlist(lapply(strsplit(split_switches, split = ""), function(x) sum(as.numeric(x))))

  # Collpase switches vector to just those of the reequired input size:
  split_switches <- split_switches[!is.na(match(partition_sizes, partitions_to_include))]

  # Subfunction to generate partition positions:
  set_partition_positions <- function(switch_sequence) {

    # How long should the output vector be?:
    vector_length <- nchar(switch_sequence) + 1

    # Turn siwtches into split after points (adding end in to complete sequence):
    split_after <- c(which(as.numeric(strsplit(switch_sequence, split = "")[[1]]) == 1), vector_length)

    # Seed start position:
    start_position <- 1

    # Set empty output list:
    output <- list()

    # For each split:
    for (i in split_after) {

      # Add to output:
      output[[(length(output) + 1)]] <- start_position:i

      # Update start position:
      start_position <- i + 1
    }

    # Return list of vectors:
    return(output)
  }

  # Convert split siwtch strings to lists of vectors of partition elements:
  partition_position_list <- lapply(as.list(split_switches), function(x) set_partition_positions(x))

  # Return list of vectors of partition elements:
  return(partition_position_list)
}
