#' Time bin partitioner
#'
#' @description
#'
#' Generates all possible contiguous partitions of N time bins.
#'
#' @param n_time_bins The number of time bins.
#' @param partition_sizes_to_include Either "all" (the default) or a vector of requested partition sizes.
#'
#' @details
#'
#' This function is designed for use with the \link{test_rates} function and generates all possible contiguous partitions of N time bins. This allows use of an information criterion like AIC to pick a "best" partition, weighing fit and partition number simultaneously.
#'
#' You can also ask for only partitions of a specific number using the \code{partition_sizes_to_include} option. For example, \code{partition_sizes_to_include = c(1, 2, 3)} will only return partitions of 1, 2, or 3 sets of elements.
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
#'
#' # Get all partitions for five time bins of size 2:
#' partition_time_bins(n_time_bins = 5, partition_sizes_to_include = 2)
#' @export partition_time_bins
partition_time_bins <- function(n_time_bins, partition_sizes_to_include = "all") {

  # Build time bins vector:
  time_bins <- 1:n_time_bins

  # Check there are multiple time bins:
  if (length(x = time_bins) < 2) stop("There must be at least two time bins.")

  # Format partition_sizes_to_include as a vector of all possible numbers if input is "all":
  if (any(partition_sizes_to_include == "all")) partition_sizes_to_include <- 1:n_time_bins

  # Get number of possible "switches", i.e., positions where a partition can (1) or cannot (0) be:
  n_switches <- length(x = time_bins) - 1

  # Calculate expeted number of partitions (if large can stop and warn user):
  n_partitions <- sum(unlist(x = lapply(X = as.list(x = partition_sizes_to_include - 1), function(x) ncol(combn(n_switches, x)))))

  # Generate starting splitsiwtches vector:
  split_switches <- as.character(0:1)

  # Generate all possible combinations of switches iteratively:
  while (nchar(x = split_switches[1]) < (length(x = time_bins) - 1)) split_switches <- apply(expand.grid(split_switches, as.character(0:1)), 1, paste, collapse = "")

  # Work out partition sizes for subsetting below:
  partition_sizes <- unlist(x = lapply(X = strsplit(split_switches, split = ""), function(x) sum(as.numeric(x)))) + 1

  # Collpase switches vector to just those of the reequired input size:
  split_switches <- split_switches[!is.na(match(partition_sizes, partition_sizes_to_include))]

  # Subfunction to generate partition positions:
  set_partition_positions <- function(switch_sequence) {

    # How long should the output vector be?:
    vector_length <- nchar(x = switch_sequence) + 1

    # Turn siwtches into split after points (adding end in to complete sequence):
    split_after <- c(which(x = as.numeric(strsplit(switch_sequence, split = "")[[1]]) == 1), vector_length)

    # Seed start position:
    start_position <- 1

    # Set empty output list:
    output <- list()

    # For each split:
    for (i in split_after) {

      # Add to output:
      output[[(length(x = output) + 1)]] <- start_position:i

      # Update start position:
      start_position <- i + 1
    }

    # Return list of vectors:
    output
  }

  # Convert split switch strings to lists of vectors of partition elements and return:
  lapply(X = as.list(x = split_switches), function(x) set_partition_positions(x))
}
