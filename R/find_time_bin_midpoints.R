#' Find time bin midpoints
#'
#' @description
#'
#' Find the midpoint values for each bin from a timeBins object
#'
#' @param time_bins A timeBins object.
#'
#' @details
#'
#' Frequently the midpoints of a series of time bins (defined by a beginning and ending) will be required, for example, when plotting binned data as a time series. Although the calculation involved is trivial (i.e., start date + end date / 2) this is a sufficiently common operation it is made into a formal function here.
#'
#' Note that this function is designed to work specifically with objects of class "timeBins" - a format specific to Claddis that looks something like this:
#'
#' \preformatted{                  fad  lad
#'     Cenomanian    99.6 93.5
#'     Turonian      93.5 89.3
#'     Coniacian     89.3 85.8
#'     Santonian     85.8 83.5
#'     Campanian     83.5 70.6
#'     Maastrichtian 70.6 65.5}
#'
#' I.e., a matrix with two columns (fad = first appearance date and lad = last appearance date) with rows corresponding to named time bins and indiviual values ages in millions of years ago (Ma). The object should also have class \code{timeBins} (see example below for hot to generate such an object). Note also that the convention here is to have time bins be ordered from oldest to youngest.
#'
#' @return A vector of time bin midpoint values.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a time bins object:
#' time_bins <- matrix(
#'   data = c(99.6, 93.5, 93.5, 89.3, 89.3, 85.8, 85.8, 83.5, 83.5, 70.6, 70.6, 65.5),
#'   ncol = 2,
#'   byrow = TRUE,
#'   dimnames = list(
#'     c("Cenomanian", "Turonian", "Coniacian", "Santonian", "Campanian", "Maastrichtian"),
#'     c("fad", "lad")
#'   )
#' )
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Return midpoints for each time bin in sequence:
#' find_time_bin_midpoints(time_bins = time_bins)
#'
#' @export find_time_bin_midpoints
find_time_bin_midpoints <- function(time_bins) {

  # Check time_bins has class timeBins and stop and warn user if not:
  if (!inherits(x = time_bins, what = "timeBins")) stop("time_bins must be an object of class \"timeBins\".")
  
  # If not a valid timeBins object then stop and provide feedback to user on what is wrong:
  if (!is.timeBins(time_bins)) stop(check_timeBins(time_bins = time_bins)[1])
  
  # Return time bin midpoints:
  apply(X = time_bins, MARGIN = 1, FUN = mean)
}
