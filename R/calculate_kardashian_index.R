#' Calculates a researcher's Kardashian Index
#'
#' @description
#'
#' Given counts of a researcher's Twitter followers and citations, returns their Kardashian Index.
#'
#' @param twitter_followers The number of twitter followers the researcher has.
#' @param total_citations The total number of citations across the researcher's publications (e.g., as garnered from a Google Scholar profile).
#'
#' @details
#'
#' This function implements the Kardashian Index of Hall (2014) and interested readers should consult that paper for more background.
#'
#' @return
#'
#' A scalar representing the ratio of expected Twitter followers (based on number of citations) to actual Twitter followers. Values greater than one indicate more Twitter followers than expected, those below one, fewer. According to Hall (2014), values above 5 are "Science Kardashians".
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Hall, N., 2014. The Kardashian index: a measure of discrepant social media profile for scientists. \emph{Genome Biology}, \bold{15}, 424.
#'
#' @examples
#'
#' # Calculate the Kardashian Index of Sam Giles (@GilesPalaeoLab)
#' # as of 10/5/21:
#' calculate_kardashian_index(
#'   twitter_followers = 6534,
#'   total_citations = 550
#' )
#'
#' # Calculate the Kardashian Index of Christopher Jackson (@seis_matters)
#' # as of 10/5/21:
#' calculate_kardashian_index(
#'   twitter_followers = 26000,
#'   total_citations = 6265
#' )
#'
#' # Calculate the Kardashian Index of Graeme T. Lloyd (@GraemeTLloyd)
#' # as of 10/5/21:
#' calculate_kardashian_index(
#'   twitter_followers = 2133,
#'   total_citations = 2780
#' )
#'
#' # Calculate the Kardashian Index of Katie Mack (@AstroKatie)
#' # as of 10/5/21:
#' calculate_kardashian_index(
#'   twitter_followers = 394900,
#'   total_citations = 1131
#' )
#'
#' @export calculate_kardashian_index
calculate_kardashian_index <- function(twitter_followers, total_citations) {
  
  # Set F_a as number of twitter followers:
  F_a <- twitter_followers
  
  # Calculate F_c using equation 1 of Hall (2014):
  F_c <- 43 * (total_citations ^ 0.32)
  
  # Retrurn Kardashian Index (equation 2 of Hall 2014):
  F_a / F_c
}
