#' Permute all ways to place n items into m bins
#'
#' @description
#'
#' Given a positive integer, n, and a number of bins (m), permutes all possible compositions.
#'
#' @param n A positive integer.
#' @param m_labels A character vector of labels for m.
#' @param allow_zero A logical indicating whether or not each bin should (\code{TRUE}) or should not (\code{FALSE}) be allowed to be zero.
#'
#' @details
#'
#' Every way that an integer (\code{n}) can be divided up into \code{m} bins can be permuted using a restricted version of the mathematical concept of compositions. In practice this function is designed to distribute the states for \code{n} tips across \code{m} states (e.g., with \link{permute_tipstates}), but many other uses are conceivable and hence this is included here as a general function.
#'
#' This algorithm reuses code from the \code{multicool} (Curran et al. 2021) and \code{partitions} (Hankin 2006) packages.
#'
#' The number of restricted compositions is given by the k-dimensional extension of triangular numbers (Baumann 2019):
#'
#' \itemize{
#'   \item{If \code{allow_zero = TRUE}, the binomial coefficient, n choose k, where n = \code{n + m} - 1 and k = \code{m}.}
#'   \item{If \code{allow_zero = FALSE}, the binomial coefficient, n choose k, where n = \code{n} - 1 and k = \code{m}.}
#' }
#'
#' @return A matrix where each row is a unique restricted composition of n and each column is a labelled bin.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Baumann, M. H., 2019. Die k-dimensionale Champagnerpyramide. Mathematische Semesterberichte, 66, 89-100.
#'
#' Curran, J., Williams, A., Kelleher, J. and Barber, D., 2021. multicool: Permutations of Multisets in Cool-Lex Order. R package version 0.1-12. https://CRAN.R-project.org/package=multicool.
#'
#' Hankin, R. K. S., 2006. Additive integer partitions in R. \emph{Journal of Statistical Software, Code Snippets}, \bold{16}, 1.
#'
#' @examples
#'
#' # Permute all the ways eight can be assigned to four bins (A, C, G, T),
#' # with each bin assigned at least one:
#' permute_restricted_compositions(
#'   n = 8,
#'   m_labels = c("A", "C", "G", "T"),
#'   allow_zero = FALSE
#' )
#'
#' @export permute_restricted_compositions
permute_restricted_compositions <- function(
  n,
  m_labels,
  allow_zero = FALSE
) {
  
  # Set m as length of m labels:
  m <- length(x = m_labels)
  
  # Special case of only one label (simple single restricted composition):
  if (m == 1) return(value = matrix(data = n, dimnames = list(c(), m_labels)))
  
  # Special case of not allowing zeroes and m and n being equal (single restricted compositon of one per label):
  if (!allow_zero && m == n) return(value = matrix(data = 1, ncol = m, dimnames = list(c(), m_labels)))
  
  # If allow_zero is FALSE but m is larger than n:
  if (!allow_zero && m > n) {
    
    # Reset allow_zero to TRUE as not possible to exclude zeroes:
    allow_zero <- TRUE
    
    # Inform user this has happened:
    cat("It is not possible to set allow_zero as FALSE as m is greater than n; allow_zero has been reset as TRUE.\n")
  }
  
  # If zeroes are not allowed reduce n by m (will add one later):
  if (!allow_zero) n <- n - m
  
  # Permute all restricted compositions:
  restricted_compositions <- do.call(
    what = rbind,
    args = apply(
      X = partitions::restrictedparts(n = n, m = m),
      MARGIN = 2,
      FUN = function(x) {
        multicool::allPerm(mcObj = multicool::initMC(x = x))
      },
      simplify = FALSE
    )
  )
  
  # If zeroes are not allowed then add one to each bin in every restricted composition:
  if (!allow_zero) restricted_compositions <- restricted_compositions + 1
  
  # Add labels to restricted compositions:
  colnames(x = restricted_compositions) <- m_labels
  
  # Output restricted compositions:
  restricted_compositions
}
