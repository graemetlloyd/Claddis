#' Get distance matrices from a cladistic matrix
#' 
#' Takes a cladistic morphological dataset and converts it into a set of pairwise distances.
#' 
#' This function is an important preliminary step in performing multivariate ordination(s) upon a cladistic dataset of discrete characters and deals with three peculiarities of such data: 1) the prevalence of missing values, 2) the use of unordered multistate characters, and 3) the presence of polymorphisms.
#' 
#' Missing data is dealt with by providing three rescaled distances as well as the uncorrected raw distances. These are: 1) the Generalised Euclidean Distance (GED) of Wills (2001) which replaces all incalculable taxon-taxon character distances with a weighted fractional mean distance, 2) Gower (1971) dissimilarity, which rescales each pairwise distance by dividing by the number of comparable characters upon which that distance is based, and 3) a set of distances rescaled against the maximum possible observable distance based on the set of comparable characters upon which that distance is based (using the difference between \code{max.vals} and \code{min.vals} supplied by the input matrix). In practice the user is not recommended to use raw distances unless the input matrix is complete.
#' 
#' The method automatically treats distances between unordered characters as zero (if the two states are identical) or one (if the two states are different). So, for example, the distances between 0 and 3 and between 2 and 3 for an unordered character are both 1.
#' 
#' Finally, polymorphisms are dealt with by using the minimum possible distance by considering all possible values implied by the polymorphism. So, for example, the distance between (01) and 3 will be considered to be 1 for an unordered character and 2 (the minimum distance, between 1 and 3) for an ordered character.
#' 
#' All metrics are rescaled according to character weightings, as supplied by the \code{weights} vector of the input matrix.
#' 
#' It is important to note that in practice not all pairwise distances can be calculated due to missing data. In these cases incomplete distance matrices will be returned, with incalculable values scored as \code{NA}. In such cases the user is advised to apply the \link{TrimMorphDistMatrix} function before ordination.
#' 
#' @param morph.matrix A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @param transform.proportional.distances Whether to transform the proportional distances (Gower and max). Options are \code{none}, \code{sqrt}, or \code{arcsine_sqrt} (the default).
#'
#' @return
#'
#' \item{raw.dist.matrix}{A distance matrix indicating the raw distances (based on the method set in \code{dist.method}) of each pairwise comparison.}
#' \item{GED.dist.matrix}{A Generalised Euclidean Distance (GED) matrix (Wills 2001).}
#' \item{gower.dist.matrix}{A distance matrix where raw distances have been rescaled using Gower (1971) dissimilarity (then rescaled from 0 to 1 and the arcsine of the square root taken).}
#' \item{max.dist.matrix}{A distance matrix where raw distances have been rescaled against the maximum possible observable distance (then the arcsine of the square root taken).}
#' \item{comp.char.matrix}{A matrix showing the number of coded characters that were used to make each pairwise comparison.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Gower, J. C., 1971. A general coefficient of similarity and some of its properties. Biometrics, 27, 857-871.
#' 
#' Wills, M. A., 2001a. Morphological disparity: a primer. In: Adrain, J. M., Edgecombe, G. D. and Lieberman, B. S. (eds.), Fossils, Phylogeny, and Form: An Analytical Approach. Kluwer Academic/Plenum Publishers, New York, p55-144.
#'
#' @keywords distance
#'
#' @examples
#' 
#' # Get morphological distances for Michaux (1989)
#' # data set:
#' distances <- MorphDistMatrix(Michaux1989)
#' 
#' # Show raw Euclidean distances:
#' distances$raw.dist.matrix
#' 
#' # Show Generailsed Euclidean Distances:
#' distances$GED.dist.matrix
#' 
#' # Show Gower rescaled distances:
#' distances$gower.dist.matrix
#' 
#' # Show maximum observable rescaled distances:
#' distances$max.dist.matrix
#' 
#' # Show number of characters that can be scored for
#' # each pairwise comparison:
#' distances$comp.char.matrix
#' 
#' @export MorphDistMatrix
MorphDistMatrix <- function(morph.matrix, transform.proportional.distances = "arcsine_sqrt") {

  # Check format of transform.proportional.distances:
  if(transform.proportional.distances != "none" && transform.proportional.distances != "sqrt" && transform.proportional.distances != "arcsine_sqrt") {

    # Give error if something other than three possible settings is given:
    stop("ERROR: transform.proportional.distances must be one of \"none\", \"sqrt\", or \"arcsine_sqrt\".")

  }

  # Isolate ordering element of morphology matrix:
  ordering <- morph.matrix$ordering

  # Isolate max values element of morphology matrix:
  max.vals <- morph.matrix$max.vals
  
  # Isolate min values element of morphology matrix:
  min.vals <- morph.matrix$min.vals
  
  # Isolate weighting element of morphology matrix:
  weights <- morph.matrix$weights
  
  # Isolate character-taxon matrix element of morphology matrix:
  morph.matrix <- morph.matrix$matrix

  # Create empty vectors to store S and W value for Wills 2001 equations (1 and 2):
  differences <- maximum.differences <- vector(mode="numeric")
    
  # Distance matrices for storing:
  comp.char.matrix <- gower.dist.matrix <- max.dist.matrix <- dist.matrix <- matrix(0, nrow=length(rownames(morph.matrix)), ncol=length(rownames(morph.matrix)))
    
  # Fill comparable characters diagonal:
  for(i in 1:length(morph.matrix[, 1])) comp.char.matrix[i,i] <- length(morph.matrix[i, ]) - length(which(is.na(morph.matrix[i, ])))

  # Set up empty matrix for storing data to calculate the Generalised Euclidean Distance of Wills (2001):
  GED.data <- matrix(nrow=0, ncol=ncol(morph.matrix))

  # Go through matrix rows:
  for(i in 1:(length(morph.matrix[, 1]) - 1)) {
        
    # Go through matrix columns:
    for(j in (i + 1):length(morph.matrix[, 1])) {
            
      # Get just the comparable characters (those coded for both taxa):
      compchar <- intersect(which(!is.na(morph.matrix[rownames(morph.matrix)[i], ])), which(!is.na(morph.matrix[rownames(morph.matrix)[j], ])))
            
      # Get comparable characters for ith taxon:
      firstrow <- morph.matrix[rownames(morph.matrix)[i], compchar]

      # Get comparable characters for jth taxon:
      secondrow <- morph.matrix[rownames(morph.matrix)[j], compchar]
            
      # Deal with polymorphic characters (if present):
      if(length(grep("&", unique(c(firstrow, secondrow)))) > 0) {
                
        # Find ampersands (polymorphisms):
        ampersand.elements <- sort(c(grep("&", firstrow), grep("&", secondrow)))
                
        # Go through each polymorphic character:
        for(k in 1:length(ampersand.elements)) {
                    
          # Find out if two codings overlap once all polymorphism resolutions are considered:
          intersection.value <- intersect(strsplit(firstrow[ampersand.elements[k]], "&")[[1]], strsplit(secondrow[ampersand.elements[k]], "&")[[1]])
                    
          # Case if polymorphic and non-polymorphic values overlap:
          if(length(intersection.value) > 0) {
                        
            # Set ith value as zero (no difference):
            firstrow[ampersand.elements[k]] <- 0

            # Set jth value as zero (no difference)
            secondrow[ampersand.elements[k]] <- 0

          }
                    
          # Case if polymorphic and non-polymorphic values do not overlap:
          if(length(intersection.value) == 0) {
                        
            # Case if character is unordered (max difference is 1):
            if(ordering[compchar[ampersand.elements[k]]] == "unord") {
                            
              # Set ith value as zero:
              firstrow[ampersand.elements[k]] <- 0

              # Set jth value as 1 (making the ij difference equal to one):
              secondrow[ampersand.elements[k]] <- 1

            }
                        
            # Case if character is ordered (max difference is > 1):
            if(ordering[compchar[ampersand.elements[k]]] == "ord") {
                            
              # Get first row value(s):
              firstrowvals <- as.numeric(strsplit(firstrow[ampersand.elements[k]], "&")[[1]])
                            
              # Get second row value(s):
              secondrowvals <- as.numeric(strsplit(secondrow[ampersand.elements[k]], "&")[[1]])
                            
              # Make mini distance matrix:
              poly.dist.mat <- matrix(0, nrow=length(firstrowvals), ncol=length(secondrowvals))
                            
              # Go through each comparison:
              for(l in 1:length(firstrowvals)) {
                                
                # Record absolute difference:
                for(m in 1:length(secondrowvals)) poly.dist.mat[l, m] <- sqrt((firstrowvals[l] - secondrowvals[m]) ^ 2)

              }
                            
              # Set first value as zero:
              firstrow[ampersand.elements[k]] <- 0
                            
              # Set second value as minimum possible difference:
              secondrow[ampersand.elements[k]] <- min(poly.dist.mat)

            }

          }

        }

      }
            
      # Get the absolute difference between the two rows:
      raw.diffs <- diffs <- abs(as.numeric(firstrow) - as.numeric(secondrow))
            
      # If there are differences greater than 1 for unordered characters then rescore as 1:
      if(length(which(diffs > 1)) > 0) diffs[which(diffs > 1)[which(ordering[compchar[which(diffs > 1)]] == "unord")]] <- 1

      # Find the incomparable characters:
      incompchar <- setdiff(1:ncol(morph.matrix), compchar)

      # Store data for GED with NAs for missing distances:
      GED.data <- rbind(GED.data, rbind(c(diffs, rep(NA, length(incompchar))), c(weights[compchar], weights[incompchar])))

      # Get weighted differences:
      diffs <- as.numeric(weights[compchar]) * diffs

      # Get raw Euclidean distance:
      raw.dist <- dist(rbind(diffs, rep(0, length(diffs))), method="euclidean")

      # Work out maximum difference (again, checked against ordering) using compchar characters only:
      raw.maxdiffs <- maxdiffs <- as.numeric(max.vals[compchar]) - as.numeric(min.vals[compchar])

      # Correct maximum possible differences for unordered characters:
      if(length(which(maxdiffs > 1)) > 0) maxdiffs[which(maxdiffs > 1)[which(ordering[compchar[which(maxdiffs > 1)]] == "unord")]] <- 1

      # Get vector of maximum differences (corrected for character weights):
      maxdiffs <- as.numeric(weights[compchar]) * maxdiffs

      # Store raw distance:
      dist.matrix[i, j] <- dist.matrix[j, i] <- raw.dist

      # Store Gower distance:
      gower.dist.matrix[i, j] <- gower.dist.matrix[j, i] <- sum(diffs) / sum(weights[compchar])

      # Store maximum-rescaled distance:
      max.dist.matrix[i, j] <- max.dist.matrix[j, i] <- sum(diffs) / sum(maxdiffs)

      # Store N comparable characters:
      comp.char.matrix[i, j] <- comp.char.matrix[j, i] <- length(compchar)

      # Add to maximum differences (S_ijk * W_ijk in equation 1 of Wills 2001):
      differences <- c(differences, diffs)

      # Add to maximum differences (S_ijk_max * W_ijk in equation 1 of Wills 2001):
      maximum.differences <- c(maximum.differences, maxdiffs)

    }

  }

  # Calculated weighted mean univariate distance for calculating GED (equation 2 in Wills 2001):
  S_ijk_bar <- sum(differences) / sum(maximum.differences)

  # Replace missing distances with S_ijk_bar (i.e., results of equation 2 in Wills 2001 into equation 1 of Wills 2001):
  GED.data[is.na(GED.data)] <- S_ijk_bar

  # Isolate the distances:
  S_ijk <- GED.data[which((1:nrow(GED.data) %% 2) == 1), ]

  # Isolate the weights:
  W_ijk <- GED.data[which((1:nrow(GED.data) %% 2) == 0), ]

  # Calculate the GED (equation 1 of Wills 2001) for each pairwise comparison (ij):
  GED_ij <- sqrt(apply(W_ijk * (S_ijk ^ 2), 1, sum))

  # Create empty GED distance matrix:
  GED.dist.matrix <- matrix(0, nrow=nrow(morph.matrix), ncol=nrow(morph.matrix))

  # Set initial value for counter:
  counter <- 1

  # Go through matrix rows:
  for(i in 1:(length(morph.matrix[, 1]) - 1)) {

    # Go through matrix columns:
    for(j in (i + 1):length(morph.matrix[, 1])) {

      # Store distance:
      GED.dist.matrix[i, j] <- GED.dist.matrix[j, i] <- GED_ij[counter]

      # Update counter:
      counter <- counter + 1

    }

  }
    
  # Set diagonals as zero:
  diag(gower.dist.matrix) <- diag(max.dist.matrix) <- 0

  # Add row and column names (taxa) to distance matrices:
  rownames(comp.char.matrix) <- colnames(comp.char.matrix) <- rownames(GED.dist.matrix) <- colnames(GED.dist.matrix) <- rownames(gower.dist.matrix) <- colnames(gower.dist.matrix) <- rownames(max.dist.matrix) <- colnames(max.dist.matrix) <- rownames(dist.matrix) <- colnames(dist.matrix) <- rownames(morph.matrix)

  # If transformation option is not "none":
  if(transform.proportional.distances != "none") {

    # If transformation option is "sqrt":
    if(transform.proportional.distances == "sqrt") {

      # Replace NaN with NA for Gower distances and take square root:
      gower.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, sqrt(gower.dist.matrix))), nrow=nrow(gower.dist.matrix), dimnames=list(rownames(gower.dist.matrix), rownames(gower.dist.matrix)))

      # Replace NaN with NA for Max distances and take square root:
      max.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, sqrt(max.dist.matrix))), nrow=nrow(max.dist.matrix), dimnames=list(rownames(max.dist.matrix), rownames(max.dist.matrix)))

    # If transformation option is "arcsine_sqrt":
    } else {

      # Establish correction factor to ensure Gower data is proportional:
      gower.correction <- max(c(max(sort(gower.dist.matrix)), 1))

      # Ensure all Gower values are on 0 to 1 scale then take arcsine of sqrt to get values that better approximate a normal distribution:
      gower.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, asin(sqrt(gower.dist.matrix / gower.correction)))), nrow=nrow(gower.dist.matrix), dimnames=list(rownames(gower.dist.matrix), rownames(gower.dist.matrix)))

      # Take arcsine square root of all MOD dist values:
      max.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, asin(sqrt(max.dist.matrix)))), nrow=nrow(max.dist.matrix), dimnames=list(rownames(max.dist.matrix), rownames(max.dist.matrix)))

    }

  # If transformation option is "none":
  } else {

    # Replace NaN with NA for Gower distances:
    gower.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, gower.dist.matrix)), nrow=nrow(gower.dist.matrix), dimnames=list(rownames(gower.dist.matrix), rownames(gower.dist.matrix)))

    # Replace NaN with NA for Max distances:
    max.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, max.dist.matrix)), nrow=nrow(max.dist.matrix), dimnames=list(rownames(max.dist.matrix), rownames(max.dist.matrix)))

  }

  # Compile results as a list:
  result <- list(dist.matrix, GED.dist.matrix, gower.dist.matrix, max.dist.matrix, comp.char.matrix)
    
  # Add names to results list:
  names(result) <- c("raw.dist.matrix", "GED.dist.matrix", "gower.dist.matrix", "max.dist.matrix", "comp.char.matrix")
    
  # Output result:
  return(result)

}
