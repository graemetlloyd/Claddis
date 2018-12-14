#' Get distance matrices from a cladistic matrix
#'
#' Takes a cladistic morphological dataset and converts it into a set of pairwise distances.
#'
#' @param morph.matrix A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @param Distance The distance metric to use. Must be one of \code{"GC"}, \code{"GED"}, \code{"RED"}, or \code{"MORD"} (the default).
#' @param TransformProportionalDistances Whether to transform the proportional distances (Gower and max). Options are \code{"none"}, \code{"sqrt"}, or \code{"arcsine_sqrt"} (the default).
#' @param PolymorphismBehaviour The distance behaviour for dealing with polymorphisms. Must be one of \code{"mean.difference"} or \code{"min.difference"} (the default.
#' @param UncertaintyBehaviour The distance behaviour for dealing with uncertainties. Must be one of \code{"mean.difference"} or \code{"min.difference"} (the default.
#'
#' @return
#'
#' \item{DistanceMetric}{The distance metric used.}
#' \item{DistanceMatrix}{The distance matrix returned.}
#' \item{ComparableCharacterMatrix}{The matrix of comparable characters used for each distance.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Thomas Guillerme \email{guillert@@tcd.ie}
#'
#' @keywords distance
#'
#' @examples
#'
#' # Get morphological distances for Day et
#' # al. (2016) data set:
#' distances <- MorphDistMatrix(Day2016)
#'
#' # Show distance metric:
#' distances$DistanceMetric
#'
#' # Show distance matrix:
#' distances$DistanceMatrix
#'
#' # Show number of characters that can be scored for
#' # each pairwise comparison:
#' distances$ComparableCharacterMatrix
#'
#' @export MorphDistMatrix
MorphDistMatrix <- function(morph.matrix, Distance = "MORD", TransformProportionalDistances = "arcsine_sqrt", PolymorphismBehaviour = "min.difference", UncertaintyBehaviour = "min.difference") {
  
  # Subfunction to find comparable characters for a pairwise taxon comparison:
  GetComparableCharacters <- function(interest.col, morph.matrix) {
    
    # Get intersection of characters that are coded for both taxa in a pair:
    output <- intersect(which(!is.na(morph.matrix[interest.col[[1]], ])), which(!is.na(morph.matrix[interest.col[[2]], ])))
    
    # Return output:
    return(list(output))
    
  }
  
  # Subfunction to get character strings for each pair of taxa:
  GetPairwiseCharacterStrings <- function(interest.col, morph.matrix) {
    
    # Get character states for first taxon in pair:
    row1 <- morph.matrix[rownames(morph.matrix)[interest.col[[1]]], ]
    
    # Get character states for second taxon in pair:
    row2 <- morph.matrix[rownames(morph.matrix)[interest.col[[2]]], ]
    
    # Return output as a list:
    return(list(row1, row2))
    
  }
  
  # Subfunction to subset pairwise comparisons by just comparable characters:
  SubsetPairwiseByComparable <- function(row.pair, comparable.characters) {
    
    # Collapse first row to just comparable characters:
    row.pair[[1]] <- row.pair[[1]][comparable.characters]
    
    # Collapse second row to just comparable characters:
    row.pair[[2]] <- row.pair[[2]][comparable.characters]
    
    # Output colapsed row pair:
    return(row.pair)
    
  }
  
  # Subfunction to edit polymorphic characters down to a single value:
  EditPolymorphisms <- function(comparisons, comparable.characters, ordering, PolymorphismBehaviour, UncertaintyBehaviour) {
    
    # Set first taxon values:
    firstrow <- comparisons[[1]]
    
    # Set second taxon values:
    secondrow <- comparisons[[2]]
    
    # Set comparable characters:
    compchar <- comparable.characters
    
    # Set ordering for comparable characters:
    charordering <- ordering[compchar]
    
    # Only if there are polymorphisms or uncertainties:
    if(length(c(grep("&", unique(c(firstrow, secondrow))), grep("/", unique(c(firstrow, secondrow))))) > 0) {
      
      # Find ampersands (polymorphisms):
      ampersand.elements <- sort(c(grep("&", firstrow), grep("&", secondrow)))
      
      # Find slashes (uncertianties):
      slash.elements <- sort(c(grep("/", firstrow), grep("/", secondrow)))
      
      # Combine to find all characters to check:
      characters.to.check <- sort(unique(c(ampersand.elements, slash.elements)))
      
      # Set behaviours as either the shared version or minimum difference if they contradict (may need to modify this later for more complex options):
      behaviour <- unlist(lapply(lapply(lapply(lapply(lapply(lapply(lapply(apply(apply(rbind(firstrow[characters.to.check], secondrow[characters.to.check]), 2, gsub, pattern = "[:0-9:]", replacement = ""), 2, list), unlist), function(x) x[nchar(x) > 0]), function(x) ifelse(nchar(x) > 0, strsplit(x, split = "")[[1]][1], x)), function(x) gsub(x, pattern = "&", replacement = PolymorphismBehaviour)), function(x) gsub(x, pattern = "/", replacement = UncertaintyBehaviour)), unique), function(x) ifelse(length(x) > 1, "min.difference", x)))
      
      # If behaviour is to find minimum differences:
      if(any(behaviour == "min.difference")) {
        
        # Set up minimum difference characters to check:
        min.characters.to.check <- characters.to.check[behaviour == "min.difference"]
        
        # Find intersecting character states for each character:
        IntersectionCharacter <- lapply(lapply(lapply(lapply(apply(rbind(firstrow[min.characters.to.check], secondrow[min.characters.to.check]), 2, strsplit, split = "&|/"), unlist), sort), rle), function(x) x$values[x$lengths > 1][1])
        
        # If at least one intersecting character state was found:
        if(any(!is.na(unlist(IntersectionCharacter)))) {
          
          # Record rows to update:
          rows.to.update <- which(!is.na(unlist(IntersectionCharacter)))
          
          # Store (first) shared state for both taxa:
          firstrow[min.characters.to.check[rows.to.update]] <- secondrow[min.characters.to.check[rows.to.update]] <- unlist(IntersectionCharacter)[rows.to.update]
          
          # Update minimum characters to check:
          min.characters.to.check <- min.characters.to.check[-rows.to.update]
          
        }
        
        # Only continue if there are still characters that need to be fixed:
        if(length(min.characters.to.check) > 0) {
          
          # Build two option matrices for every comparison:
          TwoOptionMatrices <- lapply(apply(rbind(firstrow[min.characters.to.check], secondrow[min.characters.to.check]), 2, strsplit, split = "&|/"), function(x) rbind(c(min(as.numeric(x[[1]])), max(as.numeric(x[[2]]))), c(max(as.numeric(x[[1]])), min(as.numeric(x[[2]])))))
          
          # Pick smallest difference as minimum and maximum states:
          MinMaxStates <- lapply(lapply(lapply(TwoOptionMatrices, function(x) x[which(abs(apply(x, 1, diff)) == min(abs(apply(x, 1, diff)))), ]), sort), as.character)
          
          # Set first row values(s):
          firstrow[min.characters.to.check] <- unlist(lapply(MinMaxStates, '[[', 1))
          
          # Set second row values(s):
          secondrow[min.characters.to.check] <- unlist(lapply(MinMaxStates, '[[', 2))
          
        }
        
      }
      
      # If any behaviour is to find mean differences:
      if(any(behaviour == "mean.difference")) {
        
        # Set up minimum difference characters to check:
        mean.characters.to.check <- characters.to.check[behaviour == "mean.difference"]
        
        # Build initial state matrices with column and row names as states for first and second rows:
        StateMatrices <- lapply(lapply(apply(rbind(firstrow[mean.characters.to.check], secondrow[mean.characters.to.check]), 2, list), lapply, strsplit, split = "&|/"), function(x) matrix(nrow = length(x[[1]][[1]]), ncol = length(x[[1]][[2]]), dimnames = list(x[[1]][[1]], x[[1]][[2]])))
        
        # Fill state matrices with raw differences between each state:
        StateMatrices <- lapply(StateMatrices, function(x) { for(i in 1:ncol(x)) for(j in 1:nrow(x)) x[j, i] <- abs(as.numeric(colnames(x)[i]) - as.numeric(rownames(x)[j])) ; return(x) })
        
        # If there are unordered characters present convert maximum distances to one:
        if(any(charordering[mean.characters.to.check] == "unord")) StateMatrices[which(charordering[mean.characters.to.check] == "unord")] <- lapply(StateMatrices[which(charordering[mean.characters.to.check] == "unord")], function(x) { x[x > 1] <- 1; return(x) })
        
        # Extract minimum and maximum states from each matrix with maximum being the mean distance:
        MinMaxStates <- lapply(lapply(lapply(StateMatrices, as.vector), mean), function(x) c(0, x))
        
        # Set first row values(s):
        firstrow[mean.characters.to.check] <- unlist(lapply(MinMaxStates, '[[', 1))
        
        # Set second row values(s):
        secondrow[mean.characters.to.check] <- unlist(lapply(MinMaxStates, '[[', 2))
        
      }
      
    }
    
    # Return the first and second rows either without polymorphisms or with them removed:
    return(list(firstrow, secondrow))
    
  }
  
  # Subfunction to get the absolute difference between the two rows:
  GetAbsoluteCharacterDifferences <- function(column) {
    
    # Isolate first row values:
    firstrow <- column[[1]]
    
    # Isolate second row values:
    secondrow <- column[[2]]
    
    # Get absolute differences between each pair of characters:
    return(list(abs(as.numeric(firstrow) - as.numeric(secondrow))))
    
  }
  
  # Subfunction to correct unordered distances to one::
  CorrectForUnordered <- function(differences, compchar, ordering) {
    
    # If unordered and distance greater than one replace with one:
    if(length(which(differences > 1)) > 0) differences[which(differences > 1)[which(ordering[compchar[which(differences > 1)]] == "unord")]] <- 1
    
    # Return corrected unordered distances:
    return(list(differences))
    
  }
  
  # Subfunction to find incomparable characters:
  FindIncomparableCharacters <- function(comparable.characters, morph.matrix) return(setdiff(1:ncol(morph.matrix), comparable.characters))
  
  # Subfunction to get weighted differences:
  WeightDifferences <- function(differences, comparable.characters, weights) return(list(as.numeric(weights[comparable.characters]) * differences))
  
  # Subfunction to get raw Euclidean distance:
  RawEuclideanDistance <- function(differences) return(dist(rbind(differences, rep(0, length(differences))), method = "euclidean"))
  
  # Subfunction to find maximum possible differences for the comparable characters:
  MaximumDIfferences <- function(comparable.characters, max.vals, min.vals) return(as.numeric(max.vals[comparable.characters]) - as.numeric(min.vals[comparable.characters]))
  
  # Subfunction to transform list of distances into an actual distance matrix:
  ConvertListToMatrix <- function(list, morph.matrix, diag = NULL) {
    
    # Set the number of rows:
    k <- nrow(morph.matrix)
    
    # Create the empty matrix:
    mat.out <- matrix(ncol = k, nrow = k)
    
    # Fill up the lower triangle:
    mat.out[lower.tri(mat.out)] <- unlist(list)
    
    # Make the matrix a distance matrix (both triangles have the same values):
    mat.out <- as.matrix(as.dist(mat.out))
    
    # If no diagonal is supplied:
    if(is.null(diag)) {
      
      # Set diagonal as zero:
      diag(mat.out) <- 0
      
    # If a diagonal is supplied:
    } else {
      
      # Add supplied diagonal as diagonal:
      diag(mat.out) <- diag
      
    }
    
    # Return matrix:
    return(mat.out)
    
  }
  
  # Subfunction to get count of complete characters for each taxon (diagonal in comparable characters matrix:
  CountCompleteCharacters <- function(column) return(length(column) - length(grep(TRUE, is.na(column))))
  
  # Subfunction to calculate the Gower Coefficient:
  CalculateGowerCoefficient <- function(differences, comparable.characters, weights) return(sum(differences) / sum(weights[comparable.characters]))
  
  # Subfunction to calculate MORD:
  CalculateMORD <- function(differences, maximum.differences) return(sum(differences) / sum(maximum.differences))
  
  # Subfunction for building starting GED data:
  BuildStartingGEDData <- function(differences, comparable.characters, morph.matrix, weights) return(rbind(c(differences, rep(NA, length(FindIncomparableCharacters(comparable.characters, morph.matrix)))), c(weights[comparable.characters], weights[FindIncomparableCharacters(comparable.characters, morph.matrix)])))

  # Check for step matrices and stop and warn user if found:
  if(is.list(morph.matrix$Topper$StepMatrices)) stop("Function cannot currently deal with step matrices.")
  
  # Check input of TransformProportionalDistances is valid and stop and warn if not:
  if(length(setdiff(TransformProportionalDistances, c("arcsine_sqrt", "none", "sqrt"))) > 0) stop("TransformProportionalDistances must be one of \"none\", \"sqrt\", or \"arcsine_sqrt\".")
  
  # Check input of distance is valid and stop and warn if not:
  if(length(setdiff(Distance, c("RED", "GED", "GC", "MORD"))) > 0) stop("Distance must be one or more of \"RED\", \"GED\", \"GC\", or \"MORD\".")
  
  # Check input for PolymorphismBehaviour is valid and stop and warn if not:
  if(length(setdiff(PolymorphismBehaviour, c("mean.difference", "min.difference"))) > 0) stop("PolymorphismBehaviour must be one or more of \"mean.difference\", or \"min.difference\".")
  
  # Check input for UncertaintyBehaviour is valid and stop and warn if not:
  if(length(setdiff(UncertaintyBehaviour, c("mean.difference", "min.difference"))) > 0) stop("UncertaintyBehaviour must be one or more of \"mean.difference\", or \"min.difference\".")
  
  # Isolate ordering element:
  ordering <- unlist(lapply(morph.matrix[2:length(morph.matrix)], '[[', "Ordering"))
  
  # Isolate minimum values:
  min.vals <- unlist(lapply(morph.matrix[2:length(morph.matrix)], '[[', "MinVals"))
  
  # Isolate maximum values:
  max.vals <- unlist(lapply(morph.matrix[2:length(morph.matrix)], '[[', "MaxVals"))
  
  # Isolate weights:
  weights <- unlist(lapply(morph.matrix[2:length(morph.matrix)], '[[', "Weights"))
  
  # Combine matrix blocks into a single matrix:
  morph.matrix <- do.call(cbind, lapply(morph.matrix[2:length(morph.matrix)], '[[', "Matrix"))
  
  # If there are inapplicables then convert these to NAs:
  if(any(sort(morph.matrix == ""))) morph.matrix[morph.matrix == ""] <- NA
  
  # Find all possible (symmetric) pariwise comparisons for the N taxa in the matrix (excluding self-comparisons):
  comparisons <- combn(1:nrow(morph.matrix), 2)
  
  # Find all comparable characters for each pair of taxa:
  list.of.compchar <- unlist(apply(comparisons, 2, GetComparableCharacters, morph.matrix), recursive = FALSE)
  
  # Get character states for each pairwise comparison:
  rows.pairs <- apply(comparisons, 2, GetPairwiseCharacterStrings, morph.matrix)
  
  # Subset each pairwise comparison by just the comparable characters:
  matrix.of.char.comp <- mapply(SubsetPairwiseByComparable, rows.pairs, list.of.compchar)
  
  # Deal with any polymorphisms found and collapse appropriately:
  matrix.of.char.comp <- mapply(EditPolymorphisms, unlist(apply(matrix.of.char.comp, 2, list), recursive = FALSE), list.of.compchar, MoreArgs = list(ordering, PolymorphismBehaviour, UncertaintyBehaviour))
  
  # Get the absolute differences between each comparable character for each pairwise comparison:
  raw.diffs <- diffs <- unlist(apply(matrix.of.char.comp, 2, GetAbsoluteCharacterDifferences), recursive = FALSE)
  
  # Correct distances for unordered characters where distance is greater than one:
  diffs <- mapply(CorrectForUnordered, diffs, list.of.compchar, MoreArgs = list(ordering))
  
  # Weight differences:
  diffs <- mapply(WeightDifferences, diffs, list.of.compchar, MoreArgs = list(weights))
  
  # Get raw Euclidean distance:
  raw.dist <- lapply(diffs, RawEuclideanDistance)
  
  # Only calculate the max differences for "GED" or "MORD" matrices:
  if(Distance == "GED" || Distance == "MORD") {
    
    # Find maximum possible differences for the comparable characters:
    maxdiffs <- lapply(list.of.compchar, MaximumDIfferences, max.vals, min.vals)
    
    # Correct maximum differences for unordered characters:
    maxdiffs <- mapply(WeightDifferences, mapply(CorrectForUnordered, maxdiffs, list.of.compchar, MoreArgs = list(ordering)), list.of.compchar, MoreArgs = list(weights))
    
  }
  
  # If calculating Raw Eucldiean Distances build the distance matrix:
  if(Distance == "RED") dist.matrix <- ConvertListToMatrix(raw.dist, morph.matrix)

  # If calculating the Gower Coefficient build the distance matrix:
  if(Distance == "GC") dist.matrix <- ConvertListToMatrix(as.list(mapply(CalculateGowerCoefficient, diffs, list.of.compchar, MoreArgs = list(weights))), morph.matrix)
  
  # If calculating the MORD build the distance matrix:
  if(Distance == "MORD") dist.matrix <- ConvertListToMatrix(mapply(CalculateMORD, diffs, maxdiffs), morph.matrix)
  
  # If calculating the GED:
  if(Distance == "GED") {
    
    # Build starting GED data:
    GED.data <- mapply(BuildStartingGEDData, diffs, list.of.compchar, MoreArgs = list(morph.matrix, weights), SIMPLIFY = FALSE)
    
    # Transpose matrices:
    GED.data <- lapply(GED.data, t)
    
    #TG: second, create the matrix
    GED.data <- matrix(data = (unlist(GED.data)), ncol = ncol(morph.matrix), byrow = TRUE)
    
    # Add to maximum differences (S_ijk * W_ijk in equation 1 of Wills 2001):
    differences <- unlist(diffs)
    
    # Add to maximum differences (S_ijk_max * W_ijk in equation 1 of Wills 2001):
    maximum.differences <- unlist(maxdiffs)
    
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
    
    # Create GED distance matrix:
    dist.matrix <- ConvertListToMatrix(as.list(GED_ij), morph.matrix)
    
  }
  
  # Build comparable characters matrix:
  comp.char.matrix <- ConvertListToMatrix(lapply(list.of.compchar, length), morph.matrix, diag = apply(morph.matrix, 1, CountCompleteCharacters))
  
  # Add row and column names (taxa) to distance matrices:
  rownames(dist.matrix) <- colnames(dist.matrix) <- rownames(comp.char.matrix) <- colnames(comp.char.matrix) <- rownames(morph.matrix)
  
  # If there are any NaNs replace with NAs:
  if(any(is.nan(dist.matrix))) dist.matrix[is.nan(dist.matrix)] <- NA
  
  # If transforming distance matrix by taking the square root - take the square root:
  if(TransformProportionalDistances == "sqrt") dist.matrix <- sqrt(dist.matrix)
  
  # If transforming distance matrix by taking the arcsine square root:
  if(TransformProportionalDistances == "arcsine_sqrt") {
    
    # Check for squared distances greater than 1:
    if(any(sort(sqrt(dist.matrix)) > 1)) {
      
      # Warn user that distances were rescaled:
      print("Squared distances found of greater than 1 so matrix was rescaled prior to taking arcsine.")
      
      # Take the arcsine square root of the rescaled distance matrix:
      dist.matrix <- asin(sqrt(dist.matrix) / max(sort(sqrt(dist.matrix))))
      
    # If squared distances are less than or equal to one:
    } else {
      
      # Take the arcsine square root directly:
      dist.matrix <- asin(sqrt(dist.matrix))
      
    }
    
  }
  
  # Compile results as a list:
  result <- list(Distance, dist.matrix, comp.char.matrix)
  
  # Add names to list:
  names(result) <- c("DistanceMetric", "DistanceMatrix", "ComparableCharacterMatrix")
  
  # Output result:
  return(result)
  
}
