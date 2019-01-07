#' Get distance matrices from a cladistic matrix
#'
#' There are many options to consider when generating a distance matrix from morphological data, including the metric to use, how to treat inapplicable, polymorphic (e.g., 0&1), or uncertain (e.g., 0/1) states, and whether the output should be transformed (e.g., by taking the square root so that the distances are - or approximate - Euclidean distances). Some of these issues have been discussed previously in the literature (e.g., Lloyd 2016; Hopkins and St John 2018), but all likely require further study.
#'
#' Claddis currently offers four different distance metrics: 1. Raw Euclidean Distance (\code{RED}) - this is only really applicable if there are no missing data, 2. The Gower Coefficient (\code{GC}; Gower 1971) - this rescales distances by the number of characters that can be coded for both taxa in each pairwise comparison thus correcting for missing data, 3. The Maximum Observable Rescaled Distance (\code{MORD}) - this was introduced by Lloyd (2016) as an extension of the \code{GC} designed to deal with the fact that multistate ordered characters can lead to \code{GC}s of greater than 1 and works by rescaling by the maximum possible distance that could be observed based on the number of characters codable in each pairwise comparison meaning all resulting distances are on a zero to one scale, and 4. The Generalised Euclidean Distance - this was introduced by Wills (1998) as a means of correcting for the fact that a \code{RED} metric will become increasingly non-Euclidean as the amount of missing data increases and works by filling in missing distances (for characters that are coded as missing in at least one taxon in the pairwise comparison) by using the mean pairwise dissimilarity for that taxon pair as a substitute. In effect then, \code{RED} makes no consideration of missing data, \code{GC} and \code{MORD} normalise by the available data (and are identical if there are no ordered multistate characters), and \code{GED} fills in missing distances by extrapolating from the available data.
#'
#' Note that Lloyd (2016) misidentified the substitute dissimilarity for the \code{GED} as the mean for the whole data set (Hopkins and St John 2018) and this was the way the GED implementation of Claddis operated up to version 0.2. This has now been amended (as of version 0.3) so that the function produces the \code{GED} in the form that Wills (1998) intended. However, this implementation can still be accessed as the \code{Legacy} option for \code{GEDType}, with \code{Wills} being the WIlls (1998) implementation. An advantage of this misinterpreted form of \code{GED} is that it will always return a complete pairwise distance matrix, however it is not recommended (see Lloyd 2016). Instead a third option for \code{GEDType} - (\code{Hybrid}) - offers the same outcome but only uses the mean distance from the entire matrix in the case where there are no codable characters in common in a pairwise comparison. This new hybrid option has not been used in a published study.
#'
#' Typically the resulting distance matrix will be used in an ordination procedure such as principal coordinates (effectively classical multidimensional scaling where k, the number of axes, is maximised at N - 1, where N is the number of rows (i.e., taxa) in the matrix). As such the distance should be - or approximate - Euclidean and hence a square root transformation is typically applied (\code{TransformDistances} with the \code{sqrt} option). However, if applying pre-ordination (i.e., ordination-free) disparity metrics (e.g., weighted mean pairwise distance) you may wish to avoid any transformation (\code{none} option). In particular the \code{MORD} will only fall on a zero to one scale if this is the case. However, if transforming the \code{MORD} for ordination this zero to one property may mean the arcsine square root (\code{arcsine_sqrt} option) is preferred. (Note that if using only unordered multistate or binary characters and the \code{GC} the zero to one scale will apply too.)
#'
#' An unexplored option in distance matrix construction is how to deal with polymorphisms (Lloyd 2016). Up to version 0.2 of Claddis all polymorphisms were treated the same regardless of whether they were true polymorphisms (multiple states are observed in the taxon) or uncertainties (multiple, but not all states, are posited for the taxon). Since version 0.3, however, these two forms can be distinguished by using the different #NEXUS forms (Maddison et al. 1997), i.e., (01) for polymorphisms and \{01\} for uncertainties and within Claddis these are represented as 0&1 or 0/1, respectively. Thus, since 0.3 Claddis allows these two forms to be treated separately, and hence differently (with \code{PolymorphismBehaviour} and \code{UncertaintyBehaviour}). Again, up to version 0.2 of Claddis no options for polymorphism behaviour were offered, instead only a minimum distance was employed. I.e., the distance between a taxon coded 0&1 and a taxon coded 2 would be the smaller of the comparisons 0 with 2 or 1 with 2. Since version 0.3 this is encoded in the \code{min.difference} option. Currentlly two alternatives (\code{mean.difference} and \code{random}) are offered. The first takes the mean of each possible difference and the second simply samples one of the states at random. Note this latter option makes the function stochastic and so it should be rerun multiple times (for example, with a \code{for} loop or \code{apply} function). In general this issue (and these options) are not explored in the literature and so no recommendation can be made beyond that users should think carefully about what this choice may mean for their individual data set(s) and question(s).
#'
#' A final consideration is how to deal with inapplicable characters. Up to version 0.2 Claddis treated inapplicable and missing characters the same (as NA values, i.e., missing data). However, since Claddis version 0.3 these can be imported separately, i.e., by using the "MISSING" and "GAP" states in #NEXUS format (Maddison et al. 1997) with the latter typically representing the inapplicable character. These will appear as NA and empty strings (""), respectively in Claddis format. Hopkins and St John (2018) showed how inapplicable characters - typically assumed to represent secondary characters - could be treated in generating distance matrices. These are usually hierarchical in form. E.g., a primary character might record the presence or absence of feathers and a secondary character whether those feathers are symmetric or asymmetric. The latter will generate inapplicable states for taxa without feathers and without correcting for this ranked distances can b incorrect (Hopkins and ST John 2018). Unfortunately, however, the #NEXUS format (Maddison et al. 1997) does not really allow explicit linkage between primary and secondary characters and so at present Claddis can only really operate by assuming any character with at least one taxon encoded as inapplicable as secondary and all other characters as primary.
#'
#' @param morph.matrix A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @param Distance The distance metric to use. Must be one of \code{"GC"}, \code{"GED"}, \code{"RED"}, or \code{"MORD"} (the default).
#' @param GEDType The type of GED to use. Must be one of \code{"Legacy"}, \code{"Hybrid"}, or \code{"Wills"} (the default). See details for an explanation.
#' @param TransformDistances Whether to transform the distances. Options are \code{"none"}, \code{"sqrt"}, or \code{"arcsine_sqrt"} (the default). (Note: this is only really appropriate for the proportional distances, i.e., "GC" and "MORD".)
#' @param PolymorphismBehaviour The distance behaviour for dealing with polymorphisms. Must be one of \code{"mean.difference"} or \code{"min.difference"} (the default.
#' @param UncertaintyBehaviour The distance behaviour for dealing with uncertainties. Must be one of \code{"mean.difference"} or \code{"min.difference"} (the default.
#'
#' @description
#'
#' Takes a cladistic morphological dataset and converts it into a set of pairwise distances.
#'
#' @return
#'
#' \item{DistanceMetric}{The distance metric used.}
#' \item{DistanceMatrix}{The distance matrix returned.}
#' \item{ComparableCharacterMatrix}{The matrix of characters that can be compared for each pairwise distance.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Thomas Guillerme \email{guillert@@tcd.ie}
#'
#' @references
#'
#' Gower, J. C., 1971. A general coefficient of similarity and some of its properties. Biometrics 27, 857–871.
#'
#' Hopkins, M. J. and St John, K., 2018. A new family of dissimilarity metrics for discrete character matrices that include inapplicable characters and its importance for disparity studies. Proceedings of the Royal Society of London B, 285, 20181784.
#'
#' Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. Biological Journal of the Linnean Society, 118, 131-151.
#'
#' Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. Systematic Biology, 46, 590-621.
#'
#' Wills, M. A., 1998. Crustacean disparity through the Phanerozoic: comparing morphological and stratigraphic data. Biological Journal of the Linnean Society, 65, 455-500.
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
MorphDistMatrix <- function(morph.matrix, Distance = "MORD", GEDType = "Wills", TransformDistances = "arcsine_sqrt", PolymorphismBehaviour = "min.difference", UncertaintyBehaviour = "min.difference") {
  
  # ADD HOPKINS SUGGESTION FOR FOURTH GEDTYPE WHERE MEAN DISTANCE FOR CHARACTER REPLACES MISSING VALUES.
  # CHECK POLYMORPHISM UNCERTAINTY IN GENERAL AS NOT CLEAR IT IS DOING WHAT IT SHOULD DO.
  
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
  
  # Check input of TransformDistances is valid and stop and warn if not:
  if(length(setdiff(TransformDistances, c("arcsine_sqrt", "none", "sqrt"))) > 0) stop("TransformDistances must be one of \"none\", \"sqrt\", or \"arcsine_sqrt\".")
  
  # Check input of distance is valid and stop and warn if not:
  if(length(setdiff(Distance, c("RED", "GED", "GC", "MORD"))) > 0) stop("Distance must be one or more of \"RED\", \"GED\", \"GC\", or \"MORD\".")
  
  # Check input of GED type is valid and stop and warn if not:
  if(length(setdiff(GEDType, c("Legacy", "Hybrid", "Wills"))) > 0) stop("GEDType must be one or more of \"Legacy\", \"Hybrid\", or \"Wills\".")
  
  # Check input for PolymorphismBehaviour is valid and stop and warn if not:
  if(length(setdiff(PolymorphismBehaviour, c("mean.difference", "min.difference"))) > 0) stop("PolymorphismBehaviour must be one or more of \"mean.difference\", \"min.difference\", or \"random\".")
  
  # Check input for UncertaintyBehaviour is valid and stop and warn if not:
  if(length(setdiff(UncertaintyBehaviour, c("mean.difference", "min.difference"))) > 0) stop("UncertaintyBehaviour must be one or more of \"mean.difference\", \"min.difference\", or \"random\".")
  
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
  
  # If PolymorphismBehaviour is to randomly sample one state:
  if(PolymorphismBehaviour == "random") {
    
    # Find cells with polymorphisms:
    PolymorphismCells <- grep("&", morph.matrix)
    
    # If there are polymorphisms randomly sample one value and store:
    if(length(PolymorphismCells) > 0) morph.matrix[PolymorphismCells] <- unlist(lapply(as.list(morph.matrix[PolymorphismCells]), function(x) sample(strsplit(x, split = "&")[[1]], size = 1)))
    
    # Reset behaviour as mean difference to allow it to interact correctly with UncertaintyBehaviour later:
    PolymorphismBehaviour <- "mean.difference"
    
  }
  
  # If UncertaintyBehaviour is to randomly sample one state:
  if(UncertaintyBehaviour == "random") {
    
    # Find cells with uncertainties:
    UncertaintyCells <- grep("/", morph.matrix)
    
    # If there are uncertainties randomly sample one value and store:
    if(length(UncertaintyCells) > 0) morph.matrix[UncertaintyCells] <- unlist(lapply(as.list(morph.matrix[UncertaintyCells]), function(x) sample(strsplit(x, split = "/")[[1]], size = 1)))
    
    # Reset behaviour as mean difference to allow it to interact correctly with PolymorphismBehaviour later:
    UncertaintyBehaviour <- "mean.difference"
    
  }
  
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
    
    # Now build into matrix of pairwise comparisons (odds to be compared with adjacent evens):
    GED.data <- matrix(data = (unlist(GED.data)), ncol = ncol(morph.matrix), byrow = TRUE)
    
    # Calculate single weighted mean univariate distance for calculating GED Legacy or Hybrid (after equation 2 in Wills 2001):
    if(GEDType != "Wills") NonWills_S_ijk_bar <- rep(sum(unlist(diffs)) / sum(unlist(maxdiffs)), length.out = length(diffs))
    
    # Calculate individual pairwise weighted mean univariate distance for calculating GED Hybrid or Wills (after equation 2 in Wills 2001):
    if(GEDType != "Legacy") {
      
      # Generate individual mean pairwise distance for each comparison:
      NonLegacy_S_ijk_bar <- unlist(lapply(diffs, sum)) / unlist(lapply(maxdiffs, sum))
      
      # Find NaNs (divide by zero errors for when there are no characters in common in a pairwsie comparison):
      NaNs <- which(is.nan(NonLegacy_S_ijk_bar))
      
      # If usings WIlls version replace NaNs with NA:
      if(GEDType == "Wills" && length(NaNs) > 0) NonLegacy_S_ijk_bar[NaNs] <- NA
      
      # If using Hybrid replace NaNs with single global mean distance value:
      if(GEDType == "Hybrid" && length(NaNs) > 0) NonLegacy_S_ijk_bar[NaNs] <- NonWills_S_ijk_bar[NaNs]
      
      # Set modified non-Legacy S_ijk_bar as main S_ijk_bar:
      S_ijk_bar <- NonLegacy_S_ijk_bar
    
    }
    
    # If using Legacy set NonWills_S_ijk_bar as main S_ijk_bar:
    if(GEDType == "Legacy") S_ijk_bar <- NonWills_S_ijk_bar
    
    # For each set of differences:
    for(i in seq(from = 1, to = nrow(GED.data) - 1, length.out = length(diffs))) {
      
      # Find missing distances (if any):
      MissingDistances <- which(is.na(GED.data[i, ]))
      
      # Replace missing distances with S_ijk_bar (i.e., results of equation 2 in Wills 2001 into equation 1 of Wills 2001):
      if(length(MissingDistances) > 0) GED.data[i, MissingDistances] <- S_ijk_bar[ceiling(i / 2)]
      
    }

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
  
  # If using a proportional distance:
  if(Distance == "MORD" || Distance == "GC") {
    
    # If transforming distance matrix by taking the square root - take the square root:
    if(TransformDistances == "sqrt") dist.matrix <- sqrt(dist.matrix)
    
    # If transforming distance matrix by taking the arcsine square root:
    if(TransformDistances == "arcsine_sqrt") {
      
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
    
  }
  
  # Compile results as a list:
  result <- list(Distance, dist.matrix, comp.char.matrix)
  
  # Add names to list:
  names(result) <- c("DistanceMetric", "DistanceMatrix", "ComparableCharacterMatrix")
  
  # Output result:
  return(result)
  
}








###################Hopkins and St. John 2018#####################
##function for applying proposed family of metrics. Input: character matrix in Claddis object format; table showing hierarchical relationships between characters; alpha value.
##Input required:
##matr = character matrix in Claddis object format (see "ReadMorphNexus")
##type = table with three columns with headings "Char", "Type", and "Sub"
##      "Char" is the character number
##      "Type" is the type of character, either "S" for secondary or "T" for tertiary.  In our experience, quaternary characters are rarely to never used
##      "Sub" is the character numeber that the secondary or tertiary character is contingent on
##alpha = desired alpha value.  Should vary between 0 and 1.
##      If alpha = 0, only primary characters will contribute to dissimilarity estimate
##      If alpha = 1, then shared primary characters will be weighted by the secondary characters shared; same for secondary characters with tertiary characters
##      See text for further explanation.

##note that delta{ijk} = 1 when characters are comparable, consistent with Gower's initial formulation. It may be preferable to replace the total number of
#comparable characters with the maximum character total possible (Lloyd's MORB metric) especially if there are ordered characters. This is not implemented here
#but could be done by modifying lines 90 and 212 for the new metric and Gower's, respectively
#################

alpha.coefficient <- function(matr, type, alpha) {
  
  #matr <- Day2016
  #alpha <- 0.5 # Check values make sense if incorporating in function
  #type <- matrix(c(8, "S", 7), nrow = 1, dimnames = list(c(), c("Char", "Type", "Sub")))
  
  
  #Type should more logically be a two-column matrix with first character dependent on second character then can do any level of depth (sevndary tertiary etc.)
  #Will require checking for non-circularity of character-dependency
  #Need to do operation in order from most to least "nested" character
  #Does weighting cause issue later though? I..e, clash with satrting weights?
  
  Matrix <- do.call(cbind, lapply(matr[2:length(matr)], function(x) x$Matrix))
  Ordering <- unname(unlist(lapply(matr[2:length(matr)], function(x) x$Ordering)))
  MaxVals <- unname(unlist(lapply(matr[2:length(matr)], function(x) x$MaxVals)))
  Weights <- unname(unlist(lapply(matr[2:length(matr)], function(x) x$Weights)))
  
  pairs <- combn(nrow(Matrix), 2)
  
  HSJ <- matrix(NA, nrow = nrow(Matrix), ncol = nrow(Matrix))
  
  # This does nothing????
  alpha <- alpha
  
  for(i in 1:ncol(pairs)) {
    
    sim.temp <- matrix(abs(suppressWarnings(as.numeric(Matrix[pairs[1, i], ])) - suppressWarnings(as.numeric(Matrix[pairs[2, i], ]))))
    
    #correct values for ordered multistate characters
    if(any(Ordering == 'ord')) {
      
      m.o <- which(Ordering == 'ord')
      
      sim.temp[m.o] <- abs(suppressWarnings(as.numeric(Matrix[pairs[1, i], ][m.o])) - suppressWarnings(as.numeric(Matrix[pairs[2, i], ][m.o]))) #/MaxVals[m.o]
      
    }
    
    #correct values for unordered multistate characters
    if(any(Ordering == "unord" && MaxVals > 1)) {
      
      m.u <- which(Ordering == "unord" && MaxVals > 1)
      
      sim.temp[m.u] <- replace(sim.temp[m.u], which(sim.temp[m.u] > 1), 1)
      
    }
    #correct polymorphic characters
    if(length(grep("&", unique(c(Matrix[pairs[1, i], ], Matrix[pairs[2, i], ])))) > 0) {
      
      polym.states <- sort(c(grep("&", Matrix[pairs[1, i], ]), grep("&", Matrix[pairs[2, i], ])))
      for(v in 1:length(polym.states)) {
        char.state1 <- strsplit(Matrix[pairs[1, i], ][polym.states[v]], "&")[[1]]
        char.state2 <- strsplit(Matrix[pairs[2, i], ][polym.states[v]], "&")[[1]]
        int.value <- intersect(char.state1, char.state2)
        if(length(int.value) > 0) sim.temp[polym.states[v]] <- 0
        if(length(int.value) == 0 && anyNA(c(char.state1, char.state2)) == FALSE) {
          if(Ordering[polym.states[v]] == "unord") sim.temp[polym.states[v]] <- 1
          if(Ordering[polym.states[v]] == "ord"){
            pairs.poly <- matrix(0, nrow = length(char.state1), ncol = length(char.state2))
            for(m in 1:length(char.state1)) for(n in 1:length(char.state2)) pairs.poly[m, n] <- abs(as.numeric(char.state1[m]) - as.numeric(char.state2[n]))
            sim.temp[polym.states[v]] <- min(pairs.poly)
          }
        }
      }
    }
    #account for Tertiary characters (by weighting)
    if(any(type[, "Type"] == "T")) {
      te <- which(type[, "Type"] == "T")
      for(j in 1:length(unique(type[, "Sub"][te]))) {
        te.sub <- which(type[, "Sub"] == as.character(unique(type[, "Sub"][te])[j]))
        te.sim <- sim.temp[te.sub]
        s.te <- as.numeric(as.character(unique(na.omit(type[, "Sub"][te])))[j])
        if(length(na.omit(te.sim)) > 0) {
          sim.temp[s.te] <- 1 - (alpha * (1 - (sum(na.omit(te.sim)) / length(na.omit(te.sim)))) + (1 - alpha))  #alpha is applied as if it were a similarity measure, therefore the distance is converted to a similarity by substracted from one, then the whole thing is substracted from one to convert back to dissimilarity
          sim.temp[te.sub] <- NA
        }
      }
    }
    #account for secondary characters (by weighting)
    
    # If there are secondary characters:
    if(any(type[, "Type"] == "S")) {
      
      # Rows with secondary characters:
      s <- which(type[, "Type"] == "S")
      
      # For each dependent character (i.e., unique primary character):
      for(j in 1:length(unique(type[s, "Sub"]))) {
        
        # This seems wrong? Should be the secondary characters?:
        s.sub <- which(type[, "Sub"] == as.character(unique(type[s, "Sub"])[j]))
        
        #
        s.sim <- sim.temp[s.sub]
        
        # Primary character number:
        p.s <- as.numeric(as.character(unique(na.omit(type[s, "Sub"])))[j])
        
        #
        if(length(na.omit(s.sim)) > 0) {
          
          # Set primary character as reweighted distance value:
          sim.temp[p.s] <- 1 - (alpha * (1 - (sum(na.omit(s.sim)) / length(na.omit(s.sim)))) + (1 - alpha))
          
          # Set secondary distances to NA?:
          sim.temp[s.sub] <- NA
          
        }
        
      }
      
    }
    
    #calculate total dissimilarity
    wt.comp.char <- sum(na.omit(cbind(sim.temp, Weights))[, 2])
    
    
    HSJ[pairs[1, i], pairs[2, i]] <- HSJ[pairs[2, i], pairs[1, i]] <- sum(na.omit(sim.temp * Weights)) / wt.comp.char
  
  }
  
  diag(HSJ) <- 0
  
  return(HSJ)
  
}



###############
##Implementation of Will's GED
##matr = character matrix in Claddis object format (see "ReadMorphNexus")

GED <- function(matr) {
  
  Matrix <- do.call(cbind, lapply(matr[2:length(matr)], function(x) x$Matrix))
  Ordering <- unname(unlist(lapply(matr[2:length(matr)], function(x) x$Ordering)))
  MaxVals <- unname(unlist(lapply(matr[2:length(matr)], function(x) x$MaxVals)))
  Weights <- unname(unlist(lapply(matr[2:length(matr)], function(x) x$Weights)))
  
  pairs <- combn(nrow(Matrix), 2)
  
  GED <- matrix(NA, nrow = nrow(Matrix), ncol = nrow(Matrix))
  
  for (i in 1:ncol(pairs)) {
    sim.temp <- matrix(abs(suppressWarnings(as.numeric(Matrix[pairs[1, i], ])) - suppressWarnings(as.numeric(Matrix[pairs[2, i], ]))))
    #correct values for ordered multistate characters
    if(any(Ordering == "ord")) {
      m.o <- which(Ordering == "ord")
      sim.temp[m.o] <- abs(suppressWarnings(as.numeric(Matrix[pairs[1, i], ][m.o])) - suppressWarnings(as.numeric(Matrix[pairs[2, i], ][m.o])))
    }
    #correct values for unordered multistate characters
    if(any(Ordering == "unord" && MaxVals > 1)){
      m.u <- which(Ordering == "unord" && MaxVals > 1)
      sim.temp[m.u] <- replace(sim.temp[m.u], which(sim.temp[m.u] > 1), 1)
    }
    #correct polymorphic characters
    if(length(grep("&", unique(c(Matrix[pairs[1, i], ], Matrix[pairs[2, i], ])))) > 0){
      polym.states <- sort(c(grep("&", Matrix[pairs[1, i], ]), grep("&", Matrix[pairs[2, i], ])))
      for(v in 1:length(polym.states)) {
        char.state1 <- strsplit(Matrix[pairs[1, i], ][polym.states[v]], "&")[[1]]
        char.state2 <- strsplit(Matrix[pairs[2, i], ][polym.states[v]], "&")[[1]]
        int.value <- intersect(char.state1, char.state2)
        if(length(int.value) > 0) sim.temp[polym.states[v]] <- 0
        if(length(int.value) == 0 & anyNA(c(char.state1,char.state2)) == FALSE){
          if(Ordering[polym.states[v]] == "unord") sim.temp[polym.states[v]] <- 1
          if(Ordering[polym.states[v]] == "ord") {
            pairs.poly <- matrix(0, nrow = length(char.state1), ncol = length(char.state2))
            for(m in 1:length(char.state1)) for(n in 1:length(char.state2)) pairs.poly[m, n] <- abs(as.numeric(char.state1[m]) - as.numeric(char.state2[n]))
            sim.temp[polym.states[v]] <- min(pairs.poly)
          }
        }
      }
    }
    
    #replace missing values
    Sij <- mean(na.omit(sim.temp))
    sim.temp[which(is.na(sim.temp) == TRUE)] <- Sij
    
    #calculate total dissimilarity
    #wt.comp.char<-sum(na.omit(cbind(sim.temp, Weights))[,2])
    GED[pairs[1, i], pairs[2, i]] <- GED[pairs[2, i], pairs[1, i]] <- sqrt(sum((sim.temp * Weights)^2))
  }
  diag(GED)<-0
  return(GED)
}

#DayDiscrete <- MatrixPruner(Day2016, characters2prune = 1:3)
#DayDiscrete$Matrix_2$Matrix[, 30] <- gsub("2", "1", gsub("1", "0", DayDiscrete$Matrix_2$Matrix[, 30]))
#DayDiscrete$Matrix_2$Matrix[, 5] <- gsub("", NA, DayDiscrete$Matrix_2$Matrix[, 5])
#DayDiscrete$Matrix_2$MinVals[30] <- 0
#DayDiscrete$Matrix_2$Matrix[which(DayDiscrete$Matrix_2$Matrix[, 5] == ""), 5] <- NA

#x <- as.vector(as.dist(GED(DayDiscrete)))
#y <- as.vector(as.dist(MorphDistMatrix(morph.matrix = DayDiscrete, Distance = "GED", GEDType = "Wills", TransformDistances = "none", PolymorphismBehaviour = "min.difference", UncertaintyBehaviour = "min.difference")$DistanceMatrix))
#cor.test(x,y)$estimate

#x <- as.vector(as.dist(GED(Michaux1989)))
#y <- as.vector(as.dist(MorphDistMatrix(morph.matrix = Michaux1989, Distance = "GED", GEDType = "Wills", TransformDistances = "none", PolymorphismBehaviour = "min.difference", UncertaintyBehaviour = "min.difference")$DistanceMatrix))
#cor.test(x,y)$estimate

#x <- as.vector(as.dist(GED(Gauthier1986)))
#y <- as.vector(as.dist(MorphDistMatrix(morph.matrix = Gauthier1986, Distance = "GED", GEDType = "Wills", TransformDistances = "none", PolymorphismBehaviour = "min.difference", UncertaintyBehaviour = "min.difference")$DistanceMatrix))
#cor.test(x,y)$estimate



#sim.temp[s.te] <- 1 - (alpha * (1 - (sum(na.omit(te.sim)) / length(na.omit(te.sim)))) + (1 - alpha))  #alpha is applied as if it were a similarity

