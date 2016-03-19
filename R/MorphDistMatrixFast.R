#TG: The basic idea is that now everything is done at the same time successively, i.e. in your version, every step is calculated in a serial way (get the comparable characters for i and j, then calculate the distance between i and j etc...)
#TG: as in this version, every step is calculated at the same time successively (get all the comparable characters, calculate all the distances etc...). I feel this could still be optimised but for that the best solution would be to start again from scratch.
#TG: I tried as much as possible to keep your own structure so you can find your way easily through the code.
#TG: to highlight your steps more clearly, I used the following commenting
  #~~~~~~~
  # The original description of your step
#TG: See the end of the script for performance testing

MorphDistMatrix.fast <- function(morph.matrix, distance = c("Raw", "GED", "Gower", "Max", "Comp"), transform.proportional.distances="arcsine_sqrt") {
  # Check format of transform.proportional.distances:
  if(transform.proportional.distances != "none" && transform.proportional.distances != "sqrt" && transform.proportional.distances != "arcsine_sqrt") {

    # Give error if something other than three possible settings is given:
    stop("ERROR: transform.proportional.distances must be one of \"none\", \"sqrt\", or \"arcsine_sqrt\".")

  }

  #TG: Sanitizing the distance argument
  all_distances <- c("Raw", "GED", "Gower", "Max", "Comp")
  if(all(is.na(match(distance, all_distances)))) stop("distance must be one or more of the following: ", paste(all_distances, collapse=", "), ".", sep="")
  #TG: Note that when distance != Raw, raw distance is still calculated but not exported.

  #~~~~~~~
  # Isolate ordering element of morphology matrix:
  ordering <- morph.matrix$ordering
  
  #~~~~~~~  
  # Isolate max values element of morphology matrix:
  max.vals <- morph.matrix$max.vals

  #~~~~~~~  
  # Isolate min values element of morphology matrix:
  min.vals <- morph.matrix$min.vals

  #~~~~~~~
  # Isolate weighting element of morphology matrix:
  weights <- morph.matrix$weights

  #~~~~~~~  
  # Isolate character-taxon matrix element of morphology matrix:
  morph.matrix <- morph.matrix$matrix

  #~~~~~~~
  # Create empty vectors to store S and W value for Wills 2001 equations (1 and 2):
  #differences <- maximum.differences <- vector(mode="numeric")
  #TG: This line is not necessary any more (see below sections "Add to maximum differences (S_ijk * W_ijk in equation 1 of Wills 2001)" and "Add to maximum differences (S_ijk_max * W_ijk in equation 1 of Wills 2001)")
    
  #~~~~~~~
  # Distance matrices for storing:
  #comp.char.matrix <- gower.dist.matrix <- max.dist.matrix <- dist.matrix <- matrix(0, nrow=length(rownames(morph.matrix)), ncol=length(rownames(morph.matrix)))
  #TG: This line is not necessary any more (see below when storing the matrices)
    
  #~~~~~~~
  # Fill comparable characters diagonal:
  #TG: I've moved this section way down below in the section when you are storing the results ("Store N comparable characters:")

  #~~~~~~~
  # Set up empty matrix for storing data to calculate the Generalised Euclidean Distance of Wills (2001):
  #GED.data <- matrix(nrow=0, ncol=ncol(morph.matrix))
  #TG: This line is not necessary any more (see below when storing the matrices)
  
  #TG: Creating the list of possible combinations between characters (the total number of differences)
  comparisons <- combn(1:nrow(morph.matrix), 2)
  #TG: The idea here is to deal with lists rather than a classic for-in-loop.
  #TG: However, because this functions is about distance matrices, we do not to apply all the inner functions to all the cells of the matrices, it's what you do cleverly with you for i and j loop.
  #TG: The function combn gives us the number of possible unique pairwise combinations between the matrix rows (the pairwise distances!).
  #TG: We can the base all the inner apply function on this "list"
  #TG: (that's represented as a matrix for making it easier in some functions but you can imagine it as a list where each column are the pairwise combinations and each row represent the two rows concerned per pairwise comparison).
  #TG: The idea was to replace all your i's and j's by lists to be able to work with apply class functions (mainly mapply).

  #~~~~~~~
  # Get just the comparable characters (those coded for both taxa):
  #compchar <- intersect(which(!is.na(morph.matrix[rownames(morph.matrix)[i], ])), which(!is.na(morph.matrix[rownames(morph.matrix)[j], ])))
  #TG: or in the new "version"
  #TG: where "interest.col[[1]]" is rownames(morph.matrix)[i] and "interest.col[[2]]" is rownames(morph.matrix)[j]
  get.comp.char <- function(interest.col, morph.matrix) {
    output <- intersect(which(!is.na(morph.matrix[ interest.col[[1]], ])), which(!is.na(morph.matrix[ interest.col[[2]], ])))
    return(list(output))
  }
  #TG: or, more conveniently we can feed taxa1 and taxa2 in the comparisons apply style loop
  list.of.compchar <- unlist(apply(comparisons, 2, get.comp.char, morph.matrix), recursive = FALSE)
  #TG: we now have instead directly a list of comparable characters for all differences
  #TG: Note that this list/unlist (without recursive) is just to make sure the output format is consistent
  #TG: apply loops can sometimes be funny and make their own output style depending on the size of the input.

  #TG: now a function for getting the comparable characters for the ith and jth taxa
  get.rows <- function(interest.col, morph.matrix) {
    row1 <- morph.matrix[rownames(morph.matrix)[interest.col[[1]]],]
    row2 <- morph.matrix[rownames(morph.matrix)[interest.col[[2]]],]
    return(list(row1,row2))
  }

  #TG: Get rows pairs
  rows.pairs <- apply(comparisons, 2, get.rows, morph.matrix)

  #TG: only keep the comparable characters within each rows pairs
  keep.comparable <- function(row.pair, comparable.characters) {
    row.pair[[1]] <- row.pair[[1]][comparable.characters]
    row.pair[[2]] <- row.pair[[2]][comparable.characters]
    return(row.pair)
  }

  #TG: now we can apply this function to the two lists (outputing a matrix of characters to compare)
  matrix.of.char.comp <- mapply(keep.comparable, rows.pairs, list.of.compchar)

  #~~~~~~~
  # Deal with polymorphic characters (if present):
  deal.with.poly <- function(comparisons, comparable.characters) {
    #TG: where comparisons are the columns of matrix.of.char.comp

    #TG: Setting the arguments as in the original function
    firstrow <- comparisons[[1]]
    secondrow <- comparisons[[2]]
    compchar <- comparable.characters

    if(length(grep("&", unique(c(firstrow, secondrow)))) > 0) {  
      #~~~~~~~~
      # Find ampersands (polymorphisms):
      ampersand.elements <- sort(c(grep("&", firstrow), grep("&", secondrow)))
      #~~~~~~~~
      # Go through each polymorphic character:
      go.through.poly.char <- function(amper.element, firstrow, secondrow) {
        #~~~~~~~~
        # Find out if two codings overlap once all polymorphism resolutions are considered:
        intersection.value <- intersect(strsplit(firstrow[amper.element], "&")[[1]], strsplit(secondrow[amper.element], "&")[[1]])
        #~~~~~~~~
        # Case if polymorphic and non-polymorphic values overlap:
        if(length(intersection.value) > 0) {
          #~~~~~~~~     
          # Set ith value as zero (no difference):
          firstrow[amper.element] <- 0
          #~~~~~~~~
          # Set jth value as zero (no difference)
          secondrow[amper.element] <- 0
        }
        #~~~~~~~~
        # Case if polymorphic and non-polymorphic values do not overlap:
        if(length(intersection.value) == 0) {
          #~~~~~~~~         
          # Case if character is unordered (max difference is 1):
          if(ordering[compchar[amper.element]] == "unord") {
            #~~~~~~~~
            # Set ith value as zero:
            firstrow[amper.element] <- 0
            #~~~~~~~~
            # Set jth value as 1 (making the ij difference equal to one):
            secondrow[amper.element] <- 1
          }
          #~~~~~~~~           
          # Case if character is ordered (max difference is > 1):
          if(ordering[compchar[amper.element]] == "ord") {
            #~~~~~~~~
            # Get first row value(s):
            firstrowvals <- as.numeric(strsplit(firstrow[amper.element], "&")[[1]])
            #~~~~~~~~
            # Get second row value(s):
            secondrowvals <- as.numeric(strsplit(secondrow[amper.element], "&")[[1]])
            #~~~~~~~~
            # Make mini distance matrix:
            poly.dist.mat <- matrix(0, nrow=length(firstrowvals), ncol=length(secondrowvals))
            #~~~~~~~~
            # Go through each comparison:
            for(l in 1:length(firstrowvals)) {
            #~~~~~~~~
              # Record absolute difference:
              for(m in 1:length(secondrowvals)) poly.dist.mat[l, m] <- sqrt((firstrowvals[l] - secondrowvals[m]) ^ 2)
            }
            #TG: This loop can be optimised as well but there is not many scenarios were it has to run.
            #TG: I leave it like it is for now
            #~~~~~~~~
            # Set first value as zero:
            firstrow[amper.element] <- 0
            #~~~~~~~~
            # Set second value as minimum possible difference:
            secondrow[amper.element] <- min(poly.dist.mat)
          }
        }
        #TG: note that this functions doesn't return anything per se but just modifies the firstrow and secondrow objects
      }
    
      #TG: applying this function to the ampersand elements
      silent <- lapply(as.list(ampersand.elements), go.through.poly.char, firstrow, secondrow)
    }
    return(list(firstrow, secondrow))
  }
   

  #TG: here we have to make the matrix list into a list to feed into mapply
  #unlist(apply(matrix.of.char.comp, 2, list), recursive=FALSE)
  #TG: dealing with polymorphic characters
  matrix.of.char.comp <- mapply(deal.with.poly, unlist(apply(matrix.of.char.comp, 2, list), recursive=FALSE), list.of.compchar)

  #~~~~~~~ 
  # Get the absolute difference between the two rows:
  get.differences <- function(column) {
    #TG: I'm keeping the first/second row structure here, makes it easier to read:
    firstrow <- column[[1]]
    secondrow <- column[[2]]
    return(list(abs(as.numeric(firstrow) - as.numeric(secondrow))))
  }

  #TG: again, feeding this function to the matrix of characters to compare (by columns)
  raw.diffs <- diffs <- unlist(apply(matrix.of.char.comp, 2, get.differences), recursive = FALSE)

  #~~~~~~~ 
  # If there are differences greater than 1 for unordered characters then rescore as 1:
  replace.unord <- function(differences, compchar, ordering) {
    if(length(which(differences > 1)) > 0) {
      differences[which(differences > 1)[which(ordering[compchar[which(differences > 1)]] == "unord")]] <- 1
    }
    return(list(differences))
  }

  #TG: now we can apply this function to the two lists (outputing a matrix of characters to compare)
  diffs <- mapply(replace.unord, diffs, list.of.compchar, MoreArgs=list(ordering))

  #~~~~~~~ 
  # Find the incomparable characters:
  #TG: This step will be computed internally in the other functions
  find.incomp.char <- function(comparable.characters, morph.matrix) {
    return(setdiff(1:ncol(morph.matrix), comparable.characters))
  }
  # incompchar <- lapply(list.of.compchar, find.incomp.char, morph.matrix)

  #~~~~~~~ 
  # Get weighted differences:
  weight.differences <- function(differences, comparable.characters, weights) {
    return(list(as.numeric(weights[comparable.characters]) * differences))
  }
  diffs <- mapply(weight.differences, diffs, list.of.compchar, MoreArgs=list(weights))

  #~~~~~~~ 
  # Get raw Euclidean distance:
  raw.eucl.dist <- function(differences) {
    dist(rbind(differences, rep(0, length(differences))), method="euclidean")
  }
  raw.dist <- lapply(diffs, raw.eucl.dist)

  #~~~~~~~
  # Work out maximum difference (again, checked against ordering) using compchar characters only:
  max.difference <- function(comparable.characters, max.vals, min.vals) {
    as.numeric(max.vals[comparable.characters]) - as.numeric(min.vals[comparable.characters])
  }
  raw.maxdiffs <- maxdiffs <- lapply(list.of.compchar, max.difference, max.vals, min.vals)

  #~~~~~~~
  # Correct maximum possible differences for unordered characters:
  #~~~~~~~
  # Get vector of maximum differences (corrected for character weights):
  #TG: both ones are easy since the functions are already defined above (we can even integrate them to reduce RAM demand):
  maxdiffs <- mapply(weight.differences, mapply(replace.unord, maxdiffs, list.of.compchar, MoreArgs=list(ordering)), list.of.compchar, MoreArgs=list(weights))

  #~~~~~~~
  # Store raw distance:
  #TG: here we will need a function to transform our list into the distance matrices
  list.to.matrix <- function(list, morph.matrix, diag=NULL) {
    # Set the number of rows
    k <- nrow(morph.matrix)
    # Create the empty matrix
    mat.out <- matrix(ncol=k, nrow=k)
    # Fill up the lower triangle
    mat.out[lower.tri(mat.out)] <- unlist(list)
    # Make the matrix a distance matrix (both triangles have the same values)
    mat.out<-as.matrix(as.dist(mat.out))
    # Fill up the diagonal
    if(is.null(diag)) {
      #If diag argument is null, fill the diagonal as a distance matrix (0)
      diag(mat.out) <- 0
    } else {
      #Else fill the diagonal with the argument
      diag(mat.out) <- diag
    }
    return(mat.out)
  }
  dist.matrix <- list.to.matrix(raw.dist, morph.matrix)
  #TG: note the orginal line here had a repetition (that doesn't take any computing time however, just wanted to point it out in case it's a bug) 
  # dist.matrix[i, j] <- dist.matrix[j, i] <- raw.dist
  #TG: same for the other matrices below.

  # Add row and column names (taxa) to distance matrix:
  rownames(dist.matrix) <- colnames(dist.matrix) <- rownames(morph.matrix)

  #~~~~~~~
  # Store N comparable characters:
  #TG: define the quick apply function to get the diagonal values
  count.comp.char <- function(colum) {
      length(colum) - length(grep(TRUE, is.na(colum)))
  }
  #TG: now we can feed that to the diag argument from list.to.matrix
  comp.char.matrix <- list.to.matrix(lapply(list.of.compchar, length), morph.matrix, diag=apply(morph.matrix, 1, count.comp.char))

  # Add row and column names (taxa) to distance matrix:
  rownames(comp.char.matrix) <- colnames(comp.char.matrix) <- rownames(morph.matrix)
 

  ###########
  #TG: computing the Gower distance
  ###########

  if(any(distance == "Gower")) {
    #~~~~~~~
    # Store Gower distance:
    gower.dist <- function(differences, comparable.characters, weights) {
      sum(differences) / sum(weights[comparable.characters])
    }
    gower.dist.matrix <- list.to.matrix(as.list(mapply(gower.dist, diffs, list.of.compchar, MoreArgs=list(weights))), morph.matrix)

    # Add row and column names (taxa) to distance matrix:
    rownames(gower.dist.matrix) <- colnames(gower.dist.matrix) <- rownames(morph.matrix)
 
   }


  ###########
  #TG: computing the Max distance
  ###########

  if(any(distance == "Max")) {
    #~~~~~~~
    # Store maximum-rescaled distance:
    max.distance <- function(differences, maximum.differences) {
      sum(differences) / sum(maximum.differences)
    }
    max.dist.matrix <- list.to.matrix(mapply(max.distance, diffs, maxdiffs), morph.matrix)

    # Add row and column names (taxa) to distance matrix:
    rownames(max.dist.matrix) <- colnames(max.dist.matrix) <- rownames(morph.matrix)

  }


  ###########
  #TG: computing the GED distance
  ###########

  if(any(distance == "GED")) {
    #~~~~~~~ 
    # Store data for GED with NAs for missing distances:
    #TG: this one is slightly modified, the output is now a list of matrices 2*ncol(morph.matrix) with the first row being the differences and the second row the weighting
    store.data.GED <- function(differences, comparable.characters, morph.matrix, weights) {
      return(rbind(c(differences, rep(NA, length(find.incomp.char(comparable.characters, morph.matrix)))), c(weights[comparable.characters], weights[find.incomp.char(comparable.characters, morph.matrix)])))
    }
    #TG: now with the mapply (I'll apply this scheme again further below so I'll stop writing obvious comments for this structure ;))
    GED.data <- mapply(store.data.GED, diffs, list.of.compchar, MoreArgs=list(morph.matrix, weights), SIMPLIFY = FALSE)
    #TG: output GED.data as a matrix (as in the original function)
    #TG: first let's transpose the matrices to facilitate the output
    GED.data <- lapply(GED.data, t)
    #TG: second, create the matrix
    GED.data <- matrix(data=(unlist(GED.data)), ncol=ncol(morph.matrix), byrow=TRUE)

    #~~~~~~~
    # Add to maximum differences (S_ijk * W_ijk in equation 1 of Wills 2001):
    differences <- unlist(diffs)

    #~~~~~~~
    # Add to maximum differences (S_ijk_max * W_ijk in equation 1 of Wills 2001):
    maximum.differences <- unlist(maxdiffs)
    #TG: Note that therefore the line near the start defining the differences/maximum.differences vectors is not necessary.

    #~~~~~~~
    # Calculated weighted mean univariate distance for calculating GED (equation 2 in Wills 2001):
    S_ijk_bar <- sum(differences) / sum(maximum.differences)

    #~~~~~~~
    # Replace missing distances with S_ijk_bar (i.e., results of equation 2 in Wills 2001 into equation 1 of Wills 2001):
    GED.data[is.na(GED.data)] <- S_ijk_bar

    #~~~~~~~
    # Isolate the distances:
    S_ijk <- GED.data[which((1:nrow(GED.data) %% 2) == 1), ]

    #~~~~~~~
    # Isolate the weights:
    W_ijk <- GED.data[which((1:nrow(GED.data) %% 2) == 0), ]

    #~~~~~~~
    # Calculate the GED (equation 1 of Wills 2001) for each pairwise comparison (ij):
    GED_ij <- sqrt(apply(W_ijk * (S_ijk ^ 2), 1, sum))

    #~~~~~~~
    # Create empty GED distance matrix:
    GED.dist.matrix <- list.to.matrix(as.list(GED_ij), morph.matrix)

    # Add row and column names (taxa) to distance matrix:
    rownames(GED.dist.matrix) <- colnames(GED.dist.matrix) <- rownames(morph.matrix)

  }
  

  # Add row and column names (taxa) to distance matrices:
  #rownames(comp.char.matrix) <- colnames(comp.char.matrix) <- rownames(GED.dist.matrix) <- colnames(GED.dist.matrix) <- rownames(gower.dist.matrix) <- colnames(gower.dist.matrix) <- rownames(max.dist.matrix) <- colnames(max.dist.matrix) <- rownames(dist.matrix) <- colnames(dist.matrix) <- rownames(morph.matrix)
  #TG: This step is now dealt individually for each matrix

  # If transformation option is not "none":
  if(transform.proportional.distances != "none") {

    # If transformation option is "sqrt":
    if(transform.proportional.distances == "sqrt") {

      # Replace NaN with NA for Gower distances and take square root:
      if(any(distance == "Gower")) {
        gower.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, sqrt(gower.dist.matrix))), nrow=nrow(gower.dist.matrix), dimnames=list(rownames(gower.dist.matrix), rownames(gower.dist.matrix)))
      }

      # Replace NaN with NA for Max distances and take square root:
      if(any(distance == "Max")) {
        max.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, sqrt(max.dist.matrix))), nrow=nrow(max.dist.matrix), dimnames=list(rownames(max.dist.matrix), rownames(max.dist.matrix)))
      }

    # If transformation option is "arcsine_sqrt":
    } else {

      if(any(distance == "Gower")) {
        # Establish correction factor to ensure Gower data is proportional:
        gower.correction <- max(c(max(sort(gower.dist.matrix)), 1))

        # Ensure all Gower values are on 0 to 1 scale then take arcsine of sqrt to get values that better approximate a normal distribution:
        gower.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, asin(sqrt(gower.dist.matrix / gower.correction)))), nrow=nrow(gower.dist.matrix), dimnames=list(rownames(gower.dist.matrix), rownames(gower.dist.matrix)))
      }

      if(any(distance == "Max")) {
        # Take arcsine square root of all MOD dist values:
        max.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, asin(sqrt(max.dist.matrix)))), nrow=nrow(max.dist.matrix), dimnames=list(rownames(max.dist.matrix), rownames(max.dist.matrix)))
      }

    }

  # If transformation option is "none":
  } else {

    if(any(distance == "Gower")) {
      # Replace NaN with NA for Gower distances:
      gower.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, gower.dist.matrix)), nrow=nrow(gower.dist.matrix), dimnames=list(rownames(gower.dist.matrix), rownames(gower.dist.matrix)))
    }

    if(any(distance == "Max")) {
      # Replace NaN with NA for Max distances:
      max.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, max.dist.matrix)), nrow=nrow(max.dist.matrix), dimnames=list(rownames(max.dist.matrix), rownames(max.dist.matrix)))
    }

  }

  # Compile results as a list:
  #result <- list(dist.matrix, GED.dist.matrix, gower.dist.matrix, max.dist.matrix, comp.char.matrix)

  result <- list()
  # Adding the raw distance
  if(any(distance == "Raw")) result$raw.dist.matrix <- dist.matrix
  # Adding the GED distance
  if(any(distance == "GED")) result$GED.dist.matrix <- GED.dist.matrix
  # Adding the Gower distance
  if(any(distance == "Gower")) result$gower.dist.matrix <- gower.dist.matrix
  # Adding the max distance
  if(any(distance == "Max")) result$max.dist.matrix <- max.dist.matrix
  # Adding the comparable characters
  if(any(distance == "Comp")) result$comp.char.matrix <- comp.char.matrix
    
  # Add names to results list:
  #names(result) <- c("raw.dist.matrix", "GED.dist.matrix", "gower.dist.matrix", "max.dist.matrix", "comp.char.matrix")

  # Output result:
  return(result)
}

#####################
# Performance testing
#####################

# Loading the packages
#  library(Claddis)
#  library(testthat)

# Testing with the Michaux 1989 matrix
#  time_base <- system.time(results_base <- MorphDistMatrix(Michaux1989))
#  time_fast <- system.time(results_fast <- MorphDistMatrix.fast(Michaux1989))
# What's the time difference?
#  time_base/time_fast

#Is it actually the same output?
#  test_that("Are both version giving the exact same results (no rounding nor random seeds)", {
#    expect_equal(results_base, results_fast)
#  })

# Testing with the Gauthier 1986 matrix
#  time_base <- system.time(results_base <- MorphDistMatrix(Gauthier1986))
#  time_fast <- system.time(results_fast <- MorphDistMatrix.fast(Gauthier1986))
# What's the time difference?
#  time_base/time_fast

#Is it actually the same output?
#  test_that("Are both version giving the exact same results (no rounding nor random seeds)", {
#    expect_equal(results_base, results_fast)
#  })

# Testing with a more complex (bigger) matrix
# Some random matrix
#  matrix <- matrix(sample(c("0", "1", "2", "0&1", "0&2", "1&2", "0&1&2"), 5000, replace=TRUE, prob=c(0.425, 0.42, 0.12, 0.01, 0.01, 0.01, 0.005))
#    , nrow=50, ncol=100, dimnames=list(c(1:50)))
# Adding 25% of missing characters
#  matrix[sample(1:5000, 200)] <- NA

# ordering
#  ordering <- sample(c("unord", "ord"), 100, replace=TRUE, prob=c(0.9,0.1))

# weights
#  weights <- sample(c(1:3), 100, replace=TRUE, prob=c(0.85,0.1,0.05))

# Function for making a morph.matrix like object
#  make.nexus<-function(matrix, header, ordering, weights) {
#      nexus<-list()
#      nexus$header<-header
#      nexus$matrix<-matrix
#      nexus$ordering<-ordering
#      nexus$weights<-weights
#      nexus$max.vals<-apply(matrix, 2, max, na.rm=TRUE)
#      nexus$min.vals<-apply(matrix, 2, min, na.rm=TRUE)
#      return(nexus)
#  }

#  morph.matrix <- make.nexus(matrix, header="example", ordering, weights)

# Testing with the random matrix
#  time_base <- system.time(results_base <- MorphDistMatrix(morph.matrix))
#  time_fast <- system.time(results_fast <- MorphDistMatrix.fast(morph.matrix))
# What's the time difference?
#  time_base/time_fast

#Is it actually the same output?
#  test_that("Are both version giving the exact same results (no rounding nor random seeds)", {
#    expect_equal(results_base, results_fast)
#  })

#####################
# Output testing
#####################

# results_base <- MorphDistMatrix(Michaux1989)
# results_fast <- MorphDistMatrix.fast(Michaux1989)
# results_raw <- MorphDistMatrix.fast(Michaux1989, distance = "Raw")
# results_GED <- MorphDistMatrix.fast(Michaux1989, distance = "GED")
# results_Gower <- MorphDistMatrix.fast(Michaux1989, distance = "Gower")
# results_Max <- MorphDistMatrix.fast(Michaux1989, distance = "Max")
# results_Comp <- MorphDistMatrix.fast(Michaux1989, distance = "Comp")

# # Testing if the algorithm gives the same output with default argument (all matrices)
# test_that("Are both the base and fast algorithms equal?", {
#   expect_equal(results_base, results_fast)
# })

# # Testing if each individual matrix is corrects
# test_that("Are each matrix calculated individually correct?", {
#   expect_equal(results_base$raw.dist.matrix, results_raw$raw.dist.matrix)
#   expect_equal(results_base$GED.dist.matrix, results_GED$GED.dist.matrix)
#   expect_equal(results_base$gower.dist.matrix, results_Gower$gower.dist.matrix)
#   expect_equal(results_base$max.dist.matrix, results_Max$max.dist.matrix)
#   expect_equal(results_base$comp.char.matrix, results_Comp$comp.char.matrix)
# })

