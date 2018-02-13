  # Deal with polymorphic characters (if present):
  deal.with.poly <- function(comparisons, comparable.characters, ordering) {
    #TG: where comparisons are the columns of matrix.of.char.comp

    #TG: Setting the arguments as in the original function
    firstrow <- comparisons[[1]]
    secondrow <- comparisons[[2]]
    compchar <- comparable.characters

    if(length(grep("&", unique(c(firstrow, secondrow)))) > 0) {  
      #~~~~~~~~
      # Find ampersands (polymorphisms) as a logical list:
      #TG: this is now down via a lapply and outputs a list
      grep.ampersand <- function(X) {ifelse(length(grep("&", X) > 0), TRUE, FALSE)}
      ampersand.elements <- list(sapply(firstrow, grep.ampersand), sapply(secondrow, grep.ampersand))
      #TG: select the ambiguous elements
      ambiguous_char <- ampersand.elements[[1]] | ampersand.elements[[2]]
      # #~~~~~~~~
      # # Create a list of elements to solve recursively
      # char_to_solve <- list(firstrow[ambiguous_char],secondrow[ambiguous_char])
      # #TG: select the orderings
      # orders <- ordering[ambiguous_char]


      solve.one.character <- function(first, second, order) {
        #~~~~~~~~
        # Find out if two codings overlap once all polymorphism resolutions are considered:
        first <- as.numeric(strsplit(first, "&")[[1]])
        second <- as.numeric(strsplit(second, "&")[[1]])
        intersection.value <- intersect(first, second)
        #~~~~~~~~
        # Case if polymorphic and non-polymorphic values overlap:
        if(length(intersection.value) > 0) {
          #~~~~~~~~     
          # Set ith and jth value as zero (no difference):
          first <- second <- 0
        } else {
          #~~~~~~~~  
          # Case if polymorphic and non-polymorphic values do not overlap: (length(intersection.value) == 0)
          #~~~~~~~~         
          # Case if character is unordered (max difference is 1):
          if(order == "unord") {
            #~~~~~~~~
            # Set ith value as zero:
            first <- 0
            #~~~~~~~~
            # Set jth value as 1 (making the ij difference equal to one):
            second <- 1
          } else {
            #~~~~~~~~           
            # Case if character is ordered (max difference is > 1): (ordering[compchar[increment]] == "ord")
            #~~~~~~~~
            # Make mini distance matrix:
            #TG: this is a transposed distance matrix with the states for secondrow in column and the ones for firstrow in rows
            poly.dist.mat <- matrix(0, nrow = length(first), ncol = length(second))
            #~~~~~~~~
            # Go through each comparison:
            for(l in 1:length(first)) {
                #~~~~~~~~
                # Record absolute difference:
                for(m in 1:length(second)) poly.dist.mat[l, m] <- sqrt((first[l] - second[m]) ^ 2)
            }
            #~~~~~~~~
            # Set first value as zero:
            first <- 0
            #~~~~~~~~
            # Set second value as minimum possible difference:
            secon <- min(poly.dist.mat)
          }
        }
        return(c(first, second))
      }

      #TG: Looping through the characters to solve
      for(char in which(ambiguous_char)) {
        solved_char <- solve.one.character(firstrow[char], secondrow[char], ordering[char])
        firstrow[char] <- solved_char[1]
        secondrow[char] <- solved_char[2]
      }

      #TG: return the arguments without polymorphism
      return(list(firstrow, secondrow))
    }
  }
