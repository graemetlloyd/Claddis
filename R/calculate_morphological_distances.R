#' Get distance matrices from a cladistic matrix
#'
#' @description
#'
#' Takes a cladistic morphological dataset and converts it into a set of pairwise distances.
#'
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @param distance_metric The distance metric to use. Must be one of \code{"GC"}, \code{"GED"}, \code{"RED"}, or \code{"MORD"} (the default).
#' @param ged_type The type of GED to use. Must be one of \code{"Legacy"}, \code{"Hybrid"}, or \code{"Wills"} (the default). See details for an explanation.
#' @param distance_transformation The type of distance transformation to perform. Options are \code{"none"}, \code{"sqrt"}, or \code{"arcsine_sqrt"} (the default). (Note: this is only really appropriate for the proportional distances, i.e., "GC" and "MORD".)
#' @param polymorphism_behaviour The distance behaviour for dealing with polymorphisms. Must be one of \code{"mean.difference"}, \code{"min.difference"} (the default), or \code{"random"}.
#' @param uncertainty_behaviour The distance behaviour for dealing with uncertainties. Must be one of \code{"mean.difference"}, \code{"min.difference"} (the default), or \code{"random"}.
#' @param inapplicable_behaviour The behaviour for dealing with inapplicables. Must be one of \code{"missing"} (default), or \code{"HSJ"} (Hopkins and St John 2018; see details).
#' @param character_dependencies Only relevant if using \code{inapplicable_behaviour = "HSJ"}. Must be a two-column matrix with colnames "DependentCharacter" and "IndependentCharacter" that specifies character hierarchies. See details.
#' @param alpha The alpha value (sensu Hopkins and St John 2018). Only relevant if using \code{inapplicable_behaviour = "HSJ"}. See details.
#'
#' @details
#'
#' There are many options to consider when generating a distance matrix from morphological data, including the metric to use, how to treat inapplicable, polymorphic (e.g., 0&1), or uncertain (e.g., 0/1) states, and whether the output should be transformed (e.g., by taking the square root so that the distances are - or approximate - Euclidean distances). Some of these issues have been discussed previously in the literature (e.g., Lloyd 2016; Hopkins and St John 2018), but all likely require further study.
#'
#' Claddis currently offers four different distance metrics: 1. Raw Euclidean Distance (\code{RED}) - this is only really applicable if there are no missing data, 2. The Gower Coefficient (\code{GC}; Gower 1971) - this rescales distances by the number of characters that can be coded for both taxa in each pairwise comparison thus correcting for missing data, 3. The Maximum Observable Rescaled Distance (\code{MORD}) - this was introduced by Lloyd (2016) as an extension of the \code{GC} designed to deal with the fact that multistate ordered characters can lead to \code{GC}s of greater than 1 and works by rescaling by the maximum possible distance that could be observed based on the number of characters codable in each pairwise comparison meaning all resulting distances are on a zero to one scale, and 4. The Generalised Euclidean Distance - this was introduced by Wills (1998) as a means of correcting for the fact that a \code{RED} metric will become increasingly non-Euclidean as the amount of missing data increases and works by filling in missing distances (for characters that are coded as missing in at least one taxon in the pairwise comparison) by using the mean pairwise dissimilarity for that taxon pair as a substitute. In effect then, \code{RED} makes no consideration of missing data, \code{GC} and \code{MORD} normalise by the available data (and are identical if there are no ordered multistate characters), and \code{GED} fills in missing distances by extrapolating from the available data.
#'
#' Note that Lloyd (2016) misidentified the substitute dissimilarity for the \code{GED} as the mean for the whole data set (Hopkins and St John 2018) and this was the way the GED implementation of Claddis operated up to version 0.2. This has now been amended (as of version 0.3) so that the function produces the \code{GED} in the form that Wills (1998) intended. However, this implementation can still be accessed as the \code{Legacy} option for \code{ged_type}, with \code{Wills} being the WIlls (1998) implementation. An advantage of this misinterpreted form of \code{GED} is that it will always return a complete pairwise distance matrix, however it is not recommended (see Lloyd 2016). Instead a third option for \code{ged_type} - (\code{Hybrid}) - offers the same outcome but only uses the mean distance from the entire matrix in the case where there are no codable characters in common in a pairwise comparison. This new hybrid option has not been used in a published study.
#'
#' Typically the resulting distance matrix will be used in an ordination procedure such as principal coordinates (effectively classical multidimensional scaling where k, the number of axes, is maximised at N - 1, where N is the number of rows (i.e., taxa) in the matrix). As such the distance should be - or approximate - Euclidean and hence a square root transformation is typically applied (\code{distance_transformation} with the \code{sqrt} option). However, if applying pre-ordination (i.e., ordination-free) disparity metrics (e.g., weighted mean pairwise distance) you may wish to avoid any transformation (\code{none} option). In particular the \code{MORD} will only fall on a zero to one scale if this is the case. However, if transforming the \code{MORD} for ordination this zero to one property may mean the arcsine square root (\code{arcsine_sqrt} option) is preferred. (Note that if using only unordered multistate or binary characters and the \code{GC} the zero to one scale will apply too.)
#'
#' An unexplored option in distance matrix construction is how to deal with polymorphisms (Lloyd 2016). Up to version 0.2 of Claddis all polymorphisms were treated the same regardless of whether they were true polymorphisms (multiple states are observed in the taxon) or uncertainties (multiple, but not all states, are posited for the taxon). Since version 0.3, however, these two forms can be distinguished by using the different #NEXUS forms (Maddison et al. 1997), i.e., (01) for polymorphisms and \{01\} for uncertainties and within Claddis these are represented as 0&1 or 0/1, respectively. Thus, since 0.3 Claddis allows these two forms to be treated separately, and hence differently (with \code{polymorphism_behaviour} and \code{uncertainty_behaviour}). Again, up to version 0.2 of Claddis no options for polymorphism behaviour were offered, instead only a minimum distance was employed. I.e., the distance between a taxon coded 0&1 and a taxon coded 2 would be the smaller of the comparisons 0 with 2 or 1 with 2. Since version 0.3 this is encoded in the \code{min.difference} option. Currentlly two alternatives (\code{mean.difference} and \code{random}) are offered. The first takes the mean of each possible difference and the second simply samples one of the states at random. Note this latter option makes the function stochastic and so it should be rerun multiple times (for example, with a \code{for} loop or \code{apply} function). In general this issue (and these options) are not explored in the literature and so no recommendation can be made beyond that users should think carefully about what this choice may mean for their individual data set(s) and question(s).
#'
#' A final consideration is how to deal with inapplicable characters. Up to version 0.2 Claddis treated inapplicable and missing characters the same (as NA values, i.e., missing data). However, since Claddis version 0.3 these can be imported separately, i.e., by using the "MISSING" and "GAP" states in #NEXUS format (Maddison et al. 1997), with the latter typically representing the inapplicable character. These appear as NA and empty strings (""), respectively, in Claddis format. Hopkins and St John (2018) showed how inapplicable characters - typically assumed to represent secondary characters - could be treated in generating distance matrices. These are usually hierarchical in form. E.g., a primary character might record the presence or absence of feathers and a secondary character whether those feathers are symmetric or asymmetric. The latter will generate inapplicable states for taxa without feathers and without correcting for this ranked distances can be incorrect (Hopkins and St John 2018). Unfortunately, however, the #NEXUS format (Maddison et al. 1997) does not really allow explicit linkage between primary and secondary characters and so this information must be provided separately to use the Hopkins and St John (2018) approach. This is done here with the \code{character_dependencies} option. This must be in the form of a two-column matrix with column headers of "DependentCharacter" and "IndependentCharacter". The former being secondary characters and the latter the corresponding primary character. (Note that characters are to be numbered across the whole matrix from 1 to N and do not restart with each block of the matrix.) If using \code{inapplicable_behaviour = "HSJ"} the user must also provide an \code{alpha} value between zero and one. When \code{alpha = 0} the secondary characters contribute nothing to the distance and when \code{alpha = 1} the primary character is not counted in the weight separately (see Hopkins and St John 2018). The default value (0.5) offers a compromise bteween these two extremes.
#'
#' Here the implementation of this approach differs somewhat from the code available in the supplementary materials to Hopkins and St John (2018). Specifically, this approach is incorporated (and used) regardless of the overriding distance metric (i.e., the \code{Distance} option). Additionally, the Hopkins and St John function specifically allows an extra level of dependency (secondary and tertary characters) with these being applied recursively (tertiary first then secondary). Here, though, additional levels of dependency do not need to be defined by the user as this information is already encoded in the \code{character_dependencies} option. Furthermore, because of this any level of dependency is possible (if unlikely).
#'
#' @return
#'
#' \item{distance_metric}{The distance metric used.}
#' \item{distance_matrix}{The distance matrix returned.}
#' \item{comparable_character_matrix}{The matrix of characters that can be compared for each pairwise distance.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Thomas Guillerme \email{guillert@@tcd.ie}
#'
#' @references
#'
#' Gower, J. C., 1971. A general coefficient of similarity and some of its properties. \emph{Biometrika}, \bold{27}, 857<U+2013>871.
#'
#' Hopkins, M. J. and St John, K., 2018. A new family of dissimilarity metrics for discrete character matrices that include inapplicable characters and its importance for disparity studies. \emph{Proceedings of the Royal Society of London B}, \bold{285}, 20181784.
#'
#' Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. \emph{Biological Journal of the Linnean Society}, \bold{118}, 131-151.
#'
#' Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. \emph{Systematic Biology}, \bold{46}, 590-621.
#'
#' Wills, M. A., 1998. Crustacean disparity through the Phanerozoic: comparing morphological and stratigraphic data. \emph{Biological Journal of the Linnean Society}, \bold{65}, 455-500.
#'
#' @examples
#'
#' # Get morphological distances for the Day et
#' # al. (2016) data set:
#' distances <- calculate_morphological_distances(day_2016)
#'
#' # Show distance metric:
#' distances$DistanceMetric
#'
#' # Show distance matrix:
#' distances$distance_matrix
#'
#' # Show number of characters that can be scored for
#' # each pairwise comparison:
#' distances$comparable_character_matrix
#'
#' # To repeat using the Hopkins and St John approach
#' # we first need to define the character dependency
#' # (here there is only one, character 8 is a
#' # secondary where 7 is the primary character):
#' character_dependencies <- matrix(c(8, 7),
#'   ncol = 2,
#'   byrow = TRUE, dimnames = list(
#'     c(),
#'     c(
#'       "DependentCharacter",
#'       "IndependentCharacter"
#'     )
#'   )
#' )
#'
#' # Get morphological distances for the Day et
#' # al. (2016) data set using HSJ approach:
#' distances <- calculate_morphological_distances(day_2016,
#'   inapplicable_behaviour = "HSJ",
#'   character_dependencies = character_dependencies,
#'   alpha = 0.5
#' )
#'
#' # Show distance metric:
#' distances$DistanceMetric
#'
#' # Show distance matrix:
#' distances$distance_matrix
#'
#' # Show number of characters that can be scored for
#' # each pairwise comparison:
#' distances$comparable_character_matrix
#' @export calculate_morphological_distances
calculate_morphological_distances <- function(cladistic_matrix, distance_metric = "MORD", ged_type = "Wills", distance_transformation = "arcsine_sqrt", polymorphism_behaviour = "min.difference", uncertainty_behaviour = "min.difference", inapplicable_behaviour = "missing", character_dependencies = NULL, alpha = 0.5) {

  # ADD HOPKINS SUGGESTION (VIA EMAIL) FOR FOURTH GEDTYPE WHERE MEAN DISTANCE FOR CHARACTER REPLACES MISSING VALUES.
  # CHECK POLYMORPHISM UNCERTAINTY IN GENERAL AS NOT CLEAR IT IS DOING WHAT IT SHOULD DO.
  # CHECK TRASNFORM IS APPROPRIATE AND WARN USER IF NOT
  # MAYBE ALLOW MANHATTAN TYPE DISTANCES TOO.
  # ADD LEHAMN REFERENCE!
  # ALLOW MANHATTAN DISTANCES

  # Subfunction to find comparable characters for a pairwise taxon comparison:
  find_comparable <- function(interest.col, cladistic_matrix) {

    # Get intersection of characters that are coded for both taxa in a pair:
    output <- intersect(intersect(which(x = !is.na(cladistic_matrix[interest.col[[1]], ])), which(x = cladistic_matrix[interest.col[[1]], ] != "")), intersect(which(x = !is.na(cladistic_matrix[interest.col[[2]], ])), which(x = cladistic_matrix[interest.col[[2]], ] != "")))

    # Return output:
    return(list(output))
  }

  # Subfunction to get character strings for each pair of taxa:
  get_pairwise_strings <- function(interest.col, cladistic_matrix) {

    # Get character states for first taxon in pair:
    row1 <- cladistic_matrix[rownames(x = cladistic_matrix)[interest.col[[1]]], ]

    # Get character states for second taxon in pair:
    row2 <- cladistic_matrix[rownames(x = cladistic_matrix)[interest.col[[2]]], ]

    # Return output as a list:
    return(list(row1, row2))
  }

  # Subfunction to subset pairwise comparisons by just comparable characters:
  subset_by_comparable <- function(row.pair, comparable.characters) {

    # Collapse first row to just comparable characters:
    row.pair[[1]] <- row.pair[[1]][comparable.characters]

    # Collapse second row to just comparable characters:
    row.pair[[2]] <- row.pair[[2]][comparable.characters]

    # Output colapsed row pair:
    return(row.pair)
  }

  # Subfunction to edit polymorphic characters down to a single value:
  edit_polymorphisms <- function(comparisons, comparable.characters, ordering, polymorphism_behaviour, uncertainty_behaviour) {

    # Set first taxon values:
    firstrow <- comparisons[[1]]

    # Set second taxon values:
    secondrow <- comparisons[[2]]

    # If there are any inapplicables:
    if (any(c(firstrow, secondrow) == "")) {

      # Find inapplicable positions:
      InapplicablePositions <- sort(x = unique(x = c(which(x = firstrow == ""), which(x = secondrow == ""))))

      # Find polymorphism and uncertainty positions:
      PolymorphismAndUncertaintyPositions <- sort(x = unique(x = c(grep("/|&", firstrow), grep("/|&", secondrow))))

      # If there are polymorphisms or uncertianties that match up with inapplicables:
      if (length(x = intersect(InapplicablePositions, PolymorphismAndUncertaintyPositions)) > 0) {

        # Find positions where collapsing to a single value is required:
        CollapsePositions <- intersect(InapplicablePositions, PolymorphismAndUncertaintyPositions)

        # Collapse any polymorphisms or uncertianties in first row to just first value:
        firstrow[CollapsePositions] <- unlist(x = lapply(X = strsplit(firstrow[CollapsePositions], split = "/|&"), function(x) ifelse(length(x = x) == 0, "", x[1])))

        # Collapse any polymorphisms or uncertianties in second row to just first value:
        secondrow[CollapsePositions] <- unlist(x = lapply(X = strsplit(secondrow[CollapsePositions], split = "/|&"), function(x) ifelse(length(x = x) == 0, "", x[1])))
      }
    }

    # Set comparable characters:
    compchar <- comparable.characters

    # Set ordering for comparable characters:
    charordering <- ordering[compchar]

    # Only if there are polymorphisms or uncertainties:
    if (length(x = c(grep("&", unique(x = c(firstrow, secondrow))), grep("/", unique(x = c(firstrow, secondrow))))) > 0) {

      # Find ampersands (polymorphisms):
      ampersand.elements <- sort(x = c(grep("&", firstrow), grep("&", secondrow)))

      # Find slashes (uncertianties):
      slash.elements <- sort(x = c(grep("/", firstrow), grep("/", secondrow)))

      # Combine to find all characters to check:
      characters.to.check <- sort(x = unique(x = c(ampersand.elements, slash.elements)))

      # Set behaviours as either the shared version or minimum difference if they contradict (may need to modify this later for more complex options):
      behaviour <- unlist(x = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(apply(apply(rbind(firstrow[characters.to.check], secondrow[characters.to.check]), 2, gsub, pattern = "[:0-9:]", replacement = ""), 2, list), unlist), function(x) x[nchar(x = x) > 0]), function(x) ifelse(nchar(x = x) > 0, strsplit(x, split = "")[[1]][1], x)), function(x) gsub(pattern = "&", replacement = polymorphism_behaviour, x = x)), function(x) gsub(pattern = "/", replacement = uncertainty_behaviour, x = x)), unique), function(x) ifelse(length(x) > 1, "min.difference", x)))

      # If behaviour is to find minimum differences:
      if (any(behaviour == "min.difference")) {

        # Set up minimum difference characters to check:
        min.characters.to.check <- characters.to.check[behaviour == "min.difference"]

        # Find intersecting character states for each character:
        IntersectionCharacter <- lapply(X = lapply(X = lapply(X = lapply(X = apply(rbind(firstrow[min.characters.to.check], secondrow[min.characters.to.check]), 2, strsplit, split = "&|/"), unlist), sort), rle), function(x) x$values[x$lengths > 1][1])

        # If at least one intersecting character state was found:
        if (any(!is.na(unlist(x = IntersectionCharacter)))) {

          # Record rows to update:
          rows.to.update <- which(x = !is.na(unlist(x = IntersectionCharacter)))

          # Store (first) shared state for both taxa:
          firstrow[min.characters.to.check[rows.to.update]] <- secondrow[min.characters.to.check[rows.to.update]] <- unlist(x = IntersectionCharacter)[rows.to.update]

          # Update minimum characters to check:
          min.characters.to.check <- min.characters.to.check[-rows.to.update]
        }

        # Only continue if there are still characters that need to be fixed:
        if (length(x = min.characters.to.check) > 0) {

          # Build two option matrices for every comparison:
          TwoOptionMatrices <- lapply(X = apply(rbind(firstrow[min.characters.to.check], secondrow[min.characters.to.check]), 2, strsplit, split = "&|/"), function(x) rbind(c(min(as.numeric(x[[1]])), max(as.numeric(x[[2]]))), c(max(as.numeric(x[[1]])), min(as.numeric(x[[2]])))))

          # Pick smallest difference as minimum and maximum states:
          MinMaxStates <- lapply(X = lapply(X = lapply(X = TwoOptionMatrices, function(x) x[which(x = abs(apply(x, 1, diff)) == min(abs(apply(x, 1, diff)))), ]), sort), as.character)

          # Set first row values(s):
          firstrow[min.characters.to.check] <- unlist(x = lapply(X = MinMaxStates, "[[", 1))

          # Set second row values(s):
          secondrow[min.characters.to.check] <- unlist(x = lapply(X = MinMaxStates, "[[", 2))
        }
      }

      # If any behaviour is to find mean differences:
      if (any(behaviour == "mean.difference")) {

        # Set up minimum difference characters to check:
        mean.characters.to.check <- characters.to.check[behaviour == "mean.difference"]

        # Build initial state matrices with column and row names as states for first and second rows:
        StateMatrices <- lapply(X = lapply(X = apply(rbind(firstrow[mean.characters.to.check], secondrow[mean.characters.to.check]), 2, list), lapply, strsplit, split = "&|/"), function(x) matrix(nrow = length(x = x[[1]][[1]]), ncol = length(x = x[[1]][[2]]), dimnames = list(x[[1]][[1]], x[[1]][[2]])))

        # Fill state matrices with raw differences between each state:
        StateMatrices <- lapply(X = StateMatrices, function(x) {
          for (i in 1:ncol(x)) for (j in 1:nrow(x)) x[j, i] <- abs(as.numeric(colnames(x = x)[i]) - as.numeric(rownames(x = x)[j]))
          return(x)
        })

        # If there are unordered characters present convert maximum distances to one:
        if (any(charordering[mean.characters.to.check] == "unord")) {
          StateMatrices[which(x = charordering[mean.characters.to.check] == "unord")] <- lapply(X = StateMatrices[which(x = charordering[mean.characters.to.check] == "unord")], function(x) {
            x[x > 1] <- 1
            return(x)
          })
        }

        # Extract minimum and maximum states from each matrix with maximum being the mean distance:
        MinMaxStates <- lapply(X = lapply(X = lapply(X = StateMatrices, as.vector), mean), function(x) c(0, x))

        # Set first row values(s):
        firstrow[mean.characters.to.check] <- unlist(x = lapply(X = MinMaxStates, "[[", 1))

        # Set second row values(s):
        secondrow[mean.characters.to.check] <- unlist(x = lapply(X = MinMaxStates, "[[", 2))
      }
    }

    # Return the first and second rows either without polymorphisms or with them removed:
    return(list(firstrow, secondrow))
  }

  # Subfunction to get the absolute difference between the two rows:
  calculate_absolute_difference <- function(column) {

    # Isolate first row values:
    firstrow <- column[[1]]

    # Isolate second row values:
    secondrow <- column[[2]]

    # Get absolute differences between each pair of characters:
    return(list(abs(as.numeric(firstrow) - as.numeric(secondrow))))
  }

  # Subfunction to correct unordered distances to one:
  fix_unordered <- function(differences, compchar, ordering) {

    # If unordered and distance greater than one replace with one:
    if (length(x = which(x = differences > 1)) > 0) differences[which(x = differences > 1)[which(x = ordering[compchar[which(x = differences > 1)]] == "unord")]] <- 1

    # Return corrected unordered distances:
    return(list(differences))
  }

  # Subfunction to find incomparable characters:
  find_incomparable <- function(comparable.characters, cladistic_matrix) {
    return(setdiff(x = 1:ncol(cladistic_matrix), y = comparable.characters))
  }

  # Subfunction to get weighted differences:
  weigh_differences <- function(differences, comparable.characters, character_weights) {
    return(list(as.numeric(character_weights[comparable.characters]) * differences))
  }

  # Subfunction to get raw Euclidean distance:
  calculate_red <- function(differences) {
    return(dist(rbind(differences, rep(0, length(x = differences))), method = "euclidean"))
  }

  # Subfunction to find maximum possible differences for the comparable characters:
  find_maximum_difference <- function(comparable.characters, max.vals, min.vals) {
    return(as.numeric(max.vals[comparable.characters]) - as.numeric(min.vals[comparable.characters]))
  }

  # Subfunction to transform list of distances into an actual distance matrix:
  convert_list_to_matrix <- function(list, cladistic_matrix, diag = NULL) {

    # Set the number of rows:
    k <- nrow(cladistic_matrix)

    # Create the empty matrix:
    mat.out <- matrix(ncol = k, nrow = k)

    # Fill up the lower triangle:
    mat.out[lower.tri(mat.out)] <- unlist(x = list)

    # Make the matrix a distance matrix (both triangles have the same values):
    mat.out <- as.matrix(as.dist(mat.out))

    # If no diagonal is supplied:
    if (is.null(diag)) {

      # Set diagonal as zero:
      diag(x = mat.out) <- 0

      # If a diagonal is supplied:
    } else {

      # Add supplied diagonal as diagonal:
      diag(x = mat.out) <- diag
    }

    # Return matrix:
    return(mat.out)
  }

  # Subfunction to get count of complete characters for each taxon (diagonal in comparable characters matrix:
  count_complete <- function(column) {
    return(length(x = column) - sum(is.na(column)))
  }

  # Subfunction to calculate the Gower Coefficient:
  calculate_gc <- function(differences, comparable.characters, character_weights) {
    return(sum(differences) / sum(character_weights[comparable.characters]))
  }

  # Subfunction to calculate MORD:
  calculate_mord <- function(differences, maximum.differences) {
    return(sum(differences) / sum(maximum.differences))
  }

  # Subfunction for building starting GED data:
  build_ged_data <- function(differences, comparable.characters, cladistic_matrix, character_weights) {
    return(rbind(c(differences, rep(NA, length(x = find_incomparable(comparable.characters, cladistic_matrix)))), c(character_weights[comparable.characters], character_weights[find_incomparable(comparable.characters, cladistic_matrix)])))
  }

  # Subfunction to apply Hopkins and St John (2018) Alpha weighting of inapplicables:
  weigh_inapplicable_alpha <- function(diffs, comparable.characters, ordering, character_weights, character_dependencies, charactersByLevel, alpha) {

    # Set differences:
    Differences <- diffs

    # Set comparable characters:
    Comparablecharacters <- comparable.characters

    # Set ordering for comparable characters:
    Characterordering <- ordering[Comparablecharacters]

    # Set weights for comparable characters:
    character_weights <- character_weights[Comparablecharacters]

    # Fof each character level (from most to least nested):
    for (i in length(x = charactersByLevel):2) {

      # Get independent characters for current levels dependent characters:
      Independentcharacters <- unique(x = unlist(x = lapply(X = as.list(x = charactersByLevel[[i]]), function(x) unname(character_dependencies[character_dependencies[, "DependentCharacter"] == x, "IndependentCharacter"]))))

      # For each independent character:
      for (j in Independentcharacters) {

        # Find dependent characters:
        Dependentcharacters <- unname(character_dependencies[character_dependencies[, "IndependentCharacter"] == j, "DependentCharacter"])

        # Check characters are present in current distance:
        charactersPresent <- intersect(Comparablecharacters, Dependentcharacters)

        # If characters are present:
        if (length(x = charactersPresent) > 0) {

          # Set positions of dependent characters in current differences vector:
          DependentPositions <- match(charactersPresent, Comparablecharacters)

          # Get position of independent character in current differences vector:
          IndependentPosition <- which(x = Comparablecharacters == j)

          # Stop and warn user if matrix contains an impossible coding (i.e., dependent character coded when independent character is missing):
          if (length(x = IndependentPosition) == 0) stop("Found a dependent character coded when character it depends on is missing. Check matrix codings.")

          # Overwrite independent position with alpha-weighted value:
          diffs[IndependentPosition] <- 1 - (alpha * (1 - (sum(diffs[DependentPositions] * character_weights[DependentPositions]) / sum(character_weights[DependentPositions]))) + (1 - alpha))

          # Overwrite dependent positions with NAs:
          diffs[DependentPositions] <- NA
        }
      }
    }

    # Return modified character comparisons:
    return(diffs)
  }

  # Check for step matrices and stop and warn user if found:
  if (is.list(cladistic_matrix$topper$step_matrices)) stop("Function cannot currently deal with step matrices.")

  # Check input of distance_transformation is valid and stop and warn if not:
  if (length(x = setdiff(x = distance_transformation, y = c("arcsine_sqrt", "none", "sqrt"))) > 0) stop("distance_transformation must be one of \"none\", \"sqrt\", or \"arcsine_sqrt\".")

  # Check input of distance is valid and stop and warn if not:
  if (length(x = setdiff(x = distance_metric, y = c("RED", "GED", "GC", "MORD"))) > 0) stop("distance_metric must be one or more of \"RED\", \"GED\", \"GC\", or \"MORD\".")

  # Check input of GED type is valid and stop and warn if not:
  if (length(x = setdiff(x = ged_type, y = c("Legacy", "Hybrid", "Wills"))) > 0) stop("ged_type must be one or more of \"Legacy\", \"Hybrid\", or \"Wills\".")

  # Check input for polymorphism_behaviour is valid and stop and warn if not:
  if (length(x = setdiff(x = polymorphism_behaviour, y = c("mean.difference", "min.difference", "random"))) > 0) stop("polymorphism_behaviour must be one or more of \"mean.difference\", \"min.difference\", or \"random\".")

  # Check input for uncertainty_behaviour is valid and stop and warn if not:
  if (length(x = setdiff(x = uncertainty_behaviour, y = c("mean.difference", "min.difference", "random"))) > 0) stop("uncertainty_behaviour must be one or more of \"mean.difference\", \"min.difference\", or \"random\".")

  # Check input for inapplicable_behaviour is valid and stop and warn if not:
  if (length(x = setdiff(x = inapplicable_behaviour, y = c("missing", "HSJ"))) > 0) stop("inapplicable_behaviour must be one or more of \"missing\", or \"HSJ\".")

  # Check that if using HSJ character dependencies have been specified:
  if (inapplicable_behaviour == "HSJ" && is.null(character_dependencies)) stop("If using the \"HSJ\" inapplicable_behaviour then character_dependencies must be specified.")

  # If using HSJ and character_dependencies is set (will check data are formatted correctly):
  if (inapplicable_behaviour == "HSJ" && !is.null(character_dependencies)) {

    # Check character_dependencies is a matrix and stop and warn user if not:
    if (!is.matrix(character_dependencies)) stop("character_dependencies must be in the form of a two-column matrix.")

    # Check character_dependencies has two columns and stop and warn user if not:
    if (ncol(character_dependencies) != 2) stop("character_dependencies must be in the form of a two-column matrix.")

    # Check character_dependencies column names are correct and stop and warn user if not:
    if (length(x = setdiff(x = c("DependentCharacter", "IndependentCharacter"), y = colnames(x = character_dependencies))) > 0) stop("character_dependencies column names must be exactly \"DependentCharacter\" and \"IndependentCharacter\".")

    # Check character_dependencies are numeric values and stop and warn user if not:
    if (!is.numeric(character_dependencies)) stop("character_dependencies values must be numeric.")

    # Check character_dependencies values are within range of matrix dimensions and stop and warn user if not:
    if (length(x = setdiff(x = as.vector(character_dependencies), y = 1:sum(unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], function(x) ncol(x$matrix))))))) > 0) stop("character_dependencies can only contain character numbers within the dimensions of the cladistic_matrix specified.")

    # Check character_dependencies values do not lead to duplicated parent characters and stop and warn user if not:
    if (any(duplicated(character_dependencies[, "DependentCharacter"]))) stop("character_dependencies characters can not be dependent on two or more different independent characters.")

    # Find any characters that are both dependent and independent (and hence may lead to circularity issues):
    charactersToCheckForCircularDependency <- intersect(character_dependencies[, "DependentCharacter"], character_dependencies[, "IndependentCharacter"])

    # If there is the possibility for circularity:
    if (length(x = charactersToCheckForCircularDependency) > 0) {

      # For the ith independent character:
      for (i in unique(x = character_dependencies[, "IndependentCharacter"])) {

        # Set current character as ith character:
        CurrentCharacter <- i

        # Ste starting found character as ith character:
        Foundcharacters <- i

        # Keep going until the current character is not an independent character:
        while (sum(unlist(x = lapply(X = as.list(x = CurrentCharacter), function(x) sum(character_dependencies[, "IndependentCharacter"] == x)))) > 0) {

          # Find any dependent character(s):
          DependentCharacter <- unlist(x = lapply(X = as.list(x = CurrentCharacter), function(x) unname(character_dependencies[character_dependencies[, "IndependentCharacter"] == x, "DependentCharacter"])))

          # Check character was not already found (creating a circularity) and stop and wanr user if true:
          if (length(x = intersect(DependentCharacter, Foundcharacters)) > 0) stop("Circularity found in character_dependencies.")

          # Update found characters:
          Foundcharacters <- c(Foundcharacters, DependentCharacter)

          # Update current character(s):
          CurrentCharacter <- DependentCharacter
        }
      }
    }

    # Check alpha is a value between zero and one and stop and warn user if not:
    if (alpha > 1 || alpha < 0) stop("alpha must be a value between zero and one")
  }

  # Isolate ordering element:
  ordering <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering")))

  # Isolate minimum values:
  min.vals <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values")))

  # Isolate maximum values:
  max.vals <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values")))

  # Isolate weights:
  character_weights <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "character_weights")))

  # Combine matrix blocks into a single matrix:
  cladistic_matrix <- do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"))

  # If polymorphism_behaviour is to randomly sample one state:
  if (polymorphism_behaviour == "random") {

    # Find cells with polymorphisms:
    PolymorphismCells <- grep("&", cladistic_matrix)

    # If there are polymorphisms randomly sample one value and store:
    if (length(x = PolymorphismCells) > 0) cladistic_matrix[PolymorphismCells] <- unlist(x = lapply(X = as.list(x = cladistic_matrix[PolymorphismCells]), function(x) sample(strsplit(x, split = "&")[[1]], size = 1)))

    # Reset behaviour as mean difference to allow it to interact correctly with uncertainty_behaviour later:
    polymorphism_behaviour <- "mean.difference"
  }

  # If uncertainty_behaviour is to randomly sample one state:
  if (uncertainty_behaviour == "random") {

    # Find cells with uncertainties:
    UncertaintyCells <- grep("/", cladistic_matrix)

    # If there are uncertainties randomly sample one value and store:
    if (length(x = UncertaintyCells) > 0) cladistic_matrix[UncertaintyCells] <- unlist(x = lapply(X = as.list(x = cladistic_matrix[UncertaintyCells]), function(x) sample(strsplit(x, split = "/")[[1]], size = 1)))

    # Reset behaviour as mean difference to allow it to interact correctly with polymorphism_behaviour later:
    uncertainty_behaviour <- "mean.difference"
  }

  # If there are inapplicables and using the missing option then convert these to NAs:
  if (any(sort(x = cladistic_matrix == "")) && inapplicable_behaviour == "missing") cladistic_matrix[cladistic_matrix == ""] <- NA

  # Find all possible (symmetric) pairwise comparisons for the N taxa in the matrix (excluding self-comparisons):
  comparisons <- combn(1:nrow(cladistic_matrix), 2)

  # Find all comparable characters for each pair of taxa:
  list.of.compchar <- unlist(x = apply(comparisons, 2, find_comparable, cladistic_matrix), recursive = FALSE)

  # Get character states for each pairwise comparison:
  rows.pairs <- apply(comparisons, 2, get_pairwise_strings, cladistic_matrix)

  # Subset each pairwise comparison by just the comparable characters:
  matrix.of.char.comp <- mapply(subset_by_comparable, rows.pairs, list.of.compchar)

  # Deal with any polymorphisms found and collapse appropriately:
  matrix.of.char.comp <- mapply(edit_polymorphisms, unlist(x = apply(matrix.of.char.comp, 2, list), recursive = FALSE), list.of.compchar, MoreArgs = list(ordering, polymorphism_behaviour, uncertainty_behaviour))

  # Get the absolute differences between each comparable character for each pairwise comparison:
  diffs <- unlist(x = apply(matrix.of.char.comp, 2, calculate_absolute_difference), recursive = FALSE)

  # Correct distances for unordered characters where distance is greater than one:
  diffs <- mapply(fix_unordered, diffs, list.of.compchar, MoreArgs = list(ordering))

  # If applying the Hopkins and St John alpha approach:
  if (inapplicable_behaviour == "HSJ") {

    # Set primary-level characters in a list (where secondary etc. level characters will be added in turn):
    charactersByLevel <- list(unname(setdiff(x = unique(x = character_dependencies[, "IndependentCharacter"]), y = unique(x = character_dependencies[, "DependentCharacter"]))))

    # Set starting more nested characters:
    HigherLevelcharacters <- setdiff(x = unique(x = c(character_dependencies)), y = unlist(x = charactersByLevel))

    # Whilst there are still more nested levels of characters:
    while (length(x = HigherLevelcharacters) > 0) {

      # Add next level characters to characters by level list at next level:
      charactersByLevel[[(length(x = charactersByLevel) + 1)]] <- unname(character_dependencies[unlist(x = lapply(X = as.list(x = charactersByLevel[[length(x = charactersByLevel)]]), function(x) which(x = character_dependencies[, "IndependentCharacter"] == x))), "DependentCharacter"])

      # Set new higher level characters:
      HigherLevelcharacters <- setdiff(x = unique(x = c(character_dependencies)), y = unlist(x = charactersByLevel))
    }

    # Update differences with HSJ alpha weights:
    diffs <- mapply(weigh_inapplicable_alpha, diffs, list.of.compchar, MoreArgs = list(ordering, character_weights, character_dependencies, charactersByLevel, alpha))

    # Reweight dependent characters zero:
    character_weights[unlist(x = charactersByLevel[2:length(x = charactersByLevel)])] <- 0

    # Update comparable characters by pruning out NAs:
    list.of.compchar <- mapply(function(x, y) y[!is.na(x)], x = diffs, y = list.of.compchar, SIMPLIFY = FALSE)

    # Update differences by pruning out NAs:
    diffs <- lapply(X = diffs, function(x) x[!is.na(x)])
  }

  # Weight differences:
  diffs <- mapply(weigh_differences, diffs, list.of.compchar, MoreArgs = list(character_weights))

  # Get raw Euclidean distance (if using it):
  if (distance_metric == "RED") raw.dist <- lapply(X = diffs, calculate_red)

  # Only calculate the max differences for "GED" or "MORD" matrices:
  if (distance_metric == "GED" || distance_metric == "MORD") {

    # Find maximum possible differences for the comparable characters:
    maxdiffs <- lapply(X = list.of.compchar, find_maximum_difference, max.vals, min.vals)

    # Correct maximum differences for unordered characters:
    maxdiffs <- mapply(weigh_differences, mapply(fix_unordered, maxdiffs, list.of.compchar, MoreArgs = list(ordering)), list.of.compchar, MoreArgs = list(character_weights))
  }

  # If calculating Raw Euclidean Distances build the distance matrix:
  if (distance_metric == "RED") distance_matrix <- convert_list_to_matrix(raw.dist, cladistic_matrix)

  # If calculating the Gower Coefficient build the distance matrix:
  if (distance_metric == "GC") distance_matrix <- convert_list_to_matrix(as.list(x = mapply(calculate_gc, diffs, list.of.compchar, MoreArgs = list(character_weights))), cladistic_matrix)

  # If calculating the MORD build the distance matrix:
  if (distance_metric == "MORD") distance_matrix <- convert_list_to_matrix(mapply(calculate_mord, diffs, maxdiffs), cladistic_matrix)

  # If calculating the GED:
  if (distance_metric == "GED") {

    # Build starting GED data:
    GED.data <- mapply(build_ged_data, diffs, list.of.compchar, MoreArgs = list(cladistic_matrix, character_weights), SIMPLIFY = FALSE)

    # Transpose matrices:
    GED.data <- lapply(X = GED.data, t)

    # Now build into matrix of pairwise comparisons (odds to be compared with adjacent evens):
    GED.data <- matrix(data = (unlist(x = GED.data)), ncol = ncol(cladistic_matrix), byrow = TRUE)

    # Calculate single weighted mean univariate distance for calculating GED Legacy or Hybrid (after equation 2 in Wills 2001):
    if (ged_type != "Wills") NonWills_S_ijk_bar <- rep(sum(unlist(x = diffs)) / sum(unlist(x = maxdiffs)), length.out = length(x = diffs))

    # Calculate individual pairwise weighted mean univariate distance for calculating GED Hybrid or Wills (after equation 2 in Wills 2001):
    if (ged_type != "Legacy") {

      # Generate individual mean pairwise distance for each comparison:
      NonLegacy_S_ijk_bar <- unlist(x = lapply(X = diffs, sum)) / unlist(x = lapply(X = maxdiffs, sum))

      # Find NaNs (divide by zero errors for when there are no characters in common in a pairwsie comparison):
      NaNs <- which(x = is.nan(NonLegacy_S_ijk_bar))

      # If usings WIlls version replace NaNs with NA:
      if (ged_type == "Wills" && length(x = NaNs) > 0) NonLegacy_S_ijk_bar[NaNs] <- NA

      # If using Hybrid replace NaNs with single global mean distance value:
      if (ged_type == "Hybrid" && length(x = NaNs) > 0) NonLegacy_S_ijk_bar[NaNs] <- NonWills_S_ijk_bar[NaNs]

      # Set modified non-Legacy S_ijk_bar as main S_ijk_bar:
      S_ijk_bar <- NonLegacy_S_ijk_bar
    }

    # If using Legacy set NonWills_S_ijk_bar as main S_ijk_bar:
    if (ged_type == "Legacy") S_ijk_bar <- NonWills_S_ijk_bar

    # For each set of differences:
    for (i in seq(from = 1, to = nrow(GED.data) - 1, length.out = length(x = diffs))) {

      # Find missing distances (if any):
      missingDistances <- which(x = is.na(GED.data[i, ]))

      # Replace missing distances with S_ijk_bar (i.e., results of equation 2 in Wills 2001 into equation 1 of Wills 2001):
      if (length(x = missingDistances) > 0) GED.data[i, missingDistances] <- S_ijk_bar[ceiling(i / 2)]
    }

    # Isolate the distances:
    S_ijk <- GED.data[which(x = (1:nrow(GED.data) %% 2) == 1), ]

    # Isolate the weights:
    W_ijk <- GED.data[which(x = (1:nrow(GED.data) %% 2) == 0), ]

    # Calculate the GED (equation 1 of Wills 2001) for each pairwise comparison (ij):
    GED_ij <- sqrt(apply(W_ijk * (S_ijk^2), 1, sum))

    # Create GED distance matrix:
    distance_matrix <- convert_list_to_matrix(as.list(x = GED_ij), cladistic_matrix)
  }

  # Build comparable characters matrix:
  comp.char.matrix <- convert_list_to_matrix(lapply(X = list.of.compchar, length), cladistic_matrix, diag = apply(cladistic_matrix, 1, count_complete))

  # Add row and column names (taxa) to distance matrices:
  rownames(x = distance_matrix) <- colnames(x = distance_matrix) <- rownames(x = comp.char.matrix) <- colnames(x = comp.char.matrix) <- rownames(x = cladistic_matrix)

  # If there are any NaNs replace with NAs:
  if (any(is.nan(distance_matrix))) distance_matrix[is.nan(distance_matrix)] <- NA

  # If using a proportional distance:
  if (distance_metric == "MORD" || distance_metric == "GC") {

    # If transforming distance matrix by taking the square root - take the square root:
    if (distance_transformation == "sqrt") distance_matrix <- sqrt(distance_matrix)

    # If transforming distance matrix by taking the arcsine square root:
    if (distance_transformation == "arcsine_sqrt") {

      # Check for squared distances greater than 1:
      if (any(sort(x = sqrt(distance_matrix)) > 1)) {

        # Warn user that distances were rescaled:
        print("Squared distances found of greater than 1 so matrix was rescaled prior to taking arcsine.")

        # Take the arcsine square root of the rescaled distance matrix:
        distance_matrix <- asin(sqrt(distance_matrix) / max(sort(x = sqrt(distance_matrix))))

        # If squared distances are less than or equal to one:
      } else {

        # Take the arcsine square root directly:
        distance_matrix <- asin(sqrt(distance_matrix))
      }
    }
  }

  # Compile results as a list:
  result <- list(distance_metric = distance_metric, distance_matrix = distance_matrix, comparable_character_matrix = comp.char.matrix)

  # Output result:
  return(result)
}
