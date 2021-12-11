#' Get distance matrices from a cladistic matrix
#'
#' @description
#'
#' Takes a cladistic morphological dataset and converts it into a set of pairwise distances.
#'
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @param distance_metric The distance metric to use. Must be one of \code{"gc"}, \code{"ged"}, \code{"red"}, or \code{"mord"} (the default).
#' @param ged_type The type of GED to use. Must be one of \code{"legacy"}, \code{"hybrid"}, or \code{"wills"} (the default). See details for an explanation.
#' @param distance_transformation The type of distance transformation to perform. Options are \code{"none"}, \code{"sqrt"}, or \code{"arcsine_sqrt"} (the default). (Note: this is only really appropriate for the proportional distances, i.e., "gc" and "mord".)
#' @param polymorphism_behaviour The distance behaviour for dealing with polymorphisms. Must be one of \code{"mean_difference"}, \code{"min_difference"} (the default), or \code{"random"}.
#' @param uncertainty_behaviour The distance behaviour for dealing with uncertainties. Must be one of \code{"mean_difference"}, \code{"min_difference"} (the default), or \code{"random"}.
#' @param inapplicable_behaviour The behaviour for dealing with inapplicables. Must be one of \code{"missing"} (default), or \code{"hsj"} (Hopkins and St John 2018; see details).
#' @param character_dependencies Only relevant if using \code{inapplicable_behaviour = "hsj"}. Must be a two-column matrix with colnames "dependent_character" and "independent_character" that specifies character hierarchies. See details.
#' @param alpha The alpha value (sensu Hopkins and St John 2018). Only relevant if using \code{inapplicable_behaviour = "hsj"}. See details.
#'
#' @details
#'
#' There are many options to consider when generating a distance matrix from morphological data, including the metric to use, how to treat inapplicable, polymorphic (e.g., 0&1), or uncertain (e.g., 0/1) states, and whether the output should be transformed (e.g., by taking the square root so that the distances are - or approximate - Euclidean distances). Some of these issues have been discussed previously in the literature (e.g., Lloyd 2016; Hopkins and St John 2018), but all likely require further study.
#'
#' Claddis currently offers four different distance metrics: 1. Raw Euclidean Distance (\code{"red"}) - this is only really applicable if there are no missing data, 2. The Gower Coefficient (\code{"gc"}; Gower 1971) - this rescales distances by the number of characters that can be coded for both taxa in each pairwise comparison thus correcting for missing data, 3. The Maximum Observable Rescaled Distance (\code{"mord"}) - this was introduced by Lloyd (2016) as an extension of the \code{"gc"} designed to deal with the fact that multistate ordered characters can lead to \code{"gc"}s of greater than 1 and works by rescaling by the maximum possible distance that could be observed based on the number of characters codable in each pairwise comparison meaning all resulting distances are on a zero to one scale, and 4. The Generalised Euclidean Distance - this was introduced by Wills (1998) as a means of correcting for the fact that a \code{"red"} metric will become increasingly non-Euclidean as the amount of missing data increases and works by filling in missing distances (for characters that are coded as missing in at least one taxon in the pairwise comparison) by using the mean pairwise dissimilarity for that taxon pair as a substitute. In effect then, \code{"red"} makes no consideration of missing data, \code{"gc"} and \code{"mord"} normalise by the available data (and are identical if there are no ordered multistate characters), and \code{"ged"} fills in missing distances by extrapolating from the available data.
#'
#' Note that Lloyd (2016) misidentified the substitute dissimilarity for the \code{"ged"} as the mean for the whole data set (Hopkins and St John 2018) and this was the way the GED implementation of Claddis operated up to version 0.2. This has now been amended (as of version 0.3) so that the function produces the \code{"ged"} in the form that Wills (1998) intended. However, this implementation can still be accessed as the \code{"legacy"} option for \code{ged_type}, with \code{"wills"} being the WIlls (1998) implementation. An advantage of this misinterpreted form of the GED is that it will always return a complete pairwise distance matrix, however it is not recommended (see Lloyd 2016). Instead a third option for \code{ged_type} - (\code{"hybrid"}) - offers the same outcome but only uses the mean distance from the entire matrix in the case where there are no codable characters in common in a pairwise comparison. This new hybrid option has not been used in a published study.
#'
#' Typically the resulting distance matrix will be used in an ordination procedure such as principal coordinates (effectively classical multidimensional scaling where k, the number of axes, is maximised at N - 1, where N is the number of rows (i.e., taxa) in the matrix). As such the distance should be - or approximate - Euclidean and hence a square root transformation is typically applied (\code{distance_transformation} with the \code{"sqrt"} option). However, if applying pre-ordination (i.e., ordination-free) disparity metrics (e.g., weighted mean pairwise distance) you may wish to avoid any transformation (\code{"none"} option). In particular the MORD will only fall on a zero to one scale if this is the case. However, if transforming the MORD for ordination this zero to one property may mean the arcsine square root (\code{"arcsine_sqrt"} option) is preferred. (Note that if using only unordered multistate or binary characters and the \code{"gc"} the zero to one scale will apply too.)
#'
#' An unexplored option in distance matrix construction is how to deal with polymorphisms (Lloyd 2016). Up to version 0.2 of Claddis all polymorphisms were treated the same regardless of whether they were true polymorphisms (multiple states are observed in the taxon) or uncertainties (multiple, but not all states, are posited for the taxon). Since version 0.3, however, these two forms can be distinguished by using the different #NEXUS forms (Maddison et al. 1997), i.e., (01) for polymorphisms and \{01\} for uncertainties and within Claddis these are represented as 0&1 or 0/1, respectively. Thus, since 0.3 Claddis allows these two forms to be treated separately, and hence differently (with \code{polymorphism_behaviour} and \code{uncertainty_behaviour}). Again, up to version 0.2 of Claddis no options for polymorphism behaviour were offered, instead only a minimum distance was employed. I.e., the distance between a taxon coded 0&1 and a taxon coded 2 would be the smaller of the comparisons 0 with 2 or 1 with 2. Since version 0.3 this is encoded in the \code{"min_difference"} option. Currently two alternatives (\code{"mean_difference"} and \code{"random"}) are offered. The first takes the mean of each possible difference and the second simply samples one of the states at random. Note this latter option makes the function stochastic and so it should be rerun multiple times (for example, with a \code{for} loop or \code{apply} function). In general this issue (and these options) are not explored in the literature and so no recommendation can be made beyond that users should think carefully about what this choice may mean for their individual data set(s) and question(s).
#'
#' A final consideration is how to deal with inapplicable characters. Up to version 0.2 Claddis treated inapplicable and missing characters the same (as NA values, i.e., missing data). However, since Claddis version 0.3 these can be imported separately, i.e., by using the "MISSING" and "GAP" states in #NEXUS format (Maddison et al. 1997), with the latter typically representing the inapplicable character. These appear as NA and empty strings (""), respectively, in Claddis format. Hopkins and St John (2018) showed how inapplicable characters - typically assumed to represent secondary characters - could be treated in generating distance matrices. These are usually hierarchical in form. E.g., a primary character might record the presence or absence of feathers and a secondary character whether those feathers are symmetric or asymmetric. The latter will generate inapplicable states for taxa without feathers and without correcting for this ranked distances can be incorrect (Hopkins and St John 2018). Unfortunately, however, the #NEXUS format (Maddison et al. 1997) does not really allow explicit linkage between primary and secondary characters and so this information must be provided separately to use the Hopkins and St John (2018) approach. This is done here with the \code{character_dependencies} option. This must be in the form of a two-column matrix with column headers of \code{"dependent_character"} and \code{"independent_character"}. The former being secondary characters and the latter the corresponding primary character. (Note that characters are to be numbered across the whole matrix from 1 to N and do not restart with each block of the matrix.) If using \code{inapplicable_behaviour = "hsj"} the user must also provide an \code{alpha} value between zero and one. When \code{alpha = 0} the secondary characters contribute nothing to the distance and when \code{alpha = 1} the primary character is not counted in the weight separately (see Hopkins and St John 2018). The default value (0.5) offers a compromise between these two extremes.
#'
#' Here the implementation of this approach differs somewhat from the code available in the supplementary materials to Hopkins and St John (2018). Specifically, this approach is incorporated (and used) regardless of the overriding distance metric (i.e., the \code{distance_metric} option). Additionally, the Hopkins and St John function specifically allows an extra level of dependency (secondary and tertary characters) with these being applied recursively (tertiary first then secondary). Here, though, additional levels of dependency do not need to be defined by the user as this information is already encoded in the \code{character_dependencies} option. Furthermore, because of this any level of dependency is possible (if unlikely), e.g., quarternary etc.
#'
#' @return
#'
#' \item{distance_metric}{The distance metric used.}
#' \item{distance_matrix}{The pairwise distance matrix generated.}
#' \item{comparable_character_matrix}{The matrix of characters that can be compared for each pairwise distance.}
#' \item{comparable_weights_matrix}{The matrix of weights for each pairwise distance.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Thomas Guillerme \email{guillert@@tcd.ie}
#'
#' @references
#'
#' Gower, J. C., 1971. A general coefficient of similarity and some of its properties. \emph{Biometrika}, \bold{27}, 857-871.
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
#' # Get morphological distances for the Day et al. (2016) data set:
#' distances <- calculate_morphological_distances(cladistic_matrix = day_2016)
#'
#' # Show distance metric:
#' distances$distance_metric
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
#'       "dependent_character",
#'       "independent_character"
#'     )
#'   )
#' )
#'
#' # Get morphological distances for the Day et
#' # al. (2016) data set using HSJ approach:
#' distances <- calculate_morphological_distances(
#'   cladistic_matrix = day_2016,
#'   inapplicable_behaviour = "hsj",
#'   character_dependencies = character_dependencies,
#'   alpha = 0.5
#' )
#'
#' # Show distance metric:
#' distances$distance_metric
#'
#' # Show distance matrix:
#' distances$distance_matrix
#'
#' # Show number of characters that can be scored for
#' # each pairwise comparison:
#' distances$comparable_character_matrix
#'
#' # Show weighting of each calculable pairwise distance:
#' distances$comparable_weights_matrix
#'
#' @export calculate_morphological_distances
calculate_morphological_distances <- function(cladistic_matrix, distance_metric = "mord", ged_type = "wills", distance_transformation = "arcsine_sqrt", polymorphism_behaviour = "min_difference", uncertainty_behaviour = "min_difference", inapplicable_behaviour = "missing", character_dependencies = NULL, alpha = 0.5) {

  # ADD HOPKINS SUGGESTION (VIA EMAIL) FOR FOURTH GEDTYPE WHERE MEAN DISTANCE FOR CHARACTER REPLACES MISSING VALUES.
  # CHECK POLYMORPHISM UNCERTAINTY IN GENERAL AS NOT CLEAR IT IS DOING WHAT IT SHOULD DO.
  # CHECK TRANSFORM IS APPROPRIATE AND WARN USER IF NOT
  # MAYBE ALLOW MANHATTAN TYPE DISTANCES TOO.
  # ADD LEHMANN REFERENCE!
  # ALLOW MANHATTAN DISTANCES
  # CONSIDER DISTANCES FOR POLYMORPHISMS IN SAME WAY PHYTOOLS DOES WITH POLYMK
  # ADD TRANSORMATION USED TO OUTPUT (AS MAY CHANGE IF OPTIONS COLLIDE)
  # ALLOW WMPD SOMEHOW? MAYBE A SEPARATE FUNCTION WITH GROUPS (WHICH NEEDS TO BE IMPLEMENTED ACROSS THE PACKAGE FOR DISPARITY PLOTS ETC.)
  # RETOOL AROUND STOCHASTIC CHARACTER MAPS IF DOING PHYLOGENY TOO
  # CHECK HSJ MAKES SENSE IN CONTEXT OF COMPARABLE WEIGHTS - MIGHT NEED A THEORETICAL SOLUTION BEFORE AN IMPLEMENTATION ONE.

  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")

  # Subfunction to find comparable characters for a pairwise taxon comparison:
  find_comparable <- function(taxon_pair, cladistic_matrix) {

    # Get intersection of characters that are coded for both taxa in a pair:
    output <- intersect(intersect(which(x = !is.na(cladistic_matrix[taxon_pair[[1]], ])), which(x = cladistic_matrix[taxon_pair[[1]], ] != "")), intersect(which(x = !is.na(cladistic_matrix[taxon_pair[[2]], ])), which(x = cladistic_matrix[taxon_pair[[2]], ] != "")))

    # Return output:
    list(output)
  }

  # Subfunction to get character strings for each pair of taxa:
  get_pairwise_strings <- function(taxon_pair, cladistic_matrix) {

    # Get character states for first taxon in pair:
    row1 <- cladistic_matrix[rownames(x = cladistic_matrix)[taxon_pair[[1]]], ]

    # Get character states for second taxon in pair:
    row2 <- cladistic_matrix[rownames(x = cladistic_matrix)[taxon_pair[[2]]], ]

    # Return output as a list:
    list(row1, row2)
  }

  # Subfunction to subset pairwise comparisons by just comparable characters:
  subset_by_comparable <- function(row_pair, comparable_characters) {

    # Collapse first row to just comparable characters:
    row_pair[[1]] <- row_pair[[1]][comparable_characters]

    # Collapse second row to just comparable characters:
    row_pair[[2]] <- row_pair[[2]][comparable_characters]

    # Output colapsed row pair:
    row_pair
  }

  # Subfunction to edit polymorphic characters down to a single value:
  edit_polymorphisms <- function(comparisons, comparable_characters, ordering, polymorphism_behaviour, uncertainty_behaviour) {

    # Set first taxon values:
    first_row <- comparisons[[1]]

    # Set second taxon values:
    second_row <- comparisons[[2]]

    # If there are any inapplicables:
    if (any(c(first_row, second_row) == "")) {

      # Find inapplicable positions:
      inapplicable_positions <- sort(x = unique(x = c(which(x = first_row == ""), which(x = second_row == ""))))

      # Find polymorphism and uncertainty positions:
      polymorphism_and_uncertainty_positions <- sort(x = unique(x = c(grep("/|&", first_row), grep("/|&", second_row))))

      # If there are polymorphisms or uncertianties that match up with inapplicables:
      if (length(x = intersect(inapplicable_positions, polymorphism_and_uncertainty_positions)) > 0) {

        # Find positions where collapsing to a single value is required:
        collapse_positions <- intersect(inapplicable_positions, polymorphism_and_uncertainty_positions)

        # Collapse any polymorphisms or uncertianties in first row to just first value:
        first_row[collapse_positions] <- unlist(x = lapply(X = strsplit(first_row[collapse_positions], split = "/|&"), function(x) ifelse(length(x = x) == 0, "", x[1])))

        # Collapse any polymorphisms or uncertianties in second row to just first value:
        second_row[collapse_positions] <- unlist(x = lapply(X = strsplit(second_row[collapse_positions], split = "/|&"), function(x) ifelse(length(x = x) == 0, "", x[1])))
      }
    }

    # Set ordering for comparable characters:
    character_ordering <- ordering[comparable_characters]

    # Only if there are polymorphisms or uncertainties:
    if (length(x = c(grep("&", unique(x = c(first_row, second_row))), grep("/", unique(x = c(first_row, second_row))))) > 0) {

      # Find ampersands (polymorphisms):
      ampersand_elements <- sort(x = c(grep("&", first_row), grep("&", second_row)))

      # Find slashes (uncertianties):
      slash_elements <- sort(x = c(grep("/", first_row), grep("/", second_row)))

      # Combine to find all characters to check:
      characters_to_check <- sort(x = unique(x = c(ampersand_elements, slash_elements)))

      # Set behaviours as either the shared version or minimum difference if they contradict (may need to modify this later for more complex options):
      behaviour <- unlist(x = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(apply(apply(rbind(first_row[characters_to_check], second_row[characters_to_check]), 2, gsub, pattern = "[:0-9:]", replacement = ""), 2, list), unlist), function(x) x[nchar(x = x) > 0]), function(x) ifelse(nchar(x = x) > 0, strsplit(x, split = "")[[1]][1], x)), function(x) gsub(pattern = "&", replacement = polymorphism_behaviour, x = x)), function(x) gsub(pattern = "/", replacement = uncertainty_behaviour, x = x)), unique), function(x) ifelse(length(x) > 1, "min_difference", x)))

      # If behaviour is to find minimum differences:
      if (any(behaviour == "min_difference")) {

        # Set up minimum difference characters to check:
        min_characters_to_check <- characters_to_check[behaviour == "min_difference"]

        # Find intersecting character states for each character:
        intersection_character <- lapply(X = lapply(X = lapply(X = lapply(X = apply(rbind(first_row[min_characters_to_check], second_row[min_characters_to_check]), 2, strsplit, split = "&|/"), unlist), sort), rle), function(x) x$values[x$lengths > 1][1])

        # If at least one intersecting character state was found:
        if (any(!is.na(unlist(x = intersection_character)))) {

          # Record rows to update:
          rows_to_update <- which(x = !is.na(unlist(x = intersection_character)))

          # Store (first) shared state for both taxa:
          first_row[min_characters_to_check[rows_to_update]] <- second_row[min_characters_to_check[rows_to_update]] <- unlist(x = intersection_character)[rows_to_update]

          # Update minimum characters to check:
          min_characters_to_check <- min_characters_to_check[-rows_to_update]
        }

        # Only continue if there are still characters that need to be fixed:
        if (length(x = min_characters_to_check) > 0) {

          # Build two option matrices for every comparison:
          two_option_matrices <- lapply(X = apply(rbind(first_row[min_characters_to_check], second_row[min_characters_to_check]), 2, strsplit, split = "&|/"), function(x) rbind(c(min(as.numeric(x[[1]])), max(as.numeric(x[[2]]))), c(max(as.numeric(x[[1]])), min(as.numeric(x[[2]])))))

          # Pick smallest difference as minimum and maximum states:
          min_max_states <- lapply(X = lapply(X = lapply(X = two_option_matrices, function(x) x[which(x = abs(apply(x, 1, diff)) == min(abs(apply(x, 1, diff)))), ]), sort), as.character)

          # Set first row values(s):
          first_row[min_characters_to_check] <- unlist(x = lapply(X = min_max_states, "[[", 1))

          # Set second row values(s):
          second_row[min_characters_to_check] <- unlist(x = lapply(X = min_max_states, "[[", 2))
        }
      }

      # If any behaviour is to find mean differences:
      if (any(behaviour == "mean_difference")) {

        # Set up minimum difference characters to check:
        mean_characters_to_check <- characters_to_check[behaviour == "mean_difference"]

        # Build initial state matrices with column and row names as states for first and second rows:
        state_matrices <- lapply(X = lapply(X = apply(rbind(first_row[mean_characters_to_check], second_row[mean_characters_to_check]), 2, list), lapply, strsplit, split = "&|/"), function(x) matrix(nrow = length(x = x[[1]][[1]]), ncol = length(x = x[[1]][[2]]), dimnames = list(x[[1]][[1]], x[[1]][[2]])))

        # Fill state matrices with raw differences between each state:
        state_matrices <- lapply(X = state_matrices, function(x) {
          for (i in 1:ncol(x)) for (j in 1:nrow(x)) x[j, i] <- abs(as.numeric(colnames(x = x)[i]) - as.numeric(rownames(x = x)[j]))
          return(x)
        })

        # If there are unordered characters present convert maximum distances to one:
        if (any(character_ordering[mean_characters_to_check] == "unordered")) {
          state_matrices[which(x = character_ordering[mean_characters_to_check] == "unordered")] <- lapply(X = state_matrices[which(x = character_ordering[mean_characters_to_check] == "unordered")], function(x) {
            x[x > 1] <- 1
            return(x)
          })
        }

        # Extract minimum and maximum states from each matrix with maximum being the mean distance:
        min_max_states <- lapply(X = lapply(X = lapply(X = state_matrices, as.vector), mean), function(x) c(0, x))

        # Set first row values(s):
        first_row[mean_characters_to_check] <- unlist(x = lapply(X = min_max_states, "[[", 1))

        # Set second row values(s):
        second_row[mean_characters_to_check] <- unlist(x = lapply(X = min_max_states, "[[", 2))
      }
    }

    # Return the first and second rows either without polymorphisms or with them removed:
    return(list(first_row, second_row))
  }

  # Subfunction to get the absolute difference between the two rows:
  calculate_absolute_difference <- function(column) {

    # Isolate first row values:
    first_row <- column[[1]]

    # Isolate second row values:
    second_row <- column[[2]]

    # Get absolute differences between each pair of characters:
    return(list(abs(as.numeric(first_row) - as.numeric(second_row))))
  }

  # Subfunction to correct unordered distances to one:
  fix_unordered <- function(differences, comparable_characters, ordering) {

    # If unordered and distance greater than one replace with one:
    if (length(x = which(x = differences > 1)) > 0) differences[which(x = differences > 1)[which(x = ordering[comparable_characters[which(x = differences > 1)]] == "unordered")]] <- 1

    # Return corrected unordered distances:
    return(list(differences))
  }

  # Subfunction to find incomparable characters:
  find_incomparable <- function(comparable_characters, cladistic_matrix) {
    setdiff(x = 1:ncol(cladistic_matrix), y = comparable_characters)
  }

  # Subfunction to get weighted differences:
  weigh_differences <- function(differences, comparable_characters, character_weights) {
    list(as.numeric(character_weights[comparable_characters]) * differences)
  }

  # Subfunction to get raw Euclidean distance:
  calculate_red <- function(differences) {
    dist(rbind(differences, rep(0, length(x = differences))), method = "euclidean")
  }

  # Subfunction to find maximum possible differences for the comparable characters:
  find_maximum_difference <- function(comparable_characters, maximum_values, minimum_values) {
    as.numeric(maximum_values[comparable_characters]) - as.numeric(minimum_values[comparable_characters])
  }

  # Subfunction to transform list of distances into an actual distance matrix:
  convert_list_to_matrix <- function(list, cladistic_matrix, diag = NULL) {

    # Set the number of rows:
    k <- nrow(cladistic_matrix)

    # Create the empty matrix:
    matrix_output <- matrix(ncol = k, nrow = k)

    # Fill up the lower triangle:
    matrix_output[lower.tri(matrix_output)] <- unlist(x = list)

    # Make the matrix a distance matrix (both triangles have the same values):
    matrix_output <- as.matrix(as.dist(matrix_output))

    # If no diagonal is supplied:
    if (is.null(diag)) {

      # Set diagonal as zero:
      diag(x = matrix_output) <- 0

      # If a diagonal is supplied:
    } else {

      # Add supplied diagonal as diagonal:
      diag(x = matrix_output) <- diag
    }

    # Return matrix:
    matrix_output
  }

  # Subfunction to get count of complete characters for each taxon (diagonal in comparable characters matrix:
  count_complete <- function(column) {
    length(x = column) - sum(is.na(column))
  }

  # Subfunction to calculate the Gower Coefficient:
  calculate_gc <- function(differences, comparable_characters, character_weights) {
    sum(differences) / sum(character_weights[comparable_characters])
  }

  # Subfunction to calculate MORD:
  calculate_mord <- function(differences, maximum_differences) {
    sum(differences) / sum(maximum_differences)
  }

  # Subfunction for building starting GED data:
  build_ged_data <- function(differences, comparable_characters, cladistic_matrix, character_weights) {
    rbind(c(differences, rep(NA, length(x = find_incomparable(comparable_characters, cladistic_matrix)))), c(character_weights[comparable_characters], character_weights[find_incomparable(comparable_characters, cladistic_matrix)]))
  }

  # Subfunction to apply Hopkins and St John (2018) Alpha weighting of inapplicables:
  weigh_inapplicable_alpha <- function(differences, comparable_characters, ordering, character_weights, character_dependencies, characters_by_level, alpha) {

    # Set ordering for comparable characters:
    character_ordering <- ordering[comparable_characters]

    # Set weights for comparable characters:
    character_weights <- character_weights[comparable_characters]

    # Fof each character level (from most to least nested):
    for (i in length(x = characters_by_level):2) {

      # Get independent characters for current levels dependent characters:
      independent_characters <- unique(x = unlist(x = lapply(X = as.list(x = characters_by_level[[i]]), function(x) unname(character_dependencies[character_dependencies[, "dependent_character"] == x, "independent_character"]))))

      # For each independent character:
      for (j in independent_characters) {

        # Find dependent characters:
        dependent_characters <- unname(character_dependencies[character_dependencies[, "independent_character"] == j, "dependent_character"])

        # Check characters are present in current distance:
        present_characters <- intersect(comparable_characters, dependent_characters)

        # If characters are present:
        if (length(x = present_characters) > 0) {

          # Set positions of dependent characters in current differences vector:
          dependent_positions <- match(present_characters, comparable_characters)

          # Get position of independent character in current differences vector:
          independent_position <- which(x = comparable_characters == j)

          # Stop and warn user if matrix contains an impossible coding (i.e., dependent character coded when independent character is missing):
          if (length(x = independent_position) == 0) stop("Found a dependent character coded when character it depends on is missing. Check matrix codings.")

          # Overwrite independent position with alpha-weighted value:
          differences[independent_position] <- 1 - (alpha * (1 - (sum(differences[dependent_positions] * character_weights[dependent_positions]) / sum(character_weights[dependent_positions]))) + (1 - alpha))

          # Overwrite dependent positions with NAs:
          differences[dependent_positions] <- NA
        }
      }
    }

    # Return modified character comparisons:
    differences
  }

  # Check for costmatrices and stop and warn user if found:
  if (is.list(cladistic_matrix$topper$costmatrices)) stop("Function cannot currently deal with costmatrices.")

  # Check input of distance_transformation is valid and stop and warn if not:
  if (length(x = setdiff(x = distance_transformation, y = c("arcsine_sqrt", "none", "sqrt"))) > 0) stop("distance_transformation must be one of \"none\", \"sqrt\", or \"arcsine_sqrt\".")

  # Check input of distance is valid and stop and warn if not:
  if (length(x = setdiff(x = distance_metric, y = c("red", "ged", "gc", "mord"))) > 0) stop("distance_metric must be one or more of \"red\", \"ged\", \"gc\", or \"mord\".")

  # Check input of GED type is valid and stop and warn if not:
  if (length(x = setdiff(x = ged_type, y = c("legacy", "hybrid", "wills"))) > 0) stop("ged_type must be one or more of \"legacy\", \"hybrid\", or \"wills\".")

  # Check input for polymorphism_behaviour is valid and stop and warn if not:
  if (length(x = setdiff(x = polymorphism_behaviour, y = c("mean_difference", "min_difference", "random"))) > 0) stop("polymorphism_behaviour must be one or more of \"mean_difference\", \"min_difference\", or \"random\".")

  # Check input for uncertainty_behaviour is valid and stop and warn if not:
  if (length(x = setdiff(x = uncertainty_behaviour, y = c("mean_difference", "min_difference", "random"))) > 0) stop("uncertainty_behaviour must be one or more of \"mean_difference\", \"min_difference\", or \"random\".")

  # Check input for inapplicable_behaviour is valid and stop and warn if not:
  if (length(x = setdiff(x = inapplicable_behaviour, y = c("missing", "hsj"))) > 0) stop("inapplicable_behaviour must be one or more of \"missing\", or \"hsj\".")

  # Check that if using HSJ character dependencies have been specified:
  if (inapplicable_behaviour == "hsj" && is.null(character_dependencies)) stop("If using the \"hsj\" inapplicable_behaviour then character_dependencies must be specified.")

  # If using HSJ and character_dependencies is set (will check data are formatted correctly):
  if (inapplicable_behaviour == "hsj" && !is.null(character_dependencies)) {

    # Check character_dependencies is a matrix and stop and warn user if not:
    if (!is.matrix(character_dependencies)) stop("character_dependencies must be in the form of a two-column matrix.")

    # Check character_dependencies has two columns and stop and warn user if not:
    if (ncol(character_dependencies) != 2) stop("character_dependencies must be in the form of a two-column matrix.")

    # Check character_dependencies column names are correct and stop and warn user if not:
    if (length(x = setdiff(x = c("dependent_character", "independent_character"), y = colnames(x = character_dependencies))) > 0) stop("character_dependencies column names must be exactly \"dependent_character\" and \"independent_character\".")

    # Check character_dependencies are numeric values and stop and warn user if not:
    if (!is.numeric(character_dependencies)) stop("character_dependencies values must be numeric.")

    # Check character_dependencies values are within range of matrix dimensions and stop and warn user if not:
    if (length(x = setdiff(x = as.vector(character_dependencies), y = 1:sum(unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], function(x) ncol(x$matrix))))))) > 0) stop("character_dependencies can only contain character numbers within the dimensions of the cladistic_matrix specified.")

    # Check character_dependencies values do not lead to duplicated parent characters and stop and warn user if not:
    if (any(duplicated(character_dependencies[, "dependent_character"]))) stop("character_dependencies characters can not be dependent on two or more different independent characters.")

    # Find any characters that are both dependent and independent (and hence may lead to circularity issues):
    potentially_circular_dependencies <- intersect(character_dependencies[, "dependent_character"], character_dependencies[, "independent_character"])

    # If there is the possibility for circularity:
    if (length(x = potentially_circular_dependencies) > 0) {

      # For the ith independent character:
      for (i in unique(x = character_dependencies[, "independent_character"])) {

        # Set current character as ith character:
        current_character <- i

        # Ste starting found character as ith character:
        found_characters <- i

        # Keep going until the current character is not an independent character:
        while (sum(unlist(x = lapply(X = as.list(x = current_character), function(x) sum(character_dependencies[, "independent_character"] == x)))) > 0) {

          # Find any dependent character(s):
          dependent_character <- unlist(x = lapply(X = as.list(x = current_character), function(x) unname(character_dependencies[character_dependencies[, "independent_character"] == x, "dependent_character"])))

          # Check character was not already found (creating a circularity) and stop and wanr user if true:
          if (length(x = intersect(dependent_character, found_characters)) > 0) stop("Circularity found in character_dependencies. Fix and try again.")

          # Update found characters:
          found_characters <- c(found_characters, dependent_character)

          # Update current character(s):
          current_character <- dependent_character
        }
      }
    }

    # Check alpha is a value between zero and one and stop and warn user if not:
    if (alpha > 1 || alpha < 0) stop("alpha must be a value between zero and one")
  }

  # Isolate ordering element:
  ordering <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering")))

  # Isolate minimum values:
  minimum_values <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values")))

  # Isolate maximum values:
  maximum_values <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values")))

  # Isolate weights:
  character_weights <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "character_weights")))

  # Combine matrix blocks into a single matrix:
  cladistic_matrix <- do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"))

  # If polymorphism_behaviour is to randomly sample one state:
  if (polymorphism_behaviour == "random") {

    # Find cells with polymorphisms:
    polymorphism_cells <- grep("&", cladistic_matrix)

    # If there are polymorphisms randomly sample one value and store:
    if (length(x = polymorphism_cells) > 0) cladistic_matrix[polymorphism_cells] <- unlist(x = lapply(X = as.list(x = cladistic_matrix[polymorphism_cells]), function(x) sample(strsplit(x, split = "&")[[1]], size = 1)))

    # Reset behaviour as mean difference to allow it to interact correctly with uncertainty_behaviour later:
    polymorphism_behaviour <- "mean_difference"
  }

  # If uncertainty_behaviour is to randomly sample one state:
  if (uncertainty_behaviour == "random") {

    # Find cells with uncertainties:
    uncertainty_cells <- grep("/", cladistic_matrix)

    # If there are uncertainties randomly sample one value and store:
    if (length(x = uncertainty_cells) > 0) cladistic_matrix[uncertainty_cells] <- unlist(x = lapply(X = as.list(x = cladistic_matrix[uncertainty_cells]), function(x) sample(strsplit(x, split = "/")[[1]], size = 1)))

    # Reset behaviour as mean difference to allow it to interact correctly with polymorphism_behaviour later:
    uncertainty_behaviour <- "mean_difference"
  }

  # If there are inapplicables and using the missing option then convert these to NAs:
  if (any(sort(x = cladistic_matrix == "")) && inapplicable_behaviour == "missing") cladistic_matrix[cladistic_matrix == ""] <- NA

  # Find all possible (symmetric) pairwise comparisons for the N taxa in the matrix (excluding self-comparisons):
  comparisons <- combn(1:nrow(cladistic_matrix), 2)

  # Find all comparable characters for each pair of taxa:
  comparable_character_list <- unlist(x = apply(comparisons, 2, find_comparable, cladistic_matrix), recursive = FALSE)

  # Get character states for each pairwise comparison:
  rows_pairs <- apply(comparisons, 2, get_pairwise_strings, cladistic_matrix)

  # Subset each pairwise comparison by just the comparable characters:
  matrix_of_comparable_characters <- mapply(subset_by_comparable, rows_pairs, comparable_character_list)

  # Deal with any polymorphisms found and collapse appropriately:
  matrix_of_comparable_characters <- mapply(edit_polymorphisms, unlist(x = apply(matrix_of_comparable_characters, 2, list), recursive = FALSE), comparable_character_list, MoreArgs = list(ordering, polymorphism_behaviour, uncertainty_behaviour))

  # Get the absolute differences between each comparable character for each pairwise comparison:
  absolute_differences <- unlist(x = apply(matrix_of_comparable_characters, 2, calculate_absolute_difference), recursive = FALSE)

  # Correct distances for unordered characters where distance is greater than one:
  absolute_differences <- mapply(fix_unordered, absolute_differences, comparable_character_list, MoreArgs = list(ordering))

  # If applying the Hopkins and St John alpha approach:
  if (inapplicable_behaviour == "hsj") {

    # Set primary-level characters in a list (where secondary etc. level characters will be added in turn):
    characters_by_level <- list(unname(setdiff(x = unique(x = character_dependencies[, "independent_character"]), y = unique(x = character_dependencies[, "dependent_character"]))))

    # Set starting more nested characters:
    higher_level_characters <- setdiff(x = unique(x = c(character_dependencies)), y = unlist(x = characters_by_level))

    # Whilst there are still more nested levels of characters:
    while (length(x = higher_level_characters) > 0) {

      # Add next level characters to characters by level list at next level:
      characters_by_level[[(length(x = characters_by_level) + 1)]] <- unname(character_dependencies[unlist(x = lapply(X = as.list(x = characters_by_level[[length(x = characters_by_level)]]), function(x) which(x = character_dependencies[, "independent_character"] == x))), "dependent_character"])

      # Set new higher level characters:
      higher_level_characters <- setdiff(x = unique(x = c(character_dependencies)), y = unlist(x = characters_by_level))
    }

    # Update differences with HSJ alpha weights:
    absolute_differences <- mapply(weigh_inapplicable_alpha, absolute_differences, comparable_character_list, MoreArgs = list(ordering, character_weights, character_dependencies, characters_by_level, alpha))

    # Reweight dependent characters zero:
    character_weights[unlist(x = characters_by_level[2:length(x = characters_by_level)])] <- 0

    # Update comparable characters by pruning out NAs:
    comparable_character_list <- mapply(function(x, y) y[!is.na(x)], x = absolute_differences, y = comparable_character_list, SIMPLIFY = FALSE)

    # Update differences by pruning out NAs:
    absolute_differences <- lapply(X = absolute_differences, function(x) x[!is.na(x)])
  }

  # Weight differences:
  absolute_differences <- mapply(weigh_differences, absolute_differences, comparable_character_list, MoreArgs = list(character_weights))

  # Get raw Euclidean distance (if using it):
  if (distance_metric == "red") raw_distances <- lapply(X = absolute_differences, calculate_red)

  # Only calculate the max differences for "ged" or "mord" matrices:
  if (distance_metric == "ged" || distance_metric == "mord") {

    # Find maximum possible differences for the comparable characters:
    maximum_possible_differences <- lapply(X = comparable_character_list, find_maximum_difference, maximum_values, minimum_values)

    # Correct maximum differences for unordered characters:
    maximum_possible_differences <- mapply(weigh_differences, mapply(fix_unordered, maximum_possible_differences, comparable_character_list, MoreArgs = list(ordering)), comparable_character_list, MoreArgs = list(character_weights))
  }

  # If calculating Raw Euclidean Distances build the distance matrix:
  if (distance_metric == "red") distance_matrix <- convert_list_to_matrix(raw_distances, cladistic_matrix)

  # If calculating the Gower Coefficient build the distance matrix:
  if (distance_metric == "gc") distance_matrix <- convert_list_to_matrix(as.list(x = mapply(calculate_gc, absolute_differences, comparable_character_list, MoreArgs = list(character_weights))), cladistic_matrix)

  # If calculating the MORD build the distance matrix:
  if (distance_metric == "mord") distance_matrix <- convert_list_to_matrix(mapply(calculate_mord, absolute_differences, maximum_possible_differences), cladistic_matrix)

  # If calculating the GED:
  if (distance_metric == "ged") {

    # Build starting GED data:
    ged_data <- mapply(build_ged_data, absolute_differences, comparable_character_list, MoreArgs = list(cladistic_matrix, character_weights), SIMPLIFY = FALSE)

    # Transpose matrices:
    ged_data <- lapply(X = ged_data, t)

    # Now build into matrix of pairwise comparisons (odds to be compared with adjacent evens):
    ged_data <- matrix(data = (unlist(x = ged_data)), ncol = ncol(cladistic_matrix), byrow = TRUE)

    # Calculate single weighted mean univariate distance for calculating GED Legacy or Hybrid (after equation 2 in Wills 2001):
    if (ged_type != "wills") nonwills_s_ijk_bar <- rep(sum(unlist(x = absolute_differences)) / sum(unlist(x = maximum_possible_differences)), length.out = length(x = absolute_differences))

    # Calculate individual pairwise weighted mean univariate distance for calculating GED Hybrid or Wills (after equation 2 in Wills 2001):
    if (ged_type != "legacy") {

      # Generate individual mean pairwise distance for each comparison:
      nonlegacy_s_ijk_bar <- unlist(x = lapply(X = absolute_differences, sum)) / unlist(x = lapply(X = maximum_possible_differences, sum))

      # Find NaNs (divide by zero errors for when there are no characters in common in a pairwise comparison):
      not_a_number <- which(x = is.nan(nonlegacy_s_ijk_bar))

      # If using Wills version replace NaNs with NA:
      if (ged_type == "wills" && length(x = not_a_number) > 0) nonlegacy_s_ijk_bar[not_a_number] <- NA

      # If using Hybrid replace NaNs with single global mean distance value:
      if (ged_type == "hybrid" && length(x = not_a_number) > 0) nonlegacy_s_ijk_bar[not_a_number] <- nonwills_s_ijk_bar[not_a_number]

      # Set modified non-legacy s_ijk_bar as main s_ijk_bar:
      s_ijk_bar <- nonlegacy_s_ijk_bar
    }

    # If using Legacy set nonwills_s_ijk_bar as main s_ijk_bar:
    if (ged_type == "legacy") s_ijk_bar <- nonwills_s_ijk_bar

    # For each set of differences:
    for (i in seq(from = 1, to = nrow(ged_data) - 1, length.out = length(x = absolute_differences))) {

      # Find missing distances (if any):
      missing_distances <- which(x = is.na(ged_data[i, ]))

      # Replace missing distances with s_ijk_bar (i.e., results of equation 2 in Wills 2001 into equation 1 of Wills 2001):
      if (length(x = missing_distances) > 0) ged_data[i, missing_distances] <- s_ijk_bar[ceiling(i / 2)]
    }

    # Isolate the distances:
    s_ijk <- ged_data[which(x = (1:nrow(ged_data) %% 2) == 1), ]

    # Isolate the weights:
    w_ijk <- ged_data[which(x = (1:nrow(ged_data) %% 2) == 0), ]

    # Calculate the GED (equation 1 of Wills 2001) for each pairwise comparison (ij):
    ged_ij <- sqrt(apply(w_ijk * (s_ijk^2), 1, sum))

    # Create GED distance matrix:
    distance_matrix <- convert_list_to_matrix(as.list(x = ged_ij), cladistic_matrix)
  }

  # Build comparable characters matrix:
  comparable_character_matrix <- convert_list_to_matrix(lapply(X = comparable_character_list, length), cladistic_matrix, diag = apply(cladistic_matrix, 1, count_complete))
  
  # Build comparable weights matrix:
  comparable_weights_matrix <- convert_list_to_matrix(lapply(X = comparable_character_list, FUN = function(x) sum(x = character_weights[x])), cladistic_matrix, diag = unname(obj = apply(X = cladistic_matrix, MARGIN = 1, FUN = function(x) sum(x = character_weights[which(x = !is.na(x = x))]))))
  
  # Add row and column names (taxa) to distance matrices:
  rownames(x = distance_matrix) <- colnames(x = distance_matrix) <- rownames(x = comparable_character_matrix) <- colnames(x = comparable_character_matrix) <- rownames(x = comparable_weights_matrix) <- colnames(x = comparable_weights_matrix) <- rownames(x = cladistic_matrix)

  # If there are any NaNs replace with NAs:
  if (any(is.nan(distance_matrix))) distance_matrix[is.nan(distance_matrix)] <- NA

  # If using a proportional distance:
  if (distance_metric == "mord" || distance_metric == "gc") {

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

  # Return compiled output:
  list(distance_metric = distance_metric, distance_matrix = distance_matrix, comparable_character_matrix = comparable_character_matrix, comparable_weights_matrix = comparable_weights_matrix)
}
