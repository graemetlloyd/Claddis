#' Classify a costmatrix character
#'
#' @description
#'
#' Given a costmatrix, classifies it as one of twelve distinct character types.
#'
#' @param costmatrix An object of class \code{costMatrix}.
#'
#' @details
#'
#' Type I - constant (any type but must be one state)
#' Type II - symmetric binary (any type but must be two states)
#' Type III - multistate unordered (three plus states and unordered type)
#' Type IV - linear ordered symmetric (three plus states and ordered)
#' Type V - Non-linear ordered symmetric (four plus states and custom direct costs all equal)
#' Type VI - Binary irreversible (two states and irreversible type)
#' Type VII - Multistate irreversible (three states and irreversible type)
#' Type VIII - Binary Dollo (two states and Dollo type)
#' Type IX - Multistate Dollo (three plus states and Dollo type)
#' Type X - Multistate custom symmetric (three plus states, symmetric but variable direct transition costs)
#' Type XI - Binary custom asymmetric (two states, custom type costs unequal and non-infinite)
#' Type XII - Multistate custom asymmetric (three plus states, custom type, asymmetric costs and not Type VII)
#'
#' # For whole cladistic data set:
#' Stratigraphic is something else(?)
#' Type XIII - Gap weighted (or is this just a linear ordered? I guess issue is when N states > N tips?)
#' Type XIV - Continuous
#'
#' classify_character # For costMatrix
#' classify_characters # For cladisticMatrix
#'
#' @return A vector of named edges (X->Y) with their distances. The sum of this vector is the length of the minimum spanning tree.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Example of a Type I character:
#' classify_costmatrix(
#'   costmatrix = make_costmatrix(
#'     min_state = 0,
#'     max_state= 0,
#'     character_type = "unordered"
#'   )
#' )
#'
#' # Example of a Type II character:
#' classify_costmatrix(
#'   costmatrix = make_costmatrix(
#'     min_state = 0,
#'     max_state= 1,
#'     character_type = "unordered"
#'   )
#' )
#'
#' # Example of a Type III character:
#' classify_costmatrix(
#'   costmatrix = make_costmatrix(
#'     min_state = 0,
#'     max_state= 3,
#'     character_type = "unordered"
#'   )
#' )
#'
#' # Example of a Type IV character:
#' classify_costmatrix(
#'   costmatrix = make_costmatrix(
#'     min_state = 0,
#'     max_state= 3,
#'     character_type = "ordered"
#'   )
#' )
#'
#' # Example of a Type V character:
#' classify_costmatrix(
#'   costmatrix = convert_adjacency_matrix_to_costmatrix(
#'     adjacency_matrix = matrix(
#'       data = c(
#'         0, 1, 0, 0,
#'         1, 0, 1, 1,
#'         0, 1, 0, 0,
#'         0, 1, 0, 0
#'       ),
#'       nrow = 4,
#'       dimnames = list(0:3, 0:3)
#'     )
#'   )
#' )
#'
#' # Example of a Type VI character:
#' classify_costmatrix(
#'   costmatrix = make_costmatrix(
#'     min_state = 0,
#'     max_state= 1,
#'     character_type = "stratigraphy",
#'     state_ages = c(103.7, 99.6)
#'   )
#' )
#'
#' # Example of a Type VII character:
#' classify_costmatrix(
#'   costmatrix = make_costmatrix(
#'     min_state = 0,
#'     max_state= 3,
#'     character_type = "irreversible"
#'   )
#' )
#'
#' # Example of a Type VIII character:
#' classify_costmatrix(
#'   costmatrix = make_costmatrix(
#'     min_state = 0,
#'     max_state= 1,
#'     character_type = "dollo"
#'   )
#' )
#'
#' # Example of a Type IX character:
#' classify_costmatrix(
#'   costmatrix = make_costmatrix(
#'     min_state = 0,
#'     max_state= 3,
#'     character_type = "dollo"
#'   )
#' )
#'
#' # Example of a Type X character:
#'
#' # Example of a Type XI character:
#'
#' # Example of a Type XII character:
#' classify_costmatrix(
#'   costmatrix = make_costmatrix(
#'     min_state = 0,
#'     max_state= 3,
#'     character_type = "stratigraphy",
#'     state_ages = c(103.7, 99.6, 91.0, 77.2)
#'   )
#' )
#'
#' @export classify_costmatrix
classify_costmatrix <- function(costmatrix) {
  
  # If a constant character return Type I:
  if (costmatrix$n_states == 1) return("Type I")
  
  # If a bianry character:
  if (costmatrix$n_states == 2) {
    
    # If symmetric return Type II:
    if (costmatrix$symmetry == "Symmetric") return("Type II")
    
    # If irreversible return Type VI:
    if (costmatrix$type == "irreversible" || costmatrix$type == "stratigraphy") return("Type VI")
    
    # If Dollo then return Type VIII:
    if (costmatrix$type == "dollo") return("Type VIII")
    
    # Extra check not a custom entered irreversible, return Type VI:
    if (any(costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states] == Inf)) return("Type VI")
    
    # Only other option is Type XI so return that:
    return("Type XI")
  }
  
  # If a multistate character:
  if (costmatrix$n_states > 2) {
    
    # If trasnition costs are all the same then return Type III:
    if (length(x = unique(x = costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states][costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states] > 0])) == 1) return("Type III")

    # If ordered then return Type VIII:
    if (costmatrix$type == "ordered") return("Type IV")

    # If multistate irreversible then return Type VIII:
    if (costmatrix$type == "irreversible") return("Type VII")
    
    # If Dollo then return Type IX:
    if (costmatrix$type == "dollo") return("Type IX")
    
    # If still not identifed then convert to stategraph to allow further checks:
    stategraph <- convert_costmatrix_to_stategraph(costmatrix = costmatrix)
    vertices <- stategraph$vertices[match(costmatrix$single_states, stategraph$vertices$label), ]
    arcs <- stategraph$arcs[which(
      x = apply(
        X = cbind(
          do.call(
            what = cbind,
            args = lapply(
              X = as.list(x = costmatrix$single_states),
              FUN = function(state) stategraph$arcs$from == state
            )
          ),
          do.call(
            what = cbind,
            args = lapply(
              X = as.list(x = costmatrix$single_states),
              FUN = function(state) stategraph$arcs$to == state
            )
          )
        ),
        MARGIN = 1,
        FUN = sum
      ) == 2
    ), ]
    
    # Check for a custom ordered character and return Type IV if meets criteria:
    if (all(x = c(vertices$in_degree, vertices$out_degree) <= 2) && length(x = unique(x = arcs$weight)) == 1 && costmatrix$symmetry == "Symmetric") return("Type IV")
    
    # Check for a custom irreversible character and return Type VII if meets criteria:
    if (all(x = !c(duplicated(x = arcs$from), duplicated(x = arcs$to))) && length(x = unique(x = arcs$weight)) == 1) return("Type VII")
    
    # Check for non-linear ordered and return Type V if found:
    if (any(x = c(vertices$in_degree, vertices$out_degree) > 2) && length(x = unique(x = arcs$weight)) == 1 && costmatrix$symmetry == "Symmetric") return("Type V")

    # Check for custom symmetric and return Type X if found:
    if (length(x = unique(x = arcs$weight)) > 1 && costmatrix$symmetry == "Symmetric") return("Type X")
    
    # Anything remaining with asymmetric costs must be Type XII:
    if (costmatrix$symmetry == "Asymmetric") return("Type XII")
    
    # This should not happen but is left here as an error check:
    stop("Clasification has inexplicably failed!")
  }
}
