#' Check a stateGraph object for errors
#'
#' @description
#'
#' Internal function to check a stateGraph object for errors.
#'
#' @param stategraph A stateGraph object.
#'
#' @details
#'
#' Stategraph objects are more complex than what will typically be shown to the user. This function checks this hidden structure and reports any errors it finds.
#'
#' These checks include rules 1-7 from Hoyal Cuthill and Lloyd (i prep.).
#'
#' @return An error message or empty vector if no errors found.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Make an unordered costmatrix:
#' costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered"
#' )
#'
#' # Convert costmatrix to stategraph:
#' stategraph <- convert_costmatrix_to_stategraph(costmatrix = costmatrix)
#'
#' # Check that this is a valid stateGraph object (should return empty vector):
#' check_stateGraph(stategraph = stategraph)
#'
#' @export check_stateGraph
check_stateGraph <- function(stategraph) {
  
  # Could add more checks to arcs, vertices or adjacency_matrix
  
  # Check stategraph has class stateGraph and add error message to output if true:
  if (!inherits(x = stategraph, what = "stateGraph")) return("stategraph must be an object of class \"stateGraph\".")
  
  # Check stategraph is in form of list and add error message to output if false:
  if (!is.list(x = stategraph)) return("stategraph must be in the form of a list.")
  
  # Check length of list is at least 20 and add error message to output if false:
  if (length(x = stategraph) < 20) return("stategraph must be a list with 20 items (n_vertices, n_arcs, n_states, single_states, type, arcs, vertices, radius, diameter, adjacency_matrix, directed, includes_polymorphisms, polymorphism_costs, polymorphism_geometry, polymorphism_distance, includes_uncertainties, pruned, dollo_penalty, base_age, and weight).")
  
  # Check names are correct and in order add error message to output if false:
  if (!all(x = names(x = stategraph) == c("n_vertices", "n_arcs", "n_states", "single_states", "type", "arcs", "vertices", "radius", "diameter", "adjacency_matrix", "directed", "includes_polymorphisms", "polymorphism_costs", "polymorphism_geometry", "polymorphism_distance", "includes_uncertainties", "pruned", "dollo_penalty", "base_age", "weight"))) return("Elements of stategraph must be \"n_vertices\", \"n_arcs\", \"n_states\", \"single_states\", \"type\", \"arcs\", \"vertices\", \"radius\", \"diameter\", \"adjacency_matrix\", \"directed\", \"includes_polymorphisms\", \"polymorphism_costs\", \"polymorphism_geometry\", \"polymorphism_distance\", \"includes_uncertainties\", \"pruned\", \"dollo_penalty\", \"base_age\", and \"weight\" in that order.")
  
  # Check n_vertices is formatted correctly and add error message to output if false:
  if (!is.numeric(x = stategraph$n_vertices) || length(x = stategraph$n_vertices) != 1) return("stategraph$n_vertices should be a single numeric value indicating the number of vertices in the stategraph.")
  
  # Check n_arcs is formatted correctly and add error message to output if false:
  if (!is.numeric(x = stategraph$n_arcs) || length(x = stategraph$n_arcs) != 1) return("stategraph$n_arcs should be a single numeric value indicating the number of arcs in the stategraph.")

  # Check n_states is formatted correctly and add error message to output if false:
  if (!is.numeric(x = stategraph$n_states) || length(x = stategraph$n_vertices) != 1) return("stategraph$n_states should be a single numeric value indicating the number of singles states (i.e., excluding polymorphic or uncertain values).")

  # Check single_states are formatted correctly and match other elements and add error message to output if false:
  if (!is.character(x = stategraph$single_states) || !is.vector(x = stategraph$single_states) || length(x = stategraph$single_states) != stategraph$n_states || any(x = is.na(x = match(x = stategraph$single_states, table = stategraph$vertices[, "label"]))) || length(x = grep(pattern = "/|&", x = stategraph$single_states)) > 0) return("stategraph$single_states must be a single character vector value equal in length to stategraph$n_states, contain no polymorphic or uncertain values and match the states used in stategraph$vertices.")
  
  # Check type is formatted correctly and add error message to output if false:
  if (!is.character(x = stategraph$type) || length(x = stategraph$type) != 1) return("stategraph$type should be a single character value indicating the type of stategraph.")
  
  # Check type is from the limited available list and add error message to output if false:
  if (length(x = setdiff(x = stategraph$type, y = c("ordered", "unordered", "dollo", "irreversible", "stratigraphy", "custom"))) > 0) return("stategraph$type must be one of: \"ordered\", \"unordered\", \"dollo\", \"irreversible\", \"stratigraphy\", \"custom\".")

  # Check arcs is a data frame:
  if (!is.data.frame(x = stategraph$arcs)) return("stategraph$arcs should be a data frame.")
  
  # Check names are correct and in order add error message to output if false:
  if (!all(x = colnames(x = stategraph$arcs) == c("from", "to", "weight"))) return("Elements of stategraph$arcs must be \"from\", \"to\", and \"weight\" in that order.")

  # Check vertices is a data frame:
  if (!is.data.frame(x = stategraph$vertices)) return("stategraph$vertices should be a data frame.")
  
  # Check names are correct and in order add error message to output if false:
  if (!all(x = colnames(x = stategraph$vertices) == c("label", "in_degree", "out_degree", "eccentricity", "periphery", "centre"))) return("Elements of stategraph$arcs must be \"label\", \"in_degree\", \"out_degree\", \"eccentricity\", \"periphery\", and \"centre\" in that order.")
  
  # Check radius is formatted correctly and add error message to output if false:
  if (!is.numeric(x = stategraph$radius) || length(x = stategraph$radius) != 1) return("stategraph$radius should be a single numeric value indicating the number of vertices in the stategraph.")
  
  # Check radius matches vertex eccentricity:
  if (stategraph$radius != min(x = stategraph$vertices[, "eccentricity"])) return("stategraph$radius should be the minimum eccentrcity of the vertices of the stategraph.")

  # Check diameter is formatted correctly and add error message to output if false:
  if (!is.numeric(x = stategraph$diameter) || length(x = stategraph$diameter) != 1) return("stategraph$diameter should be a single numeric value indicating the number of vertices in the stategraph.")
  
  # Check diameter matches vertex eccentricity:
  if (stategraph$diameter != max(x = stategraph$vertices[, "eccentricity"])) return("stategraph$diameter should be the maximum eccentrcity of the vertices of the stategraph.")

  # Check adjacency_matrix is a matrix:
  if (!is.matrix(x = stategraph$adjacency_matrix)) return("stategraph$adjacency_matrix should be a matrix.")

  # Check directed is formatted correctly and add error message to output if false:
  if (!is.logical(x = stategraph$directed) || length(x = stategraph$directed) != 1) return("stategraph$directed should be a single logical value indicating whether graph is directed or not.")

  # Is graph actually directed?:
  is_directed <- ifelse(test = length(x = setdiff(x = paste(stategraph$arcs[, "from"], stategraph$arcs[, "to"], sep = "_to_"), y = paste(stategraph$arcs[, "to"], stategraph$arcs[, "from"], sep = "_to_"))) == 0, yes = ifelse(test = any((rle(x = stategraph$arcs[order(x = apply(X = stategraph$arcs[, c("from", "to"), drop = FALSE], MARGIN = 1, FUN = function(i) paste(sort(x = i), collapse = ""))), "weight"])$lengths %% 2) > 0), yes = TRUE, no = FALSE), no = TRUE)

  # Check stategraph$directed actually matches whether graph is directed:
  if (stategraph$directed && !is_directed) return("stategraph$directed is set to TRUE, but graph is not directed.")
  if (!stategraph$directed && is_directed) return("stategraph$directed is set to FALSE, but graph is directed.")

  # Check includes_polymorphisms is formatted correctly and add error message to output if false:
  if (!is.logical(x = stategraph$includes_polymorphisms) || length(x = stategraph$includes_polymorphisms) != 1) return("stategraph$includes_polymorphisms should be a single logical value indicating whether polymorphisms are included or not.")
  
  # Check that if includes_polymorphisms is TRUE there actually are polymorphisms and add error message to output if not:
  if (stategraph$includes_polymorphisms && length(x = grep(pattern = "&", x = stategraph$vertices[, "label"])) == 0) return("stategraph$includes_polymorphisms is set to TRUE, but there are no actual polymorphisms in stategraph$vertices.")
  
  # Check that if includes_polymorphisms is FALSE there actually are no polymorphisms and add error message to output if not:
  if (!stategraph$includes_polymorphisms && length(x = grep(pattern = "&", x = stategraph$vertices[, "label"])) > 0) return("stategraph$includes_polymorphisms is set to FALSE, but there are actual polymorphisms in stategraph$vertices.")
  
  # Check polymorphism_costs is a valid value and stop and add error message to output if not:
  if (length(x = stategraph$polymorphism_costs) != 1 || !is.character(x = stategraph$polymorphism_costs) || length(x = setdiff(x = stategraph$polymorphism_costs, y = c("additive", "geometric", "maddison", "stratigraphic"))) > 0) return("stategraph$polymorphism_costs must be a single character value and one of \"additive\", \"geometric\", \"maddison\", \"stratigraphic\".")
  
  # Check polymorphism_geometry is a valid value and stop and add error message to output if not:
  if (length(x = stategraph$polymorphism_geometry) != 1 || !is.character(x = stategraph$polymorphism_geometry) || length(x = setdiff(x = stategraph$polymorphism_geometry, y = c("hypercube", "hypersphere", "simplex"))) > 0) stop("stategraph$polymorphism_geometry must be a single character value and one of \"hypercube\", \"hypersphere\", \"simplex\".")
  
  # Check polymorphism_distance is a valid value and stop and add error message to output if not:
  if (length(x = stategraph$polymorphism_distance) != 1 || !is.character(x = stategraph$polymorphism_distance) || length(x = setdiff(x = stategraph$polymorphism_distance, y = c("manhattan", "euclidean", "great_circle"))) > 0) stop("stategraph$polymorphism_distance must be a single character value and one of \"manhattan\", \"euclidean\", \"great_circle\".")
  
  # Check includes_polymorphisms is formatted correctly and add error message to output if false:
  if (!is.logical(x = stategraph$includes_uncertainties) || length(x = stategraph$includes_uncertainties) != 1) return("stategraph$includes_uncertainties should be a single logical value indicating whether uncertainties are included or not.")
  
  # Check that if includes_uncertainties is TRUE there actually are uncertainties and add error message to output if not:
  if (stategraph$includes_uncertainties && length(x = grep(pattern = "/", x = stategraph$vertices[, "label"])) == 0) return("stategraph$includes_uncertainties is set to TRUE, but there are no actual uncertainties in stategraph$vertices.")
  
  # Check that if includes_uncertainties is FALSE there actually are no uncertainties and add error message to output if not:
  if (!stategraph$includes_uncertainties && length(x = grep(pattern = "/", x = stategraph$vertices[, "label"])) > 0) return("stategraph$includes_uncertainties is set to FALSE, but there are actual uncertainties in stategraph$vertices.")
  
  # Check pruned is formatted correctly and add error message to output if false:
  if (!is.logical(x = stategraph$pruned) || length(x = stategraph$pruned) != 1) return("stategraph$pruned should be a single logical value indicating whether the graph is pruned (not all states are included) or not.")
  
  # Check dollo_penalty is formatted correctly and add error message to output if false:
  if (!is.numeric(x = stategraph$dollo_penalty) || length(x = stategraph$dollo_penalty) != 1 || stategraph$dollo_penalty <= 0 || stategraph$dollo_penalty == Inf) return("stategraph$dollo_penalty should be a single finite positive numeric value that indicates the penalty to be used for a Dollo character.")
  
  # Check base_age is formatted correctly and add error message to output if false:
  if (!is.numeric(x = stategraph$base_age) || length(x = stategraph$base_age) != 1 || stategraph$base_age <= 0 || stategraph$base_age == Inf) return("stategraph$base_age should be a single finite positive numeric value that indicates the oldest age for a stratigrphic character.")
  
  # Check weight is formatted correctly and add error message to output if false:
  if (!is.numeric(x = stategraph$weight) || length(x = stategraph$weight) != 1 || stategraph$weight < 0 || stategraph$weight == Inf) return("stategraph$weight should be a single finite non-negative numeric value that indictaes the weight the character has.")
  
  # Check stategraph makes sense (but only if custom):
  if (stategraph$type == "custom") {
    
    # Convert stategraph to costmatrix:
    costmatrix <- convert_stategraph_to_costmatrix(stategraph = stategraph)
    
    # Fix costmatrix so paths are shortest:
    shortest_path_costmatrix <- fix_costmatrix(costmatrix = costmatrix, message = FALSE)
    
    # RULE SEVEN: If costMatrix is not all shortest paths add error message to output and warn user:
    if (!all(x = shortest_path_costmatrix$costmatrix == costmatrix$costmatrix)) return("stategraph is not self-consistent (at least one path is shorter - lower cost - than stated). Fix using convert_stategraph_to_costmatrix, fix_costmatrix, and convert_costmatrix_to_stategraph and then try again.")
  }

  # Return empty vector:
  return(vector(mode = "character"))
}
