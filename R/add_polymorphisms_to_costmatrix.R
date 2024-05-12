#' Adds polymorphisms to a costmatrix
#'
#' @description
#'
#' Given a costmatrix and choice of method, adds polymorphic transitions to a costmatrix.
#'
#' @param costmatrix An object of class "costMatrix".
#' @param polymorphism_costs The method to use to assign costs to the transitions to and from polymorphic states. Must be one of \code{"additive"} (the default), \code{"geometric"}, \code{"maddison"}, or \code{"stratigraphic"}. See details.
#' @param polymorphism_geometry If using \code{polymorphism_costs = "geometric"}, the type of geometry to use. Must be one of  \code{"hypercube"}, \code{"hypersphere"}, or  \code{"simplex"} (the default). See details.
#' @param polymorphism_distance If using \code{polymorphism_costs = "geometric"}, the type of distance to use. Must be one of \code{"euclidean"} (the default), \code{"great_circle"}, or \code{"manhattan"}. See details.
#' @param message Logical indicating whether (\code{TRUE}, the default) or not (\code{FALSE}) to provide messages about the process of the function.
#'
#' @details
#'
#' Polymorphisms - the presence of two or more discrete character states in an operational taxonomic unit - represent perhaps the most complex problem for phylogenetic inference and ancestral state estimation. This is at least in part due to the varied nature their meaning can have (Swofford and Maddison 1992). For example, they may represent true variance within a species or composite variance in some higher taxon represented by a single OTU (something to be avoided if at all possible). Consequently, they cannot - and probably should not - be treated the same by all software in all situations.
#'
#' One solution to the problem of polymorphisms is to pretend they do not exist. This is the approach used by common software such as TNT (Goloboff et al. 2008; Goloboff and Catalano 2016), which treat polymorphisms as uncertainties instead. I.e., they assume only one state is the "true" value and use parsimony to (implicitly) estimate what this state is. This approach can be adopted in Claddis too by simply selecting options that treat polymorphisms as uncertainties instead. However, misrepresenting true polymorphisms as uncertainties can leads to incorrect estimates of the amount of evolution (tree length under parsimony; Nixon and Davis 1991) and can disallow, potentially inappropriately, polymorphic values as ancestral state estimates (reconstructions). Claddis therefore offers options to treat polymorphisms as "real".
#'
#' There is no single agreed method on dealing with "real" polymorphisms, although many have been proposed (Wiens 1999). Due to the constraints of how Claddis works the only options implemented here - aside from simple approaches like treating them as missing data - are those that work in the costmatrix setting. These are a mixture of novel approaches, extensions of other approaches, or taken directly from Maddison and Maddison (1987).
#'
#' Generally speaking the costmatrix setting makes explicit the cost of individual transitions (e.g., from the single state 0 to the polymorphic state 1 \emph{and} 2) and hence the only decision to make is how these costs should be assigned. Here three approaches are offered: 1) the generalised additive approach, 2) the geometric approach (for unordered characters), and 3) the Maddison approach (for ordered characters). These are described in more detail below.
#'
#' \bold{The additive approach}
#'
#' One way of conceiving of polymorphic tip states is to conceptually "split" a terminal branch into multiple terminals, each with a single state from the polymorphism at its' tip. For example, a branch leading to the state "0&1" would become two branches, one leading to state 0 and one to state 1. Character costs can then be calculated for each branch separately and then, in order to "merge" these into a single branch, the sum of these costs can be taken. This approach was first suggested in Hoyal Cuthill and Lloyd (in prep) and is here termed the "additive" approach (\code{polymorphism_costs = "additive"}). A major advantage to this approach is that it's rules translate across the broadest range of character types. However, there is no logical interpretation of polymorphic ancestral states. (I.e., splitting an ancestral node into two or more nodes is not possible in the same way.) Consequently, transitions "from" polymorphic states are assigned a cost of infinity, precluding them from consideration as ancestral values.
#'
#' An example of the additive approach, applied to an unordered character with three states (0, 1, and 2) would thus appear as the following costmatrix:
#'
#' \preformatted{        ---------------------------------------------
#'         |  0  |  1  |  2  | 0&1 | 0&2 | 1&2 | 0&1&2 |
#' -----------------------------------------------------
#' |   0   |  0  |  1  |  1  |  1  |  1  |  2  |   2   |
#' -----------------------------------------------------
#' |   1   |  1  |  0  |  1  |  1  |  2  |  1  |   2   |
#' -----------------------------------------------------
#' |   2   |  1  |  1  |  0  |  2  |  1  |  1  |   2   |
#' -----------------------------------------------------
#' |  0&1  | Inf | Inf | Inf |  0  | Inf | Inf |  Inf  |
#' -----------------------------------------------------
#' |  0&2  | Inf | Inf | Inf | Inf |  0  | Inf |  Inf  |
#' -----------------------------------------------------
#' |  1&2  | Inf | Inf | Inf | Inf | Inf |  0  |  Inf  |
#' -----------------------------------------------------
#' | 0&1&2 | Inf | Inf | Inf | Inf | Inf | Inf |   0   |
#' -----------------------------------------------------}
#'
#' \bold{The geometric approach}
#'
#' An alternative approach based on an idea first suggested by Maddison and Maddison (1987) does allow polymorphisms to be treated as valid ancestral states. Here the conceptual idea is that transitions involve turning "on" or "off" particular states. Thus to go from state 1 to state 0&2 would involve turning "on" states 0 and 2 turning "off" state 1, giving a total cost of three. This approach is only suitable where costs between single states are equal and symmetric - i.e., the character is unordered.
#'
#' An example of this approach, applied to an unordered character with three states (0, 1, and 2) would thus appear as the following costmatrix:
#'
#' \preformatted{        ---------------------------------------------
#'         |  0  |  1  |  2  | 0&1 | 0&2 | 1&2 | 0&1&2 |
#' -----------------------------------------------------
#' |   0   |  0  |  2  |  2  |  1  |  1  |  3  |   2   |
#' -----------------------------------------------------
#' |   1   |  2  |  0  |  2  |  1  |  3  |  1  |   2   |
#' -----------------------------------------------------
#' |   2   |  2  |  2  |  0  |  3  |  1  |  1  |   2   |
#' -----------------------------------------------------
#' |  0&1  |  1  |  1  |  3  |  0  |  2  |  2  |   1   |
#' -----------------------------------------------------
#' |  0&2  |  1  |  3  |  1  |  2  |  0  |  2  |   1   |
#' -----------------------------------------------------
#' |  1&2  |  3  |  1  |  1  |  2  |  2  |  0  |   1   |
#' -----------------------------------------------------
#' | 0&1&2 |  2  |  2  |  2  |  1  |  1  |  1  |   0   |
#' -----------------------------------------------------}
#'
#' Note: to maintain a base cost of one the original single-state-to-single-state transitions are doubled from one to two. In order to avoid this affecting tree lengths the character weight is halved.
#'
#' This Maddison and Maddison (1987) approach is extended here by drawing on an analogy, namely that this cost assignment is identical to considering a hypercube whose vertices lie at one of two (presence or absence) positions on orthogonal axes representing individual states. Thus the polymorphism 0&2 is at the coordinates (0 = present, 1 = absent, 2 = present). Here costs can be assigned by taking the minimum Manhattan distance between each pair of vertices. NB: here the origin (absence of all states) is unoccupied, but this does not affect the assignment of costs.
#'
#' This analogy is why the approach is here termed the "geometric" approach (\code{polymorphism_costs = "geometric"}). Importantly, although the example above only represents the regular cube (three states = three dimensions) the approach translates to any number of states (i.e., any number of dimensions) without loss of generality, and hence this particular option is properly termed a hypercube (i.e., \code{polymorphism_geometry = "hypercube"} and \code{polymorphism_distance = "manhattan"}).
#'
#' In of itself this does not change what was proposed by Maddison and Maddison. However, the geometric setting allows us to consider both different shapes (topologies) and different distance measures, as long as the same generality to higher dimensions holds. Two additional shapes are offered here. These are the hypersphere (the N-dimensional circle) and the simplex (the N-dimensional equilateral triangle). In addition, two distances are also offered. These are Euclidean (the straight line distance between two points) and Great circle (the shortest distance across the "surface" of a hypersphere) distances.
#'
#' These distances are intended to match respective straights (hypercube and Manhattan, hypersphere and Great Circle, simplex and Euclidean), but any combination is permitted. Furthermore, and as above, all costs are rescaled such that the base transition cost is always one with the character weight modified accordingly.
#'
#' \bold{The Maddison approach}
#'
#' Although the above approach is based on Maddison and Maddison (1987) the Maddison name (\code{polymorphism_costs = "maddison"}) is here reserved for an approach proposed by the same authors (Maddison and Maddison 2000) for use with ordered characters. The concept of this approach shares the idea of "switching" states on and off but here adds the idea of adjacency of switches. In other words, intermediate states (the essential element of an ordered character) are additional switches that must be turned "on" and "off" again in order to "reach" other states. Thus to go from state 0 to state 0 and 2 for the linear ordered character 0-1-2 we would have to turn "on" state to reach state "2", which is also turned on. Then we would turn off state "1", meaning  total of three switches were flipped which would be the cost. Note: that here we are not turning "off" state 0 as it is present at both ends of the transition.
#'
#' Extending our three state linear ordered example to all possible transitions we get the costmatrix:
#'
#' \preformatted{        ---------------------------------------------
#'         |  0  |  1  |  2  | 0&1 | 0&2 | 1&2 | 0&1&2 |
#' -----------------------------------------------------
#' |   0   |  0  |  2  |  4  |  1  |  3  |  3  |   2   |
#' -----------------------------------------------------
#' |   1   |  2  |  0  |  2  |  1  |  3  |  1  |   2   |
#' -----------------------------------------------------
#' |   2   |  4  |  2  |  0  |  3  |  3  |  1  |   2   |
#' -----------------------------------------------------
#' |  0&1  |  1  |  1  |  3  |  0  |  2  |  2  |   1   |
#' -----------------------------------------------------
#' |  0&2  |  3  |  3  |  3  |  2  |  0  |  2  |   1   |
#' -----------------------------------------------------
#' |  1&2  |  3  |  1  |  1  |  2  |  2  |  0  |   1   |
#' -----------------------------------------------------
#' | 0&1&2 |  2  |  2  |  2  |  1  |  1  |  1  |   0   |
#' -----------------------------------------------------}
#'
#' NB: Again the single-state-to-single-state cost is modified (doubled) meaning the character weight is halved.
#'
#' \bold{The stratigraphic approach}
#'
#' The stratigraphic case seems like it ought to follow the same rule as irreversible characters (see below), but here "polymorphisms" have a different logical meaning. Specifically they encode the presence of a taxon in multiple time units. Thus, transitions \emph{to} a stratigraphic polymorphism are based on the oldest state. By contrast, polymorphisms can never be ancestral states as there would be no basis to choose between (e.g.) state 0 and state 0 and 1.
#'
#' The option \code{polymorphism_costs = "stratigraphic"} is only available where \code{costmatrix$type} is \code{"stratigraphic"} and, where \code{costmatrix$type} is \code{"stratigraphic"}, \code{polymorphism_costs = "stratigraphic"} is the only option.
#'
#' \bold{Additional considerations}
#'
#' Aside from these three options there are some further considerations the user should make with respect to some more complex character types. These are discussed below under appropriate headers.
#'
#' \emph{Dollo and irreversible polymorphisms}
#'
#' At first blush Dollo and irreversible (Camin-Sokal) characters appear to confound costmatrix construction. More specifically, for a Dollo character, the state 0 and 1 would suggest state 1 was acquired previously, but lost in some member(s) of the OTU. For an irreversible character the opposite would be true, and we would assume that state 1 was acquired independently by some member(s) of the OTU and state 0 was acquired previously. It is not possible to code either of these outcomes under the geometric or Maddison approaches. However, the additive approach perfectly captures this nuance. Thus tree length calculation and ancestral state estimation for Dollo and Camin-Sokal characters is possible under the additive approach and only the possibility of polymorphic ancestral states is excluded.
#'
#' \emph{Transitions between uncertainties and polymorphisms}
#'
#' Because Claddis also permits users to include uncertainties in costmatrices transitions between these and polymorphisms must also be considered. Here this is automated by using the "minimum rule". In practice this means a transition from, for example, state 0 and 1 to state 1 or 2 will take the lowest cost of the two options (from state 0 and 1 to state 1 \emph{or} from state 0 and 1 to state 2). In the case of the additive approach these costs will always be infinite. In all cases transitions from uncertainties are assigned infinite cost as these are not considered valid ancestral states.
#'
#' \emph{Polymorphism limits and gap-weighted characters}
#'
#' Even where polymorphisms may be appropriate (i.e., for some ordered and unordered characters) they still represent a major increase to costmatrix size that can cause memory limits to be hit. Specifically, and as shown in the hypercube example above, there will be as many possible states as N^2 - 1, where N is the number of single (i.e., non-polymorphic) states and the minus one indicates the exclusion of the unrepresentable origin value. Thus the size of costmatrix required grows very quickly:
#'
#' \preformatted{------------------------------
#' | N states | Costmatrix size |
#' ------------------------------
#' |     2    |      3 x 3      |
#' |     3    |      7 x 7      |
#' |     4    |     15 x 15     |
#' |     5    |     31 x 31     |
#' |     6    |     63 x 63     |
#' |     7    |    127 x 127    |
#' |     8    |    255 x 255    |
#' |     9    |    511 x 511    |
#' |    10    |   1023 x 1023   |
#' |    11    |   2047 x 2047   |
#' |    12    |   4095 x 4095   |
#' |    13    |   8191 x 8191   |
#' |    14    |  16383 x 16383  |
#' ------------------------------}
#'
#' Because of this, the function will become extremely slow for these higher values and here is hard-capped at no more than fourteen states. Consequently any gap-weighted characters are also most likely inappropriate for polymorphism use (as well as being too memory intensive).
#'
#' @return
#'
#' An object of class "costMatrix".
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics/ \emph{Cladistics}, \bold{32}. 221-238
#'
#' Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. \emph{Cladistics}, \bold{24}, 774-786.
#'
#' Hoyal Cuthill and Lloyd, in prep.
#'
#' Maddison, W. P. and Maddison, D. R., 1987. MacClade 2.1, computer program and manual. Cambridge, Massachusetts.
#'
#' Nixon, K. C. and Davis, J. I., 1991. Polymorphic taxa, missing values and cladistic analysis. \emph{Cladistics}, \bold{7}, 233-241.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. In R. L. Mayden (ed.), \emph{Systematics, Historical Ecology, and North American Freshwater Fishes}. Stanford University Press, Stanford. pp187-223.
#'
#' Wiens, J. J., 1999. Polymorphism in systematics and comparative biology. \emph{Annual Review of Ecology and Systematics}, \bold{30}, 327-362.
#'
#' @examples
#'
#' # Generate an example three-state unordered character costmatrix:
#' unordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered"
#' )
#'
#' # Generate an example three-state ordered character costmatrix:
#' ordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "ordered"
#' )
#'
#' # Generate an example three-state ordered character costmatrix with uncertainties already included:
#' ordered_uncertainty_costmatrix <- list(
#'   size = 7,
#'   n_states = 3,
#'   single_states = c("0", "1", "2"),
#'   type = "ordered",
#'   costmatrix = matrix(
#'     data = c(
#'       0, 1, 2, 0, 0, 1, 0,
#'       1, 0, 1, 0, 1, 0, 0,
#'       2, 1, 0, 1, 0, 0, 0,
#'       Inf, Inf, Inf, 0, Inf, Inf, Inf,
#'       Inf, Inf, Inf, Inf, 0, Inf, Inf,
#'       Inf, Inf, Inf, Inf, Inf, 0, Inf,
#'       Inf, Inf, Inf, Inf, Inf, Inf, 0
#'     ),
#'     nrow = 7,
#'     byrow = TRUE,
#'     dimnames = list(
#'       c(as.character(x = 0:2), "0/1", "0/2", "1/2", "0/1/2"),
#'       c(as.character(x = 0:2), "0/1", "0/2", "1/2", "0/1/2")
#'     )
#'   ),
#'   symmetry = "Symmetric",
#'   includes_polymorphisms = FALSE,
#'   polymorphism_costs = "additive",
#'   polymorphism_geometry = "simplex",
#'   polymorphism_distance = "euclidean",
#'   includes_uncertainties = TRUE,
#'   pruned = FALSE,
#'   dollo_penalty = 999,
#'   base_age = 1,
#'   weight = 1
#' )
#' class(ordered_uncertainty_costmatrix) <- "costMatrix"
#'
#' # Generate an example five-state stratigraphic character costmatrix:
#' stratigraphic_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 4,
#'   character_type = "stratigraphy",
#'   state_ages = c(0, 1.3, 5.3, 8.0, 11.1)
#' )
#'
#' # Add polymorphisms to unordered costmatrix using additive method:
#' unordered_costmatrix_additive_polymorphisms <- add_polymorphisms_to_costmatrix(
#'   costmatrix = unordered_costmatrix,
#'   polymorphism_costs = "additive"
#' )
#'
#' # Show unordered costmatrix using additive method:
#' unordered_costmatrix_additive_polymorphisms$costmatrix
#'
#' # Add polymorphisms to unordered costmatrix using geometric simplex method:
#' unordered_costmatrix_simplex_polymorphisms <- add_polymorphisms_to_costmatrix(
#'   costmatrix = unordered_costmatrix,
#'   polymorphism_costs = "geometric",
#'   polymorphism_geometry = "simplex",
#'   polymorphism_distance = "euclidean"
#' )
#'
#' # Show unordered costmatrix using geometric simplex method:
#' unordered_costmatrix_simplex_polymorphisms$costmatrix
#'
#' # Add polymorphisms to unordered costmatrix using geometric hypercube method:
#' unordered_costmatrix_hypercube_polymorphisms <- add_polymorphisms_to_costmatrix(
#'   costmatrix = unordered_costmatrix,
#'   polymorphism_costs = "geometric",
#'   polymorphism_geometry = "hypercube",
#'   polymorphism_distance = "manhattan"
#' )
#'
#' # Show unordered costmatrix using geometric hypercube method:
#' unordered_costmatrix_hypercube_polymorphisms$costmatrix
#'
#' # Add polymorphisms to ordered costmatrix using additive method:
#' ordered_costmatrix_additive_polymorphisms <- add_polymorphisms_to_costmatrix(
#'   costmatrix = ordered_costmatrix,
#'   polymorphism_costs = "additive"
#' )
#'
#' # Show ordered costmatrix using additive method:
#' ordered_costmatrix_additive_polymorphisms$costmatrix
#'
#' # Add polymorphisms to ordered costmatrix using maddison method:
#' ordered_costmatrix_maddison_polymorphisms <- add_polymorphisms_to_costmatrix(
#'   costmatrix = ordered_costmatrix,
#'   polymorphism_costs = "maddison"
#' )
#'
#' # Show ordered costmatrix using Maddison method:
#' ordered_costmatrix_maddison_polymorphisms$costmatrix
#'
#' # Add polymorphisms to ordered uncertainty costmatrix using additive method:
#' ordered_uncertainty_costmatrix_additive_polymorphisms <- add_polymorphisms_to_costmatrix(
#'   costmatrix = ordered_uncertainty_costmatrix,
#'   polymorphism_costs = "additive"
#' )
#'
#' # Show costmatrix, with polymorphism-to-uncertainty transitions as infinities:
#' ordered_uncertainty_costmatrix_additive_polymorphisms$costmatrix
#'
#' # Add polymorphisms to ordered uncertainty costmatrix using Maddison method:
#' ordered_uncertainty_costmatrix_maddison_polymorphisms <- add_polymorphisms_to_costmatrix(
#'   costmatrix = ordered_uncertainty_costmatrix,
#'   polymorphism_costs = "maddison"
#' )
#'
#' # Show costmatrix, with polymorphism-to-uncertainty transitions
#' # interpolated using minimum cost rule:
#' ordered_uncertainty_costmatrix_maddison_polymorphisms$costmatrix
#'
#' # Add polymorphisms to stratigraphic costmatrix using stratigraphic method:
#' stratigraphic_costmatrix_polymorphisms <- add_polymorphisms_to_costmatrix(
#'   costmatrix = stratigraphic_costmatrix,
#'   polymorphism_costs = "stratigraphic"
#' )
#'
#' # Show stratigraphic costmatrix using stratigraphic method:
#' stratigraphic_costmatrix_polymorphisms$costmatrix
#'
#' @export add_polymorphisms_to_costmatrix
add_polymorphisms_to_costmatrix <- function(
  costmatrix,
  polymorphism_costs = "additive",
  polymorphism_geometry = "simplex",
  polymorphism_distance = "euclidean",
  message = TRUE
) {
  
  # TO DO:
  # - UPDATE convert_adjacency_matrix_to_costmatrix AND OTHER FUNCTIONS FOR NEW FORMAT!
  # - ALLOW SYMMETRY TO BE ON FULL MATRIX EVEN IF SINGLE STATE METRIC MAY BE DIFFERENT? OR OPPOSITE!!!!!
  # - CHANGE SOME SYMMETRY CHECKS IF WHAT THEY REALLY CARE ABOUT IS SINGLE STATE TO SINGLE STATE TRANSITIONS ONLY!!!!
  
  ### FIX COSTMATRIX WILL HAVE TO OPERATE ON NORMAL MATRIX (POSSIBLY RESCALED FIRST TO CHECK) THEN RE-ADD POLYMORPHISMS AND UNCERTAINTIE USING SAME RULES AS ORIGINALLY APPLIED MEANING WILL ALSO NEED TO STORE THESE IN COSTMATRIX FORMAT!!!!!
  
  ### NEED TO ADD GEOMETRY OF STATE GRAPH FOR CHECKING WHAT IS ALLOWABLE? DOES THAT MEAN NEED COORDINATES OF VERTICES IN SPACE TOO?
  
  ### CHECK costmatrix is of class costmatrix! (Need to rewrite that function for new structure)
  
  # Check if costmatrix is size 1 (no possibility for polymorphisms):
  if (costmatrix$size == 1) {
    
    # Message user about this if requested:
    if (message) print("Costmatrix is size 1 meaning no polymorphisms are possible. Returning original costmatrix.")
    
    # Return original costmatrix:
    return(value = costmatrix)
  }

  # Check polymorphism_costs is a valid value and stop and warn user if not:
  if (length(x = setdiff(x = polymorphism_costs, y = c("additive", "geometric", "maddison", "stratigraphic"))) > 0) stop("polymorphism_costs must be one of \"additive\", \"geometric\", \"maddison\", \"stratigraphic\".")
  
  # Check polymorphism_geometry is a valid value and stop and warn user if not:
  if (length(x = setdiff(x = polymorphism_geometry, y = c("hypercube", "hypersphere", "simplex"))) > 0) stop("polymorphism_geometry must be one of \"hypercube\", \"hypersphere\", \"simplex\".")
  
  # Check polymorphism_distance is a valid value and stop and warn user if not:
  if (length(x = setdiff(x = polymorphism_distance, y = c("manhattan", "euclidean", "great_circle"))) > 0) stop("polymorphism_distance must be one of \"manhattan\", \"euclidean\", \"great_circle\".")
  
  # If character is stratigraphic but polymorphism costs are not:
  if (polymorphism_costs != "stratigraphic" && costmatrix$type == "stratigraphic") {
    
    # Change polymorphism_costs to "stratigraphic":
    polymorphism_costs <- "stratigraphic"
    
    # If message is TRUE inform user of this change:
    if (message) print("polymorphism_costs must be \"stratigraphic\" if character type is \"stratigraphic\". polymorphism_costs changed to \"stratigraphic\" instead.")
  }

  # If polymorphism_costs is geometric but character type is not unordered:
  if (polymorphism_costs == "geometric" && costmatrix$type != "unordered") {
    
    # Reset polymorphism_costs as additive:
    polymorphism_costs <- "additive"
    
    # If message is TRUE inform user of this change:
    if (message) print("polymorphism_costs cannot be \"geometric\" if character type is not \"unordered\". polymorphism_costs changed to \"additive\" instead.")
  }
  
  # Create variable for whether character is possible ordered (Type I, Type II, Type IV or Type V):
  character_is_ordered <- ifelse(test = length(x = intersect(x = c("Type I", "Type II", "Type IV", "Type V"), y = classify_costmatrix(costmatrix = costmatrix))) == 1, yes = TRUE, no = FALSE)
  
  # If polymorphism_costs is maddison but character type is not linear or non-linear ordered:
  if (polymorphism_costs == "maddison" && !character_is_ordered) {
    
    # Reset polymorphism_costs as additive:
    polymorphism_costs <- "additive"
    
    # If message is TRUE inform user of this change:
    if (message) print("polymorphism_costs cannot be \"maddison\" if character type is not ordered (i.e., linear ordered or state tree). polymorphism_costs changed to \"additive\ instead.")
  }
  
  # Check if polymorphisms are already present:
  if (costmatrix$includes_polymorphisms) {
    
    # Message user about this if requested:
    if (message) print("costmatrix already includes polymorphisms. Returning original costmatrix.")
    
    # Return original costmatrix:
    return(value = costmatrix)
  }
  
  # Store original state labels (before adding polymorphisms to them):
  original_state_labels <- colnames(x = costmatrix$costmatrix)
  
  # Permute new polymophic states:
  polymorphic_states <- permute_all_polymorphisms(single_states = costmatrix$single_states)
  
  # Calculate number of new states:
  n_new_states <- length(x = polymorphic_states)
  
  # Check data are not too big (>= 2^14 states) and stop and warn user if so:
  if (n_new_states >= 16384) stop("Costmatrix would be too large. Use fewer states.")

  # Store original costmatrix for later use:
  original_costmatrix <- costmatrix$costmatrix
  
  # Update costmatrix size with new states:
  costmatrix$size <- costmatrix$size + n_new_states
  
  # Initialise new costmatrix with infinities:
  costmatrix$costmatrix <- matrix(
    data = Inf,
    nrow = costmatrix$size,
    ncol = costmatrix$size,
    dimnames = list(
      c(original_state_labels, polymorphic_states),
      c(original_state_labels, polymorphic_states)
    )
  )

  # Populate portion representing original costmatrix with original values:
  costmatrix$costmatrix[original_state_labels, original_state_labels] <- original_costmatrix
  
  # Set diagonal to zero (always true):
  diag(x = costmatrix$costmatrix) <- 0
  
  # If using stratigraphic approach:
  if (polymorphism_costs == "stratigraphic") {
    
    # Store stratigraphic costs as transition to oldest states in polymorphism:
    costmatrix$costmatrix[costmatrix$single_states, polymorphic_states] <- do.call(
      what = cbind,
      args = lapply(
        X = as.list(polymorphic_states),
        FUN = function(polymorphism) {
          oldest_state <- strsplit(x = polymorphism, split = "&")[[1]][1]
          costmatrix$costmatrix[costmatrix$single_states, oldest_state]
        }
      )
    )
  }
  
  # If using geometric approach:
  if (polymorphism_costs == "geometric") {
    
    # Generate and store all possible states:
    all_states <- c(costmatrix$single_states, polymorphic_states)
    
    # Create coordinate matrix and initialise with zeroes:
    state_presence_matrix <- matrix(
      data = 0,
      nrow = length(x = costmatrix$single_states),
      ncol = length(x = all_states),
      dimnames = list(costmatrix$single_states, all_states)
    )
    
    # If using the hypercube coordinate space assign coordinates accordingly:
    if (polymorphism_geometry == "hypercube") {
      for(i in 1:ncol(x = state_presence_matrix)) {
        state_presence_matrix[strsplit(x = colnames(x = state_presence_matrix)[i], split = "&")[[1]], i] <- 1
      }
    }
    
    # If using the hypersphere or simplex coordinate space:
    if (polymorphism_geometry == "hypersphere" || polymorphism_geometry == "simplex") {
      
      # For each coding:
      for(i in 1:ncol(x = state_presence_matrix)) {
        
        # Isolate components of polymorphism:
        components <- strsplit(x = colnames(x = state_presence_matrix)[i], split = "&")[[1]]
        
        # If using the hypersphere coordinate space then apply coordinates accordingly (the square root of 1/N states in polymorphism on each axis state is present):
        if (polymorphism_geometry == "hypersphere") {
          state_presence_matrix[components, i] <- sqrt(x = 1 / length(x = components))
        }
        
        # If using the simplex coordinate space then apply coordinates accordingly (1/N states in polymorphism on each axis state is present):
        if (polymorphism_geometry == "simplex") {
          state_presence_matrix[components, i] <- 1 / length(x = components)
        }
      }
    }
    
    # If using a manhattan or euclidean distance, calculate distance directly from coordinate-space:
    if (polymorphism_distance == "euclidean" || polymorphism_distance == "manhattan") {
      new_costmatrix <- as.matrix(
        x = dist(
          x = t(x = state_presence_matrix),
          method = polymorphism_distance,
          diag = TRUE,
          upper = TRUE
        )
      )
    }
    
    # If using a great circle distance:
    if (polymorphism_distance == "great_circle") {
      
      # Start by calculating the euclidean distances between each point:
      new_costmatrix <- as.matrix(
        x = dist(
          x = t(x = state_presence_matrix),
          method = "euclidean",
          diag = TRUE,
          upper = TRUE
        )
      )
      
      # Treating distances as the chord length of a circle of radius one transform them to the arc length for the same:
      new_costmatrix <- 2 * asin(x = new_costmatrix / 2)
    }
    
    # Rescale costmatrix so that base cost is one:
    new_costmatrix <- new_costmatrix / min(x = new_costmatrix[lower.tri(x = new_costmatrix)])
    
    # Update character weight based on new single state to single state base cost:
    costmatrix$weight <- costmatrix$weight / new_costmatrix[costmatrix$single_states, costmatrix$single_states][lower.tri(x = new_costmatrix[costmatrix$single_states, costmatrix$single_states])][1]

    # Insert geometric costmatrix into costmatrix:
    costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states] <- new_costmatrix[costmatrix$single_states, costmatrix$single_states]
    costmatrix$costmatrix[costmatrix$single_states, polymorphic_states] <- new_costmatrix[costmatrix$single_states, polymorphic_states]
    costmatrix$costmatrix[polymorphic_states, costmatrix$single_states] <- new_costmatrix[polymorphic_states, costmatrix$single_states]
    costmatrix$costmatrix[polymorphic_states, polymorphic_states] <- new_costmatrix[polymorphic_states, polymorphic_states]
  }

  # If using Maddison approach:
  if (polymorphism_costs == "maddison") {
    
    # Subfunction to calculate a Maddison distance:
    calculate_maddison_distance <- function(costmatrix, state_1, state_2) {
      
      # Get vectors of end states:
      states_1 <- strsplit(x = state_1, split = "&")[[1]]
      states_2 <- strsplit(x = state_2, split = "&")[[1]]

      # Collect unique end states (to help identify those that can be pruned):
      end_states <- unique(x = c(states_1, states_2))
      
      # Initialise prunable states vector with just potential prunable states (if bridge vertices will have to be retained - checked below):
      potential_prunable_states <- setdiff(x = costmatrix$single_states, y = end_states)

      # If there are potential prunable states (pendant vertices):
      if (length(potential_prunable_states) > 0) {
        
        # Collect actually prunable states:
        prunable_states <- unlist(
          x = lapply(
            X = as.list(potential_prunable_states),
            FUN = function(pruned_state) {
              
              # Make pruned costmatrix by removing pruned state:
              pruned_costmatrix <- costmatrix$costmatrix[end_states, setdiff(x = costmatrix$single_states, y = pruned_state)]
              
              # Get new minimum costs with pruned state removed:
              min_pruned_cost <- apply(
                X = pruned_costmatrix,
                MARGIN = 1,
                FUN = function(row) min(x = row[row > 0])
              )
              
              # If state is prunable (all minimum costs remains 1):
              if (all(min_pruned_cost == 1)) {
                
                # Return prunable state:
                return(value = pruned_state)
                
              # If state is not prunable:
              } else {
                
                # Return empty vector:
                return(value = c())
              }
            }
          )
        )
        
        # If there are prunable state(s):
        if (length(x = prunable_states) > 0) {
          
          # Identify retained states:
          costmatrix$single_states <- setdiff(x = costmatrix$single_states, y = prunable_states)

          # Prune costmatrix down to retained states:
          costmatrix$costmatrix <- costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states]
        }
      }
      
      # Get number of single switches (states that only needed to be turned on/off once):
      single_switches <- length(
        x = setdiff(
          x = union(x = states_1, y = states_2),
          y = intersect(x = states_1, y = states_2)
        )
      )
      
      # Get number of double switches (i.e., states that must be turned off then on again to reach other states):
      double_switches <- length(x = setdiff(x = costmatrix$single_states, y = end_states))

      # Return Maddison distance:
      (double_switches * 2) + single_switches
    }

    # Generate and store all possible states:
    all_states <- c(costmatrix$single_states, polymorphic_states)
    
    # Calculate costmatrix size:
    costmatrix_size <- length(x = all_states)
    
    # Initialise maddison costmatrix with zeroes:
    new_costmatrix <- matrix(
      data = 0,
      nrow = costmatrix_size,
      ncol = costmatrix_size,
      dimnames = list(all_states, all_states)
    )
    
    # For each costmatrix row:
    for(i in 1:(costmatrix_size - 1)) {
      
      # For each costmatrix column:
      for(j in (i + 1):costmatrix_size) {
        
        # Calculate maddison distance and store:
        new_costmatrix[all_states[i], all_states[j]] <- new_costmatrix[all_states[j], all_states[i]] <- calculate_maddison_distance(
          costmatrix = costmatrix,
          state_1 = all_states[i],
          state_2 = all_states[j]
        )
      }
    }
    
    # Insert maddison costmatrix into costmatrix:
    costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states] <- new_costmatrix[costmatrix$single_states, costmatrix$single_states]
    costmatrix$costmatrix[costmatrix$single_states, polymorphic_states] <- new_costmatrix[costmatrix$single_states, polymorphic_states]
    costmatrix$costmatrix[polymorphic_states, costmatrix$single_states] <- new_costmatrix[polymorphic_states, costmatrix$single_states]
    costmatrix$costmatrix[polymorphic_states, polymorphic_states] <- new_costmatrix[polymorphic_states, polymorphic_states]
    
    # Halve character weight on account of new minima being 2 for a single state transition:
    costmatrix$weight <- costmatrix$weight / 2
  }
  
  # If uncertainties were included (will need to calculate costs for polymorphism to uncertainty transitions):
  if (costmatrix$includes_uncertainties) {
    
    # First establish which uncertainties are present:
    uncertainty_states <- colnames(costmatrix$costmatrix)[grep(pattern = "/", x = colnames(costmatrix$costmatrix))]
    
    # Inset cost for minimum polymorphic to single state transition for each uncertainty:
    costmatrix$costmatrix[polymorphic_states, uncertainty_states] <- do.call(
      what = cbind,
      args = lapply(
        X = as.list(x = uncertainty_states),
        FUN = function(uncertainty) {
          uncertainty_components <- strsplit(x = uncertainty, split = "/")[[1]]
          apply(
            X = costmatrix$costmatrix[polymorphic_states, uncertainty_components],
            MARGIN = 1,
            FUN = min
          )
        }
      )
    )
  }

  # If using additive approach (no need to update uncertainties for this):
  if (polymorphism_costs == "additive") {
    
    # Generate and store all possible states:
    all_states <- c(costmatrix$single_states, polymorphic_states)
    
    # Calculate costmatrix size:
    costmatrix_size <- length(x = all_states)
    
    # Initialise new costmatrix with infinite costs:
    new_costmatrix <- matrix(
      data = Inf,
      nrow = costmatrix_size,
      ncol = costmatrix_size,
      dimnames = list(all_states, all_states)
    )
    
    # Populate single states part of matrix by copying from starting costmatrix:
    new_costmatrix[costmatrix$single_states, costmatrix$single_states] <- costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states]

    # Set diagonal as zeroes (as must be):
    diag(x = new_costmatrix) <- 0
    
    # Calculate and store additive costmatrix:
    for(ancestor in costmatrix$single_states) {
      for(descendant in polymorphic_states) {
        new_costmatrix[ancestor, descendant] <- sum(
          x = costmatrix$costmatrix[ancestor, strsplit(x = descendant, split = "&")[[1]]]
        )
      }
    }
    
    # Insert additive costmatrix into costmatrix:
    costmatrix$costmatrix[costmatrix$single_states, costmatrix$single_states] <- new_costmatrix[costmatrix$single_states, costmatrix$single_states]
    costmatrix$costmatrix[costmatrix$single_states, polymorphic_states] <- new_costmatrix[costmatrix$single_states, polymorphic_states]
    costmatrix$costmatrix[polymorphic_states, costmatrix$single_states] <- new_costmatrix[polymorphic_states, costmatrix$single_states]
    costmatrix$costmatrix[polymorphic_states, polymorphic_states] <- new_costmatrix[polymorphic_states, polymorphic_states]
  }
  
  # Update includes_polymorphisms to TRUE:
  costmatrix$includes_polymorphisms <- TRUE
  
  # Return costmatrix with polymorphisms added:
  return(value = costmatrix)
}
