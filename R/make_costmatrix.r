#' Make a costmatrix for a given set of states
#'
#' @description
#'
#' Given a set of discrete states and a character type will make the approriate costmatrix for them.
#'
#' @param min_state The minimum character state (defaults to \code{0}).
#' @param max_state The maximum character state. Must be \code{1} or greater.
#' @param character_type The type of character desired. Must be one of: \code{"ordered"}, \code{"unordered"}, \code{"dollo"}, \code{"irreversible"}, or \code{"stratigraphy"}.
#' @param include_polymorphisms Logical indicating whether or not to include polymorphic state combinations (defaults to \code{FALSE}).
#' @param include_uncertainties Logical indicating whether or not to include uncertainty state combinations (defaults to \code{FALSE}).
#' @param polymorphism_costs Only used if \code{include_polymorphisms = TRUE}. See \link{add_polymorphisms_to_costmatrix}.
#' @param polymorphism_geometry Only used if \code{include_polymorphisms = TRUE}. See \link{add_polymorphisms_to_costmatrix}.
#' @param polymorphism_distance Only used if \code{include_polymorphisms = TRUE}. See \link{add_polymorphisms_to_costmatrix}.
#' @param state_ages A vector of ages assigned to each state. Only triggered if \code{character_type = "stratigraphy"}.
#' @param dollo_penalty The size of the cost penalty for the acquisition of a Dollo character (defaults to \code{999}). Only triggered if \code{character_type = "dollo"}. Note: this should always be a positive real value greater than one, and never infinity (\code{Inf}), as at least one acquisition is expected.
#' @param message Logical indicating whether (\code{TRUE}, the default) or not (\code{FALSE}) to provide messages about the process of the function.
#'
#' @details
#'
#' Costmatrices encode the parsimony cost (typically the number of evolutionary "steps") for each possible state-to-state transition. They can be used to estimate total lengths and make estimates (often also called reconstructions) for ancestral states at branching points on a given topology (tree). This function automates the generation of some common costmatrix types for use elsewhere in Claddis and hence is intended primarily as an internal function. However, as costmatrices are fundamental to various key analyses their operation is detailed extensively below.
#'
#' \bold{Costmatrix basics}
#'
#' Although their usage is rare as explicit statements (see Hooker 2014 for an example), almost any character type (e.g., ordered, unordered) can be expressed as a costmatrix. These are always square, with rows and columns representing the same set of states that will usually represent, but are not restricted to, the range of values observed in the data. Individual values in a costmatrix represent the cost of a particular transition, with the row value being the "from" state and the column value being the "to" state. By convention, the cost of going from any character to itself is zero and hence so is the diagonal of the matrix. An example costmatrix for a three state unordered character would thus look like this:
#'
#' \preformatted{    -------------
#'     | 0 | 1 | 2 |
#' -----------------
#' | 0 | 0 | 1 | 1 |
#' -----------------
#' | 1 | 1 | 0 | 1 |
#' -----------------
#' | 2 | 1 | 1 | 0 |
#' -----------------}
#'
#' Hence going from state 2 (third row) to state 1 (second column) has a cost of 1. Indeed, as this is an unordered character, \emph{any} off-diagonal value has the same cost of 1. By contrast an ordered matrix might look like this:
#'
#' \preformatted{    -------------
#'     | 0 | 1 | 2 |
#' -----------------
#' | 0 | 0 | 1 | 2 |
#' -----------------
#' | 1 | 1 | 0 | 1 |
#' -----------------
#' | 2 | 2 | 1 | 0 |
#' -----------------}
#'
#' Now going from state 0 (first row) to state 2 (third column) costs two, as by implication this transition must pass through the intermediate state (1).
#'
#' So far these examples are symmetric - i.e., you can imagine the diagonal as a line of reflection, or, alternatively, that going from state X to state Y \emph{always} costs the same as going from state Y to state X. However, asymmetric costmatrices are also possible and there are multiple such character-types that may be relevant here. (Specifically, Dollo, irreversible, also known as Camin-Sokal, and stratigraphic.)
#'
#' \bold{Character types}
#'
#' This function will generate a costmatrix for every possible transition (between \code{min_state} and \code{max_state}, inclusive) for a requested \code{character_type}, each of which is detailed further below.
#'
#' \emph{Ordered} - \code{character_type = "ordered"}
#'
#' An ordered character (really a linear ordered character), also known as a Wagner character (Wagner 1961) and formalised by Farris (1970), is one where the order of states matters. State-to-state transitions must occur through a linear series. For example, for a three state character the states must be in order of transition, 0 <-> 1 <-> 2. Additionally, all transitions are symmetric. 0 -> 1 = 1 -> 0.
#'
#' \emph{Unordered} - \code{character_type = "unordered"}
#'
#' An unordered character, also known as a Fitch character (Fitch 1971), is one where the order of states does not matter, and that state-to-state transitions are all direct and symmetric. For example, for a five-state character the transition 0 -> 4 is direct and costs the same as 3 -> 1, or any other transition, excepting those from any state to itself.
#'
#' \emph{Dollo} - \code{character_type = "dollo"}
#'
#' Under a Dollo assumption the acquisition of a \emph{derived} character state is considered to be sufficiently complex (biologically difficult) that in all probability it only occurred once. Even if lost later in evolution, it is considered that it cannot be reacquired. An example of this is the frequent loss of teeth in dinosaurs (Brocklehurst and Field 2021). Teeth were never reacquired in this group, with tooth-like serrations forming in many birds instead. Modelling such characters with a costmatrix is challenging (indeed, it cannot truly been done with \emph{just} a costmatrix - see below). The way it is typically done (Swofford and Olsen 1990) is to form the costmatrix like this:
#'
#' \preformatted{    -------------------
#'     |  0  |  1  |  2  |
#' -----------------------
#' | 0 |  0  |  1D |  2D |
#' -----------------------
#' | 1 |  1  |  0  |  1D |
#' -----------------------
#' | 2 |  2  |  1  |  0  |
#' -----------------------}
#'
#' Where \code{D} is some arbitrary large number, set here using the \code{dollo_penalty} option. Note that \code{D} should still be finite as the purpose is to weight acquisitions such that they are sufficiently expensive that any most parsimonious reconstruction will favour only one transition, but not make them so expensive that they are not favoured at all. Furthermore, in this example there are two derived states (1 and 2) and hence logically the weighting should be similar to an ordered character. Specifically, that there should be a single acquisition of state 1 (from state 0) and that this should precede a single acquistion of state 2 (from state 1). Most of the time, though, a Dollo character will be binary (e.g., see \link{map_dollo_changes}).
#'
#' Importantly, and as stated above, a Dollo costmatrix, unlike the ordered and unordered costmatrices, cannot be used without further analytical restrictions. Specifically, because it assumes an asymmetric acquisition of \emph{derived} states the root or "primitive" value must be 0. (It would not be logical to apply such a costmatrix where the root state is 1 or 2.) Thus use of a Dollo character requires the additional assumption that the primitive (root) state is known \emph{a priori}. As always, it is up to the user to know that this is a valid assumption for their data, and to set a value for the \code{dollo_penalty} accordingly.
#'
#' Note, that another way to conceive of a Dollo character is that the assumption is being applied that any homoplasy is in the form of reversals only and hence can be though of as making the opposite assumption to an irreversible character.
#'
#' \emph{Irreversible (Camin-Sokal)} - \code{character_type = "irreversible"}
#'
#' Sometimes confused with a Dollo character, irreversible or Camin-Sokal (Camin and Sokal 1965) characters, do not allow for any reversals (returns to some prior state) once the derived state has been acquired. An irreversible costmatrix might look like this:
#'
#' \preformatted{    -------------------
#'     |  0  |  1  |  2  |
#' -----------------------
#' | 0 |  0  |  1  |  2  |
#' -----------------------
#' | 1 | Inf |  0  |  1  |
#' -----------------------
#' | 2 | Inf | Inf |  0  |
#' -----------------------}
#'
#' The upper triangle (representing gains) is the same as for an ordered character, but the lower triangle (representing losses, or reversals) will always be comprised of infinite costs (\code{Inf}), precluding their possibility in a most parsimonious reconstruction. Thus, in practice, any number of gains is allowed, but losses never are.
#'
#' Like a Dollo character, this approach assumes the direction of evolution is known \emph{a priori} and that 0 is always the primitive, or root, state. Although here, and unlike a Dollo character, the root state is forced by the costmatrix alone.
#'
#' \emph{Stratigraphic} - \code{character_type = "stratigraphy"}
#'
#' Stratigraphic characters, where states represent geologic time units in which an OTU is present, are also logically irreversible - as a younger unit cannot exist before (be "ancestral" to) an older one. Thus here, too, the lower triangle of the costmatrix is always populated by infinities to preclude their possibility. However, this time the upper triangle requires additional information to populate its' values, i.e., the gap between states in units of time (typically millions of years).
#'
#' This is best illustrated with an example. Let's assume we have three states and they represent three consecutive geologic stages: 0 (Santonian, 85.8-83.5 Ma), 1 (Campanian, 83.5-70.6 Ma), and 2 (Maastrichtian, 70.6-65.5 Ma). Now the cost of a transition between states should be expressed in some difference between the ages of each state. In this example we will use the midpoints as our point estimate: 0 (84.65 Ma), 1 (77.05 Ma), and 2 (68.05 Ma). These must also be supplied to the function using the \code{state_ages} option, i.e.: \code{state_ages = c(84.65, 77.05, 68.05)}. The resulting costmatrix would look like this:
#'
#' \preformatted{    --------------------
#'     |  0  |  1  |   2  |
#' ------------------------
#' | 0 |  0  | 7.6 | 16.6 |
#' ------------------------
#' | 1 | Inf |  0  |   9  |
#' ------------------------
#' | 2 | Inf | Inf |   0  |
#' ------------------------}
#'
#' Note that the user needn't use the midpoint (a start or end date may be preferable). Indeeed, a simple unit count could be applied instead but this would then be identical to an irreversible character (except, perhaps, in the case of polymorphisms - see \link{add_polymorphisms_to_costmatrix}).
#'
#' The use of stratigraphy as a character is probably not generally recommended, but is offered here for those that wish to explore stratocladistic methods that aim to incorporate stratigraphic information in phylogenetic inference and identify ancestor-descendant relationships (see Fisher 1994).
#'
#' As with the other asymmetric characters it is an additional assumption that the root represent 0 (the oldest sampled value), but as with irreverisble characters, the costmatrix alone forces this assumption.
#'
#' \emph{Other character types}
#'
#' Although no other character types are offered here, other forms of costmatrix can of course exist (e.g., see those in Hooker 2014). These can still be used elsewhere in Claddis, but users must either generate them themselves or have them be specified within a NEXUS file imported into Claddis. Note that general advice is to always have the diagonal be zero and the off-diagonal values be positive, otherwise downstream use will typically be confounded.
#'
#' \bold{Character weights}
#'
#' As discusssed in Hoyal Cuthill and Lloyd (in prep), the costs in a costmatrix are more properly considered as ratios. Consider the following two costmatrices:
#'
#' \preformatted{    -------------
#'     | 0 | 1 | 2 |
#' -----------------
#' | 0 | 0 | 1 | 2 |
#' -----------------
#' | 1 | 1 | 0 | 1 |
#' -----------------
#' | 2 | 2 | 1 | 0 |
#' -----------------}
#'
#' \preformatted{    -------------
#'     | 0 | 1 | 2 |
#' -----------------
#' | 0 | 0 | 2 | 4 |
#' -----------------
#' | 1 | 2 | 0 | 2 |
#' -----------------
#' | 2 | 4 | 2 | 0 |
#' -----------------}
#'
#' Although the specific costs differ the ratio \emph{between} costs within each costmatrix is identical. I.e., the cost of a 0 to 1 transition is half the cost of a 0 to 2 transition in both costmatrices. Or, to put it another way, the second costmatrix can be derived from the first by multiplying every cost by the same "weight", here two. As such, these two costmatrices can be considered identical for analytical purposes, save for this weight term. So, for example, the ancestral state estimates would be identical for a given tree and set of tip values. And the tree length (including the minimum and maximum length) would be the same save for the weight modifier.
#'
#' This observation is largely irrelevant to this function. However, it is noted here so the user understands why the costmatrx structure (outlined below) includes a weight term. Further, it is important to state why this information should be encoded "outside" the costmatrix rather than encodng "inside" the costmatrix by simply multiplying every cost by the desired character weight (e.g., as in the 2-4-2 example above). The reason is that this can lead to confounded results in some cases. For example, if weighting a character zero to exclude it from an analysis would result in a costmatrix where all off diagonal costs were zero (the same cost as no change). In practice this would mean any ancestral state estimation would become equally likely. If encoding the weight "outside" the costmatrix this situation can be handled correctly, both within Claddis and if exporting data for use elsewhere.
#'
#' \bold{Polymorphisms and uncertainties}
#'
#' Where a taxon (terminal) has two or more states (or a single state cannot be identified) it represents either a polymorphism or an uncertainty. For some kinds of analyses, such as calculating tree lengths or estimating ancesral states, there is no need to explicitly address uncertainties in a costmatrix setting (Swofford and Maddison 1992). However, the same cannot be said of polymorphisms (Maddison and Maddison 1987). In any case, this function offers the chance to include either or both polymorphic or uncertain states into the costmatrix by calling \link{add_polymorphisms_to_costmatrix} or \link{add_uncertainties_to_costmatrix}, as appropriate. Interested users should consult the helpfiles for those functions for more details.
#'
#' \bold{Relation to Q-matrices}
#'
#' Q-matrices (e.g., see Swofford and Olsen 1990) are used in likelihood and Bayesian approaches and also deal with state-to-state transitions, but are fundamentally distinct from costmatrices in multiple ways. For example, Q-matrices encode transition rates, not costs, and their values are typically parameters that are estimated from the data not prior statements about the data. However, both can be rendered as state graph or Markov chain diagrams. Put simply though, they are not directly equivalent or interchangable.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Brocklehurst, N. and Field, D. J., 2021. Macroevolutionary dynamics of dentition in Mesozoic birds reveal no long-term selection towards tooth loss. \emph{iScience}, \bold{24}, 102243.
#'
#' Camin, J. H. and Sokal, R. R., 1965. A method for deducing branching sequences in phylogeny. \emph{Evolution}, \bold{19}, 311-26.
#'
#' Farris, J. S., 1970. Methods for computing Wagner trees. \emph{Systematic Zoology}, \bold{19}, 83-92.
#'
#' Fisher, D. C., 1994. Stratocladistics: morphological and temporal patterns and their relation to phylogenetic process. In L. Grande and O. Rieppel (eds.), \emph{Interpreting the Hierarchy of Nature}. Academic Press, San Diego. pp133–171.
#'
#' Fitch, W. M., 1971. Towards defining the course of evolution: minimum change for a specified tree topology. \emph{Systematic Zoology}, \bold{20}, 406-416.
#'
#' Hooker, J. J., 2014. New postcranial bones of the extinct mammalian family Nyctitheriidae (Paleogene, UK): primitive euarchontans with scansorial locomotion. \emph{Palaeontologia Electronica}, \bold{17.3.47A}, 1-82.
#'
#' Maddison, W. P. and Maddison, D. R., 1987. MacClade 2.1, computer program and manual. Cambridge, Massachusetts.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. In R. L. Mayden (ed.), \emph{Systematics, Historical Ecology, and North American Freshwater Fishes}. Stanford University Press, Stanford. pp187-223.
#'
#' Swofford, D. L. and Olsen, G. J., 1990. Phylogeny reconstruction. In D. M. Hillis and C. Moritz (eds.), \emph{Molecular Systematics}. Sinauer Associates, Sunderland. pp411-501.
#'
#' Wagner, W. H., 1961. Problems in the classification of ferns. \emph{Recent Advances in Botany}, \bold{1}, 841-844.
#'
#' @return
#'
#' An object of class \code{costMatrix} containing the following elements:
#'
#' \itemize{
#'   \item{\code{size} The number of columns (or rows) in the costmatrix.}
#'   \item{\code{n_states} The number of (single) states in the costmatrix.}
#'   \item{\code{single_states} The labels of the (single) states in the costmatrix.}
#'   \item{\code{type} The type of costmatrix. One of: \code{"ordered"}, \code{"unordered"}, \code{"dollo"}, \code{"irreversible"}, \code{"stratigraphy"}, or \code{"custom"}.}
#'   \item{\code{costmatrix} The costmatrix itself.}
#'   \item{\code{symmetry} The symmetry of the costmatrix. NB: This refers to the single states part only.}
#'   \item{\code{includes_polymorphisms} A logical indicating whether or not polymorphic states are included in the costmatrix.}
#'   \item{\code{polymorphism_costs} The means by which costs were assigned to polymorphisms. Must be one of: \code{"additive"}, \code{"geometric"}, \code{"maddison"}, or \code{"stratigraphic"}.}
#'   \item{\code{polymorphism_geometry} The geometry by which costs were assigned to polymorphisms. Must be one of: \code{"hypercube"}, \code{"hypersphere"}, or \code{"simplex"}.}
#'   \item{\code{polymorphism_distance} The distance metric by which costs were assigned to polymorphisms. Must be one of: \code{"euclidean"}, \code{"great_circle"}, or \code{"manhattan"}.}
#'   \item{\code{includes_uncertainties} A logical indicating whether or not uncertain states are included in the costmatrix.}
#'   \item{\code{pruned} A logical indicating whether or not the costmatrix represents a pruned version.}
#'   \item{\code{dollo_penalty} A numeric value indicating the penalty used for a Dollo character.}
#'   \item{\code{base_age} A numeric value indicating the base (oldest) age used for a stratigraphic character.}
#'   \item{\code{weight} A numeric value indicating the character weight.}
#' }
#'
#' @seealso
#'
#' \link{permute_all_polymorphisms}
#'
#' @examples
#'
#' # Make an unordered three-state costmatrix:
#' make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "unordered"
#' )
#'
#' # Make an ordered three-state costmatrix:
#' make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "ordered"
#' )
#'
#' # Make a three-state Dollo costmatrix:
#' make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "dollo",
#'   dollo_penalty = 100
#' )
#'
#' # Make a three-state irreversible costmatrix:
#' make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "irreversible",
#' )
#'
#' # Make a three-state stratigraphic costmatrix:
#' make_costmatrix(
#'   min_state = 0,
#'   max_state = 2,
#'   character_type = "stratigraphy",
#'   state_ages = c(52, 34, 12)
#' )
#'
#' @export make_costmatrix
make_costmatrix <- function(
  min_state = 0,
  max_state,
  character_type,
  include_polymorphisms = FALSE,
  include_uncertainties = FALSE,
  polymorphism_costs = "additive",
  polymorphism_geometry = "simplex",
  polymorphism_distance = "euclidean",
  state_ages,
  dollo_penalty = 999,
  message = TRUE
) {
  
  # OPERATE ON A VECTOR OF UNCERTAINTIES/POLYMORPHISMS INSTEAD OF GENERATING ALL OF THEM? (MAYBE ALL IS BETTER FOR ANCESTRAL STATE RECONSTRUCTION?
  # - Need more incoming checks probably?
  # - Option to disallow polymorphisms as ancestors (still allows costs for transitions to polymorphic states to be taken into account)
  
  # Check character_type is a valid value and stop and warn user if not:
  if (length(x = setdiff(x = character_type, y = c("ordered", "unordered", "dollo", "irreversible", "stratigraphy"))) > 0) stop("character_type must be one of \"ordered\", \"unordered\", \"dollo\", \"irreversible\", or \"stratigraphy\".")
  
  # Get single states:
  single_states <- as.character(x = min_state:max_state)
  
  # Case if character is ordered:
  if (character_type == "ordered") {
    
    # Store symmetry of costmatrix as symmetric:
    symmetry <- "Symmetric"
    
    # Calculate costmatrix size:
    costmatrix_size <- length(x = single_states)
     
    # Initialise costmatrix with all values set to zero:
    costmatrix <- matrix(data = 0, nrow = costmatrix_size, ncol = costmatrix_size, dimnames = list(single_states, single_states))
      
    # Store ordered values in upper and lower triangles:
    costmatrix[upper.tri(x = costmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(costmatrix_size - 1)), FUN = function(cost) cost:1))
    costmatrix[lower.tri(x = costmatrix)] <- unlist(x = lapply(X = as.list(x = (costmatrix_size - 1):1), FUN = function(cost) 1:cost))
  }

  # Case if character is unordered:
  if (character_type == "unordered") {
    
    # Store symmetry of costmatrix as symmetric:
    symmetry <- "Symmetric"

    # Initialise costmatrix with all values set to one:
    costmatrix <- matrix(
      data = 1,
      nrow = length(x = single_states),
      ncol = length(x = single_states),
      dimnames = list(single_states, single_states)
    )
      
    # Make diagional of costmatrix zero:
    diag(x = costmatrix) <- 0
  }
  
  # Case if character is Dollo:
  if (character_type == "dollo") {
    
    # Store symmetry of costmatrix as asymmetric:
    symmetry <- "Asymmetric"

    # Calculate costmatrix size:
    costmatrix_size <- length(x = single_states)
      
    # Initialise costmatrix with all values set to zero:
    costmatrix <- matrix(
      data = 0,
      nrow = costmatrix_size,
      ncol = costmatrix_size,
      dimnames = list(single_states, single_states)
    )
      
    # Set upper triangle as an ordered character times the dollo_penalty
    costmatrix[upper.tri(x = costmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(costmatrix_size - 1)), FUN = function(cost) cost:1)) * dollo_penalty
      
    # Store regular ordered values in lower triangles:
    costmatrix[lower.tri(x = costmatrix)] <- unlist(x = lapply(X = as.list(x = (costmatrix_size - 1):1), FUN = function(cost) 1:cost))
  }
  
  # Case if character is irreversible:
  if (character_type == "irreversible") {
    
    # Store symmetry of costmatrix as asymmetric:
    symmetry <- "Asymmetric"

    # Calculate costmatrix size:
    costmatrix_size <- length(x = single_states)
      
    # Initialise costmatrix with all values set to zero:
    costmatrix <- matrix(data = 0, nrow = costmatrix_size, ncol = costmatrix_size, dimnames = list(single_states, single_states))
      
    # Exclude losses by setting lower triangle values to infinity:
    costmatrix[lower.tri(x = costmatrix)] <- Inf
      
    # Set upper triangle as an ordered character:
    costmatrix[upper.tri(x = costmatrix)] <- unlist(x = lapply(X = as.list(x = 1:(costmatrix_size - 1)), FUN = function(cost) cost:1))
  }
  
  # Set base age arbitrarily at 1:
  base_age <- 1
  
  # Case if character is stratigraphy:
  if (character_type == "stratigraphy") {
    
    # Store symmetry of costmatrix as asymmetric:
    symmetry <- "Asymmetric"

    # Add state names to ages:
    names(x = state_ages) <- min_state:max_state
    
    # Generate initial costmatrix using temporal distances:
    costmatrix <- as.matrix(x = dist(x = state_ages, diag = TRUE, upper = TRUE))
    
    # Set infinite cost to all reversals:
    costmatrix[lower.tri(x = costmatrix)] <- Inf
    
    # Replace base age with actual maximum stratigraphic age:
    base_age = max(x = state_ages)
  }
  
  # Create full costmatrix object:
  costmatrix <- list(
    size = ncol(x = costmatrix),
    n_states = ncol(x = costmatrix),
    single_states = single_states,
    type = character_type,
    costmatrix = costmatrix,
    symmetry = symmetry,
    includes_polymorphisms = FALSE,
    polymorphism_costs = polymorphism_costs,
    polymorphism_geometry = polymorphism_geometry,
    polymorphism_distance = polymorphism_distance,
    includes_uncertainties = FALSE,
    pruned = FALSE,
    dollo_penalty = dollo_penalty,
    base_age = base_age,
    weight = 1
  )
  
  # Set class of output as costMatrix:
  class(costmatrix) <- "costMatrix"

  # If polymorphisms were requested:
  if (include_polymorphisms) {
  
    # Add polymorphisms to costmatrix using options supplied:
    costmatrix <- add_polymorphisms_to_costmatrix(
      costmatrix = costmatrix,
      polymorphism_costs = polymorphism_costs,
      polymorphism_geometry = polymorphism_geometry,
      polymorphism_distance = polymorphism_distance,
      message = message
    )
  }
  
  # If uncertainties were requested:
  if (include_uncertainties) {
    
    # Add uncertainties to costmatrix using options supplied
    costmatrix <- add_uncertainties_to_costmatrix(
      costmatrix = costmatrix,
      message = message
    )
  }
  
  # Return costmatrix:
  costmatrix
}


