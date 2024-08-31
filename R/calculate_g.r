#' Calculate the maximum tree length, g, under parsimony
#'
#' @description
#'
#' Given a costmatrix and set of tip states returns the longest possible tree length under maximum parsimony.
#'
#' @param costmatrix An object of class \code{costMatrix}.
#' @param tip_states A character vector of tip states, with polymorphic states separated by \code{&}, uncertainties by \code{/}, missing values as \code{NA}, and inapplicables as empty strings \code{""}.
#' @param polymorphism_behaviour One of either "missing", "uncertainty", "polymorphism", or "random". See details.
#' @param uncertainty_behaviour One of either "missing", "uncertainty", "polymorphism", or "random". See details.
#'
#' @details
#'
#' The maximum cost a character could have on any tree under maximum parsimony, termed \emph{g}, depends on both the individual state-to-state transition costs (captured by a costmatrix) and the sampled states (i.e., the \code{tip_states} input). In practice this is the maximum parsimony length on the star tree. This length cannot be exceeded by any other tree (Hoyal Cuthill and Lloyd, in review). Note: this is standard practice in phylogenetics software and is also how both PAUP* (Swofford 2003) and TNT (Goloboff et al. 2008; Goloboff and Catalano 2016) calculate maximum cost.
#'
#' \bold{Special cases}
#'
#' A number of special cases apply to calculating \emph{g} and are discussed further below.
#'
#' \emph{Polymorphisms}
#'
#' Polymorphisms remain a complex problem in phylogenetics and here multiple options are provided to deal with them. These include: 1. \code{"missing"} - where they are simply replaced by a missing value (see below), 2. \code{"uncertainty"} - where they are treated as uncertainties instead (see below), 3. \code{"polymorphism"} - where they are treated as genuinely polymorphic, and 4. \code{"random"} - where one of the tip states is selected at random.
#'
#' Options 1, 2, and 4 can be seen as undercounting the true amount of evolution that has occurred. However, how to \emph{correctly} count this amount is unclear. If option 3 is chosen then polymorphic states must be present in \code{costmatrix} and users should refer to the \link{add_polymorphisms_to_costmatrix} function for details on available options.
#'
#' \emph{Uncertainties}
#'
#' Uncertainties are much simpler to deal with than polymorphisms and a means to incorporate them into length counts was laid out in Swofford and Maddison (1992). Indeed, popular software such as PAUP* (Swofford 2003) and TNT (Goloboff et al. 2008; Goloboff and Catalano 2016) simply treat polymorphisms as uncertainties perhaps because of this. There is still a concern of undercounting evolutionary change for uncertainties in the maximum parsimony context as the cheapest possible state will in effect be used everytime, whereas future study that removes uncertainty may reveal s higher cost state to be the true value. As such the same options are offered for uncertainties as polymorphisms, including to treat them as polymorphisms although this should probably only be done where they were miscoded in the first place.
#'
#' Again, if using uncertainties as uncertainties, these must be included in the costmatrix and this can be done by using the \link{add_uncertainties_to_costmatrix} function.
#'
#' \emph{Missing values}
#'
#' In practice missing values (\code{NA}) may exist amongst \code{tip_states}. These are permitted, and in practice are mathematically and practically equivalent to a statement that a tip could be any state present in the costmatrix (i.e., a special case of an uncertainty where no state can be ruled out). However, it should be considered in interpretation that \emph{g} will typically become smaller as the number of missing values increases.
#'
#' \emph{Inapplicable values}
#'
#' Inapplicable values (\code{""}) may also exist amongst \code{tip_states}. These are conceptually different to missing values as there is no possibility that they can ever be (re)coded. Currently these are treated exactly the same as missing values, but again the user should apply caution in interpreting \emph{g} in such cases as again it will be smaller than otherwise identical characters with fewer inapplicable values.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Jen Hoyal Cuthill \email{j.hoyal-cuthill@@essex.ac.uk}
#'
#' @references
#'
#' Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics/ \emph{Cladistics}, \bold{32}. 221-238
#'
#' Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. \emph{Cladistics}, \bold{24}, 774-786.
#'
#' Hoyal Cuthill, J. F. and Lloyd, G. T., in press. Measuring homoplasy I: comprehensive measures of maximum and minimum cost under parsimony across discrete cost matrix character types. \emph{Cladistics}, bold{}, .
#'
#' Swofford, D. L., 2003. \emph{PAUP*. Phylogenetic Analysis Using Parsimony (*and Other Methods). Version 4}. Sinauer Associates, Sunderland, Massachusetts.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. In R. L. Mayden (ed.), \emph{Systematics, Historical Ecology, and North American Freshwater Fishes}. Stanford University Press, Stanford. pp187-223.
#'
#' @return
#'
#' A single value indicating the maximum length, \emph{g}. Note: this is not modified by \code{costmatrix$weight}.
#'
#' @seealso
#'
#' \link{calculate_gmax}
#'
#' @examples
#'
#' # Create a Type I character costmatrix:
#' constant_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 0,
#'   character_type = "unordered"
#' )
#'
#' # Calculate g for the case of five state 0s:
#' calculate_g(
#'   costmatrix = constant_costmatrix,
#'   tip_states = c("0", "0", "0", "0", "0")
#' )
#'
#' # Create a Type II character costmatrix:
#' binary_symmetric_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state = 1,
#'   character_type = "unordered"
#' )
#'
#' # Calculate g for the case of two state 0s and three state 1s:
#' calculate_g(
#'   costmatrix = binary_symmetric_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1")
#' )
#'
#' # Create a Type III character costmatrix:
#' unordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 2,
#'   character_type = "unordered"
#' )
#'
#' # Calculate g for the case of two state 0s and three state 1s and two state 2s:
#' calculate_g(
#'   costmatrix = unordered_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1", "2", "2")
#' )
#'
#' # Create a Type IV character costmatrix:
#' linear_ordered_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 2,
#'   character_type = "ordered"
#' )
#'
#' # Calculate g for the case of two state 0s and three state 1s and two state 2s:
#' calculate_g(
#'   costmatrix = linear_ordered_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1", "2", "2")
#' )
#'
#' # Create a Type V character costmatrix:
#' nonlinear_ordered_costmatrix <- convert_adjacency_matrix_to_costmatrix(
#'   adjacency_matrix = matrix(
#'     data = c(
#'       0, 1, 0, 0,
#'       1, 0, 1, 1,
#'       0, 1, 0, 0,
#'       0, 1, 0, 0
#'     ),
#'     nrow = 4,
#'     dimnames = list(0:3, 0:3)
#'   )
#' )
#'
#' # Calculate g for the case of two state 0s, three state 1s, two state 2s and one state 3:
#' calculate_g(
#'   costmatrix = nonlinear_ordered_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1", "2", "2", "3")
#' )
#'
#' # Create a Type VI character costmatrix:
#' binary_irreversible_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 1,
#'   character_type = "irreversible"
#' )
#'
#' # Calculate g for the case of two state 0s and three state 1s:
#' calculate_g(
#'   costmatrix = binary_irreversible_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1")
#' )
#'
#' # Create a Type VII character costmatrix:
#' multistate_irreversible_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 2,
#'   character_type = "irreversible"
#' )
#'
#' # Calculate g for the case of two state 0s and three state 1s and two state 2s:
#' calculate_g(
#'   costmatrix = multistate_irreversible_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1", "2", "2")
#' )
#'
#' # Create a Type VIII character costmatrix:
#' binary_dollo_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 1,
#'   character_type = "dollo"
#' )
#'
#' # Calculate g for the case of two state 0s and three state 1s:
#' calculate_g(
#'   costmatrix = binary_dollo_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1")
#' )
#'
#' # Create a Type IX character costmatrix:
#' multistate_dollo_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 2,
#'   character_type = "dollo"
#' )
#'
#' # Calculate g for the case of two state 0s and three state 1s and two state 2s:
#' calculate_g(
#'   costmatrix = multistate_dollo_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1", "2", "2")
#' )
#'
#' # Create a Type X character costmatrix:
#' multistate_symmetric_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 5,
#'   character_type = "ordered"
#' )
#' multistate_symmetric_costmatrix$type <- "custom"
#' multistate_symmetric_costmatrix$costmatrix <- matrix(
#'   data = c(
#'     0, 1, 2, 3, 2, 3,
#'     1, 0, 3, 2, 1, 2,
#'     2, 3, 0, 3, 2, 1,
#'     3, 2, 3, 0, 1, 2,
#'     2, 1, 2, 1, 0, 1,
#'     3, 2, 1, 2, 1, 0
#'   ),
#'   nrow = multistate_symmetric_costmatrix$size,
#'   ncol = multistate_symmetric_costmatrix$size,
#'   byrow = TRUE,
#'   dimnames = list(
#'     multistate_symmetric_costmatrix$single_states,
#'     multistate_symmetric_costmatrix$single_states
#'   )
#' )
#'
#' # Calculate g for the case of two state 0s, three state 1s, two state 2s,
#' # one state 3, three state 4s and two state 5s:
#' calculate_g(
#'   costmatrix = multistate_symmetric_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1", "2", "2", "3", "4", "4", "4", "5", "5")
#' )
#'
#' # Create a Type XI character costmatrix:
#' binary_asymmetric_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 1,
#'   character_type = "ordered"
#' )
#' binary_asymmetric_costmatrix$type <- "custom"
#' binary_asymmetric_costmatrix$costmatrix <- matrix(
#'   data = c(
#'     0, 1,
#'     10, 0
#'   ),
#'   nrow = binary_asymmetric_costmatrix$size,
#'   ncol = binary_asymmetric_costmatrix$size,
#'   byrow = TRUE,
#'   dimnames = list(
#'     binary_asymmetric_costmatrix$single_states,
#'     binary_asymmetric_costmatrix$single_states
#'   )
#' )
#' binary_asymmetric_costmatrix$symmetry <- "Asymmetric"
#'
#' # Calculate g for the case of two state 0s and three state 1s:
#' calculate_g(
#'   costmatrix = binary_asymmetric_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1")
#' )
#'
#' # Create a Type XII character costmatrix:
#' multistate_asymmetric_costmatrix <- make_costmatrix(
#'   min_state = 0,
#'   max_state= 2,
#'   character_type = "ordered"
#' )
#' multistate_asymmetric_costmatrix$type <- "custom"
#' multistate_asymmetric_costmatrix$costmatrix <- matrix(
#'   data = c(
#'     0, 1, 1,
#'     1, 0, 1,
#'     10, 10, 0
#'   ),
#'   nrow = multistate_asymmetric_costmatrix$size,
#'   ncol = multistate_asymmetric_costmatrix$size,
#'   byrow = TRUE,
#'   dimnames = list(
#'     multistate_asymmetric_costmatrix$single_states,
#'     multistate_asymmetric_costmatrix$single_states
#'   )
#' )
#' multistate_asymmetric_costmatrix$symmetry <- "Asymmetric"
#'
#' # Calculate g for the case of two state 0s and three state 1s and two state 2s:
#' calculate_g(
#'   costmatrix = multistate_asymmetric_costmatrix,
#'   tip_states = c("0", "0", "1", "1", "1", "2", "2")
#' )
#'
#' @export calculate_g
calculate_g <- function(
  costmatrix,
  tip_states,
  polymorphism_behaviour = "polymorphism",
  uncertainty_behaviour = "uncertainty"
) {

  # ADD INPUT CHECKS - E.G., CHECK TIP STATE LABELS MATCH COSTMATRIX LABELS! (CURRENTLY UST DROPS TIP STATES NOT FOUND IN COSTMARIX)
  # CONSIDER SOME OTHER EDGE CASES LIKE NO TIPS (EXCEPT MISSING/INAPPLICABLE OR ONLY A SINGLE CODED STATE).
  # WHAT ABOUT POLYMORPHISMS WITH DOLLO?
  # ADD SOME POLYMORPHIC/UNCERTAINTY EXAMPLES TO THE EXAMPLES!
  # ADD POSSIBILITY TO FORCE ROOT STATE
  
  # Initialise character_is_dollo as FALSE:
  character_is_dollo <- FALSE
  
  # If character actually is Dollo:
  if (costmatrix$type == "dollo") {
    
    # Overwrite character_is_dollo with TRUE:
    character_is_dollo <- TRUE
    
    # Convert costmatrix into a linear ordered character to maintain correct costs below:
    new_costmatrix <- make_costmatrix(
      min_state = 0,
      max_state = costmatrix$n_states - 1,
      character_type = "ordered",
      include_polymorphisms = costmatrix$includes_polymorphisms,
      include_uncertainties = costmatrix$includes_uncertainties,
      polymorphism_costs = costmatrix$polymorphism_costs,
      polymorphism_geometry = costmatrix$polymorphism_geometry,
      polymorphism_distance = costmatrix$polymorphism_distance,
      dollo_penalty = costmatrix$dollo_penalty
    )
    
    # Update weight term:
    new_costmatrix$weight <- costmatrix$weight
    
    # Overwrite Dollo costmatrix with new costmatrix:
    costmatrix <- new_costmatrix
  }
  
  # Find positions of any polymorphsims or uncertainties (need to ascertain before making any changes):
  polymorphic_tips <- grep(pattern = "&", x = tip_states)
  uncertain_tips <- grep(pattern = "/", x = tip_states)
  
  # If polymorphisms are present amongst tips:
  if (length(x = polymorphic_tips) > 0) {
    
    # Replace polymorphisms with a missing state if requested:
    if (polymorphism_behaviour == "missing") tip_states[polymorphic_tips] <- NA
    
    # Replace polymorphisms with uncertainties if requested:
    if (polymorphism_behaviour == "uncertainty") tip_states[polymorphic_tips] <- gsub(pattern = "&", replacement = "/", x = tip_states[polymorphic_tips])
    
    # Randomly sample a single state of each polymorphism if requested:
    if (polymorphism_behaviour == "random") tip_states[polymorphic_tips] <- unlist(x = lapply(X = as.list(x = tip_states[polymorphic_tips]), FUN = function(x) sample(x = strsplit(x = x, split = "&")[[1]], size = 1)))
  }
  
  # If uncertainties are present amongst tips:
  if (length(x = uncertain_tips) > 0) {
    
    # Replace uncertainties with a missing state if requested:
    if (uncertainty_behaviour == "missing") tip_states[uncertain_tips] <- NA
    
    #  Replace uncertainties with polymorphisms if requested:
    if (uncertainty_behaviour == "polymorphism") tip_states[uncertain_tips] <- gsub(pattern = "/", replacement = "&", x = tip_states[uncertain_tips])
    
    # Randomly sample a single state of each uncertainty if requested:
    if (uncertainty_behaviour == "random") tip_states[uncertain_tips] <- unlist(x = lapply(X = as.list(x = tip_states[uncertain_tips]), FUN = function(x) sample(x = strsplit(x = x, split = "/")[[1]], size = 1)))
  }

  # Error check that polymorphisms are included in costmatrix if present amongst tips:
  if (length(x = grep(pattern = "&", x = tip_states)) > 0 && !costmatrix$includes_polymorphisms) stop("Polymorphisms not included in costmatrix.")

  # Error check that uncertaintes are included in costmatrix if present amongst tips:
  if (length(x = grep(pattern = "/", x = tip_states)) > 0 && !costmatrix$includes_uncertainties) stop("Uncertainties not included in costmatrix.")

  # Isolate unique sampled states:
  unique_sampled_states <- intersect(x = colnames(x = costmatrix$costmatrix), y = unique(x = tip_states))
  
  # Get frequencies that each state are sampled:
  sampled_state_frequencies <- unlist(
    x = lapply(
      X = as.list(x = unique_sampled_states),
      FUN = function(state) {
        x <- length(x = which(x = tip_states == state))
        names(x = x) <- state
        x
      }
    )
  )
  
  # Build these into a sampled state frequency matrix (to multiply against cost matrix):
  state_frequency_matrix <- costmatrix$costmatrix[, unique_sampled_states, drop = FALSE]
  state_frequency_matrix[] <- matrix(data = rep(x = sampled_state_frequencies, times = costmatrix$size), nrow = costmatrix$size, byrow = TRUE)
  
  # Get total cost for each possible ancestral state on the star tree:
  total_cost_for_each_ancestral_state <- apply(
    X = state_frequency_matrix * costmatrix$costmatrix[, unique_sampled_states],
    MARGIN = 1,
    FUN = function(x) sum(x = as.numeric(x = gsub(pattern = NaN, replacement = 0, x = x)))
  )
  
  # Special case of a Dollo character:
  if (character_is_dollo) {
    
    # Dollo root must be forced to be second most derived tip states (allows gains but not duplicate gains):
    dollo_root_state <- as.character(x = sort(x = as.numeric(x = tip_states), decreasing = TRUE)[2])
    
    # Restrict total_cost_for_each_ancestral_state to just the value with forced Dollo root:
    total_cost_for_each_ancestral_state <- total_cost_for_each_ancestral_state[dollo_root_state]
  }
  
  # Return calculated value of g:
  unname(obj = min(x = total_cost_for_each_ancestral_state))
}
