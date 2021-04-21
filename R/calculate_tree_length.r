#' Calculates the parsimony length of a set of phylogenetic tree(s)
#'
#' @description
#'
#' Given a tree, or set of trees, and a cladistic matrix returns their parsimony length in number of steps.
#'
#' @param trees A tree (phylo object) or set of trees (multiPhylo object).
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}. These should be discrete with labels matching the tip labels of \code{tree}.
#' @param inapplicables_as_missing Logical that decides whether or not to treat inapplicables as missing (TRUE) or not (FALSE, the default and recommended option).
#' @param polymorphism_behaviour One of either "missing", "uncertainty", "polymorphism", or "random". See details.
#' @param uncertainty_behaviour One of either "missing", "uncertainty", "polymorphism", or "random". See details.
#' @param polymorphism_shape Argument passed to \link{make_stepmatrix}.
#' @param polymorphism_distance Argument passed to \link{make_stepmatrix}.
#' @param state_ages Argument passed to \link{make_stepmatrix}.
#' @param dollo_penalty Argument passed to \link{make_stepmatrix}.
#'
#' @details
#'
#' Under the maximum parsimony criterion, a phylogenetic hypothesis is considered optimal if it requires the fewest number of evolutionary "steps". In order to evalulate this criterion we must therefore be able to calculate a tree's "length" (total steps assuming the fewest possible steps for every character used). Given a set of phylogenetic hypothes(es) and a cladistic matrix this function calculates the minimum length for each tree.
#'
#' \bold{Input data format}
#'
#' This function operates on a phylogenetic tree, or trees (in \code{ape} format), and a cladistic matrix (in \link{cladisticMatrix} format [UPDATE]). However, the algorithm used is based on the generalised stepmatrix approach of Swofford and Maddison (1992) and hence stepmatrices need to be defined for each character (this is done internally by calling \link{make_stepmatrix}), and some of the options are merely passed to this function.
#'
#' \bold{Algorithm}
#'
#' Technically the Swofford and Maddison (1992) algorithm is designed for ancestral state reconstruction, but as its' first pass of the tree assigns lengths for each possible state at each node the minimum value of these options at the root is also the tree length for that character and hence by skipping the later steps this can be used as a tree length algorithm by simply adding these values for each character. The choice of the Swofford and Maddison algorithm, rather than the Wagner or Fitch algorithms (for ordered and unordered characters, respectively) is to generalize to the broadest range of character types, including asymmetric characters (Camin-Sokal, Dollo, stratigraphic), custom character types (specified using stepmatrices or character state trees), as well as to any resolution of tree (i.e., including multifurcating trees - important for establishing maximum step counts for homoplasy indices). The only restriction here is that the tree must be rooted such that time's arrow is explicitly present. This is essential, as the root defines the lengths across the whole tree, but also for asymmetric characters directionality must be explicit. The two obvious downsides to this algorithm is that it can be slower and that it is not appropriate for unrooted trees.
#'
#' \bold{Stepmatrices and stepmatrix options}
#'
#' Stepmatrices are described in detail in the \link{make_stepmatrix} manual, as are the options that are passed from this function to that one.




# Polymorphism options ("missing", "uncertainty", "polymorphism", or "random")
# - Treat as uncertainty (undercounts changes, may affect ASR)
# - Treat as missing (undercounts changes, may affect ASR)
# - Treat as states (requires stepmatrix to cover possibilities)




#' Can be used to infer ancestral states or get tree lengths (could theoretically build a tree inference procedure around this by proposing different topologies, but likely to be much slower than compiled software). Also a first step to generating character maps that can be used for plotting changes on a tree, performing rate analysis (see \link{test_rates}), or generating phylomorphospaces (not implemented yet).
#'
#' \emph{Explicit options}
#'
#' - inapplicables_as_missing = FALSE, what to do with inapplicable characters. Currently treat as missing is only real option, but can be differet state that is also inferred as an ancestor.
#' - More for uncertianties and polymorphisms (see below).
#' - Add option to predict missing tip values (but definitely not recommended! - see Lloyd 2018)
#' - Collpasing to a single output with ACCTRAN etc.?
#'
#' \bold{Implementation}
#'
#' Algorithms for reconstructing ancestral states under parsimony can be specific (e.g., for ordered characters only) or general (allowing for any character type, including highly specific custom-formats). The advantage of the former is usually that they are faster as multiple assumptions are already built in, but the latter offer more flexibility and represent the version applied here.
#'
#' More specifically this function is based on the approach described by Swofford and Maddison (1992) that assumes a rooted tree and uses a stepmatrix to define the costs (in number of steps) of specific transitions, i.e., from state X to state Y. The rooting is key here as it sets the directionality of evolution (i.e., time) and hence allows the use of asymmetric characters (where the cost of going from state X to state Y is different than going from stete Y to state X) as well as other direction-based optimisations (e.g., ACCTRAN and DELTRAN).
#'
#'
#' \bold{Special character types}
#' \emph{Missing values}
#' \emph{Inapplicables}
#' \emph{Uncertainties}
#'
#' In Claddis, uncertainties are coded by separating states with the slash character, e.g., 0/1. The Swofford and Maddison (1992) algorithm - and the implementation used here - treats uncertainties (automatically) by allowing any of the possible states at the tip (here 0 or 1, but not, say, 2 or 3). In other words, it will opt for solutions that effectively prefer whichever state is easiest (fewest steps) to reach from a proposed ancestral value. This could therefore be used to predict what the "true" value may be, at least under the assumption of parsimony. However, such phylogenetic "prediction" [OPTION?] is generally advised against (Lloyd 2018).
#'
#' \emph{Polymorphisms}
#'
#' TEXT?
#'
#' \bold{Polytomies}
#'
#' Another advantage of the Swofford and Maddison (1992) approach is that it does not require fully-bifurcating topologies, hence these are "allowed" as input. However, it is important to note that the implementation here will treat any such polytomies (also known as multifurcations) as "hard". I.e., that they genuinely represent one species splitting into three or more descendant species. As always, the Ian Malcolm caveat applies - "your scientists were so preoccupied with whether or not they sould, they didn’t stop to think if they should."


#'
#'
#' # Option to fix root state (polarize) somehow
#' # From is row, to is column (directional hence requires/assumes rooted trees)
#' # Record in manual that weighted parsimony is folly, but could use to chose between MPRs (maybe?)

#' \bold{Polymorphic characters}
#'
#' Random is same as PERDA approach of Watanabe [ADD REF]



#' # MANUAL: generalized (stepmatrix) solution as accounts for more possibilities
#' # MANUAL: really only for rooted trees (unrooted possible using same methods (see Swofford and Maddison (1992), but unclear why you would want to do this with trees without direction of time (i.e., a root).
#' # MANUAL: works with multifurcations and doesn't require a fully bifurcating tree, but does assume polytomies are therefore "hard".
#' # MANUAL: Swofford and Maddison (1992; p210): "Both [ACCTRAN and DELTRAN] start by either assuming that the ancestral state is known or by choosing a state from the set of optimal assignments to the root node" - i.e., always an arbitrainess problem at root!
#' MISSING DATA (NA) IS DEALT WITH BY SETTING ALL TIP VALUES FOR A MISSING STATE TO ZERO (AS PER SWOFFORD AND MADDISON 1992).
#' UNCERTAIN STATES (E.G., 0/1) CONSTRAIN ANCESTRAL STATES SLIGHTLY MORE THAN NA (ASSUMING THEY EXCLUDE SOME STATES. IMPLEMENTATION SAME AS SWOFFORD AND MADDISON (1992).
#' CALCULATES TREE LENGTHS!
#'
#'
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Lloyd, G. T., 2018. Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. \emph{Palaeontology}, \bold{61}, 637-645.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. \emph{In} R. L. Mayden (ed.) Systematics, Historical Ecology, and North American Freshwater Fishes. Stanford University Press, Stanford, p187-223.
#'
#' @return
#'
#' A list with multiple components, including:
#'
#' \item{input_trees}{The tree(s) used as input.}
#' \item{input_matrix}{The raw (unmodified) \code{cladistic_matrix} input.}
#' \item{input_options}{The various input options used. Output for use by downstream functions, such as ancestral state estimation and stochastic character mapping.}
#' \item{stepmatrices}{The stepmatrices (one for each character) used. These are typically generated automatically by the funcion, but are output here for later use in ancestral state estimation and stochastic character mapping functions.}
#' \item{character_matrix}{The single character matrix object used. Essentially the \code{input_matrix} modified by the \link{input_options}.}
#' \item{character_lengths}{A matrix of characters (rows) and trees (columns) with values indicating the number of steps. The column sums of this matrix are the \code{tree_lengths} values. This output can also be used for homoplasy metrics.}
#' \item{character_weights}{A vector of the character weights used.}
#' \item{tree_lengths}{The primary output - the length for each input tree in total steps.}
#' \item{node_values}{The values (lengths for each state) for each node acrss trees and characters. This is used by \link{reconstruct_ancestral_states} for ancestral state reconstruction.}
#'
#' @seealso
#'
#' \link{make_stepmatrix}
#'
#' @examples
#'
#' # Calculate maximum parsimony tree length for Gauthier 1986:
#' cladistic_matrix <- Claddis::gauthier_1986
#' tree <- ape::read.tree(
#'   text = "(Outgroup,(Ornithischia,(Sauropodomorpha,(Ceratosauria,Procompsognathus,
#'   Liliensternus,(Carnosauria,(Ornithmimidae,Saurornitholestes,Hulsanpes,(Coelurus,
#'   Elmisauridae,(Compsognathus,(Ornitholestes,Microvenator,Caenagnathidae,
#'   (Deinonychosauria,Avialae))))))))));"
#' )
#' calculate_tree_length(
#'   tree = tree,
#'   cladistic_matrix = cladistic_matrix,
#'   inapplicables_as_missing = FALSE,
#'   polymorphism_behaviour = "uncertainty",
#'   uncertainty_behaviour = "uncertainty",
#'   polymorphism_shape = "hypersphere",
#'   polymorphism_distance = "great_circle",
#'   state_ages = c(),
#'   dollo_penalty = 999
#' )
#'
#' @export calculate_tree_length
calculate_tree_length <- function(trees, cladistic_matrix, inapplicables_as_missing = FALSE, polymorphism_behaviour, uncertainty_behaviour, polymorphism_shape, polymorphism_distance, state_ages, dollo_penalty) {
  
  # ADD CHECKS!
  # - is.tree/is.cladisticmatrix
  # - Trees match cladistic_matrix names.
  
  # TO DO
  # - Conditional for invariant character plus autapomorphic ones (fast solution that escapes estimates of other values)
  # - Allow for true polymorphisms (maybe conditionals for counting changes?) This is a hard problem!
  # - How to allow for missing or uncertainty? This bit should be easy as set all states to zero, or all uncertain states to zero at tips.
  # - All (format?)/ACCTRAN/DELTRAN/(Random?)/Branch lengths as weights to help collapse?/Uncertainity (e.g., 0/1)
  # - Missing and inapplicables shpuld probably be assigned to nodes on a third pass "down" the tree (tips to roots) such that any all-descendants set gets assigned an NA/"" too. Otherwise sets false certainity by "bleeding" an ancestral value "up" the tree.
  # - Tip states as entered need preserving separately as some modification will need to happen to the variable if certain options are chosen.
  # - Make options match estimate_ancestral_states:
  #   - @param estimate_all_nodes Logical that allows the user to make estimates for all ancestral values. The default (\code{FALSE}) will only make estimates for nodes that link coded terminals (recommended).
  #   - @param estimate_tip_values Logical that allows the user to make estimates for tip values. The default (\code{FALSE}) will only makes estimates for internal nodes (recommended).
  #   - @param polymorphism_behaviour One of either "equalp" or "treatasmissing".
  #   - @param uncertainty_behaviour One of either "equalp" or "treatasmissing".
  
  # CHECKS TO WRITE
  # - states are discrete
  # - states in tip states are all available in stepmatrix and vice versa
  # - stepmatrix values are zero or positive (don't have to be integers?)
  # - Check there are enough states to do something with!
  # - Put checks at higher level to avoid repeating for every character in a matrix (slow?)
  
  
  
  ### GONNA HAVE TO WORK OUT HOW TO SET THIS AUTOMATICALY?
  ### IT IS HARD THOUGH AS NODE NUMBERS VARY ACROSS TREES! HMMMM.... MAYBE DEEP SIX THIS FOR NOW AND REYURN TO IT LATER?
  ### ALSO NEEDS TO BE ADDED TO OUTPUT LATER SOMEHOW
  ### MAYBE MAKE THIS A CLASS AND/OR "APPEND" IT TO A "phylo" CLASS OBJECT
  node_constraints = NULL

  # If trees is a single topology then reformat as a list of length one:
  if (class(x = trees) == "phylo") trees <- list(trees)
  
  # Set include_polymorphisms for make_stepmatrix by polymorphism_behaviour or uncertainty_behaviour choice:
  include_polymorphisms <- ifelse(test = polymorphism_behaviour == "polymorphism" || uncertainty_behaviour == "polymorphism", yes = TRUE, no = FALSE)
  
  # Combine all blocks into a single input matrix:
  single_input_matrix <- list(
    matrix = do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$matrix)),
    ordering = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$ordering))),
    weights = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$character_weights))),
    minima = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$minimum_values))),
    maxima = unname(obj = do.call(what = c, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], FUN = function(x) x$maximum_values)))
  )
  
  # If inapplicables are to be treated as missing:
  if (inapplicables_as_missing) {
    
    # As long as inapplicables are found then convert them to NA:
    if (any(x = single_input_matrix$matrix == "", na.rm = TRUE)) single_input_matrix$matrix[single_input_matrix$matrix == ""] <- NA
    
  }
  
  # If polymorphisms are present in the matrix:
  if (length(x = grep(pattern = "&", x = single_input_matrix$matrix)) > 0) {
    
    if (polymorphism_behaviour == "missing") single_input_matrix$matrix[grep(pattern = "&", x = single_input_matrix$matrix)] <- NA
    if (polymorphism_behaviour == "uncertainty") single_input_matrix$matrix[grep(pattern = "&", x = single_input_matrix$matrix)] <- gsub(pattern = "&", replacement = "/", x = single_input_matrix$matrix[grep(pattern = "&", x = single_input_matrix$matrix)])
    if (polymorphism_behaviour == "random") single_input_matrix$matrix[grep(pattern = "&", x = single_input_matrix$matrix)] <- unlist(x = lapply(X = as.list(x = single_input_matrix$matrix[grep(pattern = "&", x = single_input_matrix$matrix)]), FUN = function(x) sample(x = strsplit(x = x, split = "&")[[1]], size = 1)))

  }
  
  # If uncertainties are present in the matrix:
  if (length(x = grep(pattern = "/", x = single_input_matrix$matrix)) > 0) {
    
    if (uncertainty_behaviour == "missing") single_input_matrix$matrix[grep(pattern = "/", x = single_input_matrix$matrix)] <- NA
    if (uncertainty_behaviour == "polymorphism") single_input_matrix$matrix[grep(pattern = "/", x = single_input_matrix$matrix)] <- gsub(pattern = "/", replacement = "&", x = single_input_matrix$matrix[grep(pattern = "/", x = single_input_matrix$matrix)])
    if (uncertainty_behaviour == "random") single_input_matrix$matrix[grep(pattern = "/", x = single_input_matrix$matrix)] <- unlist(x = lapply(X = as.list(x = single_input_matrix$matrix[grep(pattern = "/", x = single_input_matrix$matrix)]), FUN = function(x) sample(x = strsplit(x = x, split = "/")[[1]], size = 1)))
    
  }

  # Build character list:
  character_list <- lapply(
    X = as.list(x = 1:length(x = single_input_matrix$ordering)),
    FUN = function(x) list(
      tip_states = single_input_matrix$matrix[, x],
      stepmatrix = if (length(x = grep(pattern = "step", x = single_input_matrix$ordering[x])) == 1) {
        cladistic_matrix$topper$step_matrices[[single_input_matrix$ordering[x]]]
      } else {
        make_stepmatrix(
          min_state = single_input_matrix$minima[x],
          max_state = single_input_matrix$maxima[x],
          character_type = single_input_matrix$ordering[x],
          include_polymorphisms = include_polymorphisms,
          polymorphism_shape = polymorphism_shape,
          polymorphism_distance = polymorphism_distance,
          state_ages = state_ages,
          dollo_penalty = dollo_penalty
        )
      },
      weight = single_input_matrix$weights[x],
      inapplicables_as_missing = inapplicables_as_missing
    )
  )
  
  # deal with dollo and irreversible root forcing? (and polymorphisms)
  # conditional to recode polymorphisms if stratigraphic character


  
  
  
  
  
  
  # NEED DIFFERENT FUNCTION FOR CONTINUOUS CHARACTERS (AS STEPMATRIX MAKES LITTLE SENSE!).
  # CAN USE asr_squared_change_parsimony in castor but does mean adding a dependency
  
  first_pass <- function(tree, tip_states, stepmatrix, weight, node_constraints = NULL) {
    
    # Reorder tip states (1 to N) and store an unmodified version:
    pristine_tip_states <- tip_states <- tip_states[tree$tip.label]
    
    # Store pristine input tree:
    input_tree <- tree
    
    # If there are no branch durations set these as all one:
    if (is.null(x = tree$edge.length[1])) tree$edge.length <- rep(x = 1, length.out = nrow(x = tree$edge))
    
    # Reformat tip states as character vector:
    tip_states <- as.character(x = tip_states)
    
    # Re-add names to tip states:
    names(tip_states) <- tree$tip.label
    
    # Establish number of tips:
    tip_count <- ape::Ntip(phy = tree)
    
    # Establish total number of nodes (terminal and internal):
    node_count <- tip_count + tree$Nnode
    
    # Initialise node_values matrix:
    node_values <- matrix(data = 0, nrow = node_count, ncol = stepmatrix$size, dimnames = list(c(), colnames(stepmatrix$stepmatrix)))
    
    # Begin by inserting Inf as default tip value:
    node_values[1:tip_count, ] <- Inf
    
    # Populate tip values (zeroes for assigned state(s):
    for(i in 1:tip_count) node_values[i, match(x = strsplit(x = tip_states[i], split = "/")[[1]], table = colnames(x = node_values))] <- 0
    
    # If any tips are missing then set all states for those as zero:
    if (any(x = is.na(x = tip_states))) node_values[which(x = is.na(x = tip_states)), ] <- 0
    
    # If any node_constraints are supplied:
    if (length(x = node_constraints) > 0) {
      
      # Apply node constraints to node_values:
      node_values[as.numeric(x = names(x = node_constraints)), ] <- do.call(
        what = rbind,
        args = mapply(
          FUN = function(x, node_number) x[as.numeric(node_number), , drop = FALSE],
          x = lapply(
            X = as.list(x = node_constraints),
            FUN = function(x) {
              all_inf_node_values <- node_values
              all_inf_node_values[] <- Inf
              allowable_states <- strsplit(x = x, split = "/")[[1]]
              disallowable_states <- setdiff(x = colnames(x = node_values), y = allowable_states)
              cbind(x = all_inf_node_values[, disallowable_states, drop = FALSE], y = node_values[, allowable_states, drop = FALSE])[, colnames(x = node_values), drop = FALSE]
            }
          ),
          node_number = names(node_constraints),
          SIMPLIFY = FALSE
        )
      )
    }
    
    # First pass (traverse tree from tips to root):
    for(needle in (tip_count + node_count - tip_count):(tip_count + 1)) {
      
      # Find decsendants of current node:
      descendants <- tree$edge[tree$edge[, 1] == needle, 2]
      
      # Calculate and store new node values:
      node_values[needle, ] <- apply(
        X = rbind(
          node_values[needle, ],
          unlist(x = lapply(
            X = as.list(x = colnames(x = node_values)),
            FUN = function(fromstate) {
              sum(x = unlist(x = lapply(
                X = as.list(x = descendants),
                FUN = function(descendant) min(x = node_values[descendant, ] + stepmatrix$stepmatrix[fromstate, ]))))
            })
          )
        ),
        MARGIN = 2,
        FUN = max
      )
    }
    
    # Check there are no all Inf rows that confound further usage and stop and warn user:
    if (any(x = apply(X = node_values == Inf, MARGIN = 1, FUN = all))) stop("Options have lead to an all Inf row in node_values. Consider removing or reducing node_constraints and try again.")
    
    # Store tree length:
    tree_length <- min(x = node_values[tip_count + 1, ]) * weight
    
    # Build input value list:
    input_values <- list(tree = tree, tip_states = pristine_tip_states, stepmatrix = stepmatrix, weight = weight, inapplicables_as_missing = inapplicables_as_missing, node_constraints = node_constraints)
    
    # Return output:
    list(length = tree_length, node_values = node_values)
    
  }


  #
  trees_list <- lapply(
    X = trees,
    FUN = function(y) {
      first_pass_values = lapply(
        X = character_list,
        function(x) first_pass(
          tree = y,
          tip_states = x$tip_states,
          stepmatrix = x$stepmatrix,
          weight = x$weight,
          node_constraints = node_constraints
        )
      )
    }
  )
  
  # Build matrix of individual character tree lengths for output:
  character_lengths <- do.call(what = cbind, args = lapply(X = trees_list, FUN = function(y) unlist(x = lapply(X = y, FUN = function(x) x$length))))
  
  # Return compiled output:
  list(
    input_trees = trees,
    input_matrix = cladistic_matrix,
    input_options = list(
      inapplicables_as_missing = inapplicables_as_missing,
      polymorphism_behaviour = polymorphism_behaviour,
      uncertainty_behaviour = uncertainty_behaviour,
      polymorphism_shape = polymorphism_shape,
      polymorphism_distance = polymorphism_distance,
      state_ages = state_ages,
      dollo_penalty = dollo_penalty
    ),
    stepmatrices = lapply(X = character_list, FUN = function(x) x$stepmatrix),
    character_matrix = single_input_matrix,
    character_lengths = character_lengths,
    character_weights = single_input_matrix$weights,
    tree_lengths = apply(X = character_lengths, MARGIN = 2, FUN = sum),
    node_values = lapply(X = trees_list, FUN = function(y) lapply(X = y, FUN = function(x) x$node_values))
  )

}

#trees = ape::read.tree("http://www.graemetlloyd.com/mpts/Hooker_2014a.tre")
#cladistic_matrix = Claddis::read_nexus_matrix("http://www.graemetlloyd.com/nexus/Hooker_2014a.nex")
#inapplicables_as_missing = FALSE
#polymorphism_behaviour = "uncertainty"
#uncertainty_behaviour = "uncertainty"
#polymorphism_shape = "hypersphere"
#polymorphism_distance = "great_circle"
#state_ages = c()
#dollo_penalty = 999

#calculate_tree_length(trees = trees, cladistic_matrix = cladistic_matrix, inapplicables_as_missing = inapplicables_as_missing, polymorphism_behaviour = "uncertainty", uncertainty_behaviour = "uncertainty", polymorphism_shape = polymorphism_shape, polymorphism_distance = polymorphism_distance, state_ages = state_ages, dollo_penalty = dollo_penalty)$tree_lengths



# FIXING NODE VALUES RELATES TO ANCESTRY STATEMENTS TOO (EFFECTIVELY CAN BE USED TO ASSUME ANCESTOR - LINKS TO STRATOCLADISTICS)
























# NEED TO COUNT DOLLO CHARACTERS DIFFERENTLY FOR PARSIMONY OR WILL OVERWHELM BROADER CHARACTER TRENDS.
# ALSO NEEDS TO BE RECORDED AS SUCH FOR STOCHASTIC CHARACTER MAPS AND HOMOPLASY METRICS.

# ADD SWOFFORD REF TO DESCRIPTION FILE (PROLLY NO DOI!)








# OUTPUT:
# - Most of this needs to be put in different functions!
# - Ancestral states
# - Ambiguities?
# - Stepmatrix used
# - Tree used
# - Tip states used
# - Tree lengths
# - How root value was chosen (arbitrary forced, unambiguous?)
# - Algorithm (parsimony, ML, etc.) - for later integration into single ASR function with other options (likelihood, rerooting etc.)
# - CI? (Characters as vector then can do summed for matrix)
# - RI? (Need to know max value which is just fit on star tree)
# - N reversals (Will depend on root and "direction")
# - N parallelisms (Will depend on root and "direction")
# - Character map? (need to think carefully what this format should be and hence how stuff like ASE needs to inform it.)
# - MAKE THIS A CLASS TOO!



# WRITE SOME KIND OF CHECK FOR SELF-CONSISTENCY IN STEP MATRICES. I.E., THE DIRECT PATH IS ALWAYS THE SAME LENGTH OR SHORTER THAN AN INDIRECT ONE. E.G., If 0 -> 2 is 4 steps, but 0->1->2 is only 3 steps. (Citation is MacClade 4 manual!)
# OPTION TO CONSTRAIN PARTICULAR NODES AND/OR MAKE TAXA ANCESTORS (MIGHT BE TRICKY AND A DOWN THE LINE ADDITION!) HMMM. OR HAVE CONDITIONAL THAT SETS INFINITE VALUES TO THE NON PRESENT STATES AT THIS NODE?
# MORE OUTPUT OR SUMMARY OPTIONS COULD BE N CHANGES ON EACH BRANCH ACROSS MPRS, MEAN OF SAME, MIN AND MAX OF SAME. SIMILAR FOR CHANGE TYPES I.E., STEP MATRIX BUT ACTUALLY FREQUENCY OF CHANGES FOR EACH TRANSITION (THIS IS REALLY AN SCM OUTPUT).
# AS RATE TENDS TOWARD INFINITY THEN RECONSTRUCTION AT A NODE TENDS TOWARDS EQUAL FREQUENCY OF EACH STATE (I.E., MAXIMAL UNCERTAINTY).
# ALLOW TREES N STEPS LONGER SOMEHOW? OR RATHER MPRs ONE STEP LONGER, OR SOME OTHER VALUE OF USERS CHOOSING.
# "This illustrates the fact that ACCTRAN and DELTRAN do not always choose a single one of the most parsimonious reconstructions." MacClade 4 manual, page 99).
# "Note that MacClade, unlike PAUP*, does not choose the lowest-valued state at the root to begin these processes. Thus MacClade's ACCTRAN and DELTRAN may not fully resolve ambiguity in the ancestral state reconstruction."
# INTERMEDIATES ARE GONNA MATTER IF MOVING TO CHARACTER MAPS. E.G., IF ONLY 0 AND 2 ARE SAMPLED BUT A CHARACTER IS ORDERED THEN THERE ARE TWO CHANGES ALONG THE BRANCH NOT ONE TWO-STEP CHANGE.
# POLYMORPHISM APPROACH SHOULDN'T CONFOUND THE "NORMAL" BEHAVIOUR OF AN ORDERED OR UNORDERED CHARACTER.
# TEST WEIGHTING OF STRATOCLADISTICS BY USING A STRATIGRAPHIC CHARACTER AND WEIGHTING IT MULTIPLE WAYS (AS A SLIDER) AND SHOW HOW IT EFFECTS PARSIMONY VS STRAT CONGRUENCE TREE LANDSCAPES
# MONOFURCATIONS IDEA IN CASTOR IS INTERESTING AS ALLOWS POINT ESTIMATES ALONG BRANCHES (E.G., FOR SCM OF A CONTINUOUS CHARACTER)







### MPR function below here

#' Determine maximum parsimony ancestral state reconstruction(s)
#'
#' @description
#'
#' Given a tree, or set of trees, and a cladistic matrix returns all most parsimonious reconstruction(s) for every character.
#'
#' @param trees A tree (phylo object) or set of trees (multiPhylo object).
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}. These should be discrete with labels matching the tip labels of \code{tree}.
#' @param estimate_all_nodes Logical that allows the user to make estimates for all ancestral values. The default (\code{FALSE}) will only make estimates for nodes that link coded terminals (recommended).
#' @param estimate_tip_values Logical that allows the user to make estimates for tip values. The default (\code{FALSE}) will only makes estimates for internal nodes (recommended).
#' @param inapplicables_as_missing Logical that decides whether or not to treat inapplicables as missing (TRUE) or not (FALSE, the default and recommended option).
#' @param polymorphism_behaviour One of either "missing", "uncertainty", "polymorphism", or "random". See details.
#' @param uncertainty_behaviour One of either "missing", "uncertainty", "polymorphism", or "random". See details.
#' @param polymorphism_shape Argument passed to \link{make_stepmatrix}.
#' @param polymorphism_distance Argument passed to \link{make_stepmatrix}.
#' @param state_ages Argument passed to \link{make_stepmatrix}.
#' @param dollo_penalty Argument passed to \link{make_stepmatrix}.
#'
#' @details
#'
#' Users may also wish to refer to the more complex whole-matrix, likelihood-based function \link{estimate_ancestral_states}. Although, note that eventually parsimony will simply be an option to that function.
#'


# REFERENCE castor asr_max_parsimony and phangorn acctran as alternates

#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. \emph{In} R. L. Mayden (ed.) Systematics, Historical Ecology, and North American Freshwater Fishes. Stanford University Press, Stanford, p187-223.
#'
#' @return
#'
#' A list with multiple components, including:
#'
#' \item{ITEM}{ITEM EXPLANATION.}
#'
#' @seealso
#'
#' \link{estimate_ancestral_states}
#'
#' @examples
#'
#' # TO DO
#'
#' @export reconstruct_ancestral_states
reconstruct_ancestral_states <- function(trees, cladistic_matrix, estimate_all_nodes = FALSE, estimate_tip_values = FALSE, inapplicables_as_missing = FALSE, polymorphism_behaviour = "uncertainty", uncertainty_behaviour = "uncertainty", polymorphism_shape, polymorphism_distance, state_ages, dollo_penalty) {
  
  # COLLAPSE OPTIONS (MAYBE FOR ANOTHER FUNCTION?):
  # - ACCTRAN
  # - DELTRAN
  # - MINF (Swofford and Maddison 1987) - I THINK NOT!
  # - MINSTATE *shrug emoji*
  # - MAXSTATE *shrug emoji*
  # - WEIGHTED BY BRANCH DURATION? (E.G., ((0:1,1:4)); SHOULD FAVOUR A 0 ROOT STATE). NEED WAY TO "SCORE" THESE SUCH THAT CAN FIND MPR WITH BEST SCORE.
  
  # TO ADD INTO NODE ESTIMATE VERSION OF FUNCTION:
  #estimate_all_nodes <- FALSE
  #estimate_tip_values <- FALSE # Should not happen for inapplicables when they are treated as inapplicables

  # WILL NEED TO MODIFY BELOW TO DEAL WITH UNCERTAINTIES AND POLYMORPHISMS AT TIPS
  
  
  
  
  tree_length_output <- calculate_tree_length(
    trees = trees,
    cladistic_matrix = cladistic_matrix,
    inapplicables_as_missing = inapplicables_as_missing,
    polymorphism_behaviour = polymorphism_behaviour,
    uncertainty_behaviour = uncertainty_behaviour,
    polymorphism_shape = polymorphism_shape,
    polymorphism_distance = polymorphism_distance,
    state_ages = state_ages,
    dollo_penalty = dollo_penalty
  )
  
  
  # THEN:
  # second_pass # All MPRs (OUTPUTS: all_mprs)
  # third_pass # Missing/inapplicables applied to internal nodes (OUTPUTS: modified_mprs)


  
  # Build node estimates from single state results:
  node_estimates <- matrix(
    data = apply(
      X = node_values,
      MARGIN = 1,
      FUN = function(x) {
        x <- names(x = x[x == min(x = x)])
        ifelse(
          test = length(x = x) == 1,
          yes = x,
          no = NA
        )
      }),
    ncol = 1,
    nrow = node_count
  )
  
  # Update tip states with their original values:
  node_estimates[1:tip_count] <- tip_states
  
  ### ABOVE NEEDS TO BE INPUT TIP STATES NOT MODIFIED VERSION! BUT LATER FOR OUTPUT? NEED TO RECORD SOMEWHERE WHAT HAS HAPPENED THOUGH SO CAN PASS TO A CHARACTER MAP FUNCTION AND TREAT EVERYTHING CORRECTLY.
  ### ACTUALLY SHOULD PROBABLY BE
  
  # Only continue if there is any ambiguity at internal nodes:
  if (any(x = is.na(x = node_estimates[(tip_count + 1):node_count]))) {
    
    # Find all possible root states:
    possible_root_states <- colnames(x = node_values)[node_values[tip_count + 1, ] == min(x = node_values[tip_count + 1, ])]
    
    # Make new node estimates with every possible root state:
    node_estimates <- do.call(what = cbind, args = lapply(X = as.list(x = possible_root_states), FUN = function(x) {y <- node_estimates; y[tip_count + 1, ] <- x; y}))
    
    # For each internal node above the root to the tips:
    for(needle in (tip_count + 2):(tip_count + node_count - tip_count)) {
      
      # Establish ancestor of current node:
      ancestor_node <- tree$edge[tree$edge[, 2] == needle, 1]
      
      # Reformat node_estimates as a list to enable lapply:
      node_estimates <- split(x = node_estimates, f = col(x = node_estimates))
      
      # Permute all possible values and reformat as matrix:
      node_estimates <- do.call(what = cbind, args = lapply(X = node_estimates, FUN = function(x) {
        
        # Get updated tree lengths for current node as per Swofford and Maddison 1992 second pass:
        updated_tree_lengths <- stepmatrix$stepmatrix[x[ancestor_node], ] + node_values[needle, ]
        
        # Store all possible most parsimonious state(s) for current node:
        possible_states <- names(x = updated_tree_lengths[updated_tree_lengths == min(x = updated_tree_lengths)])
        
        # Store total number of possible states:
        n_states <- length(x = possible_states)
        
        # Create new node estimates to store possible state(s):
        new_estimates <- matrix(data = rep(x = x, times = n_states), ncol = n_states)
        
        # Add possible state(s)
        new_estimates[needle, ] <- possible_states
        
        # Return new estimates only:
        new_estimates
        
      }))
      
    }
    
  }
  
  # Return compiled output:
  list(length = tree_length, most_parsimonious_reconstructions = node_estimates, input_tree = input_tree)
  
}

# HOMOPLASY FUNCTION
# - Option to define outgroup (root state?) as will affect min steps value
# - Otherwise need to do path both ways, e.g., 0->1->2 and 2->1->0 in case of asymmetric characters
# - Keep check for variance as no changes at all if invariant.
# - How to handle uncertainties?

