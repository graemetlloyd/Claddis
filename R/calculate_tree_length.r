#' Calculates the parsimony length of a set of phylogenetic tree(s)
#'
#' @description
#'
#' Given a tree, or set of trees, and a cladistic matrix returns their parsimony length in number of steps.
#'
#' @param trees A tree (\code{phylo} object) or set of trees (\code{multiPhylo} object).
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}. These should be discrete with rownames (ytaxon labels) matching the tip labels of \code{trees}.
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
#' This function operates on a phylogenetic tree, or trees (in \code{ape} format), and a cladistic matrix (in \code{cladisticMatrix} format). However, the algorithm used is based on the generalised stepmatrix approach of Swofford and Maddison (1992) and hence stepmatrices need to be defined for each character (this is done internally by calling \link{make_stepmatrix}), and some of the options are merely passed to this function.
#'
#' \bold{Algorithm}
#'
#' Technically the Swofford and Maddison (1992) algorithm is designed for ancestral state reconstruction, but as its' first pass of the tree assigns lengths for each possible state at each node the minimum value of these options at the root is also the tree length for that character and hence by skipping the later steps this can be used as a tree length algorithm by simply summing the values across each character. The choice of the Swofford and Maddison algorithm, rather than the Wagner or Fitch algorithms (for ordered and unordered characters, respectively) is to generalize to the broadest range of character types, including asymmetric characters (Camin-Sokal, Dollo, stratigraphic), custom character types (specified using stepmatrices or character state trees), as well as to any resolution of tree (i.e., including multifurcating trees - important for establishing maximum step counts for homoplasy indices). The only restriction here is that the tree must be rooted such that time's arrow is explicitly present. This is essential, as the root defines the lengths across the whole tree, but also for asymmetric characters directionality must be explicit, as well as some downstream appraoches (such as ACCTRAN and DELTRAN). The two obvious drawbacks to this algorithm are that it can be slower and that it is not appropriate for unrooted trees.
#'
#' \bold{Stepmatrices and stepmatrix options}
#'
#' Stepmatrices are described in detail in the \link{make_stepmatrix} manual, as are the options that are passed from this function to that one. Thus, the user is directed there for a more in-depth discussion of options.
#'
#' \bold{Inapplicable and missing characters}
#'
#' In practice these two character types are treated the same for length calculations - in effect these are "free" characters that do not constrain the tree length calculation in the same way that a coded character would (because a coded character's transition cost must be accounted for; Swofford and Maddison 1992). Note that there \emph{are} reasons to take differences into account in phylogenetic inference itself (see papers by Brazeau et al. 2019 and Goloboff et al. in press). The option to treat them differently here is therefore only important in terms of downstream analyses, such as ancestral state reconstruction (see \link{reconstruct_ancestral_states} for details).
#'
#' \bold{Polymorphisms and uncertainties}
#'
#' Polymorphisms (coded with empersands between states) and uncertainties (coded with slashes between states) can be interpreted in different ways, including those that affect estimates of tree length. Hence four options are provided to the user here:
#'
#'  \enumerate{
#'    \item{Missing (\code{polymorphism_behaviour = "missing"} or \code{uncertainty_behaviour = "missing"}). }{Here polymorphisms are simply replaced by the missing character (\code{NA}). This removes polymorphisms and uncertainties from the calculation process completely (likely leading to undercounts), and hence is not generally recommended.}
#'    \item{Uncertainty (\code{polymorphism_behaviour = "uncertainty"} or \code{uncertainty_behaviour = "uncertainty"}). }{This is the intended use of uncertain codings (e.g., \code{0/1}) and constrains the tree length calculation to having to explain the \emph{least} costly transition of those in the uncertainty. This is recommended for uncertain characters (although note that it biases the result towards the shortest possible length), but not truly polymorphic characters (as one or more state acquisitions are being missed, see Nixon and Davis 1991 and \link{make_stepmatrix} for discussion of this). This is also - to the best of my knowledge - the approach used by most parsimony software, such as PAUP* (Swofford 2003) and TNT (Goloboff et al. 2008; Goloboff and Catalano 2016).}
#'    \item{Polymorphism (\code{polymorphism_behaviour = "polymorphism"} or \code{uncertainty_behaviour = "polymorphism"}). }{If polymorphisms are real then some means of accounting for the changes that produce them seems appropriate, albeit difficult (see Nixon and Davis 1991 and Swofford and Maddison 1992 for discussions). If this option is applied it triggers the downstream options in \emph{make_stepmatrix} (by internally setting \code{include_polymorphisms = TRUE}), and the user should look there for more information. This is tentatively recommended for true polymorphisms (but note that it complicates interpretation), but not uncertainties.}
#'    \item{Random (\code{polymorphism_behaviour = "random"} or \code{uncertainty_behaviour = "random"}). }{Another means of dealing with multiple-state characters is simply to sample a single state at random for each one, for example as Watanabe (2016) did with their PERDA algorithm. This simplifies the process, but also logically requires running the function multiple times to quantify uncertainty. This is not recommended for true polymorphisms (as interpretation is confounded), but may be appropriate for a less downwards biased tree count than \code{"uncertainty"}.}
#' }
#'
#' These choices can also effect ancestral state estimation (see \link{reconstruct_ancestral_states}).
#'
#' \bold{Polytomies}
#'
#' Polytomies are explicitly allowed by the function, but will always be treated as "hard" (i.e., literal multifurcations). Note that typically these will lead to higher tree lengths than fully bifurcating trees and indeed that the maximum number of steps is typically calculated from the star tree (single multifurcation).
#'
#' \bold{Further constraints}
#'
#' In future the function will allow restrictions to be placed on the state at particular internal nodes. This can have multiple applications, including (for example) treating some taxa as ancestral such that their states are directly tied to specific nodes, e.g., in stratocladistics (Fisher 1994; Marcot and Fox 2008).
#'
#' \bold{Character weights}
#'
#' Tree lengths output already include corrections for character weights as supplied in the \code{cladistic_matrix} input. So, for example, if a binary character has two steps on the tree, but is weighted five then it will contribute a value of ten steps to the result.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Brazeau, M. D., Guillerme, T. and Smith, M. R., 2019. An algorithm for morphological phylogenetic analysis with inapplicable data. \emph{Systematic Biology}, \bold{68}, 619-631.
#'
#' Fisher, D. C., 1994. Stratocladistics: morphological and temporal patterns and their relation to phylogenetic process. In L. Grande and O. Rieppel (eds.), \emph{Interpreting the Hierarchy of Nature}. Academic Press, San Diego. pp133–171.
#'
#' Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics/ \emph{Cladistics}, \bold{32}. 221-238
#'
#' Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. \emph{Cladistics}, \bold{24}, 774-786.
#'
#' Goloboff, P. A., De Laet, J., Rios-Tamayo, D. and Szumik, C. A., in press. A reconsideration of inapplicable characters, and an approximation with step‐matrix recoding. \emph{Cladistics}.
#'
#' Marcot, J. D. and Fox, D. L., 2008. StrataPhy: a new computer program for stratocladistic analysis. \emph{Palaeontologia Electronica}, \bold{11}, 5A.
#'
#' Nixon, K. C. and Davis, J. I., 1991. Polymorphic taxa, missing values and cladistic analysis. \emph{Cladistics}, \bold{7}, 233-241.
#'
#' Swofford, D. L., 2003. \emph{PAUP*. Phylogenetic Analysis Using Parsimony (*and Other Methods). Version 4}. Sinauer Associates, Sunderland, Massachusetts.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. \emph{In} R. L. Mayden (ed.) Systematics, Historical Ecology, and North American Freshwater Fishes. Stanford University Press, Stanford, p187-223.
#'
#' Watanabe, A., 2016. The impact of poor sampling of polymorphism on cladistic analysis. \emph{Cladistics}, \bold{32}, 317-334.
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
#' \link{make_stepmatrix}, \link{reconstruct_ancestral_states}
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
  
  # CHECKS TO WRITE
  # - states are discrete
  # - states in tip states are all available in stepmatrix and vice versa
  # - stepmatrix values are zero or positive (don't have to be integers?)
  # - Check there are enough states to do something with!
  # - Put checks at higher level to avoid repeating for every character in a matrix (slow?)
  
  # deal with dollo and irreversible root forcing? (and polymorphisms)
  # conditional to recode polymorphisms if stratigraphic character
  # NEED DIFFERENT FUNCTION FOR CONTINUOUS CHARACTERS (AS STEPMATRIX MAKES LITTLE SENSE!).
  # CAN USE asr_squared_change_parsimony in castor but does mean adding a dependency
  
  # NEED TO COUNT DOLLO CHARACTERS DIFFERENTLY FOR PARSIMONY OR WILL OVERWHELM BROADER CHARACTER TRENDS.

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
  
  # Subfunction to do Swofford and Maddison 1992 first pass of tree:
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
#' Lloyd, G. T., 2018. Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. \emph{Palaeontology}, \bold{61}, 637-645.
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
  
  # ALLOW RECONSTRUCTIONS N STEPS LONGER SOMEHOW?
  
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

# WHAT IF ALL CHARACTERS ARE UNCERTAIN!!!!!!
# MANUAL: Swofford and Maddison (1992; p210): "Both [ACCTRAN and DELTRAN] start by either assuming that the ancestral state is known or by choosing a state from the set of optimal assignments to the root node" - i.e., always an arbitrainess problem at root!


# ALSO NEEDS TO BE RECORDED AS SUCH FOR STOCHASTIC CHARACTER MAPS AND HOMOPLASY METRICS.

# ADD SWOFFORD REF TO DESCRIPTION FILE (PROLLY NO DOI!)



# Record in manual that weighted parsimony is folly, but could use to chose between MPRs (maybe?)




# OUTPUT:
# - Ambiguities?
# - Tip states used
# - How root value was chosen (arbitrary forced, unambiguous?)
# - Algorithm (parsimony, ML, etc.) - for later integration into single ASR function with other options (likelihood, rerooting etc.)
# - Character map? (need to think carefully what this format should be and hence how stuff like ASE needs to inform it.)
# - MAKE THIS A CLASS TOO!



# MORE OUTPUT OR SUMMARY OPTIONS COULD BE N CHANGES ON EACH BRANCH ACROSS MPRS, MEAN OF SAME, MIN AND MAX OF SAME. SIMILAR FOR CHANGE TYPES I.E., STEP MATRIX BUT ACTUALLY FREQUENCY OF CHANGES FOR EACH TRANSITION (THIS IS REALLY AN SCM OUTPUT).
# AS RATE TENDS TOWARD INFINITY THEN RECONSTRUCTION AT A NODE TENDS TOWARDS EQUAL FREQUENCY OF EACH STATE (I.E., MAXIMAL UNCERTAINTY).
# "This illustrates the fact that ACCTRAN and DELTRAN do not always choose a single one of the most parsimonious reconstructions." MacClade 4 manual, page 99).
# "Note that MacClade, unlike PAUP*, does not choose the lowest-valued state at the root to begin these processes. Thus MacClade's ACCTRAN and DELTRAN may not fully resolve ambiguity in the ancestral state reconstruction."
# INTERMEDIATES ARE GONNA MATTER IF MOVING TO CHARACTER MAPS. E.G., IF ONLY 0 AND 2 ARE SAMPLED BUT A CHARACTER IS ORDERED THEN THERE ARE TWO CHANGES ALONG THE BRANCH NOT ONE TWO-STEP CHANGE.
# TEST WEIGHTING OF STRATOCLADISTICS BY USING A STRATIGRAPHIC CHARACTER AND WEIGHTING IT MULTIPLE WAYS (AS A SLIDER) AND SHOW HOW IT EFFECTS PARSIMONY VS STRAT CONGRUENCE TREE LANDSCAPES
# MONOFURCATIONS IDEA IN CASTOR IS INTERESTING AS ALLOWS POINT ESTIMATES ALONG BRANCHES (E.G., FOR SCM OF A CONTINUOUS CHARACTER)


# HOMOPLASY FUNCTION:
# - Option to define outgroup (root state?) as will affect min steps value
# - Otherwise need to do path both ways, e.g., 0->1->2 and 2->1->0 in case of asymmetric characters
# - Keep check for variance as no changes at all if invariant.
# - How to handle uncertainties?
# TO INCLUDE:
# - CI? (Characters as vector then can do summed for matrix)
# - RI? (Need to know max value which is just fit on star tree)
# - N reversals (Will depend on root and "direction")
# - N parallelisms (Will depend on root and "direction")
# - Record character types too (Dollo important as really want to count single steps but will affect min and max values)

