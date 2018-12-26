#' Discrete character rates across trees, time, and character types
#'
#' \bold{Introduction}
#'
#' Morphological change can be captured by discrete characters and their evolution modelled as occurring along the branches of a phylogenetic tree. This function takes as primary input a character-taxon matrix of discrete characters (in the format imported by \link{ReadMorphNexus}) and a time-scaled phylogenetic tree (in the format of \pkg{paleotree} or \pkg{strap}) and begins by inferring ancestral states at the tree's internal nodes using the \link{AncStateEstMatrix} function. From here changes along individual branches can be estimated (only the minimum number of changes are inferred; see \link{GetAllStateChanges} for an alternative but unfininshed approach) and hence rates can be calculated.
#'
#' A discrete character rate can be expressed as the mean number of changes per character per million years and can be calculated for a branch (edge) of the tree, a clade (a mean rate for the edges descended from a single node), a character partition (the mean rate for a subset of the characters across all edges), or, most complex, the mean rate across the edges (or parts of edges) present in a time bin (defined by two values denoting the beginning and end of the time bin). In an ideal scenario these rates could be compared at face value, but that would require a large number of characters and very minimal (or zero) missing data. I.e., at an extreme of missing data if only one character can be observed along a branch it will either change (the maximum possible rate of evolution) or it will not (the minimum possible rate of evolution). In such cases it would be unwise to consider either outcome as being a significant departure from the mean rate.
#'
#' Because of these complications Lloyd et al. (2012) devised tests by which the significance of an edge (or other paritioning of the data, i.e., a clade, time bin etc.) could be considered to be significantly high or low in comparison to the mean rate for the whole tree (i.e., whether a two-rate model could be considered more likely than a one-rate model). This is achieved through a likelihood ratio test:
#'
#' \deqn{LR = value of likehood function under the null (one-rate) hypothesis / maximum possible value of likehood function under the alternative (two-rate) hypotheses}
#'
#' Typically we might expect the two hypotheses to be well defined a priori. E.g., an expectation that a specific branch of the tree might have a higher or lower rate than background due to some evolutionary shift. However, Lloyd et al. (2012) instead provided an exploratory approach whereby every possible one edge value was compared with the rate for the rest of the tree (and the equivalent with clades adn time bins). This was the default in Claddis up to version 0.2, but this has now been replaced (since version 0.3) with a more customisable set of options that allows different types of hypotheses (e.g., partitioning the data by character), as well as more complex hypotheses (e.g., a 3-rate model), to be tested.
#'
#' \bold{The four types of rate hypothesis}
#'
#' Following Cloutier (1991), Lloyd (2016) extended the four main types of rate hypotheses to:
#'
#' \enumerate{
#'   \item A branch rate (available here with the \code{BranchPartitionsToTest} option).
#'   \item A clade rate (available here with the \code{CladePartitionsToTest} option).
#'   \item A time bin rate (available here with the \code{TimeBinPartitionsToTest} option).
#'   \item A character partition rate (available here with the \code{CharacterPartitionsToTest} option).
#' }
#'
#' In Claddis (>=0.3) these partitions are defined as a list of lists where only the first N - 1 partitions need be defined. E.g., if comparing the first edge value to the rest of the tree then the user only needs to define the value "1" and the function will automatically add a second partition containing all other edges. This can be set with the option \code{BranchPartitionsToTest = lapply(as.list(1), as.list)}. Similarly, to do what Lloyd et al. (2012) did and repeat the test for every edge in the tree (and assuming this variable is already named "tree") you could use, \code{BranchPartitionsToTest = lapply(as.list(1:nrow(tree$edge)), as.list)}.
#'
#' Because of the flexibility of this function the user can define any set of edges as well. For example, they could test whether terminal branches have a different rate from internal branches with \code{BranchPartitionsToTest = list(list(match(1:Ntip(tree), tree$edge[, 2])))}. The \code{CladePartitionsToTest} is really just a special subset of this type of hypothesis, but with edges being defined as descending from a specific internal node in the tree. Once again, an explorartory approach like that of Lloyd et al. (2012) can be used with: \code{CladePartitionsToTest = lapply(as.list(Ntip(tree) + (2:Nnode(tree))), as.list)}. Note that this excludes the root node as this would define a single partition and hence would represent the null hypothesis (a single rate model for the whole tree). More generally clades must be defined by the node numbers they correspond to. In R an easy way to identify these is with: \code{plot(tree); nodelabels()}.
#'
#' Time bin paritions are defined in a similar way, but are numbered 1:N beginning with the oldest time bin. So if wanting to do an exploratory test of single bin partitions (and only four time bins were specified) you could use: \code{TimeBinPartitionsToTest = lapply(as.list(1:4), as.list)}. Bins can be combined too, just as edges are above. For example, time bins 1 and 2 could form a single partition with: \code{TimeBinPartitionsToTest = list(list(1:2))}. Or if looking to test a model where each bin has its' own rate value you could use: \code{TimeBinPartitionsToTest = list(as.list(1:3))}. Note, as before we do not need to specify the fourth bin as this will be automatically done by the function, however, \code{TimeBinPartitionsToTest = list(as.list(1:4))} will also work. Some caution needs to be applied with N-rate models (where N is three or larger) as a result favouring such models does not necessarily endorse N-separate rates. I.e., it could simply be that one bin has such a larger excursion that overall the N-rate model fits better than the 1-rate model, but some bins could be combined. It is up to the user to check this themselves by exploring smaller combinations of bins. For example, if the four rate model was considered significant the user could explore all three- and two-rate models as well to check that four-rates are really the optimal explanantion.
#'
#' Finally, character partitions allow the user to explore whether rates vary across different character types, e.g., skeletal characters versus soft tissue characters, or cranial characters versus postcranial characters. Here characters are simply numbered 1:N, but here single character partitions are less likely to be of interest. As an example of use lets say the first ten characters are what we are interested in as a partition (the second partition being the remaining characters), we could use: \code{CharacterPartitionsToTest = list(list(1:10))}.
#'
#' Note that the list structure is critical to defining partitions as it allows partitions to be of different sizes and numbers. For example, one partition of three and another of six, or one set of two partitions and another set of four partitions. However, it may not be intuitive to some users so the examples it is recommended that the user refers to the examples above as a guide.
#'
#' Additionally it should be noted that the user can test multiple types of hypotheses simultaneously with the function. For example, performing several branch tests whilst performing clade tests. However, they needn't perform all types simultaneously (and unused partition types can be set to NULL, the default in each case).
#'
#' \bold{Other options}
#'
#' Since Claddis version 0.3 this function ahs allowed the user to set many more options than were offered previously and these should be considered carefully before running any tests.
#'
#' Firstly, the user can pick an option for \code{ChangeTimes} which sets the times character changes are inferred to occur. This is only relevant when the user is performing time bin partition test(s) as this requires some inference to be made about when changes occur on branches that may span multiple time bins. The current options are: \code{"midpoint"} (all changes are inferred to occur midway along the branch, effectively mimicking the approach of Ruta et al. 2006), \code{"spaced"} (all changes are inferred to occur equally spaced a long the branch, with changes occurring in character number order), or \code{"random"} (changes are assigned a random time by drawing from a uniform distribution between the beginning and end of each branch). The first of these is likely to lead to unrealistically "clumped" changes and by extension implies completely correlated character change that would violate the assumptions of the Poisson distribution that underlies the significance tests here (Lloyd et al. 2012). At the other extreme, the equally spaced option will likely unrealistically smooth out time series and potnentially make it harder to reject the single-rate null. For these reasons, the random option is recommended and is set as the default. However, because it is random this makes the function stochastic (the answer can vary each time it is run) and so the user should therefore run the function multiple times if using this option (i.e., by using a for loop) and aggregating the results at the end (e.g., as was done by previous authors; Lloyd et al. 2012; Close et al. 2015).
#'
#' Secondly, the \code{Alpha} value sets the significance threshold by which the likelihood ratio test's resulting p-value is compared. Following lloyd et al. (2012) this is set lower (0.01) than the standard 0.05 value by default as those authors found rates to be highly heterogenous in their data set (fossil lungfish). However, this should not be adopted as a "standard" value without question. Note that the function also corrects for multiple comparisons to avoid Type I errors (false positives). It does so (following Lloyd et al. 2012) using the Benjamini-Hochberg (Benjamini and Hochberg 1995) False Discovery Rate approach (see Lloyd et al. 2012 for a discussion of why).
#'
#' Thirdly, polymorphisms and uncertainities create complications for assessing character changes along branches. These can occur at the tips (true polymorphisms or uncertainties in sampled taxa) and internal nodes (uncertainty over the estimated ancestral state). There are two options presented here, and applicable to both \code{PolymorphismState} and \code{UncertaintyState} (allowing these to be set separately). These are to convert such values to missing (NA) or to pick one of the possible states at random. Using missing values will increase overall uncertainty and potentially lead to Type II errors (false negatives), but represents a conservative solution. The random option is an attempt to avoid Type II errors, but can be considered unrealistic, especially if there are true polymorphisms. Additionally, the random option will again make the function stochastic meaning the user should run it multiple times amd aggregate the results. Note that if there are no polymorphisms or uncertianties in the character-taxon matrix the latter can still occur with ancestral state estimates, espcially if the threshold value is set at a high value (see \link{AncStateEstMatrix} for details).
#'
#' Fourthly, inapplicable characters can additionally complicate matters as they are not quite the same as missing data. I.e., they can mean that change in a particular character is not even possible along a branch. However, there is no easy way to deal with such characters at present so the user is not presented with a true option here - currently all inapplicable states are simply converted to missing values by the function. In future, though, other options may be available here. For now it is simply noted that users should be careful in making inferences if there are inapplcaible characters in their data and should perhaps consider removing them with \link{MatrixPruner} to gauge their effect.
#'
#' Fifthly, there are currenty two further options for assessing rates across time bins. As noted above a complication here is that character changes (the rate numerator) and character completeness (part of the rate denominator) are typically assessed on branches. However, branches will typically span time bin boundaries and hence many bins will contain only some portion of particular branches. The exact portion can be easily calculated for branch durations (the other part of the rate denominator) and the \code{ChangeTimes} option above is used to set the rate numerator, however, the completeness remains complex to deal with. The first attempt to deal with this was made by Close et al. (2015) who simply did w eighted mean completeness by using the proportion of a branch in each bin as the weight and multiplying this by each branches completeness (the \code{"Close"} option here). However, this may lead to unrealistic "smoothing" of the data and perhaps more importantly makes no account of which characters are known in a bin. Lloyd (2016) proposed an alternative "subtree" approach which assesses completeness by considering each character to be represented by a subtree where only branches that are complete are retained then branch durations in each bin are summed across subtrees such that the duration term automatically includes completenes (the \code{"Lloyd"} option here). Here the latter is recommended, but note that no true comparison has yet been made in the literature.
#'
#' Sixthly, all character changes are weighted according to the weights provided in the input character-taxon matrix. In many cases these will simply all be one, although see the equalise weights option in \link{ReadMorphNexus}. However, when weights vary they can create some issues for the function. Specifcally, changes are expected to be in the same (integer) units, but if weights vary then they have to be modelled accordingly. I.e., a character twice the weight of another may lead to a single change being counted as two changes. This is most problematic when the user has continuous characters which are automatically converted to gap-weighted (Thiele 1993) characters. However, this conversion creates drastically down-weighted characters nd hence the user may wish to set the \code{EnsureAllWeightsAreIntegers} option to TRUE. Note that reweighting will affect the results and hence shifting the weights of characters up or down will necessarily lead to shifts in the relative Type I and II errors. This is an unexplored aspect of such approaches, but is something the user should be aware of. More broadly it is recommended that continuous (or gap-weighted) characters be avoided when using these approaches.
#'
#' Finally, the ramaining options (\code{EstimateAllNodes}, \code{EstimateTipValues}, \code{InapplicablesAsMissing}, \code{PolymorphismBehaviour}, \code{UncertaintyBehaviour}, and \code{Threshold}) are all simply passed directly to \link{AncStateEstMatrix} for estimating the ancestral states and users should consult the help file for that function for further details.
#'
#' Note that currently the function cannot deal with step matrices.
#'
#' Note that the terminal versus internal option from Brusatte et al. (2014) is yet to be implemented by default.
#'
#' @description
#'
#' Given a tree and a cladistic-type matrix uses likelihood ratio tests to compare N-rate and 1-rate models across branches, clades, time bins, or character partitions.
#'
#' @param tree A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{clad.matrix}.
#' @param clad.matrix A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @param TimeBins A vector of ages (in millions of years) indicating the boundaries of a series of time bins in order from oldest to youngest.
#' @param BranchPartitionsToTest A list of branch(es) (edge numbers) to test for a 2-rate parameter model (i.e., one rate for the edge and another for the rest of the tree). If NULL (the default) then no partition test(s) will be made.
#' @param CharacterPartitionsToTest A list of character partition(s) (character numbers) to test for a 2-rate parameter model (i.e., one rate for the partition and another for the remaining characters). If NULL (the default) then no partition test(s) will be made.
#' @param CladePartitionsToTest A list of clade partition(s) (node numbers) to test for a 2-rate parameter model (i.e., one rate for the clade and another for the rest of the tree). If NULL (the default) then no partition test(s) will be made.
#' @param TimeBinPartitionsToTest A list of time bin partition(s) (numbered 1 to N) to test for a 2-rate parameter model (i.e., one rate for the time bin(s) and another for the remaining time bins). If NULL (the default) then no partition test(s) will be made.
#' @param ChangeTimes The time at which to record the character changes. One of \code{"midpoint"} (changes occur at the midpoint of the branch), \code{"spaced"} (changes equally spaced along branch), or \code{"random"} (change times drawn at random from a uniform distribution; the default and recommended option). Note: this is only meaningful if testing for time bin partitions.
#' @param Alpha The alpha value to be used for the significance tests. The default is 0.01.
#' @param PolymorphismState One of \code{"missing"} (converts polymorphic values to NA; the default) or \code{"random"} (picks one of the possible polymorphic states at random).
#' @param UncertaintyState One of \code{"missing"} (converts uncertain values to NA; the default) or \code{"random"} (picks one of the possible uncertain states at random).
#' @param InapplicableState The only current option is \code{"missing"} (converts value to NA).
#' @param TimeBinApproach One of \code{"Close"} or \code{"Lloyd"} (the default).
#' @param EnsureAllWeightsAreIntegers Logical for whether (\code{TRUE}) to reweight non-integer weights until all weights are integers or to leave them as they are (\code{FALSE}; the default).
#' @param EstimateAllNodes Option passed to internal use of \link{AncStateEstMatrix}.
#' @param EstimateTipValues Option passed to internal use of \link{AncStateEstMatrix}.
#' @param InapplicablesAsMissing Option passed to internal use of \link{AncStateEstMatrix}.
#' @param PolymorphismBehaviour Option passed to internal use of \link{AncStateEstMatrix}.
#' @param UncertaintyBehaviour Option passed to internal use of \link{AncStateEstMatrix}.
#' @param Threshold Option passed to internal use of \link{AncStateEstMatrix}.
#'
#' @return
#'
#' \item{character.changes}{A matrix containing information on all the character changes reconstructed and used to measure rates.}
#' \item{node.results}{A table displaying the results of the per-clade rate tests.}
#' \item{branch.results}{A table displaying the results of the per-branch rate tests.}
#' \item{per.bin.rates}{Per time-bin rates (use with caution).}
#' \item{per.bin.rates.tb}{Per time-bin rates for terminal branches (use with caution).}
#' \item{per.bin.rates.ib}{Per time-bin rates for internal branches (use with caution).}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Steve C. Wang \email{scwang@@swarthmore.edu}
#'
#' @references
#'
#' Benjamini, Y. and Hochberg, Y., 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society, Series B, 57, 289-300.
#'
#' Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. Current Biology, 24, 2386-2392.
#'
#' Close, R. A., Friedman, M., Lloyd, G. T. and Benson, R. B. J., 2015. Evidence for a mid-Jurassic adaptive radiation in mammals. Current Biology, 25, 2137-2142.
#'
#' Cloutier, R., 1991. Patterns, trends, and rates of evolution within the Actinistia. Environmental Biology of Fishes, 32, 23–58.
#'
#' Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. Biological Journal of the Linnean Society, 118, 131-151.
#'
#' Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying heterogeneity in rates of morphological evolution: discrete character change in the evolution of lungfish (Sarcopterygii; Dipnoi). Evolution, 66, 330-348.
#'
#' Ruta, M., Wagner, P. J. and Coates, M. I., 2006. Evolutionary patterns in early tetrapods. I. Rapid initial diversification followed by decrease in rates of character change. Proceedinsg of the Royal Society of London B, 273, 2107–2111.
#'
#' Thiele, K.. 1993. The Holy Grail of the perfect character: the cladistic treatment of morphometric data. Cladistics, 9, 275-304.
#'
#' @keywords evolution,rates
#'
#' @examples
#' 
#' # Set random seed:
#' set.seed(17)
#' 
#' # Generate a random tree for the Michaux data set:
#' tree <- rtree(nrow(Michaux1989$Matrix_1$Matrix))
#' 
#' # Update taxon names to match those in the data matrix:
#' tree$tip.label <- rownames(Michaux1989$Matrix_1$Matrix)
#' 
#' # Set root time by making youngest taxon extant:
#' tree$root.time <- max(diag(vcv(tree)))
#' 
#' # Get discrete character rates:
#' x <- DiscreteCharacterRate(tree = tree, clad.matrix =
#'   Michaux1989, TimeBins = seq(from = tree$root.time,
#'   to = 0, length.out = 5), BranchPartitionsToTest =
#'   lapply(as.list(1:nrow(tree$edge)), as.list),
#'   CharacterPartitionsToTest = lapply(as.list(1:3),
#'   as.list), CladePartitionsToTest =
#'   lapply(as.list(Ntip(tree) + (2:Nnode(tree))),
#'   as.list), TimeBinPartitionsToTest =
#'   lapply(as.list(1:4), as.list), ChangeTimes =
#'   "random", Alpha = 0.01, PolymorphismState =
#'   "missing", UncertaintyState = "missing",
#'   InapplicableState = "missing", TimeBinApproach =
#'   "Lloyd", EnsureAllWeightsAreIntegers = FALSE,
#'   EstimateAllNodes = FALSE, EstimateTipValues =
#'   FALSE, InapplicablesAsMissing = FALSE,
#'   PolymorphismBehaviour = "equalp",
#'   UncertaintyBehaviour = "equalp", Threshold = 0.01)
#'
#' @export DiscreteCharacterRate

# OPTION FOR MULTIPLE COMPARISON CORRECTION APPROACH TO APPLY (ALTHOUGH PROLLY ONLY USE BH)
# LIST OF LISTS FOR N-RATE PARAMETER TESTS? YEP. AND THEN CHECK THESE SO CAN SET "OTHER" RATE SECTIONS. NEEDS CHECK THAT DEFINED PARTITIONS DO NOT CONTAIN
# WRITE SEARCH VERSION FOR FINDING RATE SHIFTS? SHOULD THIS EVEN BE AN OPTION?
# MAYBE MAKE ANCESTRAL STATE UNCERTIANTY DIFFERENT FOR TIPS THAN NODES?
# CHANGE TIMES CANNOT BE COMPLETELY RANDOM AS MULTISTEP (I.E. ORDERED MULTISTEP) CHANHES MUST BE IN A SEQUENCE. ALTHOUGH CAN BE MODELLED AS SUCH?
# TEST FOR OVERLAP BETWEEN CLADES IS MORE COMPLEX, MAYBE CONVERT THESE TO BRANCH PARTITIONS IN FUNCTION INSTEAD AS THIS CAN THEN USE EXISTING CODE.
# ADD TERMINAL VERSUS INTERNAL OPTION SOMEHOW/SOMEWHERE.
# ALLOW OPTION TO IGNORE SOME PARTS OF THE TREE FOR PARTITON TESTS? MAKES CALCULATING THE MAN RATE TRICKIER BUT MIGHT MKAE SENSE E.G. FOR INGROUP ONLY TESTS.
# THE MEAN RATE SHOULD BE THE MEAN RATE FOR ANY PARTITION TYPE? THEN THE EXPECTED NUMBER OF CHANGES IS SIMPLY A PRODUCT OF TIME AND COMPLETENESS?
# ADD MEAN RATE TO OUTPUT?
# NEED EXTRA FUNCTION(S) TO VISUALISE RESULTS MOST LIKELY
# ALLOW REWEIGHTING OF INAPPLICABLES ZERO AS AN OPTION FOR THEM?

DiscreteCharacterRate <- function(tree, clad.matrix, TimeBins, BranchPartitionsToTest = NULL, CharacterPartitionsToTest = NULL, CladePartitionsToTest = NULL, TimeBinPartitionsToTest = NULL, ChangeTimes = "random", Alpha = 0.01, PolymorphismState = "missing", UncertaintyState = "missing", InapplicableState = "missing", TimeBinApproach = "Lloyd", EnsureAllWeightsAreIntegers = FALSE, EstimateAllNodes = FALSE, EstimateTipValues = FALSE, InapplicablesAsMissing = FALSE, PolymorphismBehaviour = "equalp", UncertaintyBehaviour = "equalp", Threshold = 0.01) {
  
  # Check for step matrices and stop and warn user if found:
  if(is.list(morph.matrix$Topper$StepMatrices)) stop("Function cannot currently deal with step matrices.")

  # Check tree has branch lengths:
  if(is.null(tree$edge.length)) stop("Tree does not have branch lengths (durations). Try timescaling the tree, e.g., with DatePhylo.")
  
  # Check tree has root age:
  if(is.null(tree$root.time)) stop("Tree is missing $root.time. Try setting this before continuing, e.g., tree$root.time <- 104.2.")
  
  # Check ChangeTimes is correctly formatted or stop and warn user:
  if(length(setdiff(ChangeTimes, c("midpoint", "spaced", "random"))) > 0) stop("ChangeTimes must be one of \"midpoint\", \"spaced\", or \"random\".")
  
  # Check PolymorphismState is correctly formatted or stop and warn user:
  if(length(setdiff(PolymorphismState, c("missing", "random"))) > 0) stop("PolymorphismState must be one of \"missing\" or \"random\".")
  
  # Check UncertaintyState is correctly formatted or stop and warn user:
  if(length(setdiff(UncertaintyState, c("missing", "random"))) > 0) stop("UncertaintyState must be one of \"missing\" or \"random\".")
  
  # Check InapplicableState is correctly formatted or stop and warn user:
  if(length(setdiff(InapplicableState, c("missing"))) > 0) stop("InapplicableState must be \"missing\".")
  
  # Check TimeBinApproach is correctly formatted or stop and warn user:
  if(length(setdiff(TimeBinApproach, c("Close", "Lloyd"))) > 0) stop("TimeBinApproach must be one of \"Close\" or \"Lloyd\".")
  
  # Check partitions are not all NULL values:
  if(is.null(BranchPartitionsToTest) && is.null(CharacterPartitionsToTest) && is.null(CladePartitionsToTest) && is.null(TimeBinPartitionsToTest)) stop("No partitions are requested. Set at least one of BranchPartitionsToTest, CharacterPartitionsToTest, CladePartitionsToTest, or TimeBinPartitionsToTest to a list of appropriate values. Type \"?DiscreteCharacterRate\" for help.")
  
  # Get internal node numbers:
  InternalNodeNumbers <- 1:ape::Nnode(tree) + ape::Ntip(tree)
  
  # Get edge numbers:
  EdgeNumbers <- 1:nrow(tree$edge)
  
  # Get character numbers:
  CharacterNumbers <- 1:sum(unlist(lapply(lapply(clad.matrix[2:length(clad.matrix)], '[[', "Matrix"), ncol)))
  
  # Ensure time bins are in correct order:
  TimeBins <- sort(unique(TimeBins), decreasing = TRUE)
  
  # Find the Time bin midpoints:
  TimeBinMidpoints <- (TimeBins[2:length(TimeBins)] + TimeBins[1:(length(TimeBins) - 1)]) / 2
  
  # Get the numbers for each time bins:
  TimeBinNumbers <- 1:length(TimeBinMidpoints)
  
  # Subfunction to ensure partitions are formatted correctly:
  PartitionFormatter <- function(PartitionsToTest, ValidValues, PartitionName) {
    
    # Check partitions are in the form of a list of lists:
    if(!all(c(all(unlist(lapply(PartitionsToTest, is.list))), is.list(PartitionsToTest)))) stop(paste(PartitionName, " must be in the form of a list of lists.", sep = ""))
    
    # Get a vector of any non-present valid values:
    NonPresentValues <- setdiff(unique(unlist(PartitionsToTest)), ValidValues)
    
    # Check valid values have been used and if not stop and warn user:
    if(length(NonPresentValues) > 0) stop(paste(PartitionName, "Partitions to test must be defined using the valid range of values (", paste(range(ValidValues), collapse = " to "), ") only.", sep = ""))
    
    # Check partitions never overlap and stop and warn user if they do:
    Check <- lapply(PartitionsToTest, function(x) if(any(duplicated(sort(unlist(x))))) stop(paste("Each partition of ", PartitionName, " must not contain overlapping values (e.g., can not have 1:3 and 3:5 as both contain 3).", sep = "")))
    
    # Subfunction to ad the missing partition (if exists):
    AddMissingPartitions <- function(x, ValidValues) {
      
      # Define any missing values:
      MissingValues <- setdiff(ValidValues, unlist(x))
      
      # If there are missing values add them to list at end:
      if(length(MissingValues) > 0) x[[(length(x) + 1)]] <- MissingValues
      
      # Return x:
      return(x)
      
    }
    
    # Add in missing partitions (if any):
    PartitionsToTest <- lapply(PartitionsToTest, AddMissingPartitions, ValidValues = ValidValues)
    
    # Check partitions are all at least two in size or else no comparison can be made:
    if(any(unlist(lapply(PartitionsToTest, length)) == 1)) stop("Partitions must divide the available data into at least two parts.")

    # Return formatted partitions to test:
    return(PartitionsToTest)

  }

  # If performing branch partition test(s) check and reformat branch partitions:
  if(!is.null(BranchPartitionsToTest)) BranchPartitionsToTest <- PartitionFormatter(PartitionsToTest = BranchPartitionsToTest, ValidValues = EdgeNumbers, PartitionName = "BranchPartitionsToTest")
  
  # If performing character partition test(s) check and reformat character partitions:
  if(!is.null(CharacterPartitionsToTest)) CharacterPartitionsToTest <- PartitionFormatter(PartitionsToTest = CharacterPartitionsToTest, ValidValues = CharacterNumbers, PartitionName = "CharacterPartitionsToTest")
  
  # If performing clade partition test(s)
  if(!is.null(CladePartitionsToTest)) {
    
    # Convert clade partitions to edge partitions:
    CladePartitionsToTest <- lapply(CladePartitionsToTest, lapply, GetDescendantEdges, tree = tree)
    
    # Check and reformat clade partitions:
    CladePartitionsToTest <- PartitionFormatter(PartitionsToTest = CladePartitionsToTest, ValidValues = EdgeNumbers, PartitionName = "CladePartitionsToTest")
    
  }

  # If performing time bin partition test(s) check and reformat time bin partitions:
  if(!is.null(TimeBinPartitionsToTest)) TimeBinPartitionsToTest <- PartitionFormatter(PartitionsToTest = TimeBinPartitionsToTest, ValidValues = TimeBinNumbers, PartitionName = "TimeBinPartitionsToTest")
  
  # Subfunction to calculate maximum likelihood p value:
  GetMaximumLikelihoodPValue <- function(MeanRate, SampledRates, SampledChanges, SampledCompleteness, SampledTime) {
    
    # Set maximum likelihood numerator:
    MaximumLikelihoodNumerator <- MeanRate
    
    # Set maximum likelihood denominator:
    MaximumLikelihoodDenominator <- SampledRates
    
    # Get log numerator:
    LogNumerator <- sum(log(dpois(round(SampledChanges), MaximumLikelihoodNumerator * SampledCompleteness * SampledTime)))
    
    # Get log denominator:
    LogDenominator <- sum(log(dpois(round(SampledChanges), MaximumLikelihoodDenominator * SampledCompleteness * SampledTime)))
    
    # Get test statistic:
    TestStatistic <- -2 * (LogNumerator - LogDenominator)
    
    # Calculate position of test statistic in chi-square distribution to get probability:
    PValue <- pchisq(TestStatistic, length(SampledRates) - 1, lower.tail = FALSE)
    
    # Output probability for later alpha comparison:
    return(PValue)
    
  }
  
  # Subfunction to perform branch partition tests:
  PerformBranchPartitionTests <- function(x) {
    
    # Subfunction to get edge values for likelihood test(s):
    GetEdgeValuesForLikelihood <- function(x) {
      
      # Isolate edge changes:
      EdgeChanges <- sum(unlist(lapply(EdgeList[x], function(x) sum(x$CharacterChanges[, "Weight"]))))
      
      # Isolate edge completeness (includes time term):
      EdgeCompleteness <- sum(unlist(lapply(EdgeList[x], function(x) sum(Weights[x$ComparableCharacters] / sum(Weights)))) * unlist(lapply(EdgeList[x], function(x) x$BranchDuration)))
      
      # Set edge duration as one as already included in completeness term:
      EdgeDuration <- 1
      
      # Set edge rate:
      EdgeRate <- EdgeChanges / (EdgeCompleteness * EdgeDuration)
      
      # COmbine values into vector for output:
      output <- c(EdgeChanges, EdgeCompleteness, EdgeDuration, EdgeRate)
      
      # Add naems to vector:
      names(output) <- c("EdgeChanges", "EdgeCompleteness", "EdgeDuration", "EdgeRate")
      
      # Return values as output:
      return(output)
      
    }
    
    # Get edge values for each partition:
    EdgeTestValues <- do.call(rbind, lapply(x, GetEdgeValuesForLikelihood))
    
    # Get p value for likelihood ratio test:
    TestPValue <- GetMaximumLikelihoodPValue(MeanRate = MeanEdgeRate, SampledRates = EdgeTestValues[, "EdgeRate"], SampledChanges = EdgeTestValues[, "EdgeChanges"], SampledCompleteness = EdgeTestValues[, "EdgeCompleteness"], SampledTime = EdgeTestValues[, "EdgeDuration"])
    
    # Combine output into a list:
    Output <- list(EdgeTestValues[, "EdgeRate"], TestPValue)
    
    # Add names to output:
    names(Output) <- c("BranchRates", "PValue")
    
    # Return output:
    return(Output)
    
  }
  
  # Subfunction to perform time bin partition tests:
  PerformTimeBinPartitionTests <- function(x) {
    
    # Subfunction to get time bin values for likelihood test(s):
    GetTimeBinValuesForLikelihood <- function(x) {
      
      # Isolate time bin changes:
      TimeBinChanges <- sum(unlist(lapply(WeightedChangesByBin[x], sum)))
      
      # If using Close approach isolate time bin completeness (includes time term):
      if(TimeBinApproach == "Close") TimeBinCompleteness <- sum(unlist(lapply(ProportionalCompletenessByBin[x], sum)))
      
      # If using Lloyd approach isolate time bin completeness (includes time term):
      if(TimeBinApproach == "Lloyd") TimeBinCompleteness <- sum(unlist(lapply(DurationCompletenessByBin[x], sum)))
      
      # Set time bin duration as one as already included in completeness term:
      TimeBinDuration <- 1
      
      # Set time bin rate:
      TimeBinRate <- TimeBinChanges / (TimeBinCompleteness * TimeBinDuration)
      
      # COmbine values into vector for output:
      output <- c(TimeBinChanges, TimeBinCompleteness, TimeBinDuration, TimeBinRate)
      
      # Add names to vector:
      names(output) <- c("TimeBinChanges", "TimeBinCompleteness", "TimeBinDuration", "TimeBinRate")
      
      # Return values as output:
      return(output)
      
    }
    
    # Get edge values for each partition:
    TimeBinTestValues <- do.call(rbind, lapply(x, GetTimeBinValuesForLikelihood))
    
    # Get p value for likelihood ratio test:
    TestPValue <- GetMaximumLikelihoodPValue(MeanRate = MeanTimeBinRate, SampledRates = TimeBinTestValues[, "TimeBinRate"], SampledChanges = TimeBinTestValues[, "TimeBinChanges"], SampledCompleteness = TimeBinTestValues[, "TimeBinCompleteness"], SampledTime = TimeBinTestValues[, "TimeBinDuration"])
    
    # Combine output into a list:
    Output <- list(TimeBinTestValues[, "TimeBinRate"], TestPValue)
    
    # Add names to output:
    names(Output) <- c("TimeBinRates", "PValue")
    
    # Return output:
    return(Output)
    
  }

  # Get ages for each (tip and internal) node:
  NodeAges <- GetNodeAges(tree)

  # Get branch ages (from and to):
  BranchAges <- unname(cbind(NodeAges[as.character(tree$edge[, 1])], NodeAges[as.character(tree$edge[, 2])]))

  # Build edge list from node numbers (from-to) for each branch:
  EdgeList <- lapply(apply(tree$edge, 1, list), function(x) {names(x) <- "NodeNumberFromTo"; return(x)})
  
  # Add node ages to edge list:
  for(i in 1:length(EdgeList)) EdgeList[[i]]$NodeAgeFromTo <- BranchAges[i, ]
  
  # Add node ages (from-to) to each edge in list:
  EdgeList <- lapply(EdgeList, function(x) {x$BranchDuration <- x$NodeAgeFromTo[1] - x$NodeAgeFromTo[2]; return(x)})
  
  # Get vector of branch types:
  BranchTypes <- gsub("0", "Internal", gsub("1", "Terminal", as.numeric(tree$edge[, 2] <= ape::Ntip(tree))))
  
  # Add branch type to edge list:
  for(i in 1:length(EdgeList)) EdgeList[[i]]$BranchType <- BranchTypes[i]
  
  # Find descendant edges for each internal node:
  DescendantEdgesForEachInternalNode <- lapply(as.list(InternalNodeNumbers), GetDescendantEdges, tree = tree)
  
  # Use these to assign clade memberships for each branch:
  for(i in 1:length(EdgeList)) EdgeList[[i]]$CladeMembership <- InternalNodeNumbers[unlist(lapply(lapply(DescendantEdgesForEachInternalNode, '==', i), sum)) == 1]
  
  # Add opposite (clade opposition) to each branch:
  EdgeList <- lapply(EdgeList, function(x) {x$CladeOpposition <- setdiff(InternalNodeNumbers, x$CladeMembership); return(x)})
  
  # Get ancestral character states:
  AncestralStates <- AncStateEstMatrix(InputMatrix = clad.matrix, Tree = tree, EstimateAllNodes = EstimateAllNodes, EstimateTipValues = EstimateTipValues, InapplicablesAsMissing = InapplicablesAsMissing, PolymorphismBehaviour = PolymorphismBehaviour, UncertaintyBehaviour = UncertaintyBehaviour, Threshold = Threshold)
  
  # Build single matrix of all states in tip label then node number order:
  AllStates <- do.call(cbind, lapply(lapply(AncestralStates[2:length(AncestralStates)], '[[', "Matrix"), function(x) x[c(tree$tip.label, 1:ape::Nnode(tree) + ape::Ntip(tree)), , drop = FALSE]))
  
  # Make vector of ordering of characters:
  Ordering <- unname(unlist(lapply(clad.matrix[2:length(clad.matrix)], '[[', "Ordering")))
  
  # Make vector of weights of characters:
  Weights <- unname(unlist(lapply(clad.matrix[2:length(clad.matrix)], '[[', "Weights")))
  
  # Make vector of minimum values:
  MinVals <- unname(unlist(lapply(clad.matrix[2:length(clad.matrix)], '[[', "MinVals")))
  
  # Make vector of maximum values:
  MaxVals <- unname(unlist(lapply(clad.matrix[2:length(clad.matrix)], '[[', "MaxVals")))
  
  # Find positions in matrix with polymorphisms:
  PolymorphismPositions <- grep("&", AllStates)
  
  # Find positions in matrix with uncertainties:
  UncertaintyPositions <- grep("/", AllStates)
  
  # Find positions in matrix with inapplicables:
  InapplicablePositions <- which(AllStates == "")

  # If polymorphisms were found:
  if(length(PolymorphismPositions) > 0) {
    
    # If replacing polymorphsims with missing do so:
    if(PolymorphismState == "missing") AllStates[PolymorphismPositions] <- NA
    
    # If replacing polymorphisms with random values draw and replace:
    if(PolymorphismState == "random") AllStates[PolymorphismPositions] <- unlist(lapply(strsplit(AllStates[PolymorphismPositions], "&"), sample, size = 1))
    
  }
  
  # If uncertainties were found:
  if(length(UncertaintyPositions) > 0) {
    
    # If replacing uncertainties with missing do so:
    if(UncertaintyState == "missing") AllStates[UncertaintyPositions] <- NA
    
    # If replacing uncertainties with random values draw and replace:
    if(UncertaintyState == "random") AllStates[UncertaintyPositions] <- unlist(lapply(strsplit(AllStates[UncertaintyPositions], "&"), sample, size = 1))
    
  }
  
  # If inapplicable states were found:
  if(length(InapplicablePositions) > 0) {
    
    # If replacing inapplicables with missing do so:
    if(InapplicableState == "missing") AllStates[InapplicablePositions] <- NA
    
  }
  
  # Check for continuous characters as these will need to be modified for modelling to discrete characters:
  if(any(Ordering == "cont")) {
    
    # Tell user this is happening:
    cat("Continuous characters found. Converting to gap-weighted discrete characters.\n")
    
    # Find out which characters are continuous:
    ContinuousCharactersFound <- which(Ordering == "cont")
    
    # Rescale continous characters as zero to one values:
    ListOfContinuousValuesRescaledZeroToOne <- lapply(lapply(lapply(apply(AllStates[, ContinuousCharactersFound, drop = FALSE], 2, list), unlist), as.numeric), function(x) {x <- x - min(sort(x)); x <- x / max(sort(x)); return(x)})
    
    # Now discretize and store these characters (0 to 31 scale):
    AllStates[, ContinuousCharactersFound] <- do.call(cbind, lapply(lapply(lapply(ListOfContinuousValuesRescaledZeroToOne, function(x) as.list(x)), lapply, function(x) ifelse(is.na(x), NA, max(which(x >= (0:31) / 31)) - 1)), unlist))
    
    # Convert character type to ordered:
    Ordering[ContinuousCharactersFound] <- "ord"
    
    # Convert weights to 1/31:
    Weights[ContinuousCharactersFound] <- 1 / 31
    
    # Set minimum value to zero:
    MinVals[ContinuousCharactersFound] <- 0
    
    # Set maximum value to 31:
    MaxVals[ContinuousCharactersFound] <- 31
    
  }
  
  # If EnsureAllWeightsAreIntegers is TRUE rescale weights until they are all integers so can model appropriately with Poisson later:
  if(EnsureAllWeightsAreIntegers) while(is.character(all.equal(sum(Weights %% 1), 0))) Weights <- (1 / (Weights %% 1)[(Weights %% 1) > 0])[1] * Weights

  # Add from-to node states for each character to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$CharacterStatesFromTo <- matrix(AllStates[x$NodeNumberFromTo, , drop = FALSE], nrow = 2, dimnames = list(c("From", "To"))); return(x)})
  
  # Subfunction to define character changes:
  BuildChangesMatrix <- function(x) {
    
    # Find only comparable characters (those scored for both from and to states):
    ComparableCharacters <- which(apply(!apply(x$CharacterStatesFromTo, 2, is.na), 2, all))
    
    # Isolate comparable ordering:
    ComparableOrdering <- Ordering[ComparableCharacters]
    
    # Isolate comparable weights:
    ComparableWeights <- Weights[ComparableCharacters]
    
    # Isolate only characters that actually differ (change):
    CharacterDifferences <- which(x$CharacterStatesFromTo["From", ComparableCharacters] != x$CharacterStatesFromTo["To", ComparableCharacters])
    
    # Build character change matrix:
    CharacterChanges <- matrix(nrow = 0, ncol = 5, dimnames = list(c(), c("Character", "From", "To", "Steps", "Weight")))
    
    # If characters change then make a matrix from them:
    if(length(CharacterDifferences) > 0) CharacterChanges <- rbind(CharacterChanges, cbind(as.numeric(ComparableCharacters[CharacterDifferences]), as.numeric(x$CharacterStatesFromTo["From", ComparableCharacters[CharacterDifferences]]), as.numeric(x$CharacterStatesFromTo["To", ComparableCharacters[CharacterDifferences]]), ifelse(ComparableOrdering[CharacterDifferences] == "unord", 1, abs(as.numeric(x$CharacterStatesFromTo["To", ComparableCharacters[CharacterDifferences]]) - as.numeric(x$CharacterStatesFromTo["From", ComparableCharacters[CharacterDifferences]]))), ComparableWeights[CharacterDifferences]))
    
    # Store character changes as new sublist for x:
    x$CharacterChanges <- CharacterChanges
    
    # Store comparable characters as new sublist of x:
    x$ComparableCharacters <- ComparableCharacters
    
    # Return x:
    return(x)
    
  }
  
  # Get character changes and comparable characters and add to edge list:
  EdgeList <- lapply(EdgeList, BuildChangesMatrix)
  
  # Check whether time bins are being compared (otherwise no need to assign character changes):
  if(!is.null(TimeBinPartitionsToTest)) {
    
    # Subfunction to add change times to character changes:
    AddChangeTimes <- function(x, ChangeTimes) {
      
      # Isolate character changes:
      CharacterChanges <- x$CharacterChanges
      
      # If any changes involve two or more steps (requiring replacement with multiple changes):
      if(any(CharacterChanges[, "Steps"] > 1)) {
        
        # Get multistep character changes:
        MultiStepCharacters <- which(CharacterChanges[, "Steps"] > 1)
        
        # For each multistep character change:
        for(i in rev(MultiStepCharacters)) {
          
          # Isolate other rows:
          OtherRowNumbers <- setdiff(1:nrow(CharacterChanges), i)
          
          # Get unpacked changes (X:Y, e.g., 0:2 would become 0 1 2):
          UnpackedChanges <- CharacterChanges[i, "From"]:CharacterChanges[i, "To"]
          
          # Update character changes with multistep changes unpacked:
          CharacterChanges <- rbind(CharacterChanges[OtherRowNumbers, ], unname(cbind(rep(CharacterChanges[i, "Character"], length.out = length(UnpackedChanges) - 1), UnpackedChanges[1:(length(UnpackedChanges) - 1)], UnpackedChanges[2:length(UnpackedChanges)], rep(1, length.out = length(UnpackedChanges) - 1), rep(CharacterChanges[i, "Weight"], length.out = length(UnpackedChanges) - 1))))
          
        }
        
        # Resort character changes by character number:
        CharacterChanges <- CharacterChanges[order(CharacterChanges[, "Character"]), ]
        
      }
      
      # If using midpoint option set character change times as midpoint of branch:
      if(ChangeTimes == "midpoint") CharacterChanges <- cbind(CharacterChanges, rep(x$NodeAgeFromTo[1] - (x$BranchDuration / 2), length.out = nrow(CharacterChanges)))
      
      # If using spaced then set character change times as equally spaced along branch:
      if(ChangeTimes == "spaced") CharacterChanges <- cbind(CharacterChanges, x$NodeAgeFromTo[1] - (seq(from = 0, to = x$BranchDuration, length.out = nrow(CharacterChanges) + 1)[1:nrow(CharacterChanges)] + (diff(seq(from = 0, to = x$BranchDuration, length.out = nrow(CharacterChanges) + 1))[1] / 2)))
      
      # If using random then set character change times as random draws from a uniform distribution:
      if(ChangeTimes == "random") CharacterChanges <- cbind(CharacterChanges, x$NodeAgeFromTo[1] - runif(n = nrow(CharacterChanges), min = 0, max = x$BranchDuration))
      
      # Add column name to change time column:
      colnames(CharacterChanges)[ncol(CharacterChanges)] <- "Time"
      
      # Add bin for character change as las column:
      CharacterChanges <- cbind(CharacterChanges, unlist(lapply(as.list(CharacterChanges[, "Time"]), function(x) max(which(x <= TimeBins)))))
      
      # Add column name to change time column:
      colnames(CharacterChanges)[ncol(CharacterChanges)] <- "Bin"
      
      # Overwrite character changes with new version with changes added:
      x$CharacterChanges <- CharacterChanges
      
      # Return x:
      return(x)
    
    }
    
    # Add character change times to edge list:
    EdgeList <- lapply(EdgeList, AddChangeTimes, ChangeTimes = ChangeTimes)
    
  }
  
  # Subfunction to get edge sections in time bins:
  EdgeSectionsInBins <- function(x, TimeBins = TimeBins) {
    
    # Set first appearance datum of edge:
    FAD <- x$NodeAgeFromTo[1]
    
    # Set last appearance datum of edge:
    LAD <- x$NodeAgeFromTo[2]
    
    # Get any time bin boundaries crossed (an be empty if none are):
    BoundariesCrossed <- TimeBins[2:(length(TimeBins) - 1)][intersect(which(TimeBins[2:(length(TimeBins) - 1)] > LAD), which(TimeBins[2:(length(TimeBins) - 1)] < FAD))]
    
    # If boundaries are crossed:
    if(length(BoundariesCrossed) > 0) {
      
      # Break up branch into binned sections as vector of FADs:
      FAD <- c(FAD, BoundariesCrossed)
      
      # Break up branch into binned sections as vector of LADs:
      LAD <- c(BoundariesCrossed, LAD)
      
    }
    
    # Build matrix of branch sections with FADs and LADs:
    BranchSections <- rbind(FAD, LAD)
    
    # Add bin number present in to column names:
    colnames(BranchSections) <- unlist(lapply(lapply(lapply(as.list(BranchSections["FAD", ]), '<=', TimeBins), which), max))
    
    # Add new list section for branch (edge) sections binned by time:
    x$BinnedEdgeSections <- BranchSections
    
    # Return output:
    return(x)
    
  }
  
  # Get edge sections in time bins:
  EdgeList <- lapply(EdgeList, EdgeSectionsInBins, TimeBins = TimeBins)
  
  # Add binned branch durations to edge list:
  EdgeList <- lapply(EdgeList, function(x) {BranchDurations <- rep(0, length(TimeBins) - 1); BranchDurations[as.numeric(colnames(x$BinnedEdgeSections))] <- abs(apply(x$BinnedEdgeSections, 2, diff)); x$BinnedBranchDurations <- BranchDurations; return(x)})
  
  # Add proportional binned branch lengths to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$ProportionalBinnedEdgeDurations <- x$BinnedBranchDurations / sum(x$BinnedBranchDurations); return(x)})
  
  # Add time bin membership to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$TimeBinMembership <- which(x$ProportionalBinnedEdgeDurations > 0); return(x)})
  
  # Add time bin opposition to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$TimeBinOpposition <- setdiff(1:(length(TimeBins) - 1), x$TimeBinMembership); return(x)})
  
  # If performing branch partition tests:
  if(!is.null(BranchPartitionsToTest)) {
    
    # Set mean edge rate:
    MeanEdgeRate <- sum(unlist(lapply(EdgeList, function(x) x$CharacterChanges[, "Weight"]))) / sum(unlist(lapply(EdgeList, function(x) x$BranchDuration)) * (unlist(lapply(EdgeList, function(x) sum(Weights[x$ComparableCharacters]))) / sum(Weights)))
    
    # Perform branch partition tests:
    BranchPartitionTestResults <- lapply(BranchPartitionsToTest, PerformBranchPartitionTests)
    
  # If not performing branch partition tests:
  } else {
    
    # Make empty branch partition result output:
    BranchPartitionTestResults <- NULL
    
  }
  
  
  
  CharacterPartitionTestResults <- NULL #CHARACTER HERE EVENTUALLY!
  
  
  
  
  # If performing clade partition tests:
  if(!is.null(CladePartitionsToTest)) {
    
    # Set mean edge rate:
    MeanEdgeRate <- sum(unlist(lapply(EdgeList, function(x) x$CharacterChanges[, "Weight"]))) / sum(unlist(lapply(EdgeList, function(x) x$BranchDuration)) * (unlist(lapply(EdgeList, function(x) sum(Weights[x$ComparableCharacters]))) / sum(Weights)))
    
    # Perform clade partition tests:
    CladePartitionTestResults <- lapply(CladePartitionsToTest, PerformBranchPartitionTests)
    
    # Update branch rates to clade rates in output names:
    CladePartitionTestResults <- lapply(CladePartitionTestResults, function(x) {names(x) <- gsub("BranchRates", "CladeRates", names(x)); return(x)})
    
  # If not performing clade partition tests:
  } else {
    
    # Make empty clade partition result output:
    CladePartitionTestResults <- NULL
    
  }
  
  # If performing time bin partition tests:
  if(!is.null(TimeBinPartitionsToTest)) {
    
    # Build all changes matrix:
    AllCharacterChanges <- do.call(rbind, lapply(EdgeList, function(x) x$CharacterChanges))
    
    # Get weighted changes from each bin:
    WeightedChangesByBin <- unlist(lapply(as.list(1:(length(TimeBins) - 1)), function(x) sum(AllCharacterChanges[which(AllCharacterChanges[, "Bin"] == x), "Weight"])))
    
    # If using Close approach set completeness (including time component) for each time bin:
    if(TimeBinApproach == "Close") ProportionalCompletenessByBin <- apply(apply(do.call(rbind, lapply(EdgeList, function(x) x$ProportionalBinnedEdgeDurations)), 2, '*', unlist(lapply(EdgeList, function(x) sum(Weights[x$ComparableCharacters])))), 2, sum)
    
    # If using Lloyd approach set completeness (including time component) for each time bin:
    if(TimeBinApproach == "Lloyd") DurationCompletenessByBin <- apply(apply(do.call(rbind, lapply(EdgeList, function(x) x$BinnedBranchDurations)), 2, '*', unlist(lapply(EdgeList, function(x) sum(Weights[x$ComparableCharacters])))), 2, sum)
    
    # If using Close approach set mean time bin rate:
    if(TimeBinApproach == "Close") MeanTimeBinRate <- sum(WeightedChangesByBin) / sum(ProportionalCompletenessByBin)
    
    # If using Lloyd approach set mean time bin rate:
    if(TimeBinApproach == "Lloyd") MeanTimeBinRate <- sum(WeightedChangesByBin) / sum(DurationCompletenessByBin)
    
    # Perform time bin partition tests:
    TimeBinTestResults <- lapply(TimeBinPartitionsToTest, PerformTimeBinPartitionTests)
    
  # If not performing time bin partition tests:
  } else {
    
    # Make empty time bin partition result output:
    TimeBinTestResultsPartitionTestResults <- NULL
    
  }
  
  # Compile output:
  Output <- list(BranchPartitionTestResults, CharacterPartitionTestResults, CladePartitionTestResults, TimeBinTestResults)
  
  # Add naems to output:
  names(Output) <- c("BranchPartitionResults", "CharacterPartitionResults", "CladePartitionResults", "TimeBinResults")
  
  # Return output:
  return(Output)

}

#set.seed(17)
#tree <- rtree(nrow(Day2016$Matrix_1$Matrix))
#tree$tip.label <- rownames(Day2016$Matrix_1$Matrix)
#tree$root.time <- max(diag(vcv(tree)))
#TimeBins <- seq(from = tree$root.time, to = 0, length.out = 5)

#x <- DiscreteCharacterRate(tree = tree, clad.matrix = Day2016, TimeBins = TimeBins, BranchPartitionsToTest = lapply(as.list(1:nrow(tree$edge)), as.list), CharacterPartitionsToTest = lapply(as.list(1:3), as.list), CladePartitionsToTest = lapply(as.list(Ntip(tree) + (2:Nnode(tree))), as.list), TimeBinPartitionsToTest = lapply(as.list(1:(length(TimeBins) - 1)), as.list), ChangeTimes = "random", Alpha = 0.01, PolymorphismState = "missing", UncertaintyState = "missing", InapplicableState = "missing", TimeBinApproach = "Lloyd", EnsureAllWeightsAreIntegers = FALSE, EstimateAllNodes = FALSE, EstimateTipValues = FALSE, InapplicablesAsMissing = FALSE, PolymorphismBehaviour = "equalp", UncertaintyBehaviour = "equalp", Threshold = 0.01)

#x <- DiscreteCharacterRate(tree = tree, clad.matrix = Day2016, TimeBins = TimeBins, BranchPartitionsToTest = list(list(match(1:Ntip(tree), tree$edge[, 2]))), CharacterPartitionsToTest = lapply(as.list(1:3), as.list), CladePartitionsToTest = lapply(as.list(Ntip(tree) + (2:Nnode(tree))), as.list), TimeBinPartitionsToTest = lapply(as.list(1:(length(TimeBins) - 1)), as.list), ChangeTimes = "random", Alpha = 0.01, PolymorphismState = "missing", UncertaintyState = "missing", InapplicableState = "missing", TimeBinApproach = "Lloyd", EnsureAllWeightsAreIntegers = FALSE, EstimateAllNodes = FALSE, EstimateTipValues = FALSE, InapplicablesAsMissing = FALSE, PolymorphismBehaviour = "equalp", UncertaintyBehaviour = "equalp", Threshold = 0.01)

#x



# CHECK BELOW IS TRUE

    # Number of zero-duration branches:
#nzero.tb <- nzero.ib <- nzero <- 0

    # Do not count zero-duration branches in sample:
#m <- nrow(bin.results) - nzero

    # Do not count zero-duration branches in sample:
#m.tb <- nrow(bin.results.tb) - nzero.tb

    # Do not count zero-duration branches in sample:
#m.ib <- nrow(bin.results.ib) - nzero.ib

    # Calculate cutoffs (significance thresholds) based on Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR):
#cutoffs <- c((1:m) / m * Alpha, rep(0, nzero))

    # Calculate cutoffs (significance thresholds) based on Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR):
#cutoffs.tb <- c((1:m.tb) / m.tb * Alpha, rep(0, nzero.tb))

    # Calculate cutoffs (significance thresholds) based on Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR):
#cutoffs.ib <- c((1:m.ib) / m.ib * Alpha, rep(0, nzero.ib))

    # Get indices ready for identifying significant p values:
#ifelse(length(which(sort(bin.pvals) <= cutoffs)) > 0, indices <- 1:max((1:m)[sort(bin.pvals) <= cutoffs]), indices <- c(1)[-1])

    # Get indices ready for identifying significant p values:
#ifelse(length(which(sort(bin.pvals.tb) <= cutoffs.tb)) > 0, indices.tb <- 1:max((1:m.tb)[sort(bin.pvals.tb) <= cutoffs.tb]), indices.tb <- c(1)[-1])

    # Get indices ready for identifying significant p values:
#ifelse(length(which(sort(bin.pvals.ib) <= cutoffs.ib)) > 0, indices.ib <- 1:max((1:m.ib)[sort(bin.pvals.ib) <= cutoffs.ib]), indices.ib <- c(1)[-1])

    # Isolate significant p values:
#signif <- order(bin.pvals)[indices]

    # Isolate significant p values:
#signif.tb <- order(bin.pvals.tb)[indices.tb]

    # Isolate significant p values:
#signif.ib <- order(bin.pvals.ib)[indices.ib]

    # Add 1 in significance column for significant p-values after FDR correction:
#bin.results[signif, "ml.signif"] <- 1

    # Add 1 in significance column for significant p-values after FDR correction:
#bin.results.tb[signif.tb, "ml.signif"] <- 1

    # Add 1 in significance column for significant p-values after FDR correction:
#bin.results.ib[signif.ib, "ml.signif"] <- 1

    # Indicate significantly high rate bins:
#bin.results[(bin.results[, "ml.signif"] & bin.results[, "in.rate"] > bin.results[, "out.rate"]), "ml.signif.hi"] <- 1

    # Indicate significantly high rate bins:
#bin.results.tb[(bin.results.tb[, "ml.signif"] & bin.results.tb[, "in.rate"] > bin.results.tb[, "out.rate"]), "ml.signif.hi"] <- 1

    # Indicate significantly high rate bins:
#bin.results.ib[(bin.results.ib[, "ml.signif"] & bin.results.ib[, "in.rate"] > bin.results.ib[, "out.rate"]), "ml.signif.hi"] <- 1
  
    # Indicate significantly low rate bins:
#bin.results[(bin.results[, "ml.signif"] & bin.results[, "in.rate"] < bin.results[, "out.rate"]), "ml.signif.lo"] <- 1

    # Indicate significantly low rate bins:
#bin.results.tb[(bin.results.tb[, "ml.signif"] & bin.results.tb[, "in.rate"] < bin.results.tb[, "out.rate"]), "ml.signif.lo"] <- 1

    # Indicate significantly low rate bins:
#bin.results.ib[(bin.results.ib[, "ml.signif"] & bin.results.ib[, "in.rate"] < bin.results.ib[, "out.rate"]), "ml.signif.lo"] <- 1

    # Round chi-square values in output:
#bin.results[, "ml.chisq"] <- round(bin.results[, "ml.chisq"], 2)

    # Round chi-square values in output:
#bin.results.tb[, "ml.chisq"] <- round(bin.results.tb[, "ml.chisq"], 2)

    # Round chi-square values in output:
#bin.results.ib[, "ml.chisq"] <- round(bin.results.ib[, "ml.chisq"], 2)
  
    # Round p-values in output:
#bin.results[, "ml.pval"] <- round(bin.results[, "ml.pval"], 3)

    # Round p-values in output:
#bin.results.tb[, "ml.pval"] <- round(bin.results.tb[, "ml.pval"], 3)

    # Round p-values in output:
#bin.results.ib[, "ml.pval"] <- round(bin.results.ib[, "ml.pval"], 3)
  
# WHAT TO DO WITH ZERO VALUES IN TIME SERIES? EXCLUDE?
# TIME IS KEY THING TO CHECK (IF ZERO THEN NO CHANCE TO OBSERVE ANYTHING)

  # Case if equal rates cannot be rejected:
#} else {

    # If not then print notification and stop:
#cat(paste("H_0 - all rates equal across time bins - cannot be rejected at an Alpha of ", Alpha, " (Actual p = ", chisq.p, ").\nCalculations of per-bin rates aborted.", sep=""))
    
# NEED TO MODIFY THE BELOW AS MAY BE SIG DIFFS FOR INTERNAL OR TERMINAL BRANCHES NOT SEEN IN POOLED

    # Create NULL outputs:
#bin.results.tb <- bin.results.ib <- bin.results <- message(paste("H_0 - all rates equal across time bins - cannot be rejected at an Alpha of ", Alpha, " (Actual p = ", chisq.p, ").\nA single rate of ", mlenumer, "is preferred.", sep = ""))
    
#}

  # Get number of branches:
#n <- length(tree$edge[, 1])
  
  # Get number of tips (and terminal branches):
#ntips <- n.tb <- ape::Ntip(tree)
  
  # Get number of internal branches:
#n.ib <- n - n.tb
  
  # Get total number of character changes on tree:
#changes <- observed.tree$edge.length
  
  # Get number of character changes on terminal branches:
#changes.tb <- observed.tree$edge.length[terminal.branches]
  
  # Get number of character changes on internal branches:
#changes.ib <- observed.tree$edge.length[internal.branches]
  
  # Get list of nonterminal and nonroot nodes:
#nodes <- (ape::Ntip(tree) + 2):(ape::Ntip(tree) + ape::Nnode(tree))
  
  # Percentage of characters that are observable:
#pctcomp <- completeness.tree$edge.length / Nchar
  
  # Percentage of characters that are observable on terminal branches:
#pctcomp.tb <- completeness.tree$edge.length[terminal.branches] / Nchar
  
  # Percentage of characters that are observable on internal branches:
#pctcomp.ib <- completeness.tree$edge.length[internal.branches] / Nchar
  
  # Get branch durations:
#Time <- tree$edge.length
  
  # Get terminal branch durations:
#Time.tb <- tree$edge.length[terminal.branches]
  
  # Get internal branch durations:
#Time.ib <- tree$edge.length[internal.branches]
  
  # Get number of zero length branches:
#nzero <- length(which(Time == 0))
  
  # Get number of zero length terminal branches:
#nzero.tb <- length(which(Time.tb == 0))
  
  # Get number of zero length internal branches:
#nzero.ib <- length(which(Time.ib == 0))
  
  # Get maximum likelihood numerator:
#mlenumer <- sum(changes) / sum(Time * pctcomp)
  
  # Get maximum likelihood numerator for terminal branches:
#mlenumer.tb <- sum(changes.tb) / sum(Time.tb * pctcomp.tb)
  
  # Get maximum likelihood numerator for internal branches:
#mlenumer.ib <- sum(changes.ib) / sum(Time.ib * pctcomp.ib)
  
  # Get maximum likelihood denominator:
#mledenom <- changes / (Time * pctcomp)
  
  # Get maximum likelihood denominator for terminal branches:
#mledenom.tb <- changes.tb / (Time.tb * pctcomp.tb)
  
  # Get maximum likelihood denominator for internal branches:
#mledenom.ib <- changes.ib / (Time.ib * pctcomp.ib)
  
  # Get log numerator:
#lognumer <- sum(log(dpois(changes, mlenumer * Time * pctcomp)))
  
  # Get log numerator for terminal branches:
#lognumer.tb <- sum(log(dpois(changes.tb, mlenumer.tb * Time.tb * pctcomp.tb)))
  
  # Get log numerator for internal branches:
#lognumer.ib <- sum(log(dpois(changes.ib, mlenumer.ib * Time.ib * pctcomp.ib)))
  
  # Get log denominator with zero-length branches removed:
#logdenom <- sum(log(dpois(changes, mledenom * Time * pctcomp)), na.rm = TRUE)
  
  # Get log denominator for terminal branches with zero-length branches removed:
#logdenom.tb <- sum(log(dpois(changes.tb, mledenom.tb * Time.tb * pctcomp.tb)), na.rm = TRUE)
  
  # Get log denominator for internal branches with zero-length branches removed:
#logdenom.ib <- sum(log(dpois(changes.ib, mledenom.ib * Time.ib * pctcomp.ib)), na.rm = TRUE)
  
  # Get test statistic:
#teststat <- -2 * (lognumer - logdenom)
  
  # Get test statistic for terminal branches:
#teststat.tb <- -2 * (lognumer.tb - logdenom.tb)
  
  # Get test statistic for internal branches:
#teststat.ib <- -2 * (lognumer.ib - logdenom.ib)
  
  # Calculate position of test statistic in chi-square distribution to get probability (zero-length branches not calculated in df):
#chisq.p <- pchisq(teststat, n - 1 - nzero, lower.tail = FALSE)
  
  # Calculate position of test statistic in chi-square distribution for terminal branches to get probability (zero-length branches not calculated in df):
#chisq.p.tb <- pchisq(teststat.tb, n.tb - 1 - nzero.tb, lower.tail = FALSE)
  
  # Calculate position of test statistic in chi-square distribution for internal branches to get probability (zero-length branches not calculated in df):
#chisq.p.ib <- pchisq(teststat.ib, n.ib - 1 - nzero.ib, lower.tail = FALSE)
  
  # Check to see if null hypothesis of equal rates across the tree can be rejected:
#if(chisq.p < Alpha) {
    
    # If so then print notification and carry on:
    #cat(paste("H_0 - all rates equal across the tree - is rejected at an Alpha of ", Alpha, " (actual p = ", chisq.p, ").\nContinuing to per-branch and per-clade rate calculations.", sep = ""))
    
    # Create matrix to store branch results:
    #branch.results <- matrix(0, nrow = n, ncol = 20)
    
    # Create matrix to store terminal branch results:
    #branch.results.tb <- matrix(0, nrow = n.tb, ncol = 20)
    
    # Create matrix to store internal branch results:
    #branch.results.ib <- matrix(0, nrow = n.ib, ncol = 20)
    
    # Create matrix to store node results:
    #node.results <- matrix(0, nrow = length(nodes), ncol = 29)
    
    # Add column names for branches:
    #colnames(branch.results) <- colnames(branch.results.tb) <- colnames(branch.results.ib) <- c("branch", "from", "to", "in.rate", "out.rate", "ml.chisq", "ml.pval", "ml.signif.hi", "ml.signif.hi.ti", "ml.signif.lo", "ml.signif.lo.ti", "ml.signif", "ml.signif.ti", "rand.val", "rand.mean", "rand.sd", "rand.pval", "rand.signif.hi", "rand.signif.lo", "rand.signif")
    
    # Add column names for nodes:
    #colnames(node.results) <- c("node", "in.rate", "in.rate.tb", "in.rate.ib", "out.rate", "out.rate.tb", "out.rate.ib", "ml.chisq", "ml.chisq.tb", "ml.chisq.ib", "ml.pval", "ml.pval.tb", "ml.pval.ib", "ml.signif.hi", "ml.signif.hi.tb", "ml.signif.hi.ib", "ml.signif.lo", "ml.signif.lo.tb", "ml.signif.lo.ib", "ml.signif", "ml.signif.tb", "ml.signif.ib", "rand.val", "rand.mean", "rand.sd", "rand.pval", "rand.signif.hi", "rand.signif.lo", "rand.signif")
    
    # Number branches 1 to N:
    #branch.results[, "branch"] <- 1:n
    
    # Number terminal branches:
    #branch.results.tb[, "branch"] <- terminal.branches
    
    # Number internal branches:
    #branch.results.ib[, "branch"] <- internal.branches
    
    # Add from and to node numbers:
    #branch.results[, c("from", "to")] <- tree$edge[, 1:2]
    
    # Add from and to node numbers for terminal branches:
    #branch.results.tb[, c("from", "to")] <- tree$edge[terminal.branches, 1:2]
    
    # Add from and to node numbers for internal branches:
    #branch.results.ib[, c("from", "to")] <- tree$edge[internal.branches, 1:2]
    
    # Number nodes:
    #node.results[, "node"] <- nodes
    
    # For each branch:
    #for (i in 1:n) {
      
      # Check if branch is terminal:
      #if(length(sort(match(i, terminal.branches)))) {
        
        # If yes set as TRUE:
        #branch.is.terminal <- TRUE
        
      # Branch is terminal:
      #} else {
        
        # If no set as FALSE:
        #branch.is.terminal <- FALSE
        
        #}
      
      # Get numbers of other (not ith) branches:
      #other <- (1:n)[-i]
      
      # Get numbers of other (not ith) terminal or internal branches:
      #ifelse(branch.is.terminal, other.tb <- terminal.branches[-match(i, terminal.branches)], other.ib <- internal.branches[-match(i, internal.branches)])
      
      # Create maximum likelihood denominator variable:
      #mledenom <- rep(NA, n)
      
      # Create maximum likelihood denominator variable for terminal or internal branches:
      #ifelse(branch.is.terminal, mledenom.tb <- rep(NA, n.tb), mledenom.ib <- rep(NA, n.ib))
      
      # Calculate maximum likelihood denominator for branch:
      #mledenom[i] <- mlebranch <- changes[i] / (Time[i] * pctcomp[i])
      
      # Calculate maximum likelihood denominator for terminal or internal branch:
      #ifelse(branch.is.terminal, mledenom.tb[match(i, terminal.branches)] <- mlebranch.tb <- changes.tb[match(i, terminal.branches)] / (Time.tb[match(i, terminal.branches)] * pctcomp.tb[match(i, terminal.branches)]), mledenom.ib[match(i, internal.branches)] <- mlebranch.ib <- changes.ib[match(i, internal.branches)] / (Time.ib[match(i, internal.branches)] * pctcomp.ib[match(i, internal.branches)]))
      
      # Calculate maximum likelihood denominator for other branches:
      #mledenom[other] <- mleother <- sum(changes[other]) / sum(Time[other] * pctcomp[other])
      
      # Calculate maximum likelihood denominator for other terminal or internal branches:
      #ifelse(branch.is.terminal, mledenom.tb[-match(i, terminal.branches)] <- mleother.tb <- sum(changes[other.tb]) / sum(Time[other.tb] * pctcomp[other.tb]), mledenom.ib[-match(i, internal.branches)] <- mleother.ib <- sum(changes[other.ib]) / sum(Time[other.ib] * pctcomp[other.ib]))
      
      # Calculate log numerator:
      #lognumer <- sum(log(dpois(changes, mlenumer * Time * pctcomp)), na.rm = TRUE)
      
      # Calculate log numerator for terminal or internal branches:
      #ifelse(branch.is.terminal, lognumer.tb <- sum(log(dpois(changes.tb, mlenumer.tb * Time.tb * pctcomp.tb)), na.rm = TRUE), lognumer.ib <- sum(log(dpois(changes.ib, mlenumer.ib * Time.ib * pctcomp.ib)), na.rm = TRUE))
      
      # Calculate log denominator:
      #logdenom <- sum(log(dpois(changes, mledenom * Time * pctcomp)), na.rm = TRUE)
      
      # Calculate log denominator for terminal or internal branches:
      #ifelse(branch.is.terminal, logdenom.tb <- sum(log(dpois(changes.tb, mledenom.tb * Time.tb * pctcomp.tb)), na.rm = TRUE), logdenom.ib <- sum(log(dpois(changes.ib, mledenom.ib * Time.ib * pctcomp.ib)), na.rm = TRUE))
      
      # Store rate for branch:
      #branch.results[i, "in.rate"] <- round(mlebranch, 2)
      
      # Store rate for terminal or internal branch:
      #ifelse(branch.is.terminal, branch.results.tb[match(i, branch.results.tb[, "branch"]), "in.rate"] <- round(mlebranch.tb, 2), branch.results.ib[match(i, branch.results.ib[, "branch"]), "in.rate"] <- round(mlebranch.ib, 2))
      
      # Store rate for outside branches:
      #branch.results[i, "out.rate"] <- round(mleother, 2)
      
      # Store rate for outside terminal or internal branches:
      #ifelse(branch.is.terminal, branch.results.tb[match(i, branch.results.tb[, "branch"]), "out.rate"] <- round(mleother.tb, 2), branch.results.ib[match(i, branch.results.ib[, "branch"]), "out.rate"] <- round(mleother.ib, 2))
      
      # Store chi-squared value (will be zero for zero-duration branch):
      #branch.results[i, "ml.chisq"] <- teststat <- -2 * (lognumer - logdenom)
      
      # Store chi-squared value (will be zero for zero-duration branch) for internal or terminal branches:
      #ifelse(branch.is.terminal, branch.results.tb[match(i, branch.results.tb[, "branch"]), "ml.chisq"] <- teststat.tb <- -2 * (lognumer.tb - logdenom.tb), branch.results.ib[match(i, branch.results.ib[, "branch"]), "ml.chisq"] <- teststat.ib <- -2 * (lognumer.ib - logdenom.ib))
      
      # Store probability for branch (will be 1 for zero-duration branch):
      #branch.results[i, "ml.pval"] <- pchisq(teststat, 2 - 1, lower.tail = FALSE)
      
      # Store probability for branch (will be 1 for zero-duration branch) for terminal or internal branches:
      #ifelse(branch.is.terminal, branch.results.tb[match(i, branch.results.tb[, "branch"]), "ml.pval"] <- pchisq(teststat.tb, 2 - 1, lower.tail = FALSE), branch.results.ib[match(i, branch.results.ib[, "branch"]), "ml.pval"] <- pchisq(teststat.ib, 2 - 1, lower.tail = FALSE))
      
      #}
    
    # Set initial row number for storing data:
    #j <- 1
    
    # For each node (excluding the root):
    #for (i in nodes) {
      
      # Identify branches within clade:
      #clade <- GetDescendantEdges(i, tree)
      
      # Identify terminal branches within clade:
      #clade.tb <- clade[sort(match(terminal.branches, clade))]
      
      # Identify internal branches within clade:
      #clade.ib <- clade[sort(match(internal.branches, clade))]
      
      # Identify branches outside clade:
      #nonclade <- (1:n)[-clade]
      
      # Identify terminal branches outside clade:
      #nonclade.tb <- setdiff(terminal.branches, clade.tb)
      
      # Identify internal branches outside clade:
      #nonclade.ib <- setdiff(internal.branches, clade.ib)
      
      # Set empty maximum likelihood denominator vector:
      #mledenom <- rep(NA, n)
      
      # Set empty maximum likelihood denominator vector for terminal branches:
      #mledenom.tb <- rep(NA, n.tb)
      
      # Set empty maximum likelihood denominator vector for internal branches:
      #mledenom.ib <- rep(NA, n.ib)
      
      # Fill within clade values for maximum likelihood denominator:
      #mledenom[clade] <- claderate <- sum(changes[clade]) / sum(Time[clade] * pctcomp[clade])
      
      # Fill within clade values for maximum likelihood denominator for terminal branches:
      #mledenom.tb[match(clade.tb, terminal.branches)] <- claderate.tb <- sum(changes[clade.tb]) / sum(Time[clade.tb] * pctcomp[clade.tb])
      
      # Fill within clade values for maximum likelihood denominator for internal branches (if present):
      #if(length(clade.ib) > 0) mledenom.ib[match(clade.ib, internal.branches)] <- claderate.ib <- sum(changes[clade.ib]) / sum(Time[clade.ib] * pctcomp[clade.ib])
      
      # Fill outside clade values for maximum likelihood denominator:
      #mledenom[nonclade] <- noncladerate <- sum(changes[nonclade]) / sum(Time[nonclade] * pctcomp[nonclade])
      
      # Fill outside clade values for maximum likelihood denominator for terminal branches:
      #mledenom.tb[match(nonclade.tb, terminal.branches)] <- noncladerate.tb <- sum(changes[nonclade.tb]) / sum(Time[nonclade.tb] * pctcomp[nonclade.tb])
      
      # Fill outside clade values for maximum likelihood denominator for internal branches (if present):
      #if(length(clade.ib) > 0) mledenom.ib[match(nonclade.ib, internal.branches)] <- noncladerate.ib <- sum(changes[nonclade.ib]) / sum(Time[nonclade.ib] * pctcomp[nonclade.ib])
      
      # Set log-likelihood numerator:
      #lognumer <- sum(log(dpois(changes, mlenumer * Time * pctcomp)), na.rm = TRUE)
      
      # Set log-likelihood numerator for terminal branches:
      #lognumer.tb <- sum(log(dpois(changes.tb, mlenumer.tb * Time.tb * pctcomp.tb)), na.rm = TRUE)
      
      # Set log-likelihood numerator for internal branches (if present):
      #if(length(clade.ib) > 0) lognumer.ib <- sum(log(dpois(changes.ib, mlenumer.ib * Time.ib * pctcomp.ib)), na.rm = TRUE)
      
      # Set log-likelihood denominator:
      #logdenom <- sum(log(dpois(changes, mledenom * Time * pctcomp)), na.rm = TRUE)
      
      # Set log-likelihood denominator for terminal branches:
      #logdenom.tb <- sum(log(dpois(changes.tb, mledenom.tb * Time.tb * pctcomp.tb)), na.rm = TRUE)
      
      # Set log-likelihood denominator for internal branches (if present):
      #if(length(clade.ib) > 0) logdenom.ib <- sum(log(dpois(changes.ib, mledenom.ib * Time.ib * pctcomp.ib)), na.rm = TRUE)
      
      # Record within clade rate:
      #node.results[j, "in.rate"] <- round(claderate, 2)
      
      # Record within clade rate for terminal branches:
      #node.results[j, "in.rate.tb"] <- round(claderate.tb, 2)
      
      # Record within clade rate for internal branches (if present):
      #ifelse(length(clade.ib) > 0, node.results[j, "in.rate.ib"] <- round(claderate.ib, 2), node.results[j, "in.rate.ib"] <- NA)
      
      # Record outside clade rate:
      #node.results[j, "out.rate"] <- round(noncladerate, 2)
      
      # Record outside clade rate for terminal branches:
      #node.results[j, "out.rate.tb"] <- round(noncladerate.tb, 2)
      
      # Record outside clade rate for internal branches (if present):
      #ifelse(length(clade.ib) > 0, node.results[j, "out.rate.ib"] <- round(noncladerate.ib, 2), node.results[j, "out.rate.ib"] <- NA)
      
      # Record chi-square test statistic:
      #node.results[j, "ml.chisq"] <- teststat <- -2 * (lognumer - logdenom)
      
      # Record chi-square test statistic for terminal branches:
      #node.results[j, "ml.chisq.tb"] <- teststat.tb <- -2 * (lognumer.tb - logdenom.tb)
      
      # Record chi-square test statistic for internal branches (if present):
      #ifelse(length(clade.ib) > 0, node.results[j, "ml.chisq.ib"] <- teststat.ib <- -2 * (lognumer.ib - logdenom.ib), node.results[j, "ml.chisq.ib"] <- NA)
      
      # Record p-value for chi-squared test:
      #node.results[j, "ml.pval"] <- testresult <- pchisq(teststat, 2 - 1, lower.tail = FALSE)
      
      # Record p-value for chi-squared test for terminal branches:
      #node.results[j, "ml.pval.tb"] <- testresult.tb <- pchisq(teststat.tb, 2 - 1, lower.tail = FALSE)
      
      # Record p-value for chi-squared test for terminal branches:
      #ifelse(length(clade.ib) > 0, node.results[j, "ml.pval.ib"] <- testresult.ib <- pchisq(teststat.ib, 2 - 1, lower.tail=F), node.results[j, "ml.pval.ib"] <- NA)
      
      # Update row number:
      #j <- j + 1
      
      #}
    
    # Combine terminal and internal branch rate results and order as branch.results:
    #termandintern.results <- rbind(branch.results.tb, branch.results.ib)[order(rbind(branch.results.tb, branch.results.ib)[, "branch"]), ]
    
    # Create vector to store true (1; terminal) / false (0; internal) branch type:
    #branchisterminal <- rep(0, nrow(tree$edge))
    
    # Store 1 for terminals:
    #branchisterminal[terminal.branches] <- 1
    
    # Combine and reduce branch type and terminal and internal branch rate results (ignores randomisations):
    #termandintern.results <- cbind(branchisterminal, termandintern.results[, c("in.rate", "out.rate", "ml.chisq", "ml.pval")])
    
    # Update column names (to avoid clashes with branch.results):
    #colnames(termandintern.results) <- c("is.term", "in.rate.ti", "out.rate.ti", "ml.chisq.ti", "ml.pval.ti")
    
    # Combine all branch rates versus split terminal-internal branches:
    #branch.results <- cbind(branch.results, termandintern.results)
    
    # Reorder and cut down (ignores randomisation results):
    #branch.results <- branch.results[, c("branch", "from", "to", "is.term", "in.rate", "in.rate.ti", "out.rate", "out.rate.ti", "ml.chisq", "ml.chisq.ti", "ml.pval", "ml.pval.ti", "ml.signif.hi", "ml.signif.hi.ti", "ml.signif.lo", "ml.signif.lo.ti", "ml.signif", "ml.signif.ti")]
    
    # Get just the branch p-values:
    #branch.pvals <- branch.results[, "ml.pval"]
    
    # Get just the terminal branch p-values:
    #branch.pvals.tb <- branch.results[terminal.branches, "ml.pval.ti"]
    
    # Get just the internal branch p-values:
    #branch.pvals.ib <- branch.results[internal.branches, "ml.pval.ti"]
    
    # Do not count zero-duration branches in sample:
    #m <- n - nzero
    
    # Do not count zero-duration terminal branches in sample:
    #m.tb <- n.tb - nzero.tb
    
    # Do not count zero-duration internal branches in sample:
    #m.ib <- n.ib - nzero.ib
    
    # Calculate cutoffs (significance thresholds) based on Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR):
    #cutoffs <- c((1:m) / m * Alpha, rep(0, nzero))
    
    # Calculate cutoffs (significance thresholds) for terminal branches:
    #cutoffs.tb <- c((1:m.tb) / m.tb * Alpha, rep(0, nzero.tb))
    
    # Calculate cutoffs (significance thresholds) for internal branches:
    #cutoffs.ib <- c((1:m.ib) / m.ib * Alpha, rep(0, nzero.ib))
    
    # Get indices ready for identifying significant p values:
    #ifelse(length(which(sort(branch.pvals) <= cutoffs)) > 0, indices <- 1:max((1:m)[sort(branch.pvals) <= cutoffs]), indices <- c(1)[-1])
    
    # Get indices ready for identifying significant terminal p values:
    #ifelse(length(which(sort(branch.pvals.tb) <= cutoffs.tb)) > 0, indices.tb <- 1:max((1:m.tb)[sort(branch.pvals.tb) <= cutoffs.tb]), indices.tb <- c(1)[-1])
    
    # Get indices ready for identifying significant internal p values:
    #ifelse(length(which(sort(branch.pvals.ib) <= cutoffs.ib)) > 0, indices.ib <- 1:max((1:m.ib)[sort(branch.pvals.ib) <= cutoffs.ib]), indices.ib <- c(1)[-1])
    
    # Isolate significant p values:
    #signif <- order(branch.pvals)[indices]
    
    # Isolate significant terminal p values:
    #signif.tb <- terminal.branches[order(branch.pvals.tb)[indices.tb]]
    
    # Isolate significant internal p values:
    #signif.ib <- internal.branches[order(branch.pvals.ib)[indices.ib]]
    
    # Add 1 in significance column for significant p-values after FDR correction:
    #branch.results[signif, "ml.signif"] <- 1
    
    # Add 1 in significance column for significant terminal p-values after FDR correction:
    #branch.results[signif.tb, "ml.signif.ti"] <- 1
    
    # Add 1 in significance column for significant internal p-values after FDR correction:
    #branch.results[signif.ib, "ml.signif.ti"] <- 1
    
    # Indicate significantly high rate branches:
    #branch.results[(branch.results[, "ml.signif"] & branch.results[, "in.rate"] > branch.results[, "out.rate"]), "ml.signif.hi"] <- 1
    
    # Indicate significantly high rate terminal and internal branches:
    #branch.results[(branch.results[, "ml.signif.ti"] & branch.results[, "in.rate.ti"] > branch.results[, "out.rate.ti"]), "ml.signif.hi.ti"] <- 1
    
    # Indicate significantly low rate branches:
    #branch.results[(branch.results[, "ml.signif"] & branch.results[, "in.rate"] < branch.results[, "out.rate"]), "ml.signif.lo"] <- 1
    
    # Indicate significantly low rate terminal and internal branches:
    #branch.results[(branch.results[, "ml.signif.ti"] & branch.results[, "in.rate.ti"] < branch.results[, "out.rate.ti"]), "ml.signif.lo.ti"] <- 1
    
    # Round chi-square values in output:
    #branch.results[, "ml.chisq"] <- round(branch.results[, "ml.chisq"], 2)
    
    # Round chi-square values in output:
    #branch.results[, "ml.chisq.ti"] <- round(branch.results[, "ml.chisq.ti"], 2)
    
    # Round p-values in output:
    #branch.results[, "ml.pval"] <- round(branch.results[, "ml.pval"], 3)

    # Round p-values in output:
    #branch.results[, "ml.pval.ti"] <- round(branch.results[, "ml.pval.ti"], 3)
    
    # Isolate the node p-values:
    #node.pvals <- node.results[, "ml.pval"]
    
    # Isolate the node p-values for terminal branches:
    #node.pvals.tb <- node.results[, "ml.pval.tb"]
    
    # Isolate the node p-values for internal branches:
    #node.pvals.ib <- node.results[, "ml.pval.ib"]
    
    # Set m as number of nodes (excluding root), serves for terminal branches alone too:
    #m <- length(nodes)
    
    # Set m as number of non-cherry nodes (excluding root and NAs):
    #m.ib <- length(sort(node.pvals.ib))
    
    # Calculate cutoffs (significance thresholds) based on Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR):
    #cutoffs <- (1:m) / m * Alpha
    
    # Calculate cutoffs for internal branches:
    #cutoffs.ib <- (1:m.ib) / m.ib * Alpha
    
    # Get indices ready for identifying significant p values:
    #ifelse(length(which(sort(node.pvals) <= cutoffs)) > 0, indices <- 1:max((1:m)[sort(node.pvals) <= cutoffs]), indices <- c(1)[-1])
    
    # Get indices ready for identifying significant p values for terminal branches:
    #ifelse(length(which(sort(node.pvals.tb) <= cutoffs)) > 0, indices.tb <- 1:max((1:m)[sort(node.pvals.tb) <= cutoffs]), indices.tb <- c(1)[-1])
    
    # Get indices ready for identifying significant p values for internal branches:
    #ifelse(length(which(sort(node.pvals.ib) <= cutoffs.ib)) > 0, indices.ib <- 1:max((1:m.ib)[sort(node.pvals.ib) <= cutoffs.ib]), indices.ib <- c(1)[-1])
    
    # Isolate significant p-values:
    #signif <- order(node.pvals)[indices]
    
    # Isolate significant p-values for terminal branches:
    #signif.tb <- order(node.pvals.tb)[indices.tb]
    
    # Isolate significant p-values for internal branches:
    #signif.ib <- order(node.pvals.ib)[indices.ib]
    
    # Add NAs for terminal branch only clades (cherries):
    #node.results[is.na(node.results[, "ml.chisq.ib"]), c("ml.signif.hi.ib", "ml.signif.lo.ib", "ml.signif.ib")] <- NA
    
    # Record significant clades:
    #node.results[signif, "ml.signif"] <- 1
    
    # Record significant terminal branches clades:
    #node.results[signif.tb, "ml.signif.tb"] <- 1
    
    # Record significant internal branches clades:
    #node.results[signif.ib, "ml.signif.ib"] <- 1
    
    # Round chi-squared test statistic to 1dp:
    #node.results[, "ml.chisq"] <- round(node.results[, "ml.chisq"], 1)
    
    # Round chi-squared test statistic to 1dp:
    #node.results[, "ml.chisq.tb"] <- round(node.results[, "ml.chisq.tb"], 1)
    
    # Round chi-squared test statistic to 1dp:
    #node.results[, "ml.chisq.ib"] <- round(node.results[, "ml.chisq.ib"], 1)
    
    # Round chi-squared p-value to 3dp:
    #node.results[, "ml.pval"] <- round(node.results[, "ml.pval"], 3)
    
    # Round chi-squared p-value to 3dp:
    #node.results[, "ml.pval.tb"] <- round(node.results[, "ml.pval.tb"], 3)
    
    # Round chi-squared p-value to 3dp:
    #node.results[, "ml.pval.ib"] <- round(node.results[, "ml.pval.ib"], 3)
    
    # Record significantly high clade rates:
    #node.results[(node.results[, "ml.signif"] & node.results[, "in.rate"] > node.results[, "out.rate"]), "ml.signif.hi"] <- 1
    
    # Record significantly high clade rates:
    #node.results[(node.results[, "ml.signif.tb"] & node.results[, "in.rate.tb"] > node.results[, "out.rate.tb"]), "ml.signif.hi.tb"] <- 1
    
    # Record significantly high clade rates:
    #node.results[intersect(which(node.results[, "ml.signif.ib"] == 1), which(node.results[, "in.rate.ib"] > node.results[, "out.rate.ib"])), "ml.signif.hi.ib"] <- 1
    
    # Record signficiantly low rates:
    #node.results[(node.results[, "ml.signif"] & node.results[, "in.rate"] < node.results[, "out.rate"]), "ml.signif.lo"] <- 1
    
    # Record signficiantly low rates:
    #node.results[(node.results[, "ml.signif.tb"] & node.results[, "in.rate.tb"] < node.results[, "out.rate.tb"]), "ml.signif.lo.tb"] <- 1
    
    # Record significantly low rates:
    #node.results[intersect(which(node.results[, "ml.signif.ib"] == 1), which(node.results[, "in.rate.ib"] < node.results[, "out.rate.ib"])), "ml.signif.lo.ib"] <- 1

  # Case if equal rates cannot be rejected:
  #} else {
    
    # If not then print notification and stop:
    #cat(paste("H_0 - all rates equal across the tree - cannot be rejected at an Alpha of ", Alpha, " (Actual p = ", chisq.p, ").\nCalculations of per-branch and per-clade rates aborted.", sep = ""))
    
    # Create NULL outputs:
    #branch.results <- node.results <- NULL
    
    #}
  
# COMBINE TIME BIN RESULTS INTO SINGLE VARIABLE?

  # List output matrices:
#out <- list(character.changes, node.results, branch.results, bin.results, bin.results.tb, bin.results.ib)
  
  # Add names to them:
#names(out) <- c("character.changes", "node.results", "branch.results", "per.bin.rates", "per.bin.rates.tb", "per.bin.rates.ib")
  
  # Return results:
#return(out)
  
#}
