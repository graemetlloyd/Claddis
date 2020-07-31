#' Discrete character rates across trees, time, and character types
#'
#' @description
#'
#' Given a tree and a cladistic-type matrix uses either likelihood ratio tests or the Akaike Information Criterion to compare rate models across branches, clades, time bins, or character partitions.
#'
#' @param tree A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{cladistic.matrix}.
#' @param cladistic.matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @param TimeBins A vector of ages (in millions of years) indicating the boundaries of a series of time bins in order from oldest to youngest.
#' @param BranchPartitionsToTest A list of branch(es) (edge number) partitions to test as N-rate parameter model (where N is the total number of partitions). If NULL (the default) then no partition test(s) will be made.
#' @param CharacterPartitionsToTest A list of character partition(s) (character numbers) to test as N-rate parameter model (where N is the total number of partitions). If NULL (the default) then no partition test(s) will be made.
#' @param CladePartitionsToTest A list of clade partition(s) (node numbers) to test as N-rate parameter model (where N is the total number of partitions). If NULL (the default) then no partition test(s) will be made.
#' @param TimeBinPartitionsToTest A list of time bin partition(s) (numbered 1 to N) to test as N-rate parameter model (where N is the total number of partitions). If NULL (the default) then no partition test(s) will be made.
#' @param ChangeTimes The time at which to record the character changes. One of \code{"midpoint"} (changes occur at the midpoint of the branch), \code{"spaced"} (changes equally spaced along branch), or \code{"random"} (change times drawn at random from a uniform distribution; the default and recommended option). Note: this is only meaningful if testing for time bin partitions.
#' @param LikelihoodTest Whether to apply an Akaike Information Criterion (\code{"AIC"}; the default) or likelihood ratio test (\code{"LRT"}).
#' @param alpha The alpha value to be used for the significance tests (only relevant if using the likelihood ratio test). The default is 0.01.
#' @param MultipleComparisonCorrection One of \code{"BenjaminiHochberg"} (the Benjamini and Hochberg 1995 false discovery rate approach; default and recommended) or \code{"Bonferroni"} (the Bonferroni correction). Only relevant if using the likelihood ratio test.
#' @param PolymorphismState One of \code{"missing"} (converts polymorphic values to NA; the default) or \code{"random"} (picks one of the possible polymorphic states at random).
#' @param UncertaintyState One of \code{"missing"} (converts uncertain values to NA; the default) or \code{"random"} (picks one of the possible uncertain states at random).
#' @param InapplicableState The only current option is \code{"missing"} (converts value to NA).
#' @param TimeBinApproach One of \code{"Close"} or \code{"Lloyd"} (the default).
#' @param EnsureAllweightsAreIntegers Logical for whether (\code{TRUE}) to reweight non-integer weights until all weights are integers or to leave them as they are (\code{FALSE}; the default).
#' @param estimate.all.nodes Option passed to internal use of \link{estimate_ancestral_states}.
#' @param estimate.tip.values Option passed to internal use of \link{estimate_ancestral_states}.
#' @param inapplicables.as.missing Option passed to internal use of \link{estimate_ancestral_states}.
#' @param polymorphism.behaviour Option passed to internal use of \link{estimate_ancestral_states}.
#' @param uncertainty.behaviour Option passed to internal use of \link{estimate_ancestral_states}.
#' @param threshold Option passed to internal use of \link{estimate_ancestral_states}.
#' @param allow.all.missing Logical to allow all missing character values - see \link{estimate_ancestral_states} for details.
#'
#' @details
#'
#' \bold{Introduction}
#'
#' Morphological change can be captured by discrete characters and their evolution modelled as occurring along the branches of a phylogenetic tree. This function takes as primary input a character-taxon matrix of discrete characters (in the format imported by \link{read_nexus_matrix}) and a time-scaled phylogenetic tree (in the format of \pkg{paleotree} or \pkg{strap}) and begins by inferring ancestral states at the tree's internal nodes using the \link{estimate_ancestral_states} function. From here changes along individual branches can be estimated (only the minimum number of changes are inferred; see \link{map_stochastic_changes} for an alternative but unfinished approach) and hence rates can be calculated.
#'
#' A discrete character rate can be expressed as the mean number of changes per million years (users may wish to normalise this by the number of characters for interpretation) and can be calculated for a branch (edge) of the tree, a clade (a mean rate for the edges descended from a single node), a character partition (the mean rate for a subset of the characters across all edges), or, most complex (see Lloyd 2016), the mean rate across the edges (or parts of edges) present in a time bin (defined by two values denoting the beginning and end of the time bin). In an ideal scenario these rates could be compared at face value, but that would require a large number of characters and very minimal (or zero) missing data. I.e., at an extreme of missing data if only one character can be observed along a branch it will either change (the maximum possible inferrable rate of evolution) or it will not (the minimum possible inferrable rate of evolution). In such cases it would be unwise to consider either outcome as being a significant departure from the mean rate.
#'
#' Because of these complications Lloyd et al. (2012) introduced tests by which the significance of an edge (or other partitioning of the data, i.e., a clade, time bin etc.) could be considered to be significantly high or low in comparison to the mean rate for the whole tree (i.e., whether a two-rate model could be considered more likely than a one-rate model). This is achieved through a likelihood ratio test:
#'
#' \deqn{LR = value of likehood function under the null (one-rate) hypothesis / maximum possible value of likehood function under the alternative (two-rate) hypotheses}
#'
#' Typically we might expect the two hypotheses to be well defined a priori. E.g., an expectation that a specific branch of the tree might have a higher or lower rate than background due to some evolutionary shift. However, Lloyd et al. (2012) instead provided an exploratory approach whereby every possible one edge value was compared with the rate for the rest of the tree (and the equivalent with clades and time bins). This was the default in Claddis up to version 0.2, but this has now been replaced (since version 0.3) with a more customisable set of options that allows different types of hypotheses (e.g., partitioning the data by character), as well as more complex hypotheses (e.g., a three-rate model), to be tested. Since version 0.4 the option to replace likelihood ratio tests with the Akaike Information Criterion has also been added.
#'
#' \bold{The four types of rate hypothesis}
#'
#' Following Cloutier (1991), Lloyd (2016) extended the two main types of rate hypotheses to four:
#'
#' \enumerate{
#'   \item A branch rate (available here with the \code{BranchPartitionsToTest} option).
#'   \item A clade rate (available here with the \code{CladePartitionsToTest} option).
#'   \item A time bin rate (available here with the \code{TimeBinPartitionsToTest} option).
#'   \item A character partition rate (available here with the \code{CharacterPartitionsToTest} option).
#' }
#'
#' In Claddis (>=0.3) these partitions are defined as a list of lists of vectors where only the first N - 1 partitions need be defined. E.g., if comparing the first edge value (based on \pkg{ape} numbering, i.e., \code{plot(tree); edgelabels()}) to the rest of the tree then the user only needs to define the value "1" and the function will automatically add a second partition containing all other edges. This can be set with the option \code{BranchPartitionsToTest = list(list(1))}. Similarly, to do what Lloyd et al. (2012) did and repeat the test for every edge in the tree (and assuming this variable is already named "tree") you could use, \code{BranchPartitionsToTest = lapply(as.list(1:nrow(tree$edge)), as.list)}.
#'
#' Because of the flexibility of this function the user can define any set of edges. For example, they could test whether terminal branches have a different rate from internal branches with \code{BranchPartitionsToTest = list(list(match(1:Ntip(tree), tree$edge[, 2])))}. The \code{CladePartitionsToTest} is really just a special subset of this type of hypothesis, but with edges being defined as descending from a specific internal node in the tree. Once again, an exploratory approach like that of Lloyd et al. (2012) can be used with: \code{CladePartitionsToTest = lapply(as.list(Ntip(tree) + (2:Nnode(tree))), as.list)}. Note that this excludes the root node as this would define a single partition and hence would represent the null hypothesis (a single rate model for the whole tree). (If using \code{LikelihoodTest = "AIC"} then the user typically \emph{will} want a value for a single partition.) More generally clades must be defined by the node numbers they correspond to. In R an easy way to identify these is with: \code{plot(tree); nodelabels()}.
#'
#' Time bin partitions are defined in a similar way, but are numbered 1:N starting from the oldest time bin. So if wanting to do an exploratory test of single bin partitions (and only four time bins were specified) you could use: \code{TimeBinPartitionsToTest = lapply(as.list(1:4), as.list)}. Bins can be combined too, just as edges are above. For example, time bins 1 and 2 could form a single partition with: \code{TimeBinPartitionsToTest = list(list(1:2))}. Or if looking to test a model where each bin has its' own rate value you could use: \code{TimeBinPartitionsToTest = list(as.list(1:3))}. Note, as before we do not need to specify the fourth bin as this will be automatically done by the function, however, \code{TimeBinPartitionsToTest = list(as.list(1:4))} will also work. Some caution needs to be applied with N-rate models (where N is three or larger) and \code{LikelihoodTest = "LRT"} as a result favouring such models does not necessarily endorse N-separate rates. I.e., it could simply be that one bin has such a large excursion that overall the N-rate model fits better than the 1-rate model, but some 2-rate models might be better still. It is up to the user to check this themselves by exploring smaller combinations of bins and more genrally if exploring partitions of three or more use of the Akaike Information Criterion (\code{LikelihoodTest = "AIC"}) is recommended.
#'
#' Finally, character partitions allow the user to explore whether rates vary across different character types (numbers), e.g., skeletal characters versus soft tissue characters, or cranial characters versus postcranial characters. Here characters are simply numbered 1:N (across all blocks of a matrix), but single character partitions are less likely to be of interest. As an example of use lets say the first ten characters are what we are interested in as a partition (the second partition being the remaining characters), we could use: \code{CharacterPartitionsToTest = list(list(1:10))} to test for a two-rate model with \code{LikelihoodTest = "LRT"}.
#'
#' Note that the list of lists structure is critical to defining partitions as it allows them to be of different sizes and numbers. For example, one partition of three and another of six, or one set of two partitions and another set of four partitions - structures not easily realizable using vectors or matrices. However, it may not be intuitive to some users so it is recommended that the user refers to the examples above as a guide.
#'
#' Additionally, it should be noted that the user can test multiple types of hypotheses simultaneously with the function. For example, performing several branch tests whilst also performing clade tests. However, it is not necessary to perform all types simultaneously (as was the case up to version 0.2) and unused partition types can be set to NULL, the default in each case.
#'
#' \bold{AIC vs LRT}
#'
#' Since Claddis version 0.4 the option to use the Akaike Information Criterion (AIC) instead of likelihood ratio tests (LRTs) has been added (although it was not properly functional until version 0.4.2). Practically speaking the AIC uses something similar to the denominator term from the LRT (see equation above) and adds a penalty for the number of parameters (partitions). However, it also fundamentally alters the comparative framework applied and hence needs more careful attention from the user to be applied correctly. Specifically, the LRT is by its nature comparative, always comparing an N-rate partition with a one-rate partition. By contrast the AIC does not directly apply any comparison and so the user must logically supply multiple partitionings of the data in order for the results to be meaningful. It might be assumed that a user will always want to apply a single partition that pools all the data for each type of test, whether this is all the edges (branches), time bins, or characters. This will thus be the obvious comparator for any multiple partition supplied, ensuring that any more complex partitioning is minimally superior to this. Additionally, it is also logical to consider each possible way of joining partitions simpler than the most complex partition being considered. E.g., if considering a four-partition model then the user should also consider all possible three-partition and two-partition combinations of that four-partition model. This can obviously lead to some complexity in supplying partitions on the user's part and so some automating of this process is planned in future (but is not fully available yet). For an example, see \link{partition_time_bins}.
#'
#' Additionally, AIC values are not simple to compare as there is no direct equivalent of the alpha value from the LRT. Instead the user can modify the AIC values returned themselves to get delta-AIC or Akaike weight values (e.g., with \code{geiger::aicw}). (NB: I will not explain these here as there are better explanations online.) Furthermore, since version 0.4.2 sample-size corrected AIC (AICc) is also available in the output. Note that some caution should be used in applying this if the number of partitions is equal to the sample size or only one fewer. E.g., if you have ten time bins and ten partitions you can get a negative value due to the denominator term in the AICc calculation. Thus it is advised to use the raw AIC values as first approximation and be wary if the AICc flips the preferred model to a more complex model or models (i.e., those with more partitions) as this is the opposite of the intent of the AICc.
#'
#' \bold{High versus low rates}
#'
#' Prior to Claddis version 0.3, rate results were classified as significantly high or significntly low as part of the output. This was done simply on whether the estimated p-value fell above or below the corrected alpha level (the significance threshold) and whether the first part of the two-rate partition had a higher or lower rate than the second part (the pooling of all other values). This simple interpretation is no longer valid here as the function can consider more than two partitions (high versus low is meaningless) and the allowing of AIC values means a significance test need not be performed either. Although the same interepretation can still be applied manually when only using two-partition tests and the LRT, other partition sizes and use of the AIC complicate this interpretation (in the same way an ANOVA is more complex than a two-sample t-test). This will also affect visualisation of the data (i.e. the simple pie chart coluring of non-significant, significntly high or significantly low rates seeen since Lloyd et al. 2012 will no longer apply). Instead the user should isolate the best model(s) and attempt to visualise these, perhaps using something like a heat map, with the mean rate for each partition being represented by an appropriate colour.
#'
#' \bold{Other options}
#'
#' Since Claddis version 0.3 this function has allowed the user greater control with many more options than were offered previously and these should be considered carefully before running any tests.
#'
#' Firstly, the user can pick an option for \code{ChangeTimes} which sets the times character changes are inferred to occur. This is only relevant when the user is performing time bin partition test(s) as this requires some inference to be made about when changes occur on branches that may span multiple time bins. The current options are: \code{"midpoint"} (all changes are inferred to occur midway along the branch, effectively mimicking the approach of Ruta et al. 2006), \code{"spaced"} (all changes are inferred to occur equally spaced along the branch, with changes occurring in character number order), or \code{"random"} (changes are assigned a random time by drawing from a uniform distribution between the beginning and end of each branch). The first of these is likely to lead to unrealistically "clumped" changes and by extension implies completely correlated character change that would violate the assumptions of the Poisson distribution that underlies the significance tests here (Lloyd et al. 2012). At the other extreme, the equally spaced option will likely unrealistically smooth out time series and potentially make it harder to reject the single-rate null (leading to Type II errors). For these reasons, the random option is recommended and is set as the default. However, because it is random this makes the function stochastic (the answer can vary each time it is run) and so the user should therefore run the function multiple times if using this option (i.e., by using a for loop) and aggregating the results at the end (e.g., as was done by previous authors; Lloyd et al. 2012; Close et al. 2015).
#'
#' Secondly, the \code{alpha} value sets the significance threshold by which the likelihood ratio test's resulting p-value is compared (i.e., it is only reelevant when \code{LikelihoodTest = "LRT"}. Following Lloyd et al. (2012) this is set lower (0.01) than the standard 0.05 value by default as those authors found rates to be highly heterogenous in their data set (fossil lungfish). However, this should not be adopted as a "standard" value without question (just as 0.05 shouldn't). Note that the function also corrects for multiple comparisons (using the \code{MultipleComparisonCorrection} option) to avoid Type I errors (false positives). It does so (following Lloyd et al. 2012) using the Benjamini-Hochberg (Benjamini and Hochberg 1995) False Discovery Rate approach (see Lloyd et al. 2012 for a discussion of why), but the Bonferroni correction is also offered here (albeit not recommended).
#'
#' Thirdly, polymorphisms and uncertainities create complications for assessing character changes along branches. These can occur at the tips (true polymorphisms or uncertainties in sampled taxa) and internal nodes (uncertainty over the estimated ancestral state). There are two options presented here, and applicable to both \code{PolymorphismState} and \code{UncertaintyState} (allowing these to be set separately). These are to convert such values to missing (NA) or to pick one of the possible states at random. Using missing values will increase overall uncertainty and potentially lead to Type II errors (false negatives), but represents a conservative solution. The random option is an attempt to avoid Type II errors, but can be considered unrealistic, especially if there are true polymorphisms. Additionally, the random option will again make the function stochastic meaning the user should run it multiple times amd aggregate the results. Note that if there are no polymorphisms or uncertainties in the character-taxon matrix the latter can still occur with ancestral state estimates, especially if the threshold value is set at a high value (see \link{estimate_ancestral_states} for details).
#'
#' Fourthly, inapplicable characters can additionally complicate matters as they are not quite the same as missing data. I.e., they can mean that change in a particular character is not even possible along a branch. However, there is no easy way to deal with such characters at present so the user is not presented with a true option here - currently all inapplicable states are simply converted to missing values by the function. In future, though, other options may be available here. For now it is simply noted that users should be careful in making inferences if there are inapplicable characters in their data and should perhaps consider removing them with \link{prune_cladistic_matrix} to gauge their effect.
#'
#' Fifthly, there are currenty two further options for assessing rates across time bins. As noted above a complication here is that character changes (the rate numerator) and character completeness (part of the rate denominator) are typically assessed on branches. However, branches will typically span time bin boundaries and hence many bins will contain only some portion of particular branches. The exact portion can be easily calculated for branch durations (the other part of the rate denominator) and the \code{ChangeTimes} option above is used to set the rate numerator, however, completeness remains complex to deal with. The first attempt to deal with this was made by Close et al. (2015) who simply used weighted mean completeness by calculating the proportion of a branch in each bin as the weight and multiplying this by each branch's completeness (the \code{"Close"} option here). However, this may lead to unrealistic "smoothing" of the data and perhaps more importantly makes no account of which characters are known in a bin. Lloyd (2016) proposed an alternative "subtree" approach which assesses completeness by considering each character to be represented by a subtree where only branches that are complete are retained then branch durations in each bin are summed across subtrees such that the duration term automatically includes completeness (the \code{"Lloyd"} option here). Here the latter is strongly recommended, for example, because this will lead to the same global rate across the whole tree as the branch, clade or character partitions, whereas the Close approach will not.
#'
#' Sixthly, all character changes are weighted according to the weights provided in the input character-taxon matrix. In many cases these will simply all be one, although see the equalise weights option in \link{read_nexus_matrix}. However, when weights vary they can create some issues for the function. Specifically, changes are expected to be in the same (integer) units, but if weights vary then they have to be modelled accordingly. I.e., a character twice the weight of another may lead to a single change being counted as two changes. This is most problematic when the user has continuous characters which are automatically converted to gap-weighted (Thiele 1993) characters. However, this conversion creates drastically down-weighted characters and hence the user may wish to set the \code{EnsureAllweightsAreIntegers} option to TRUE. Note that reweighting will affect the results and hence shifting the weights of characters up or down will necessarily lead to shifts in the relative Type I and II errors. This is an unexplored aspect of such approaches, but is something the user should be aware of. More broadly it is recommended that continuous (or gap-weighted) characters be avoided when using this function.
#'
#' Finally, the remaining options (\code{estimate.all.nodes}, \code{estimate.tip.values}, \code{inapplicables.as.missing}, \code{polymorphism.behaviour}, \code{uncertainty.behaviour}, and \code{threshold}) are all simply passed directly to \link{estimate_ancestral_states} for estimating the ancestral states and users should consult the help file for that function for further details.
#'
#' Note that currently the function cannot deal with step matrices and that the terminal versus internal option from Brusatte et al. (2014) is yet to be implemented.
#'
#' \bold{Output}
#'
#' The output for each LRT test (i.e., the components of the \code{BranchPartitionResults}, \code{CharacterPartitionResults}, \code{CladePartitionResults} and \code{TimeBinResults} parts of the output) includes three main parts:
#'
#' \enumerate{
#'   \item Rates.
#'   \item PValue.
#'   \item CorrectedAlpha.
#' }
#'
#' Or for each AIC test there are:
#'
#' \enumerate{
#'   \item Rates.
#'   \item AIC.
#'   \item AICc.
#' }
#'
#' For each rate test the \code{Rates} part of the output is a vector of the absolute rate (number of changes per million years) for each partition in the test (in the order they were supplied to the function). So, for example, a branch rate for the sixth edge in a tree would be the rate for the sixth edge followed by the pooled rate for all other edges. The length of the vector is the length of the number of partitions.
#'
#' The PValue is a single value indicating the probability that the likelihood ratio (see above and Lloyd et al. 2012) is one, i.e., that the likelihoods of the one-rate and N-rate models are the same.
#'
#' The CorrectedAlpha is the alpha-value that should be used to determine the significance of the current partition test (i.e., The PValue, above). If the PValue exceeds the CorrectedAlpha then the null (single-rate) hypothesis should be accepted, if lower then the null should be rejected in favour of the N-rate hypthesis. Note that the CorrectedAlpha will not typically be the same for each partition and will also typically be different from the input \code{alpha} value due to the \code{MultipleComparisonCorrection} option used.
#'
#' The AIC is the Akaike Information Criterion, and is relatively meaningless on its own and can only really be used to compare with the AIC values for other partitions of the data. The AICc is simply the sample-size corrected version of the AIC and is preferable when sample sizes are small.
#'
#' @return
#'
#' \item{TimeBinsUsed}{The time binning used (NB: May be slightly altered from the input values).}
#' \item{InferredCharacterChanges}{Matrix of inferred character changes.}
#' \item{IntrinsicCharacterRate}{The intrinsic (global) character rate in changes per million years.}
#' \item{ContinuouscharactersConvertedToDiscrete}{Whether or not continuous characters were converted to discrete characters (important for handling the data in downstream analys(es)).}
#' \item{BranchPartitionResults}{List of branch partition results (corresponding to \code{BranchPartitionsToTest}. NULL if not requested.}
#' \item{CharacterPartitionResults}{List of character partition results (corresponding to \code{CharacterPartitionsToTest}. NULL if not requested.}
#' \item{CladePartitionResults}{List of clade partition results (corresponding to \code{CladePartitionsToTest}. NULL if not requested.}
#' \item{TimeBinResults}{List of time bin partition results (corresponding to \code{TimeBinPartitionsToTest}. NULL if not requested.}
#' \item{BranchRates}{Matrix showing calculated rates for each branch. NULL if \code{BranchPartitionsToTest} is not requested.}
#' \item{CharacterRates}{Matrix showing calculated rates for each character. NULL if \code{CharacterPartitionsToTest} is not requested.}
#' \item{CladeRates}{Matrix showing calculated rates for each clade. NULL if \code{CladePartitionsToTest} is not requested.}
#' \item{TimeRates}{Matrix showing calculated rates for each time bin. NULL if \code{TimeBinPartitionsToTest} is not requested.}
#' \item{Tree}{The time-scaled input tree used as input (provided as output for use with visualisation functions).}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Steve C. Wang \email{scwang@@swarthmore.edu}
#'
#' @references
#'
#' Benjamini, Y. and Hochberg, Y., 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society, Series B}, \bold{57}, 289-300.
#'
#' Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. \emph{Current Biology}, \bold{24}, 2386-2392.
#'
#' Close, R. A., Friedman, M., Lloyd, G. T. and Benson, R. B. J., 2015. Evidence for a mid-Jurassic adaptive radiation in mammals. \emph{Current Biology}, \bold{25}, 2137-2142.
#'
#' Cloutier, R., 1991. Patterns, trends, and rates of evolution within the Actinistia. \emph{Environmental Biology of Fishes}, \bold{32}, 23–58.
#'
#' Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. \emph{Biological Journal of the Linnean Society}, \bold{118}, 131-151.
#'
#' Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying heterogeneity in rates of morphological evolution: discrete character change in the evolution of lungfish (Sarcopterygii; Dipnoi). \emph{Evolution}, \bold{66}, 330-348.
#'
#' Ruta, M., Wagner, P. J. and Coates, M. I., 2006. Evolutionary patterns in early tetrapods. I. Rapid initial diversification followed by decrease in rates of character change. \emph{Proceedinsg of the Royal Society of London B}, \bold{273}, 2107–2111.
#'
#' Thiele, K.. 1993. The Holy Grail of the perfect character: the cladistic treatment of morphometric data. \emph{Cladistics}, \bold{9}, 275-304.
#'
#' @examples
#' 
#' # Set random seed:
#' set.seed(17)
#' 
#' # Generate a random tree for the Michaux data set:
#' tree <- rtree(nrow(Michaux1989$matrix_1$Matrix))
#' 
#' # Update taxon names to match those in the data matrix:
#' tree$tip.label <- rownames(Michaux1989$matrix_1$Matrix)
#' 
#' # Set root time by making youngest taxon extant:
#' tree$root.time <- max(diag(vcv(tree)))
#' 
#' # Get discrete character rates:
#' x <- test_rates(tree = tree, cladistic.matrix =
#'   Michaux1989, TimeBins = seq(from = tree$root.time,
#'   to = 0, length.out = 5), BranchPartitionsToTest =
#'   lapply(as.list(1:nrow(tree$edge)), as.list),
#'   CharacterPartitionsToTest = lapply(as.list(1:3),
#'   as.list), CladePartitionsToTest =
#'   lapply(as.list(Ntip(tree) + (2:Nnode(tree))),
#'   as.list), TimeBinPartitionsToTest =
#'   lapply(as.list(1:4), as.list), ChangeTimes =
#'   "random", alpha = 0.01, PolymorphismState =
#'   "missing", UncertaintyState = "missing",
#'   InapplicableState = "missing", TimeBinApproach =
#'   "Lloyd")
#'
#' @export test_rates
test_rates <- function(tree, cladistic.matrix, TimeBins, BranchPartitionsToTest = NULL, CharacterPartitionsToTest = NULL, CladePartitionsToTest = NULL, TimeBinPartitionsToTest = NULL, ChangeTimes = "random", LikelihoodTest = "AIC", alpha = 0.01, MultipleComparisonCorrection = "BenjaminiHochberg", PolymorphismState = "missing", UncertaintyState = "missing", InapplicableState = "missing", TimeBinApproach = "Lloyd", EnsureAllweightsAreIntegers = FALSE, estimate.all.nodes = FALSE, estimate.tip.values = FALSE, inapplicables.as.missing = FALSE, polymorphism.behaviour = "equalp", uncertainty.behaviour = "equalp", threshold = 0.01, allow.all.missing = FALSE) {
  
  # ADD EXAMPLES OF VISUALISED OUTPUT
  
  # AICc BREAKS IF MORE THAN N-2 PARAMETERS (INFINITY OR WRONGLY NEGATIVE OUTPUT CAN OCCUR THIS WAY)
  # GLOBAL RATE NEEDS TO GO INTO LRT OPTION IF NOT USED IN AIC
  # SOMEHOW NEED TO ALLOW MULTIPLE VERSIONS OF MAIN PIEPLINE IF DOING RANDOM ASSIGNMENTS OF CHANGE TIMES
  # NEED TO ADD PARTITIONS USED TO OUTPUT SOMEHOW...
  
  # MAYBE DO SOME KIND OF WEIGHTING FOR AIC AS SOME PARTITIONS WILL CONTAIN VERY LITTLE DATA AND BE EXTREME OUTLIERS?
  # NEED TO CHECK FOR SINGLE PARTITION WITH LRT (ALLOWED WITH AIC)
  # NEED TO CHECK FOR FULLY BIFURCATING IF IMPLEMENTING WANG STUDENTS APPROACH? OR IS THAT DONE ANYWAY?
  # IF USING AIC NEED TO CHECK FOR EACH TEST TYPE AT LEAST TWO PARTITIONS ARE SUPPLIED OR ELSE THE RESULTS ARE MEANINGLESS
  # ADD CHECK TO INCLUDE ALL LESS COMPLEX SUBPARTITIONS OF ANY 3 OR MORE SIZE PARTITOINING OF THE DATA
  
  # DESIDERATA (STUFF IT WOULD BE NICE TO ADD IN FUTURE):
  #
  # WRITE SEARCH VERSION FOR FINDING RATE SHIFTS? SHOULD THIS EVEN BE AN OPTION? DOES THIS REQUIRE MODIFYING LRT TO COMPARE E.G. 2-RATE DIRECTLY WITH 3-RATE MODEL? WOULD NEED TO PERMUTE ALL POSSIBLE COMBOS AND NOT SURE HOW LARGE THESE MIGHT GET (VERY FOR EDGES).
  # MAYBE MAKE ANCESTRAL STATE UNCERTAINTY DIFFERENT FOR TIPS THAN NODES? I.E., HOW IT GETS RESOLVED CAN BE DIFFERENT (MORE OPTIONS TO FUNCTION)
  # THESE TWO ARE RELATED: 1. ADD TERMINAL VERSUS INTERNAL OPTION SOMEHOW/SOMEWHERE (BRANCH TYPE ALREADY RECORDED ON EDGE LIST!), 2. ALLOW OPTION TO IGNORE SOME PARTS OF THE TREE FOR PARTITION TESTS? MAKES CALCULATING THE MEAN RATE TRICKIER BUT MIGHT MAKE SENSE E.G. FOR INGROUP ONLY TESTS. EXCLUDE EDGES AFTER DOING ANCESTRAL STATES? OR SET THESE TO ALL NAS TO ENSURE OTHER THINGS WORK FINE? FOR EXAMPLE, USE OUTGROUPS TO SET ANCESTOR THEN EXCLUDE THEM FROM THE PROCESS.
  # EXTRA FUNCTION(S) TO VISUALISE RESULTS MOST LIKELY. DEFO NEEDED! HEAT MAP WITH EDGE BLOCKS?
  # CHECK FOR AUTAPOMORPHIES AND INFORM USER IF FOUND?
  # ADD CONTRIVED EXAMPLES (UNIT TESTS) TO SHOW HOW FUNCTION WORKS, E.G. RATE OF ONE CHANGES PER MILLION YEARS THEN DUPLICATED BLOCK WITH CHARACTER PARTITION TEST.
  # PROBABLY NEED MORE CAREFUL CHECKS FOR ZERO VALUES GENERALLY, E.G., CHARACTER WITH ALL MISSING DATA
  # ALLOW REWEIGHTING OF INAPPLICABLES ZERO AS AN OPTION FOR THEM?
  # HOW TO FORMAT OUTPUT? GET CIS FOR EACH PARTITION FOR VISUALISATION (E.G., BARPLOT OF PARTITION VALUES WITH DASHED LINE FOR MEAN AND ERROR BARS FOR CIS)? STICK WITH LIST OR COLLAPSE TO A MATRIX SOMEHOW?
  # TIME BINS WITH NOTHING IN WILL CAUSE ISSUES AS DIVIDE BY ZEROES WILL OCCUR - ADD CHECK FOR THIS.
  # WHAT IS SIGNIFICANTLY HIGH OR LOW IF THERE ARE THREE OR MORE PARTITIONS? THIS IS NOT EVEN IN OUTPUT YET. PROLLY CANNOT DO FULL STOP NOW PARTITIONS ARE MORE COMPLEX

  # Check for step matrices and stop and warn user if found:
  if (is.list(cladistic.matrix$topper$step_matrices)) stop("Function cannot currently deal with step matrices.")

  # Check tree has branch lengths:
  if (is.null(tree$edge.length)) stop("Tree does not have branch lengths (durations). Try timescaling the tree, e.g., with DatePhylo.")
  
  # Check tree has root age:
  if (is.null(tree$root.time)) stop("Tree is missing $root.time. Try setting this before continuing, e.g., tree$root.time <- 104.2.")
  
  # Check ChangeTimes is correctly formatted or stop and warn user:
  if (length(setdiff(ChangeTimes, c("midpoint", "spaced", "random"))) > 0) stop("ChangeTimes must be one of \"midpoint\", \"spaced\", or \"random\".")
  
  # Check MultipleComparisonCorrection is correctly formatted or stop and warn user:
  if (length(setdiff(MultipleComparisonCorrection, c("BenjaminiHochberg", "Bonferroni"))) > 0) stop("MultipleComparisonCorrection must be one of \"BenjaminiHochberg\" or \"Bonferroni\".")
  
  # Check PolymorphismState is correctly formatted or stop and warn user:
  if (length(setdiff(PolymorphismState, c("missing", "random"))) > 0) stop("PolymorphismState must be one of \"missing\" or \"random\".")
  
  # Check UncertaintyState is correctly formatted or stop and warn user:
  if (length(setdiff(UncertaintyState, c("missing", "random"))) > 0) stop("UncertaintyState must be one of \"missing\" or \"random\".")
  
  # Check InapplicableState is correctly formatted or stop and warn user:
  if (length(setdiff(InapplicableState, c("missing"))) > 0) stop("InapplicableState must be \"missing\".")
  
  # Check TimeBinApproach is correctly formatted or stop and warn user:
  if (length(setdiff(TimeBinApproach, c("Close", "Lloyd"))) > 0) stop("TimeBinApproach must be one of \"Close\" or \"Lloyd\".")
  
  # Check partitions are not all NULL values:
  if (is.null(BranchPartitionsToTest) && is.null(CharacterPartitionsToTest) && is.null(CladePartitionsToTest) && is.null(TimeBinPartitionsToTest)) stop("No partitions are requested. Set at least one of BranchPartitionsToTest, CharacterPartitionsToTest, CladePartitionsToTest, or TimeBinPartitionsToTest to a list of appropriate values. Type \"?test_rates\" for help.")
  
  # Get internal node numbers:
  InternalNodeNumbers <- 1:ape::Nnode(tree) + ape::Ntip(tree)
  
  # Get edge numbers:
  EdgeNumbers <- 1:nrow(tree$edge)
  
  # Get character numbers:
  CharacterNumbers <- 1:sum(unlist(lapply(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "Matrix"), ncol)))
  
  # Ensure time bins are in correct order:
  TimeBins <- sort(unique(TimeBins), decreasing = TRUE)
  
  # Find the Time bin midpoints:
  TimeBinMidpoints <- (TimeBins[2:length(TimeBins)] + TimeBins[1:(length(TimeBins) - 1)]) / 2
  
  # Get the numbers for each time bins:
  TimeBinNumbers <- 1:length(TimeBinMidpoints)
  
  # Subfunction to ensure partitions are formatted correctly:
  format_partition <- function(PartitionsToTest, ValidValues, PartitionName) {
    
    # Check partitions are in the form of a list of lists:
    if (!all(c(all(unlist(lapply(PartitionsToTest, is.list))), is.list(PartitionsToTest)))) stop(paste(PartitionName, " must be in the form of a list of lists.", sep = ""))
    
    # Get a vector of any non-present valid values:
    NonPresentValues <- setdiff(unique(unlist(PartitionsToTest)), ValidValues)
    
    # Check valid values have been used and if not stop and warn user:
    if (length(NonPresentValues) > 0) stop(paste(PartitionName, "Partitions to test must be defined using the valid range of values (", paste(range(ValidValues), collapse = " to "), ") only.", sep = ""))
    
    # Check partitions never overlap and stop and warn user if they do:
    Check <- lapply(PartitionsToTest, function(x) if (any(duplicated(sort(unlist(x))))) stop(paste("Each partition of ", PartitionName, " must not contain overlapping values (e.g., can not have 1:3 and 3:5 as both contain 3).", sep = "")))
    
    # Subfunction to ad the missing partition (if exists):
    add_missing_partitions <- function(x, ValidValues) {
      
      # Define any missing values:
      missingValues <- setdiff(ValidValues, unlist(x))
      
      # If there are missing values add them to list at end:
      if (length(missingValues) > 0) x[[(length(x) + 1)]] <- missingValues
      
      # Return x:
      return(x)
      
    }
    
    # Add in missing partitions (if any):
    PartitionsToTest <- lapply(PartitionsToTest, add_missing_partitions, ValidValues = ValidValues)
    
    # Check partitions are all at least two in size or else no comparison can be made:
    if (any(unlist(lapply(PartitionsToTest, length)) == 1) && LikelihoodTest == "LRT") stop("Partitions must divide the available data into at least two parts if performing likelihood ratio tests.")

    # Return formatted partitions to test:
    return(PartitionsToTest)

  }
  
  # Subfunction to pack partitions to short format for output:
  pack_partitions <- function(FormattedPartitions) unlist(lapply(FormattedPartitions, function(x) paste(unlist(lapply(x, function(y) {
    
    # First make sure y is sorted:
    y <- sort(y)
    
    # Covnvert y to a list (splitting if gaps greater than 1 are found):
    y <- unname(split(y, cumsum(c(TRUE, diff(y) > 1))))
    
    # Collapse gaps of one with hyphens:
    paste0(unlist(lapply(y, function(z) {res <- as.character(z); if (length(z) > 1) {r <- rle(c(1, pmin(diff(z), 2))); res <- paste0(z[c(1, cumsum(r$lengths))], c("-", " ")[r$values], collapse = ""); res <- substr(res, 1, nchar(res) - 1)}; res})), collapse = " ")
    
  })), collapse = " | ")))
  
  # If performing branch partition test(s) check and reformat branch partitions:
  if (!is.null(BranchPartitionsToTest)) BranchPartitionsToTest <- format_partition(PartitionsToTest = BranchPartitionsToTest, ValidValues = EdgeNumbers, PartitionName = "BranchPartitionsToTest")
  
  # If performing character partition test(s) check and reformat character partitions:
  if (!is.null(CharacterPartitionsToTest)) CharacterPartitionsToTest <- format_partition(PartitionsToTest = CharacterPartitionsToTest, ValidValues = CharacterNumbers, PartitionName = "CharacterPartitionsToTest")
  
  # If performing clade partition test(s)
  if (!is.null(CladePartitionsToTest)) {
    
    # Convert clade partitions to edge partitions:
    CladePartitionsToTest <- lapply(CladePartitionsToTest, lapply, find_descendant_edges, tree = tree)
    
    # Check and reformat clade partitions:
    CladePartitionsToTest <- format_partition(PartitionsToTest = CladePartitionsToTest, ValidValues = EdgeNumbers, PartitionName = "CladePartitionsToTest")
    
  }

  # If performing time bin partition test(s) check and reformat time bin partitions:
  if (!is.null(TimeBinPartitionsToTest)) TimeBinPartitionsToTest <- format_partition(PartitionsToTest = TimeBinPartitionsToTest, ValidValues = TimeBinNumbers, PartitionName = "TimeBinPartitionsToTest")
  
  # Check LikelihoodTest is correctly formatted or stop and warn user:
  if (length(setdiff(LikelihoodTest, c("AIC", "LRT"))) > 0) stop("LikelihoodTest must be one of \"AIC\" or \"LRT\".")

  # Subfunction to calculate maximum likelihood p value:
  get_likelihood_p <- function(MeanRate, SampledRates, SampledChanges, SampledCompleteness, SampledTime) {
    
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
  
  # Subfunction to calculate AIC:
  calculate_AIC <- function(SampledRates, SampledChanges, SampledCompleteness, SampledTime) {
    
    # Get log maximum likelihood estimate:
    LogMLE <- sum(log(dpois(x = round(SampledChanges), lambda = SampledRates * SampledCompleteness * SampledTime)))
    
    # Calculate AIC:
    AIC <- (2 * length(SampledRates)) - (2 * LogMLE)
    
    # Return AIC:
    return(AIC)
    
  }
  
  # Subfunction to calculate AIC from partition (with columns labelled Partition, Rate, Completeness, Duration):
  calculate_partition_AIC <- function(Partition, AICc = FALSE) {
    
    # Get log maximum likelihood estimate:
    LogMLE <- sum(log(dpois(round(Partition[, "Changes"]), Partition[, "Rate"] * Partition[, "Completeness"] * Partition[, "Duration"])))
    
    # Get k (number of parameters) term:
    k <- max(Partition[, "Partition"])
    
    # Calculate AIC:
    AIC <- (2 * k) - (2 * LogMLE)
    
    # If AICc is desired then calculate this and overwrite AIC with it:
    if (AICc) AIC <- AIC + (((2 * (k ^ 2)) + (2 * k)) / (nrow(Partition) - k - 1))
    
    # Return AIC:
    return(AIC)
    
  }
 
  # Get ages for each (tip and internal) node:
  date_nodes <- date_nodes(tree)

  # Get branch ages (from and to):
  BranchAges <- unname(cbind(date_nodes[as.character(tree$edge[, 1])], date_nodes[as.character(tree$edge[, 2])]))

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
  find_descendant_edgesForEachInternalNode <- lapply(as.list(InternalNodeNumbers), find_descendant_edges, tree = tree)
  
  # Get ancestral character states:
  AncestralStates <- estimate_ancestral_states(cladistic.matrix = cladistic.matrix, time.tree = tree, estimate.all.nodes = estimate.all.nodes, estimate.tip.values = estimate.tip.values, inapplicables.as.missing = inapplicables.as.missing, polymorphism.behaviour = polymorphism.behaviour, uncertainty.behaviour = uncertainty.behaviour, threshold = threshold, allow.all.missing = allow.all.missing)
  
  # Build single matrix of all states in tip label then node number order:
  AllStates <- do.call(cbind, lapply(lapply(AncestralStates[2:length(AncestralStates)], '[[', "Matrix"), function(x) x[c(tree$tip.label, 1:ape::Nnode(tree) + ape::Ntip(tree)), , drop = FALSE]))
  
  # Make vector of ordering of characters:
  ordering <- unname(unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "ordering")))
  
  # Make vector of weights of characters:
  weights <- unname(unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "weights")))
  
  # Make vector of minimum values:
  MinVals <- unname(unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "MinVals")))
  
  # Make vector of maximum values:
  MaxVals <- unname(unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "MaxVals")))
  
  # Find positions in matrix with polymorphisms:
  PolymorphismPositions <- grep("&", AllStates)
  
  # Find positions in matrix with uncertainties:
  UncertaintyPositions <- grep("/", AllStates)
  
  # Find positions in matrix with inapplicables:
  InapplicablePositions <- which(AllStates == "")

  # If polymorphisms were found:
  if (length(PolymorphismPositions) > 0) {
    
    # If replacing polymorphsims with missing do so:
    if (PolymorphismState == "missing") AllStates[PolymorphismPositions] <- NA
    
    # If replacing polymorphisms with random values draw and replace:
    if (PolymorphismState == "random") AllStates[PolymorphismPositions] <- unlist(lapply(strsplit(AllStates[PolymorphismPositions], "&"), sample, size = 1))
    
  }
  
  # If uncertainties were found:
  if (length(UncertaintyPositions) > 0) {
    
    # If replacing uncertainties with missing do so:
    if (UncertaintyState == "missing") AllStates[UncertaintyPositions] <- NA
    
    # If replacing uncertainties with random values draw and replace:
    if (UncertaintyState == "random") AllStates[UncertaintyPositions] <- unlist(lapply(strsplit(AllStates[UncertaintyPositions], "/"), sample, size = 1))
    
  }
  
  # If inapplicable states were found:
  if (length(InapplicablePositions) > 0) {
    
    # If replacing inapplicables with missing do so:
    if (InapplicableState == "missing") AllStates[InapplicablePositions] <- NA
    
  }
  
  # Set default converted continuous characters to FALSE:
  ContinuouscharactersConverted <- FALSE
  
  # Check for continuous characters as these will need to be modified for modelling to discrete characters:
  if (any(ordering == "cont")) {
    
    # Tell user this is happening:
    cat("Continuous characters found. Converting to gap-weighted discrete characters.\n")
    
    # Set default converted continuous characters to TRUE:
    ContinuouscharactersConverted <- TRUE
    
    # Find out which characters are continuous:
    ContinuouscharactersFound <- which(ordering == "cont")
    
    # Rescale continous characters as zero to one values:
    ListOfContinuousValuesRescaledZeroToOne <- lapply(lapply(lapply(apply(AllStates[, ContinuouscharactersFound, drop = FALSE], 2, list), unlist), as.numeric), function(x) {x <- x - min(sort(x)); x <- x / max(sort(x)); return(x)})
    
    # Now discretize and store these characters (0 to 31 scale):
    AllStates[, ContinuouscharactersFound] <- do.call(cbind, lapply(lapply(lapply(ListOfContinuousValuesRescaledZeroToOne, function(x) as.list(x)), lapply, function(x) ifelse(is.na(x), NA, max(which(x >= (0:31) / 31)) - 1)), unlist))
    
    # Convert character type to ordered:
    ordering[ContinuouscharactersFound] <- "ord"
    
    # Convert weights to 1/31:
    weights[ContinuouscharactersFound] <- 1 / 31
    
    # Set minimum value to zero:
    MinVals[ContinuouscharactersFound] <- 0
    
    # Set maximum value to 31:
    MaxVals[ContinuouscharactersFound] <- 31
    
  }
  
  # If EnsureAllweightsAreIntegers is TRUE rescale weights until they are all integers so can model appropriately with Poisson later:
  if (EnsureAllweightsAreIntegers) while(is.character(all.equal(sum(weights %% 1), 0))) weights <- (1 / (weights %% 1)[(weights %% 1) > 0])[1] * weights

  # Add from-to node states for each character to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$CharacterStatesFromTo <- matrix(AllStates[x$NodeNumberFromTo, , drop = FALSE], nrow = 2, dimnames = list(c("From", "To"))); return(x)})
  
  # Subfunction to define character changes:
  build_changes_matrix <- function(x) {
    
    # Find only comparable characters (those scored for both from and to states):
    Comparablecharacters <- which(apply(!apply(x$CharacterStatesFromTo, 2, is.na), 2, all))
    
    # Isolate comparable ordering:
    Comparableordering <- ordering[Comparablecharacters]
    
    # Isolate comparable weights:
    Comparableweights <- weights[Comparablecharacters]
    
    # Isolate only characters that actually differ (change):
    CharacterDifferences <- which(x$CharacterStatesFromTo["From", Comparablecharacters] != x$CharacterStatesFromTo["To", Comparablecharacters])
    
    # Build character change matrix:
    CharacterChanges <- matrix(nrow = 0, ncol = 5, dimnames = list(c(), c("Character", "From", "To", "Steps", "Weight")))
    
    # If characters change then make a matrix from them:
    if (length(CharacterDifferences) > 0) CharacterChanges <- rbind(CharacterChanges, cbind(as.numeric(Comparablecharacters[CharacterDifferences]), as.numeric(x$CharacterStatesFromTo["From", Comparablecharacters[CharacterDifferences]]), as.numeric(x$CharacterStatesFromTo["To", Comparablecharacters[CharacterDifferences]]), ifelse(Comparableordering[CharacterDifferences] == "unord", 1, abs(as.numeric(x$CharacterStatesFromTo["To", Comparablecharacters[CharacterDifferences]]) - as.numeric(x$CharacterStatesFromTo["From", Comparablecharacters[CharacterDifferences]]))), Comparableweights[CharacterDifferences]))
    
    # Store character changes as new sublist for x:
    x$CharacterChanges <- CharacterChanges
    
    # Store comparable characters as new sublist of x:
    x$Comparablecharacters <- Comparablecharacters
    
    # Return x:
    return(x)
    
  }
  
  # Get character changes and comparable characters and add to edge list:
  EdgeList <- lapply(EdgeList, build_changes_matrix)
  
  # Check whether time bins are being compared (otherwise no need to assign character changes):
  if (!is.null(TimeBinPartitionsToTest)) {
    
    # Subfunction to add change times to character changes:
    add_change_times <- function(x, ChangeTimes) {
      
      # Isolate character changes:
      CharacterChanges <- x$CharacterChanges
      
      # If any changes involve two or more steps (requiring replacement with multiple changes):
      if (any(CharacterChanges[, "Steps"] > 1)) {
        
        # Get multistep character changes:
        MultiStepcharacters <- which(CharacterChanges[, "Steps"] > 1)
        
        # For each multistep character change:
        for(i in rev(MultiStepcharacters)) {
          
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
      if (ChangeTimes == "midpoint") CharacterChanges <- cbind(CharacterChanges, rep(x$NodeAgeFromTo[1] - (x$BranchDuration / 2), length.out = nrow(CharacterChanges)))
      
      # If using spaced then set character change times as equally spaced along branch:
      if (ChangeTimes == "spaced") CharacterChanges <- cbind(CharacterChanges, x$NodeAgeFromTo[1] - (seq(from = 0, to = x$BranchDuration, length.out = nrow(CharacterChanges) + 1)[1:nrow(CharacterChanges)] + (diff(seq(from = 0, to = x$BranchDuration, length.out = nrow(CharacterChanges) + 1))[1] / 2)))
      
      # If using random then set character change times as random draws from a uniform distribution:
      if (ChangeTimes == "random") CharacterChanges <- cbind(CharacterChanges, x$NodeAgeFromTo[1] - stats::runif(n = nrow(CharacterChanges), min = 0, max = x$BranchDuration))
      
      # Add column name to change time column:
      colnames(CharacterChanges)[ncol(CharacterChanges)] <- "Time"
      
      # Subfunction to re-sort character change times so they occur in correct order:
      sort_change_times <- function(CharacterChanges) {
        
        # Sort change time for each character from oldest (first) to youngest (last) and store it:
        CharacterChanges[, "Time"] <- unname(unlist(lapply(as.list(unique(CharacterChanges[, "Character"])), function(x) sort(CharacterChanges[which(CharacterChanges[, "Character"] == x), "Time"], decreasing = TRUE))))
        
        # Return sorted character changes:
        return(CharacterChanges)
        
      }
      
      # Re-sort character change times so they occur in correct order:
      CharacterChanges <- sort_change_times(CharacterChanges)
      
      # Add bin for character change as last column:
      CharacterChanges <- cbind(CharacterChanges, unlist(lapply(as.list(CharacterChanges[, "Time"]), function(x) max(which(x <= TimeBins)))))
      
      # Add column name to change time column:
      colnames(CharacterChanges)[ncol(CharacterChanges)] <- "Bin"
      
      # Overwrite character changes with new version with changes added:
      x$CharacterChanges <- CharacterChanges
      
      # Return x:
      return(x)
    
    }
    
    # Add character change times to edge list:
    EdgeList <- lapply(EdgeList, add_change_times, ChangeTimes = ChangeTimes)
    
  }
  
  # Subfunction to get edge sections in time bins:
  get_edge_sections_in_bins <- function(x, TimeBins = TimeBins) {
    
    # Set first appearance datum of edge:
    FAD <- x$NodeAgeFromTo[1]
    
    # Set last appearance datum of edge:
    LAD <- x$NodeAgeFromTo[2]
    
    # Get any time bin boundaries crossed (can be empty if none are):
    BoundariesCrossed <- TimeBins[2:(length(TimeBins) - 1)][intersect(which(TimeBins[2:(length(TimeBins) - 1)] > LAD), which(TimeBins[2:(length(TimeBins) - 1)] < FAD))]
    
    # If boundaries are crossed:
    if (length(BoundariesCrossed) > 0) {
      
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
  EdgeList <- lapply(EdgeList, get_edge_sections_in_bins, TimeBins = TimeBins)
  
  # Add binned branch durations to edge list:
  EdgeList <- lapply(EdgeList, function(x) {BranchDurations <- rep(0, length(TimeBins) - 1); BranchDurations[as.numeric(colnames(x$BinnedEdgeSections))] <- abs(apply(x$BinnedEdgeSections, 2, diff)); x$BinnedBranchDurations <- BranchDurations; return(x)})
  
  # Add proportional binned branch lengths to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$ProportionalBinnedEdgeDurations <- x$BinnedBranchDurations / sum(x$BinnedBranchDurations); return(x)})
  
  # Start to build matrix of all changes with list of character changes:
  AllChanges <- lapply(EdgeList, function(x) x$CharacterChanges)
  
  # Add edge number to each matrix of character changes:
  for(i in 1:length(AllChanges)) AllChanges[[i]] <- cbind(rep(i, times = nrow(AllChanges[[i]])), AllChanges[[i]])
  
  # Combine all changes into a single matrix:
  AllChanges <- do.call(rbind, lapply(AllChanges, function(x) {colnames(x)[1] <- "Edge"; x}))
  
  # Remove silly rownames from all changes:
  rownames(AllChanges) <- NULL
  
  # Create NULL output variables (to be overwritten if called):
  BranchRates <- CladeRates <- TimeRates <- CharacterRates <- NULL

  # If doing some kind of edge test (branch or clade):
  if (!is.null(BranchPartitionsToTest) || !is.null(CladePartitionsToTest)) {
    
    # Get (weighted) number of changes on each edge:
    EdgeChanges <- unlist(lapply(EdgeList, function(x) sum(x$CharacterChanges[, "Steps"] * x$CharacterChanges[, "Weight"])))

    # Get completeness for each edge:
    EdgeCompleteness <- unlist(lapply(EdgeList, function(x) sum(weights[x$Comparablecharacters]) / sum(weights)))
    
    # Get duration of each edge:
    EdgeDurations <- unlist(lapply(EdgeList, function(x) x$BranchDuration))
    
    # Set global rate:
    GlobalRate <- sum(EdgeChanges) / sum(EdgeCompleteness * EdgeDurations)
    
    # If performing branch partition tests:
    if (!is.null(BranchPartitionsToTest)) {
      
      # Create branch rates for output:
      BranchRates <- lapply(list(as.list(1:nrow(tree$edge))), function(x) matrix(unlist(lapply(x, function(y) c(sum(EdgeChanges[y]), sum(EdgeCompleteness[y] * EdgeDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))[[1]]
      
      # Add edge numbers and rates:
      BranchRates <- cbind(Edge = 1:nrow(tree$edge), Rate = as.numeric(gsub(NaN, 0, BranchRates[, "Changes"] / BranchRates[, "Completeness"])), BranchRates)
      
      # If using Likelihood Ratio Test:
      if (LikelihoodTest == "LRT") {
        
        # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
        PartitionedData <- lapply(BranchPartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(EdgeChanges[y]), sum(EdgeCompleteness[y] * EdgeDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
        
        # Add sampled rate to paritioned data matrices:
        PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
        
        # Get LRT p-values and combine output as edge test results:
        BranchPartitionTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], get_likelihood_p(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
        
      }
      
      # If using AIC:
      if (LikelihoodTest == "AIC") {
        
        # Build partitioned data for AIC:
        PartitionedData <- lapply(BranchPartitionsToTest, function(x) {y <- cbind(Partition = rep(NA, length(tree$edge.length)), Rate = rep(NA, length(tree$edge.length)), Changes = EdgeChanges, Completeness = EdgeCompleteness, Duration = EdgeDurations); y[, "Rate"] <- as.numeric(gsub(NaN, 0, unlist(lapply(x, function(x) rep(sum(y[x, "Changes"]) / (sum(y[x, "Completeness"]) * sum(y[x, "Duration"])), length(x))))[order(unlist(x))])); y[, "Partition"] <- rep(1:length(x), unlist(lapply(x, length)))[order(unlist(x))]; y})
        
        # Get AIC, AICc and rate results:
        BranchPartitionTestResults <- lapply(PartitionedData, function(x) list(Rates = unname(unlist(lapply(as.list(unique(x[, "Partition"])), function(y) x[x[, "Partition"] == y, "Rate"][1]))), AIC = calculate_partition_AIC(x), AICc = calculate_partition_AIC(x, AICc = TRUE)))
        
      }
      
      # Pack branch partitions to test into single strings for output:
      PackedBranchPartitions <- pack_partitions(BranchPartitionsToTest)
      
      # Add packed partitions to results:
      for(i in 1:length(BranchPartitionTestResults)) BranchPartitionTestResults[[i]]$Partition <- PackedBranchPartitions[i]
    
    # If not performing branch partition tests:
    } else {
      
      # Make empty branch partition result output:
      BranchPartitionTestResults <- NULL
      
    }
    
    # If performing clade partition tests:
    if (!is.null(CladePartitionsToTest)) {
      
      # Create clade rates for output:
      CladeRates <- do.call(rbind, lapply(lapply(as.list(ape::Ntip(tree) + c(1:tree$Nnode)), function(x) find_descendant_edges(x, tree = tree)), function(y) c(Changes = sum(EdgeChanges[y]), Completeness = sum(EdgeCompleteness[y] * EdgeDurations[y]), Duration = 1)))
      
      # Add rates and node numbers:
      CladeRates <- cbind(Node = ape::Ntip(tree) + c(1:tree$Nnode), Rate = as.numeric(gsub(NaN, 0, CladeRates[, "Changes"] / CladeRates[, "Completeness"])), CladeRates)
      
      # If using Likelihood Ratio Test:
      if (LikelihoodTest == "LRT") {
        
        # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
        PartitionedData <- lapply(CladePartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(EdgeChanges[y]), sum(EdgeCompleteness[y] * EdgeDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
        
        # Add sampled rate to paritioned data matrices:
        PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
        
        # Get LRT p-values and combine output as edge test results:
        CladePartitionTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], get_likelihood_p(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
        
      }
      
      # If using AIC:
      if (LikelihoodTest == "AIC") {
        
        # Build partitioned data for AIC:
        PartitionedData <- lapply(CladePartitionsToTest, function(x) {y <- cbind(Partition = rep(NA, length(tree$edge.length)), Rate = rep(NA, length(tree$edge.length)), Changes = EdgeChanges, Completeness = EdgeCompleteness, Duration = EdgeDurations); y[, "Rate"] <- as.numeric(gsub(NaN, 0, unlist(lapply(x, function(x) rep(sum(y[x, "Changes"]) / (sum(y[x, "Completeness"]) * sum(y[x, "Duration"])), length(x))))[order(unlist(x))])); y[, "Partition"] <- rep(1:length(x), unlist(lapply(x, length)))[order(unlist(x))]; y})
        
        # Get AIC, AICc and rate results:
        CladePartitionTestResults <- lapply(PartitionedData, function(x) list(Rates = unname(unlist(lapply(as.list(unique(x[, "Partition"])), function(y) x[x[, "Partition"] == y, "Rate"][1]))), AIC = calculate_partition_AIC(x), AICc = calculate_partition_AIC(x, AICc = TRUE)))
        
      }
      
      # Pack clade partitions to test into single strings for output:
      PackedCladePartitions <- pack_partitions(CladePartitionsToTest)
      
      # Add packed partitions to results:
      for(i in 1:length(CladePartitionTestResults)) CladePartitionTestResults[[i]]$Partition <- PackedCladePartitions[i]


    # If not performing clade partition tests:
    } else {
      
      # Make empty clade partition result output:
      CladePartitionTestResults <- NULL
      
    }
    
  # If not doing clade OR branch tests:
  } else {
    
    # Make empty branch partition result output:
    BranchPartitionTestResults <- NULL

    # Make empty clade partition result output:
    CladePartitionTestResults <- NULL

  }

  # If performing character partition tests:
  if (!is.null(CharacterPartitionsToTest)) {
    
    # Get vector of (weighted) changes for each character:
    CharacterChanges <- unlist(lapply(as.list(CharacterNumbers), function(x) {CharacterRows <- which(AllChanges[, "Character"] == x); sum(AllChanges[CharacterRows, "Steps"] * AllChanges[CharacterRows, "Weight"])}))
    
    # Get vector of weighted durations for each character:
    CharacterDurations <- (weights / sum(weights)) * sum(tree$edge.length)
    
    # Get vector of completness (opportunity to observe changes) for each character:
    CharacterCompleteness <- apply(do.call(rbind, lapply(EdgeList, function(x) {CharacterPresence <- rep(0, times = length(CharacterNumbers)); CharacterPresence[x$Comparablecharacters] <- 1; CharacterPresence * x$BranchDuration})), 2, sum) / sum(tree$edge.length)
    
    # Set global rate:
    GlobalRate <- sum(CharacterChanges) / sum(CharacterCompleteness * CharacterDurations)
    
    # Create character rates for output:
    CharacterRates <- lapply(list(as.list(1:max(CharacterNumbers))), function(x) matrix(unlist(lapply(x, function(y) c(sum(CharacterChanges[y]), sum(CharacterCompleteness[y] * CharacterDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))[[1]]
    
    # Add character numbers and rates:
    CharacterRates <- cbind(Character = 1:max(CharacterNumbers), Rate = as.numeric(gsub(NaN, 0, CharacterRates[, "Changes"] / CharacterRates[, "Completeness"])), CharacterRates)
    
    # If using likelihood ratio test:
    if (LikelihoodTest== "LRT") {
      
      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
      PartitionedData <- lapply(CharacterPartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(CharacterChanges[y]), sum(CharacterCompleteness[y] * CharacterDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
      
      # Add sampled rate to paritioned data matrices:
      PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
      
      # Get P-Values and combine output as edge test results:
      CharacterPartitionTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], get_likelihood_p(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
      
    }
    
    # If using AIC:
    if (LikelihoodTest== "AIC") {
      
      # Build partitioned data for AIC:
      PartitionedData <- lapply(CharacterPartitionsToTest, function(x) {y <- cbind(Partition = rep(NA, max(CharacterNumbers)), Rate = rep(NA, max(CharacterNumbers)), Changes = CharacterChanges, Completeness = CharacterCompleteness, Duration = CharacterDurations); y[, "Rate"] <- as.numeric(gsub(NaN, 0, unlist(lapply(x, function(x) rep(sum(y[x, "Changes"]) / (sum(y[x, "Completeness"]) * sum(y[x, "Duration"])), length(x))))[order(unlist(x))])); y[, "Partition"] <- rep(1:length(x), unlist(lapply(x, length)))[order(unlist(x))]; y})
      
      # Get AIC, AICc and rate results:
      CharacterPartitionTestResults <- lapply(PartitionedData, function(x) list(Rates = unname(unlist(lapply(as.list(unique(x[, "Partition"])), function(y) x[x[, "Partition"] == y, "Rate"][1]))), AIC = calculate_partition_AIC(x), AICc = calculate_partition_AIC(x, AICc = TRUE)))
      
    }
    
    # Pack character partitions to test into single strings for output:
    PackedCharacterPartitions <- pack_partitions(CharacterPartitionsToTest)
    
    # Add packed partitions to results:
    for(i in 1:length(CharacterPartitionTestResults)) CharacterPartitionTestResults[[i]]$Partition <- PackedCharacterPartitions[i]


  # If performing branch partition tests:
  } else {
    
    # Make empty character partition result output:
    CharacterPartitionTestResults <- NULL
    
  }

  # If performing time bin partition tests:
  if (!is.null(TimeBinPartitionsToTest)) {
    
    # Get weighted number of changes from each time bin:
    Timebin_changes <- unlist(lapply(as.list(1:(length(TimeBins) - 1)), function(x) {ChangeRows <- AllChanges[, "Bin"] == x; sum(AllChanges[ChangeRows, "Steps"] * AllChanges[ChangeRows, "Weight"])}))
    
    # If using the Close time bin completeness approach get completeness value for each time bin:
    if (TimeBinApproach == "Close") TimeBinCompleteness <- apply(do.call(rbind, lapply(EdgeList, function(x) x$ProportionalBinnedEdgeDurations * (sum(weights[x$Comparablecharacters]) / sum(weights)))), 2, sum) / apply(do.call(rbind, lapply(EdgeList, function(x) x$ProportionalBinnedEdgeDurations)), 2, sum)
    
    # If using the Lloyd time bin completeness approach get completeness value for each time bin::
    if (TimeBinApproach == "Lloyd") TimeBinCompleteness <- apply(do.call(rbind, lapply(EdgeList, function(x) apply(matrix(weights[x$Comparablecharacters], ncol = 1) %*% x$BinnedBranchDurations, 2, sum))), 2, sum) / apply(do.call(rbind, lapply(EdgeList, function(x) apply(matrix(weights, ncol = 1) %*% x$BinnedBranchDurations, 2, sum))), 2, sum)
    
    # Get durations of edges in each time bin:
    TimeBinDurations <- apply(do.call(rbind, lapply(EdgeList, function(x) x$BinnedBranchDurations)), 2, sum)
    
    # Set global rate (NB: will differ between Close and Lloyd approaches, but Lloyd approach will match edge or character global rate):
    GlobalRate <- sum(Timebin_changes) / sum(TimeBinCompleteness * TimeBinDurations)
    
    # Create time rates for output:
    TimeRates <- lapply(list(as.list(1:(length(TimeBins) - 1))), function(x) matrix(unlist(lapply(x, function(y) c(sum(Timebin_changes[y]), sum(TimeBinCompleteness[y] * TimeBinDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))[[1]]
    
    # Add bin numbers and rates:
    TimeRates <- cbind(Bin = 1:(length(TimeBins) - 1), Rate = as.numeric(gsub(NaN, 0, TimeRates[, "Changes"] / TimeRates[, "Completeness"])), TimeRates)
    
    # If using Likelihood Ratio Test to compare partitions:
    if (LikelihoodTest == "LRT") {
      
      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later) for LRT:
      PartitionedData <- lapply(TimeBinPartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(Timebin_changes[y]), sum(TimeBinCompleteness[y] * TimeBinDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
      
      # Add sampled rate to paritioned data matrices for LRT:
      PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
      
      # Get P-Values and combine output as edge test results:
      TimeBinTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], get_likelihood_p(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
      
    }
    
    # If using AIC to compare partitions:
    if (LikelihoodTest == "AIC") {
      
      # Build partitioned data for AIC:
      PartitionedData <- lapply(TimeBinPartitionsToTest, function(x) {y <- cbind(Partition = rep(NA, length(Timebin_changes)), Rate = rep(NA, length(Timebin_changes)), Changes = Timebin_changes, Completeness = TimeBinCompleteness * TimeBinDurations, Duration = rep(1, length(Timebin_changes))); y[, "Rate"] <- as.numeric(gsub(NaN, 0, unlist(lapply(x, function(x) rep(sum(y[x, "Changes"]) / sum(y[x, "Completeness"]), length(x)))))); y[, "Partition"] <- rep(1:length(x), unlist(lapply(x, length))); y})
      
      # Get AIC, AICc and rate results:
      TimeBinTestResults <- lapply(PartitionedData, function(x) list(Rates = unname(unlist(lapply(as.list(unique(x[, "Partition"])), function(y) x[x[, "Partition"] == y, "Rate"][1]))), AIC = calculate_partition_AIC(x), AICc = calculate_partition_AIC(x, AICc = TRUE)))
      
    }
    
    # Pack time bin partitions to test into single strings for output:
    PackedTimeBinPartitions <- pack_partitions(TimeBinPartitionsToTest)
    
    # Add packed partitions to results:
    for(i in 1:length(TimeBinTestResults)) TimeBinTestResults[[i]]$Partition <- PackedTimeBinPartitions[i]

  # If not performing time bin partition tests:
  } else {
    
    # Make empty time bin partition result output:
    TimeBinTestResults <- NULL
    
  }

  # Set global rate for output:
  GlobalRate <- sum(unlist(lapply(EdgeList, function(x) sum(x$CharacterChanges[, "Steps"] * x$CharacterChanges[, "Weight"])))) / sum(unlist(lapply(EdgeList, function(x) sum(weights[x$Comparablecharacters]) / sum(weights))) * unlist(lapply(EdgeList, function(x) x$BranchDuration)))
  
  # If performing Likelihood Ratio Test:
  if (LikelihoodTest == "LRT") {
    
    # Subfunction to calculate adjusted alphas for multiple comparison corrections:
    add_cutoffs <- function(TestResults, alpha, MultipleComparisonCorrection = MultipleComparisonCorrection) {
      
      # Get number of comparisons performed:
      NComparisons <- length(TestResults)
      
      # If using the Benjamini-Hochberg false discovery rate approach:
      if (MultipleComparisonCorrection == "BenjaminiHochberg") {
        
        # Set cutoff values:
        CutoffValues <- ((1:NComparisons) / NComparisons) * alpha
        
        # Get actual p-values found:
        PValues <- unlist(lapply(TestResults, '[[', "PValue"))
        
        # Order cutoffs by p-value rank:
        CutoffValues <- CutoffValues[rank(PValues, ties.method = "random")]
        
      }
      
      # If using the Bonferroni correction set cutoff values as alpha over N:
      if (MultipleComparisonCorrection == "Bonferroni") CutoffValues <- alpha / NComparisons
      
      # Add cutoffs to output:
      for(i in 1:length(TestResults)) TestResults[[i]]$CorrectedAlpha <- CutoffValues[i]
      
      # Return modified test results:
      return(TestResults)
      
    }
    
    # If doing branch partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(BranchPartitionsToTest)) BranchPartitionTestResults <- add_cutoffs(TestResults = BranchPartitionTestResults, alpha = alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
    
    # If doing character partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(CharacterPartitionsToTest)) CharacterPartitionTestResults <- add_cutoffs(TestResults = CharacterPartitionTestResults, alpha = alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
    
    # If doing clade partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(CladePartitionsToTest)) CladePartitionTestResults <- add_cutoffs(TestResults = CladePartitionTestResults, alpha = alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
    
    # If doing time bin partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(TimeBinPartitionsToTest)) TimeBinTestResults <- add_cutoffs(TestResults = TimeBinTestResults, alpha = alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
    
  }

  # Compile output:
  Output <- list(TimeBinsUsed = TimeBins, InferredCharacterChanges = AllChanges, IntrinsicCharacterRate = GlobalRate, ContinuouscharactersConvertedToDiscrete = ContinuouscharactersConverted, BranchPartitionResults = BranchPartitionTestResults, CharacterPartitionResults = CharacterPartitionTestResults, CladePartitionResults = CladePartitionTestResults, TimeBinResults = TimeBinTestResults, BranchRates = BranchRates, CharacterRates = CharacterRates, CladeRates = CladeRates, TimeRates = TimeRates, Tree = tree)
  
  # Return output:
  return(Output)

}

# Ages <- read.table("~/Documents/Packages/Claddis/LungfishTest/ages.txt", sep = ",")
# Matrix <- Claddis::read_nexus_matrix("~/Documents/Packages/Claddis/LungfishTest/Lloyd_etal_2012a.nex")
# Tree <- ape::read.tree("~/Documents/Packages/Claddis/LungfishTest/Lloyd_etal_2012a.tre")
# TimeBins <- c(443.8, 419.2, 358.9, 298.9, 251.9, 201.3, 145.0, 0.0)
# 
# Tree <- Tree[sample(1:100000, 100)]
# Tree <- lapply(Tree, function(x) strap::DatePhylo(x, Ages, rlen = 2, method = "equal"))
# class(Tree) <- "multiPhylo"
# 
# tree <- Tree[[1]]
# cladistic.matrix <- Matrix
# TimeBins <- TimeBins
# BranchPartitionsToTest <- lapply(as.list(1:nrow(tree$edge)), as.list)
# CharacterPartitionsToTest <- list(list(1:91), list(Cranial = 1:81, Postcranial = 82:91))
# CladePartitionsToTest <- lapply(as.list(Ntip(tree) + (2:Nnode(tree))), as.list)
# TimeBinPartitionsToTest <- partition_time_bins(7)
# ChangeTimes = "random"
# LikelihoodTest = "AIC"
# alpha = 0.01
# MultipleComparisonCorrection = "BenjaminiHochberg"
# PolymorphismState = "missing"
# UncertaintyState = "missing"
# InapplicableState = "missing"
# TimeBinApproach = "Lloyd"
# EnsureAllweightsAreIntegers = FALSE
# estimate.all.nodes = FALSE
# estimate.tip.values = FALSE
# inapplicables.as.missing = FALSE
# polymorphism.behaviour = "equalp"
# uncertainty.behaviour = "equalp"
# threshold = 0.01

