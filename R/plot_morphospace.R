#' Plot Morphopace
#'
#' @description
#'
#' Plots a morphospace using the output from ordinate_cladistic_matrix.
#'
#' @param pcoa_input The main input in the format output from \link{ordinate_cladistic_matrix}.
#' @param x_axis Which ordination axis to plot as the x-axis (defaults to 1).
#' @param y_axis Which ordination axis to plot as the y-axis (defaults to 2).
#' @param z_axis Which ordination axis to plot as the z-axis (defaults to NULL, i.e., is not plotted).
#' @param taxon_groups A named list of groups to which taxa are assigned (optional). This is used to plot points or convex hulls in different colours corresponding to each group. As the user names the groups these can represent any grouping of interest (e.g., taxonomic, ecological, temporal, spatial). \link{assign_taxa_to_bins} can automate temporal assignments.
#' @param plot_taxon_names Logical indicating whether to plot the names of the taxa (defaults to FALSE).
#' @param plot_convex_hulls Logical indicating whether to plot convex hulls around any taxon_groups (if used).
#' @param plot_internal_nodes Logical indicating whether to plot the internal nodes of the tree (if included in \code{pcoa_input}) (defaults to FALSE).
#' @param plot_edges Logical indicating whether to plot the branches of the tree (if included in \code{pcoa_input}) (defaults to TRUE).
#' @param plot_root Logical indicating whether to plot the root separately (defaults to FALSE).
#' @param root_colour If plotting the root separately (previous option) sets the root colour.
#' @param palette The palette to use for plotting each element of taxon_groups. See \link[grDevices]{palette}.
#' @param plot_group_legend Logical indicating whether to plot a legend for taxon_groups. (Default is TRUE.)
#' @param group_legend_position Position to plot the group legend. Must be one of \code{bottom_left}, \code{bottom_right}, \code{top_left}, or \code{top_right} (the default).
#' @param plot_z_legend Logical indicating whether to plot a legend for the z-axis. (Default is TRUE.)
#' @param z_legend_position Position to plot the group legend. Must be one of \code{bottom_left}, \code{bottom_right} (the default), \code{top_left}, or \code{top_right}.
#' @param inform Logical indicating whether to inform the user of any taxon pruning. (Default is TRUE.)
#' @param x_limits Plot limits to use for x-axis. Only intended for use by \link{plot_multi_morphospace}.
#' @param y_limits Plot limits to use for y-axis. Only intended for use by \link{plot_multi_morphospace}.
#'
#' @details
#'
#' Uses output from \link{ordinate_cladistic_matrix} to make morphospace plots.
#'
#' Allows plotting of a third axis using the technique of Wills et al. (1994; their Figures 4 and 8; Wills 1998; his Figure 4), where solid and open indicate positive and negative values respectively, and the size of points their magnitudes.
#'
#' Will automatically generate phylomorphospaces if a tree was included in the ordination.
#'
#' Can also plot groups of points - whether they represent taxonomic, ecological, temporal, or spatial groupings - in different colours as well as plot translucent convex hulls around these groups, by using the \code{taxon_groups} and \code{plot_convex_hulls = TRUE} options, respectively. Note that \code{taxon_groups} should be in the form of a named list (see example below for how these should be formatted).
#'
#' Various other options allow toggling of particular features on or off. For example, the taxon names can be shown with \code{plot_taxon_names = TRUE}.
#'
#' Note that some features will generate legends that may initially appear to disappear off the sides of the plot, but simple resizing of the plot window (or increasing the width:height ratio if outputting to a file) should fix this.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Emma Sherratt \email{emma.sherratt@@gmail.com}
#'
#' @seealso
#'
#' \link{assign_taxa_to_bins}, \link{plot_chronophylomorphospace}, \link{plot_morphospace_stack}, \link{plot_multi_morphospace}, \link{ordinate_cladistic_matrix}
#'
#' @references
#'
#' Wills, M. A., 1998. Cambrian and Recent disparity: the picture from priapulids. \emph{Paleobiology}, \bold{24}, 177-199.
#'
#' Wills, M. A., Briggs, D. E. G. and Fortey, R. A., 1994. Disparity as an evolutionary index: a comparison of Cambrian and Recent arthropods. \emph{Paleobiology}, \bold{20}, 93-130.
#'
#' @examples
#'
#' \donttest{
#' # Perform a PCoA ordination on the day_2016 data set:
#' pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = day_2016)
#'
#' # Plot this as a simple bivarate morphospace:
#' plot_morphospace(pcoa_input = pcoa_input)
#'
#' # Use the Wills technique to add a third axis (PC3):
#' plot_morphospace(pcoa_input = pcoa_input, z_axis = 3)
#'
#' # You may need to resize the plot to see the legend for the z-axis
#'
#' # Add taxon names as well:
#' plot_morphospace(pcoa_input = pcoa_input, z_axis = 3, plot_taxon_names = TRUE)
#'
#' # Define some simple taxon groups for the data as a named list:
#' taxon_groups <- list(nonBurnetiamorpha = c("Biarmosuchus_tener",
#'   "Hipposaurus_boonstrai", "Bullacephalus_jacksoni", "Pachydectes_elsi",
#'   "Ictidorhinus_martinsi", "RC_20", "Herpetoskylax_hopsoni"),
#'   Burnetiamorpha = c("Lemurosaurus_pricei", "Lobalopex_mordax",
#'   "Lophorhinus_willodenensis", "Proburnetia_viatkensis", "Lende_chiweta",
#'   "Paraburnetia_sneeubergensis", "Burnetia_mirabilis", "BP_1_7098"))
#'
#' # Plot taxon groups including convex hulls:
#' plot_morphospace(pcoa_input = pcoa_input, z_axis = 3, plot_taxon_names = TRUE,
#'   taxon_groups = taxon_groups, plot_convex_hulls = TRUE)
#'
#' # Make time-scaled first MPT for Day 2016 data set:
#' time_tree <- ape::read.tree(text = paste0("(Biarmosuchus_tener:0.5,",
#'   "(((Hipposaurus_boonstrai:3.5,(Bullacephalus_jacksoni:0.75,",
#'   "Pachydectes_elsi:0.75):0.75):0.75,(Lemurosaurus_pricei:7.166666667,",
#'   "(Lobalopex_mordax:4.333333333,((Lophorhinus_willodenensis:3.666666667,",
#'   "(Proburnetia_viatkensis:0.8333333333,(Lende_chiweta:2,",
#'   "(Paraburnetia_sneeubergensis:1,Burnetia_mirabilis:2):1):1.833333333)",
#'   ":0.8333333333):0.8333333333,(BP_1_7098:2.25,Niuksenitia_sukhonensis:",
#'   "1.25):1.25):0.8333333333):0.8333333333):3.083333333):1.95,",
#'   "(Ictidorhinus_martinsi:15.9,(RC_20:11.6,(Herpetoskylax_hopsoni:11.3,",
#'   "Lycaenodon_longiceps:0.3):0.3):0.3):0.3):0.3);"))
#'
#' # Add root age to tree:
#' time_tree$root.time <- 269.5
#'
#' # Prune incomplete taxa from tree:
#' time_tree <- ape::drop.tip(phy = time_tree, tip = c("Lycaenodon_longiceps",
#'   "Niuksenitia_sukhonensis"))
#'
#' # Prune incomplete taxa from cladistic matrix:
#' cladistic_matrix <- prune_cladistic_matrix(cladistic_matrix = day_2016,
#'   taxa2prune = c("Lycaenodon_longiceps", "Niuksenitia_sukhonensis"))
#'
#' # Note: the above pruning is simply to run this example and should not be
#' # done manually as a matter of course as the functions will automatically
#' # prune tips and nodes as required.
#'
#' # Make new ordination with tree included (enabling phylomorphospace):
#' pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = cladistic_matrix,
#'   time_tree = time_tree)
#'
#' # Plot this as a simple bivarate phylomorphospace:
#' plot_morphospace(pcoa_input = pcoa_input)
#'
#' # Use the Wills technique to add a third axis (PC3):
#' plot_morphospace(pcoa_input = pcoa_input, z_axis = 3)
#'
#' # You may need to resize the plot to see the legend for the z-axis
#'
#' # Add taxon names as well:
#' plot_morphospace(pcoa_input = pcoa_input, z_axis = 3, plot_taxon_names = TRUE)
#'
#' # Add taxon groups including convex hulls:
#' plot_morphospace(pcoa_input = pcoa_input, z_axis = 3, plot_taxon_names = TRUE,
#'   taxon_groups = taxon_groups, plot_convex_hulls = TRUE)
#' }
#' @export plot_morphospace
plot_morphospace <- function(pcoa_input, x_axis = 1, y_axis = 2, z_axis = NULL, taxon_groups = NULL, plot_taxon_names = FALSE, plot_convex_hulls = FALSE, plot_internal_nodes = FALSE, plot_edges = TRUE, plot_root = TRUE, root_colour = "red", palette = "viridis", plot_group_legend = TRUE, group_legend_position = "top_right", plot_z_legend = TRUE, z_legend_position = "bottom_right", inform = TRUE, x_limits = NULL, y_limits = NULL) {
  
  # TO DO:
  #
  # - Order points by z-value so they are actually plotted from "back" to "front".
  # - Add plot_node_names option.
  # - Check inputs.

  # Check group_legend_position is a valid value and stop and warn user if not:
  if (!group_legend_position %in% c("bottom_left", "bottom_right", "top_left", "top_right")) stop("group_legend_position must be one of \"bottom_left\", \"bottom_right\", \"top_left\", or \"top_right\".")
  
  # Check z_legend_position is a valid value and stop and warn user if not:
  if (!z_legend_position %in% c("bottom_left", "bottom_right", "top_left", "top_right")) stop("z_legend_position must be one of \"bottom_left\", \"bottom_right\", \"top_left\", or \"top_right\".")
  
  # Check that if both legends are used that they are in different positions:
  if (plot_group_legend && plot_z_legend && group_legend_position == z_legend_position) stop("plot_group_legend and plot_z_legend must be different values or they will plot on top of each other.")

  # Create logical for whether taxon groups are used or not:
  taxon_groups_used <- ifelse(test = is.null(x = taxon_groups), yes = FALSE, no = TRUE)
  
  # Create logical for whether tree is used or not:
  tree_used <- ifelse(test = is.null(x = pcoa_input$time_tree), yes = FALSE, no = TRUE)
  
  # If using taxon groups:
  if (taxon_groups_used) {
    
    # Find any taxa to prune from taxon groups:
    taxa_to_prune <- setdiff(x = unique(x = unlist(x = taxon_groups)), y = rownames(x = pcoa_input$vectors))
    
    # If taxa to prune are found:
    if (length(x = taxa_to_prune) > 0) {
      
      # Go through taxon groups:
      taxon_groups <- lapply(X = taxon_groups, function(y) {
        
        # Remove any taxa not found in pcoa data:
        if (length(x = sort(x = match(x = taxa_to_prune, table = y))) > 0) y <- y[-sort(x = match(x = taxa_to_prune, table = y))]
        
        # Return
        y
        
      })
      
      # Warn user that this has happened in case it is an error:
      if (inform) print(paste0("Warning: The following taxa were removed from taxon_groups as they do not appear in pcoa_input: ", paste(taxa_to_prune, collapse = ", "), ". You may wish to double check this makes sense (e.g., because of incomplete taxa being removed by trim_matrix) and is not due to a typographical or other error which means names are not an exact match."))
    }
    
  }
  
  # If plot limits aren't set use x and y ranges to set them:
  if (is.null(x_limits)) x_limits <- range(pcoa_input$vectors[, x_axis])
  if (is.null(y_limits)) y_limits <- range(pcoa_input$vectors[, y_axis])
  
  # Set default tip numbers as just 1 to N:
  tip_numbers <- 1:nrow(x = pcoa_input$vectors)
  
  # Case if tree supplied:
  if (tree_used) {
    
    # Set basic tree information:
    n_tips <- ape::Ntip(phy = pcoa_input$time_tree)
    tip_numbers <- c(1:n_tips)
    node_numbers <- setdiff(x = 1:nrow(pcoa_input$vectors), y = tip_numbers)
    root_number <- n_tips + 1
  }
  
  # Get vector of values that correspond to scree plot:
  scree_values <- apply(pcoa_input$vectors, 2, var) / sum(apply(pcoa_input$vectors, 2, var)) * 100
  
  # Set default solid colour to black:
  solid_colours <- "black"
  
  # Set default transparent colour to 50% black:
  translucent_colours <- grDevices::rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
  
  # Set open colour (completely transparent):
  transparent_colour <- grDevices::rgb(red = 0, green = 0, blue = 0, alpha = 0)
  
  # If using taxon groups:
  if (taxon_groups_used) {
    
    # Set colours for each group:
    solid_colours <- grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1)
    translucent_colours <- grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 0.5)
    transparent_colour <- grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 0)
    
    # If not using taxon groups:
  } else {
    
    # Create dummy taxon group of all data as one:
    taxon_groups <- list(Data = rownames(pcoa_input$vector))
  }
  
  # If using a tree make sure taxon_groups only includes tips:
  if (tree_used) taxon_groups <- lapply(X = taxon_groups, FUN = function(y) intersect(y, pcoa_input$time_tree$tip.label[tip_numbers]))

  # Set point colours (background and border):
  point_col <- point_bg <- lapply(X = as.list(x = 1:length(taxon_groups)), FUN = function(y) rep(x = solid_colours[y], length.out = length(taxon_groups[[y]])))
  if (!is.null(x = z_axis)) point_bg <- lapply(X = as.list(x = 1:length(taxon_groups)), FUN = function(y) as.vector(unlist(x = lapply(X = as.list(x = taxon_groups[[y]]), FUN = function(z) ifelse(test = pcoa_input$vectors[z, z_axis] > 0, yes = solid_colours[[y]], no = transparent_colour[[y]]))), mode = "character"))

  # Make axis labels:
  x_lab <- paste("PC", x_axis, " (", round(scree_values[x_axis], 2), "% of total variance)", sep = "")
  y_lab <- paste("PC", y_axis, " (", round(scree_values[y_axis], 2), "% of total variance)", sep = "")
  if (!is.null(x = z_axis)) z_lab <- paste("PC", z_axis, " (", round(scree_values[z_axis], 2), "% of total variance)", sep = "")

  # Make all points equal in size by default:
  point_sizes <- rep(1, nrow(pcoa_input$vectors))
  
  # Set point sizes as absolute z-value if using z-axis:
  if (!is.null(x = z_axis)) point_sizes <- abs(x = pcoa_input$vectors[, z_axis]) / max(abs(x = pcoa_input$vectors[, z_axis])) * 3
  
  # Add taxon names to point_sizes:
  names(x = point_sizes) <- rownames(x = pcoa_input$vectors)
  
  # Create the basic plot space (will be empty for now):
  graphics::plot(x = pcoa_input$vectors[, x_axis], y = pcoa_input$vectors[, y_axis], type = "n", bg = "black", xlab = x_lab, ylab = y_lab, asp = TRUE, xlim = x_limits, ylim = y_limits)
  
  # Sort vectors by node number (1:N):
  if (tree_used) pcoa_input$vectors <- pcoa_input$vectors[c(pcoa_input$time_tree$tip.label, setdiff(x = rownames(x = pcoa_input$vectors), y = pcoa_input$time_tree$tip.label)), ]
  
  # Plot branches of tree (if a tree is used plotting requested):
  if (tree_used && plot_edges) for (i in 1:nrow(pcoa_input$time_tree$edge)) lines(x = pcoa_input$vectors[pcoa_input$time_tree$edge[i, ], x_axis], y = pcoa_input$vectors[pcoa_input$time_tree$edge[i, ], y_axis], col = grDevices::rgb(red = 0.65, green = 0.65, blue = 0.65, alpha = 1))
  
  # Set node colours for plotting:
  if (tree_used) node_bg <- unlist(x = lapply(X = as.list(pcoa_input$vectors[as.character(node_numbers), z_axis]), FUN = function(y) ifelse(y > 0, grDevices::rgb(red = 0.65, green = 0.65, blue = 0.65, alpha = 1), grDevices::rgb(red = 0.65, green = 0.65, blue = 0.65, alpha = 0))))

  # Plot internal nodes, if requested:
  if (tree_used && plot_internal_nodes) graphics::points(pcoa_input$vectors[node_numbers, x_axis], pcoa_input$vectors[node_numbers, y_axis], pch = 21, bg = node_bg[as.character(x = node_numbers)], col = grDevices::rgb(red = 0.65, green = 0.65, blue = 0.65, alpha = 1), cex = point_sizes[node_numbers])

  # Plot root separetely, if requested:
  if (tree_used && plot_root) graphics::points(pcoa_input$vectors[root_number, x_axis], pcoa_input$vectors[root_number, y_axis], pch = 21, col = root_colour, bg = root_colour, cex = point_sizes[root_number])

  # If convex hulls are requested:
  if(taxon_groups_used && plot_convex_hulls) {
      
    # For each group:
    x <- lapply(X = as.list(x = 1:length(x = taxon_groups)), function(y) {
        
      # Make convex hull for data:
      convex_hull <- grDevices::chull(x = pcoa_input$vectors[taxon_groups[[y]], x_axis], y = pcoa_input$vectors[taxon_groups[[y]], y_axis])
        
      # Plot convex hull as translucent polygon:
      graphics::polygon(x = pcoa_input$vectors[taxon_groups[[y]][convex_hull], x_axis], y = pcoa_input$vectors[taxon_groups[[y]][convex_hull], y_axis], col = translucent_colours[[y]], border = NA)
    })
  }
    
  # Add points to plot:
  x <- lapply(X = as.list(x = 1:length(x = taxon_groups)), function(y) graphics::points(x = pcoa_input$vectors[taxon_groups[[y]], x_axis], y = pcoa_input$vectors[taxon_groups[[y]], y_axis], pch = 21, bg = point_bg[[y]], col = point_col[[y]], cex = point_sizes[taxon_groups[[y]]]))
    
  # If plotting taxon names:
  if (plot_taxon_names) {
    
    # First establish a default position for names (to the left of the point):
    x_positions <- rep(2, nrow(pcoa_input$vectors))
    
    # Now changes negative values to plot on the right instead:
    x_positions[which(x = pcoa_input$vectors[, x_axis] < 0)] <- 4
    
    # Plot taxon names (for tips only):
    graphics::text(x = pcoa_input$vectors[tip_numbers, x_axis], y = pcoa_input$vectors[tip_numbers, y_axis], labels = rownames(x = pcoa_input$vectors)[tip_numbers], pos = x_positions[tip_numbers], cex = 0.7)
  }

  # If plotting a group legend:
  if(taxon_groups_used && plot_group_legend) {
    
    # Add groups legend to plot:
    if(group_legend_position == "bottom_left") graphics::legend(x = min(pcoa_input$vectors[, x_axis]), y = min(pcoa_input$vectors[, y_axis]), legend = names(taxon_groups), fill = solid_colours, bg = "white", xjust = 1, yjust = 0)
    if(group_legend_position == "bottom_right") graphics::legend(x = max(pcoa_input$vectors[, x_axis]), y = min(pcoa_input$vectors[, y_axis]), legend = names(taxon_groups), fill = solid_colours, bg = "white", xjust = 0, yjust = 0)
    if(group_legend_position == "top_left") graphics::legend(x = min(pcoa_input$vectors[, x_axis]), y = max(pcoa_input$vectors[, y_axis]), legend = names(taxon_groups), fill = solid_colours, bg = "white", xjust = 1, yjust = 1)
    if(group_legend_position == "top_right") graphics::legend(x = max(pcoa_input$vectors[, x_axis]), y = max(pcoa_input$vectors[, y_axis]), legend = names(taxon_groups), fill = solid_colours, bg = "white", xjust = 0, yjust = 1)
  }
  
  # If plotting a z-axis legend:
  if (plot_z_legend && !is.null(z_axis)) {
    
    # Collapse range of z-values to a spread f six:
    z_values <- seq(from = min(x = pcoa_input$vectors[, z_axis]), to = max(x = pcoa_input$vectors[, z_axis]), length.out = 6)
    
    # Make z point sizes for legend:
    z_sizes <- abs(x = z_values) / max(abs(x = z_values)) * 3
    
    # Add z legend to plot:
    if(z_legend_position == "bottom_left") graphics::legend(x = min(pcoa_input$vectors[, x_axis]), y = min(pcoa_input$vectors[, y_axis]), legend = signif(x = z_values, digits = 4), pch = 21, bg = "white", xjust = 1, yjust = 0, pt.cex = z_sizes, pt.bg = unlist(lapply(X = as.list(x = z_values), FUN = function(y) ifelse(test = y > 0, yes = "black", no = "white"))), col = "black")
    if(z_legend_position == "bottom_right") graphics::legend(x = max(pcoa_input$vectors[, x_axis]), y = min(pcoa_input$vectors[, y_axis]), legend = signif(x = z_values, digits = 4), pch = 21, bg = "white", xjust = 0, yjust = 0, pt.cex = z_sizes, pt.bg = unlist(lapply(X = as.list(x = z_values), FUN = function(y) ifelse(test = y > 0, yes = "black", no = "white"))), col = "black")
    if(z_legend_position == "top_left") graphics::legend(x = min(pcoa_input$vectors[, x_axis]), y = max(pcoa_input$vectors[, y_axis]), legend = signif(x = z_values, digits = 4), pch = 21, bg = "white", xjust = 1, yjust = 1, pt.cex = z_sizes, pt.bg = unlist(lapply(X = as.list(x = z_values), FUN = function(y) ifelse(test = y > 0, yes = "black", no = "white"))), col = "black")
    if(z_legend_position == "top_right") graphics::legend(x = max(pcoa_input$vectors[, x_axis]), y = max(pcoa_input$vectors[, y_axis]), legend = signif(x = z_values, digits = 4), pch = 21, bg = "white", xjust = 0, yjust = 1, pt.cex = z_sizes, pt.bg = unlist(lapply(X = as.list(x = z_values), FUN = function(y) ifelse(test = y > 0, yes = "black", no = "white"))), col = "black")
  }
  
  # If using a z-axis add label as plot title:
  if (!is.null(x = z_axis)) graphics::title(main = z_lab)
}

#time_tree <- ape::read.tree(text = "(Psarolepis_romeri:1,(Diabolepis_speratus:1.85,((Dipnorhynchus_kiandrensis:7.4,(Archaeonectes_pertusus:28.3,(Uranolophus_wyomingensis:1.6,(Speonesydrion_iani:0.8,(Jarvikia_arctica:36.5173913,(((Adololopas_moyasmithae:10.775,((Adelargo_schultzei:14.05,Chirodipterus_australis:3.25):6.1,(Chirodipterus_rhenanus:1.425,(Chirodipterus_wildungensis:3.25,Dipterus_cf_valenciennesi:3.25):4.675):1.425):1.425):10.31485507,(Barwickia_downunda:10.14492754,Dipterus_valenciennesi:4.444927536):4.444927536):4.444927536,(Pillararhynchus_longi:25.35217391,(((Gogodipterus_paddyensis:24.80434783,((Tarachomylax_oepiki:1.947826087,(Amadeodipterus_kencampbelli:0.9739130435,Stomiahykus_thlaodus:10.47391304):0.9739130435):0.9739130435,(Iowadipterus_halli:17.93913043,((Delatitia_breviceps:50.17391304,(Phaneropleuron_andersoni:34.69130435,((Orlovichthys_limnatis:34.32608696,(Howidipterus_donnae:16.84347826,(((Andreyevichthys_epitomus:16.88913043,Oervigia_nordica:16.88913043):16.88913043,(Grossipterus_crassus:22.79565217,(Fleurantia_denticulata:22.61304348,((Robinsondipterus_longi:16.275,(Asthenorhynchus_meemannae:10.85,(Holodipterus_elderae:5.425,Holodipterus_gogoensis:5.425):5.425):5.425):6.155434783,((Griphognathus_minutidens:14.46666667,(Griphognathus_sculpta:7.233333333,Griphognathus_whitei:7.233333333):7.233333333):7.78115942,(Rhynchodipterus_elginensis:32.86521739,(Jessenia_concentrica:0.1826086957,Soederberghia_groenlandica:21.8826087):0.1826086957):0.1826086957):0.1826086957):0.1826086957):0.1826086957):0.1826086957):0.1826086957,(Pentlandia_macroptera:8.330434783,Scaumenacia_curta:14.83043478):8.330434783):0.1826086957):0.1826086957):0.1826086957,(Holodipterus_santacrucensis:21.07439614,((Ganopristodus_splendens:55.8057971,(Megapleuron_zangerli:86.77149758,(Sagenodus_inaequalis:50.53719807,(((Eoctenodus_microsoma:2.634299517,Tranodis_castrensis:61.53429952):2.634299517,(Ctenodus_romeri:15.68429952,Straitonia_waterstoni:29.58429952):15.68429952):2.634299517,((Parasagenodus_sibiricus:11.33429952,(Gnathorhiza_serrata:20.55,((Beltanodus_ambilobensis:53.1125,(Namatozodia_pitikanta:45.525,(Ariguna_formosa:37.9375,(((Aphelodus_anapes:11.38125,Ceratodus_formosa:11.38125):11.38125,((Asiatoceratodus_sharovi:7.5875,Gosfordia_truncata:30.5875):7.5875,(Neoceratodus_forsteri:77.0875,(Mioceratodus_gregoryi:37,(Lepidosiren_paradoxa:22.4,Protopterus_annectens:9.5):9.5):86.5875):77.0875):7.5875):7.5875,(Archaeoceratodus_avus:26.675,Tellerodus_sturi:38.175):26.675):7.5875):7.5875):7.5875):12.3875,(Microceratodus_angolensis:63.9,(Palaeophichthys_parvulus:1.6,(Ptychoceratodus_serratus:45.525,(Paraceratodus_germaini:31.65,(Arganodus_atlantis:15.175,Ferganoceratodus_jurassicus:104.975):15.175):15.175):16.775):1.6):1.6):22.15):31.88429952):11.33429952,(Ceratodus_latissimus:76.93429952,Metaceratodus_wollastoni:159.4342995):76.93429952):11.33429952):2.634299517):2.634299517):2.634299517):2.634299517,(Nielsenia_nordica:14.62004831,Conchopoma_gadiforme:77.42004831):14.62004831):2.634299517):2.634299517):0.1826086957):0.1826086957):0.1826086957,(Rhinodipterus_secans:15.37826087,Rhinodipterus_ulrichi:8.87826087):8.87826087):0.1826086957):0.1826086957):0.1826086957):0.1826086957,(Palaeodaphus_insignis:12.49347826,Sunwapta_grandiceps:23.29347826):12.49347826):0.1826086957,(Melanognathus_canadensis:1.734782609,Sorbitorhynchus_deleaskitus:1.734782609):1.734782609):0.1826086957):0.1826086957):0.1826086957):0.9826086957):0.8):0.8):0.8):0.8,(Westollrhynchus_lehmanni:2,(Ichnomylax_kurnai:3.36,(Dipnorhynchus_sussmilchi:2.52,(Chirodipterus_onawwayensis:16.88,(Dipnorhynch_cathlesae:0.84,Dipnorhynchus_kurikae:0.84):0.84):0.84):0.84):2.84):2):2.65):1.85);")
#time_tree$root.time <- 419.7
#cladistic_matrix <- read_nexus_matrix("http://www.graemetlloyd.com/nexus/Lloyd_etal_2012a.nex")
#pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = cladistic_matrix, time_tree = time_tree)
#x_axis = 1
#y_axis = 2
#taxon_groups = assign_taxa_to_bins(taxon_ages, named_time_bins) ### NEW PARAM
#plot_taxon_names = TRUE
#plot_internal_nodes = TRUE
#plot_root = TRUE
#root_colour = "red"
#palette = "viridis"
#plot_group_legend = TRUE
#group_legend_position = "top_right"
#plot_z_legend = TRUE
#z_legend_position = "bottom_right"
