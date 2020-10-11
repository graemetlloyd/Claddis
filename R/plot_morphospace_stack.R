#' Plot stacked ordination spaces
#'
#' @description
#'
#' Plots a stack of ordination spaces representing multiple time-slices.
#'
#' @param pcoa_input The main input in the format output from \link{ordinate_cladistic_matrix}.
#' @param taxon_ages A two-column matrix of the first and last apperance dates (columns; \code{"fad"} and \code{"lad"}) for the taxa (as rownames) from \code{pcoa_input}.
#' @param taxon_groups An object of class \code{taxonGroups}.
#' @param time_bins An object of class \code{timeBins}.
#' @param shear A single value (between 0 and 1) that determines the "sheared" visual appearance of the platforms.
#' @param x_axis The ordination axis to plot on the x-axis.
#' @param y_axis The ordination axis to plot nn the y-axis.
#' @param palette The palette to use for plotting each element of \code{taxon_groups}. See \link[grDevices]{palette}.
#' @param plot_cushion A number determining the "cushion" around the edge of each stack in which no data will be plotted. This should be larger than zero or points will "hang" over the edge. Additionally, if using a \code{platform_size} value in excess of one this will avoid points being hidden under overlying platforms. Note that this effectively adds empty plot space around the data, it does not remove anything.
#' @param platform_size The size of each platform as a proportion. Values of less than one will show slight gaps between platforms, whereas values in excess of one will mean platforms will appear to overlap.
#' @param plot_pillars Logical indicating whether or not to plot the pillars linking the corners of each platform.
#' @param plot_crosshair Logical indicating whether or not to plot the "crosshair" (i.e., the zero-zero lines that run through the origin of the morphospace).
#' @param plot_grid_cells Logical indicating whether or not to plot grid cells that help visualise the distorted aspect ratio of the plot. Each cell is a square in the ordination space.
#' @param plot_convex_hulls Logical indicating whether or not to plot convex hulls around the taxonomic groupings. Only relevant if \code{taxon_groups} is in use.
#' @param plot_timebin_names Logical indicating whether or not to plot the names of each time bin next to each platform. I.e., the rownames from \code{time_bins}. Note if these are long they may disappear behind overlying platforms. To avoid this try using a smaller \code{platform_size} value, a larger \code{shear} value, or simply shorter or abbreviated names.
#' @param plot_tickmarks Logical indicating whether or not to plot tickmarks next to the bottom platform.
#' @param plot_group_legend Logical indicating whether or not to plot a legend. Only relevant if using \code{taxon_groups}. Note this may obscure some points so use with caution and try picking different values for \code{group_legend_position} to avoid this.
#' @param group_legend_position The position the group legend should be plotted. Only relevant if using \code{taxon_groups} and \code{plot_group_legend = TRUE}. Options are: \code{"bottom_left"}, \code{"bottom_right"}, \code{"top_left"}, and \code{"top_right"}.
#' @param point_size The size at which the points should be plotted. Note that here points are custom polygons and hence are not editable by normal plot options, e.g., \code{pch} or \code{cex}. At present all points are plotted as circles.
#'
#' @details
#'
#' This style of plot is taken from various papers by Michael Foote (Foote 1993; his Figures 2, 4, 6, 8, 10, 12, and 14; Foote 1994; his Figure 2; Foote 1995; his Figure 3; Foote 1999; his Figure 22), and can be seen elsewhere in the literature (e.g., Friedman and Coates 2006; their Figure 2c). Here multiple ordination (or morpho-) spaces are plotted in stratigraphic order (oldest at bottom) as a stacked series of "platforms" representing named time bins.
#'
#' The user needs to supply three main pieces of information to use the function: 1) ordination data that includes rows (taxa) and columns (ordination axes), 2) the ages (first and last appearance dates) of the taxa sampled, and 3) the ages (first and last appearance dates) of the named time bins used.
#'
#' Note that since version 0.6.1 this function has been completely rewritten to better reflect the usage of these type of figures (see citations above) as well as allow additional features. This was also done in part to standardise the function to fit the style of the other major disparity plotting functions in Claddis, such as \link{plot_morphospace}. This means the input data is now assumed to come directly from \link{ordinate_cladistic_matrix}, but the user could easily still bring in data from elsewhere (the way the function worked previously) by reformatting it using something like:
#'
#' \code{pcoa_input <- list(vectors = my_imported_data)}
#'
#' Where my_imported_data has columns representing ordination axes (1 to N) and rownames corresponding to taxon names.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Emma Sherratt \email{emma.sherratt@@gmail.com}
#'
#' @seealso
#'
#' \link{assign_taxa_to_bins}, \link{plot_chronophylomorphospace}, \link{plot_morphospace}, \link{plot_multi_morphospace}, \link{ordinate_cladistic_matrix}
#'
#' @references
#'
#' Foote, M., 1993. Discordance and concordance between morphological and taxonomic diversity. \emph{Paleobiology}, \bold{19}, 185-204.
#'
#' Foote, M., 1994. Morphological disparity in Ordovician-Devonian crinoids and the early saturation of morphological space. \emph{Paleobiology}, \bold{20}, 320-344.
#'
#' Foote, M., 1995. Morphological diversification of Paleozoic crinoids. \emph{Paleobiology}, \bold{21}, 273-299.
#'
#' Foote, M., 1999. Morphological diversity in the evolutionary radiation of Paleozoic and post-Paleozoic crinoids. \emph{Paleobiology}, \bold{25}, 1-115.
#'
#' Friedman, M. and Coates, M. I., 2006. A newly recognized fossil coelacanth highlights the early morphological diversification of the clade. \emph{Proceedings of the Royal Society of London B}, \bold{273}, 245-250.

#' @examples
#'
#' # Build taxon ages matrix for Day et al 2016 data:
#' taxon_ages <- matrix(data = c(269, 267, 263, 260, 265, 265, 265, 265, 257, 255, 259, 259, 258, 258,
#'   260, 257, 257, 255, 257, 257, 255, 252, 259, 259, 260, 258, 253, 252, 257, 255, 257, 255),
#'   ncol = 2, byrow = TRUE, dimnames = list(c("Biarmosuchus_tener", "Hipposaurus_boonstrai",
#'   "Bullacephalus_jacksoni", "Pachydectes_elsi", "Lemurosaurus_pricei", "Lobalopex_mordax",
#'   "Lophorhinus_willodenensis", "Proburnetia_viatkensis", "Lende_chiweta",
#'   "Paraburnetia_sneeubergensis", "Burnetia_mirabilis", "BP_1_7098", "Niuksenitia_sukhonensis",
#'   "Ictidorhinus_martinsi", "RC_20", "Herpetoskylax_hopsoni"), c("FAD", "LAD")))
#'
#' # Ordinate Day et al 2016 data set:
#' pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = prune_cladistic_matrix(
#'   cladistic_matrix = day_2016,
#'   taxa2prune = "Lycaenodon_longiceps"))
#'
#' # Build simple taxonomic groups for Day et al 2016 daat set:
#' taxon_groups <- list(nonBurnetiamorpha = c("Biarmosuchus_tener", "Hipposaurus_boonstrai",
#'   "Bullacephalus_jacksoni", "Pachydectes_elsi", "Niuksenitia_sukhonensis", "Ictidorhinus_martinsi",
#'   "RC_20", "Herpetoskylax_hopsoni"), Burnetiamorpha = c("Lemurosaurus_pricei", "Lobalopex_mordax",
#'   "Lophorhinus_willodenensis", "Proburnetia_viatkensis", "Lende_chiweta",
#'   "Paraburnetia_sneeubergensis", "Burnetia_mirabilis", "BP_1_7098"))
#'
#' # Set class as taxonGroups:
#' class(taxon_groups) <- "taxonGroups"
#'
#' # Build a sequence of equally spaced time bins spanning Day et al. 2016 data:
#' time_sequence <- seq(from = 270, to = 252, length.out = 6)
#'
#' # Reformat this sequence into named time bin matrix:
#' time_bins <- matrix(
#'   data = c(time_sequence[1:(length(x = time_sequence) - 1)],
#'   time_sequence[2:length(x = time_sequence)]),
#'   ncol = 2,
#'   dimnames = list(c("Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5"), c("fad", "lad"))
#' )
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Plot morphospace stack using named time bins:
#' plot_morphospace_stack(
#'   pcoa_input = pcoa_input,
#'   taxon_ages = taxon_ages,
#'   taxon_groups = taxon_groups,
#'   time_bins = time_bins,
#' )
#'
#' @export plot_morphospace_stack
plot_morphospace_stack <- function(pcoa_input, taxon_ages, taxon_groups, time_bins, shear = 0.2, x_axis = 1, y_axis = 2, palette = "viridis", plot_cushion = 0.3, platform_size = 0.95, plot_pillars = TRUE, plot_crosshair = TRUE, plot_grid_cells = TRUE, plot_convex_hulls = TRUE, plot_timebin_names = TRUE, plot_tickmarks = TRUE, plot_group_legend = TRUE, group_legend_position = "bottom_right", point_size = 1.5) {

  # TO DO:
  #
  # - Needs plenty of checks on input data.
  # - Add more option plots to example.
  # - PC2 axis label could be in a more sensible position, preferably aligned with rightnand side of bottom platform (dunno how to make this work on a rescalable plot) - Looks like a separate vertical axis currently (actually time!)
  # - How to allow a phylomorphospace options?
  # - Allow Wills style z-ais maybe?
  # - Other options for treatment of age data than just ranges?
  # - Check axes asked for exist in data (top level conditional).
  # - Add geologic time at left with geoscale at some point?
  # - Am assuming PC axes, but could be Relative Warp...
  # - Maybe move this all to plot3D in future (or rgl or an optional choice between them)
  # - Option to plot isometric style view so bottom corner Vs away from viewer? (Will require mainly an edit to translation function, but also y_addition)
  # - Do named time bins (matrix with rownames rather than vector) instead and push this throughout the rest of the code.
  
  # Check time_bins is in a valid format and stop and warn user if not:
  if (!is.timeBins(x = time_bins)) stop(check_timeBins(time_bins = time_bins))
  
  # If not a valid taxonGroups object then stop and provide feedback to user on what is wrong:
  if (!is.taxonGroups(x = taxon_groups)) stop(check_taxonGroups(taxon_groups = x)[1])
  
  # Subfunction to translate input ordination coordinates to stack plotting coordinates:
  translate_to_stack_coordinates <- function(x, y, x_range, y_range, shear, n_stacks, platform_size) {
    if (length(x = x) > 0) {
      x <- ((x - x_range[1]) / diff(x = x_range) * 100) * (1 - shear)
      y <- ((y - y_range[1]) / diff(x = y_range) * ((100 / n_stacks) * platform_size))
      x <- x + (y / ((100 / n_stacks) * platform_size)) * shear * 100 # Second part of shearing (smearing x-axis of points rightwards)
      translated_ccordinates <- cbind(x = x, y = y)
    } else {
      translated_ccordinates <- matrix(nrow = 0, ncol = 2, dimnames = list(c(), c("x", "y")))
    }
    translated_ccordinates
  }
  
  # Subfunction to make circular polygon around a point:
  make_circle_polygon <- function(x, y, radius, resolution = 100) {
    cbind(x = unlist(x = lapply(X = as.list(1:resolution), FUN = function(z) x + (radius * cos(2 * pi * z / resolution)))), y = unlist(x = lapply(X = as.list(1:resolution), FUN = function(z) y + (radius * sin(2 * pi * z / resolution)))))
  }
  
  # If using taxon_groups set to TRUE, else FALSE:
  taxon_groups_used <- ifelse(test = methods::hasArg(name = "taxon_groups"), yes = TRUE, no = FALSE)
  
  # If using taxon_groups:
  if (taxon_groups_used) {

    # Find any taxa duplicated across taxon_group:
    duplicated_taxa <- sort(x = unique(x = unname(obj = unlist(x = taxon_groups))[duplicated(x = unname(obj = unlist(x = taxon_groups)))]))

    # If these exist stop and warn user:
    if (length(x = duplicated_taxa) > 0) paste0("The following taxa are duplicted in taxon_groups: ", paste(duplicated_taxa, collapse = ", "), ". Taxa can only be in one group.")

  }

  # Get vector of values that correspond to scree plot:
  scree_values <- apply(pcoa_input$vectors, 2, var) / sum(apply(pcoa_input$vectors, 2, var)) * 100

  # Make axis labels:
  x_lab <- paste("PC", x_axis, " (", round(scree_values[x_axis], 2), "% of total variance)", sep = "")
  y_lab <- paste("PC", y_axis, " (", round(scree_values[y_axis], 2), "% of total variance)", sep = "")
  
  # Find the taxa assigned to each time bin:
  taxa_assigned_to_bins <- assign_taxa_to_bins(taxon_ages = taxon_ages, time_bins = time_bins)
  
  # Set default solid colour for each taxon to black:
  solid_colours <- rep(x = "black", length.out = nrow(x = pcoa_input$vectors))
  names(solid_colours) <- rownames(x = pcoa_input$vectors)
  
  # If using taxon groups (need different colours for each one):
  if (taxon_groups_used) {
    
    # Build new solid colours with different colours for each group:
    solid_colours <- unlist(x = lapply(X = as.list(1:length(x = taxon_groups)), function(x) {
      
      # Build vector of group colour:
      group_solid_colours <- rep(x = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1)[x], length.out = length(x = taxon_groups[[x]]))
      
      # Add taxon names to colours:
      names(group_solid_colours) <- taxon_groups[[x]]
      
      # Return new group colours vector:
      group_solid_colours
    }))
    
    # Set transparent colours (for use with plotting convex hulls):
    transparent_colours <- grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 0.5)
  
  # If no taxon_groups were used:
  } else {
    
    # Make single vector of everything so plotting functions below work fine:
    taxon_groups <- list(All = rownames(x = pcoa_input$vectors))
  }
  
  # Record the number of stacks to plot:
  n_stacks <- nrow(x = time_bins)
  
  # Set default ranges:
  x_range <- range(pcoa_input$vectors[, x_axis])
  y_range <- range(pcoa_input$vectors[, y_axis])
  
  # Add plot cushion to ranges:
  if (plot_cushion > 0) {
    x_range[1] <- x_range[1] - (diff(x = x_range) * plot_cushion) / 2
    x_range[2] <- x_range[2] + (diff(x = x_range) * plot_cushion) / 2
    y_range[1] <- y_range[1] - (diff(x = y_range) * plot_cushion) / 2
    y_range[2] <- y_range[2] + (diff(x = y_range) * plot_cushion) / 2
  }
  
  # Set some margins before plotting:
  graphics::par(mar = c(0, 0, 0, 1), oma = c(1, 1, 1, 1))
  
  # Create basic (empty) plot)
  graphics::plot(0:100, 0:100, type = "n", axes = FALSE, xlab = "", ylab = "")
  
  # Calsulate platform_addition value (used to define position of "front" edge of each platform):
  platform_addition <- (100 - ((100 / n_stacks) * platform_size)) / (n_stacks - 1)
  
  # If plotting pillars:
  if (plot_pillars) {
    
    # Set pillar coordintaes (base points to start with):
    pillar_coordinates <- list(
      front_left = translate_to_stack_coordinates(x = c(x_range[1], x_range[1]), y = c(y_range[1], y_range[1]), x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size),
      front_right = translate_to_stack_coordinates(x = c(x_range[2], x_range[2]), y = c(y_range[1], y_range[1]), x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size),
      back_left = translate_to_stack_coordinates(x = c(x_range[1], x_range[1]), y = c(y_range[2], y_range[2]), x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size),
      back_right = translate_to_stack_coordinates(x = c(x_range[2], x_range[2]), y = c(y_range[2], y_range[2]), x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)
    )
    
    # Update second point to join to top platform:
    pillar_coordinates <- lapply(X = pillar_coordinates, FUN = function(z) {
      z[2, "y"] <- z[2, "y"] + 100 - ((100 / n_stacks) * platform_size)
      z
    })
    
    # Plot back pillars:
    graphics::lines(x = pillar_coordinates$back_left[, "x"], y = pillar_coordinates$back_left[, "y"])
    graphics::lines(x = pillar_coordinates$back_right[, "x"], y = pillar_coordinates$back_right[, "y"])
  }
  
  # Set the point_size (uses 100th the size of the space times a user chosen factor):
  point_size <- (max(diff(x = x_range), diff(x = y_range)) / 100) * point_size
  
  # Define grid size:
  grid_size <- signif(x = min(diff(x = x_range), diff(x = y_range)) / 5, digits = 1)
  
  # Generate base grid sequence:
  x_grid_sequence <- seq(from = floor(x_range[1]), to = ceiling(x_range[2]), by = grid_size)
  y_grid_sequence <- seq(from = floor(y_range[1]), to = ceiling(y_range[2]), by = grid_size)
  
  # Shift grids to include zero:
  x_grid_sequence <- x_grid_sequence - x_grid_sequence[abs(x = x_grid_sequence) == min(abs(x = x_grid_sequence))]
  y_grid_sequence <- y_grid_sequence - y_grid_sequence[abs(x = y_grid_sequence) == min(abs(x = y_grid_sequence))]
  
  # Prune grid cell edges that fall outside plot area:
  x_grid_sequence <- x_grid_sequence[intersect(x = which(x_grid_sequence > x_range[1]), y = which(x_grid_sequence < x_range[2]))]
  y_grid_sequence <- y_grid_sequence[intersect(x = which(y_grid_sequence > y_range[1]), y = which(y_grid_sequence < y_range[2]))]

  # Main data plotting part:
  x <- lapply(X = as.list(x = 1:length(x = taxa_assigned_to_bins)), FUN = function(y) {
    
    # Build list of plot parameters (to plot successuvely with lapply later):
    plot_parameters <- list(
    
      # Store stack number (to use later for y offset for each platform):
      number = y,
      
      # Store time bin name for potential plotting later:
      bin_name = names(x = taxa_assigned_to_bins)[y],
      
      # Build coordinates for platform base:
      platform_coordinates = translate_to_stack_coordinates(x = c(x_range[1], x_range[2], x_range[2], x_range[1]), y = c(y_range[1], y_range[1], y_range[2], y_range[2]), x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size),
      
      # Build coordinates for crosshairs:
      vertical_crosshair = translate_to_stack_coordinates(x = c(0, 0), y = y_range, x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size),
      horizontal_crosshair = translate_to_stack_coordinates(x = x_range, y = c(0, 0), x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size),
      
      # Build gridlines:
      vertical_gridlines = lapply(X = as.list(x_grid_sequence), FUN = function(z) translate_to_stack_coordinates(x = c(z, z), y = y_range, x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)),
      horizontal_gridlines = lapply(X = as.list(y_grid_sequence), FUN = function(z) translate_to_stack_coordinates(x = x_range, y = c(z, z), x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)),

      # Build point coordinates (really the polygon for each point):
      point_coordinates = lapply(X = taxon_groups, FUN = function(z) {
        
        # Find the taxa present in the time bin:
        taxa_in_bin <- as.list(x = intersect(x = z, y = taxa_assigned_to_bins[[y]]))
        
        # Build coordinates for each point in the bin:
        output <- lapply(X = taxa_in_bin, FUN = function(q) {translate_to_stack_coordinates(x = make_circle_polygon(x = pcoa_input$vectors[q, x_axis], y = pcoa_input$vectors[q, y_axis], radius = point_size)[, "x"], y = make_circle_polygon(x = pcoa_input$vectors[q, x_axis], y = pcoa_input$vectors[q, y_axis], radius = point_size)[, "y"], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)})
        
        # Add taxon names to output:
        names(output) <- taxa_in_bin
        
        # Return output to list:
        output
      }),
      
      # Build point coordinates (really the polygon for each point):
      group_polygon_coordinates = lapply(X = taxon_groups, FUN = function(z) {
        
        # Find the taxa present in the time bin:
        taxa_in_bin <- intersect(x = z, y = taxa_assigned_to_bins[[y]])
        
        # Build coordinates for each point in the bin:
        translate_to_stack_coordinates(x = pcoa_input$vectors[taxa_in_bin, x_axis], y = pcoa_input$vectors[taxa_in_bin, y_axis], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)
      }),

      # Build coordinates to plot time bin names:
      text_coordinates = translate_to_stack_coordinates(x = x_range[2], y = y_range[2], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)
    )
    
    # Set up y_addition value so platforms are stacked at correct height:
    y_addition = (plot_parameters$number - 1) * platform_addition
    
    # Plot platform base:
    graphics::polygon(x = plot_parameters$platform_coordinates[, "x"], y = plot_parameters$platform_coordinates[, "y"] + y_addition, col = "white")
    
    # If plotting crosshair add 0,0 lines to platform:
    if (plot_crosshair) graphics::lines(x = plot_parameters$vertical_crosshair[, "x"], y = plot_parameters$vertical_crosshair[, "y"] + y_addition, lwd = 0.5)
    if (plot_crosshair) graphics::lines(x = plot_parameters$horizontal_crosshair[, "x"], y = plot_parameters$horizontal_crosshair[, "y"] + y_addition, lwd = 0.5)
    
    # Plot grid cells if user requested them:
    if (plot_grid_cells) plotted_gridlines <- lapply(X = plot_parameters$vertical_gridlines, FUN = function(z) graphics::lines(x = z[, "x"], y = z[, "y"] + y_addition, lwd = 0.5, lty = 3))
    if (plot_grid_cells) plotted_gridlines <- lapply(X = plot_parameters$horizontal_gridlines, FUN = function(z) graphics::lines(x = z[, "x"], y = z[, "y"] + y_addition, lwd = 0.5, lty = 3))
    
    # If user has requested convex hulls (and their are groups that make these logical):
    if (plot_convex_hulls && taxon_groups_used) {
      
      # Go through group polygons:
      lapply(X = as.list(1:length(x = plot_parameters$group_polygon_coordinates)), FUN = function(z) {
     
        # Work out which points form the hull:
        convex_hull_points <- grDevices::chull(x = plot_parameters$group_polygon_coordinates[[z]][, "x"], y = plot_parameters$group_polygon_coordinates[[z]][, "y"])
     
        # Plot convex hulls in transparent colour:
        graphics::polygon(x = plot_parameters$group_polygon_coordinates[[z]][convex_hull_points, "x"], y = plot_parameters$group_polygon_coordinates[[z]][convex_hull_points, "y"] + y_addition, col = transparent_colours[z], border = transparent_colours[z])
      })
    }

    # Plot points as coloured polygons:
    plotted_points <- lapply(X = plot_parameters$point_coordinates, FUN = function(z) lapply(X = z, FUN = function(q) graphics::polygon(x = q[, "x"], y = q[, "y"] + y_addition, border = "black", lwd = 0.5, col = solid_colours[names(z)])))
    
    # If plotting time bin names, add these at back right of each platform:
    if (plot_timebin_names) graphics::text(x = plot_parameters$text_coordinates[, "x"] - 1, y = plot_parameters$text_coordinates[, "y"] + y_addition + 2, labels = plot_parameters$bin_name, adj = 1, xpd = NA)
  })
  
  # If plotting pillars:
  if (plot_pillars) {
    
    # Plot front pillars:
    graphics::lines(x = pillar_coordinates$front_left[, "x"], y = pillar_coordinates$front_left[, "y"])
    graphics::lines(x = pillar_coordinates$front_right[, "x"], y = pillar_coordinates$front_right[, "y"])
  }
  
  # If user requested tickmarks:
  if (plot_tickmarks) {
    
    # Make starting x tickmarks:
    x_tick_min <- translate_to_stack_coordinates(x = c(x_range[1], x_range[1]), y = y_range, x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)
    x_tick_max <- translate_to_stack_coordinates(x = c(x_range[2], x_range[2]), y = y_range, x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)
    
    # Determine rescale factor required:
    x_mark_rescale_factor <- 1.5 / sqrt(x = diff(x = x_tick_max[, "x"]) ^ 2 + diff(x = x_tick_max[, "y"]) ^ 2)
    
    # Dtermine how to shift values for x tick makrs:
    x_mark_shift_values <- (x_mark_rescale_factor * x_tick_min[2, ] - x_tick_min[1, ])
    
    # Shift valus for x tick marks:
    x_tick_min[2, ] <- x_tick_min[1, ] - x_mark_shift_values
    x_tick_max[2, ] <- x_tick_max[1, ] - x_mark_shift_values
    
    # Determine y tickmarks:
    y_tick_min <- rbind(translate_to_stack_coordinates(x = x_range[2], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size), c(translate_to_stack_coordinates(x = x_range[2], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "x"] + 1.5, translate_to_stack_coordinates(x = x_range[2], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "y"]))
    y_tick_max <- rbind(translate_to_stack_coordinates(x = x_range[2], y = y_range[2], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size), c(translate_to_stack_coordinates(x = x_range[2], y = y_range[2], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "x"] + 1.5, translate_to_stack_coordinates(x = x_range[2], y = y_range[2], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "y"]))
    
    # Plot tickmark lines:
    graphics::lines(x = x_tick_min[, "x"], y = x_tick_min[, "y"])
    graphics::lines(x = x_tick_max[, "x"], y = x_tick_max[, "y"])
    graphics::lines(x = y_tick_min[, "x"], y = y_tick_min[, "y"])
    graphics::lines(x = y_tick_max[, "x"], y = y_tick_max[, "y"])
    
    # Plot tickmark values:
    graphics::text(x = x_tick_min[2, "x"], y = x_tick_min[2, "y"], pos = 1, labels = as.character(x = signif(x = x_range[1], digits = 2)), xpd = NA)
    graphics::text(x = x_tick_max[2, "x"], y = x_tick_max[2, "y"], pos = 1, labels = as.character(x = signif(x = x_range[2], digits = 2)), xpd = NA)
    graphics::text(x = y_tick_min[2, "x"], y = y_tick_min[2, "y"], pos = 4, labels = as.character(x = signif(x = y_range[1], digits = 2)), xpd = NA)
    graphics::text(x = y_tick_max[2, "x"], y = y_tick_max[2, "y"], pos = 4, labels = as.character(x = signif(x = y_range[2], digits = 2)), xpd = NA)
  }
  
  # If plotting a group legend:
  if (taxon_groups_used && plot_group_legend) {
    
    # Amount to nudge legend away from croners of plot:
    legend_adjustment <- 1
    
    # Add group legend to plot:
    if(group_legend_position == "bottom_left") graphics::legend(x = translate_to_stack_coordinates(x = x_range[1], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "x"] + legend_adjustment, y = translate_to_stack_coordinates(x = x_range[1], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "y"] + legend_adjustment, legend = names(taxon_groups), fill = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1), bg = "white", xjust = 0, yjust = 0)
    if(group_legend_position == "bottom_right") graphics::legend(x = translate_to_stack_coordinates(x = x_range[2], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "x"] - legend_adjustment, y = translate_to_stack_coordinates(x = x_range[2], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "y"] + legend_adjustment, legend = names(taxon_groups), fill = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1), bg = "white", xjust = 1, yjust = 0)
    if(group_legend_position == "top_left") graphics::legend(x = translate_to_stack_coordinates(x = x_range[1], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "x"] + legend_adjustment, y = translate_to_stack_coordinates(x = x_range[1], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "y"] + (platform_addition * (n_stacks - 1)) - legend_adjustment, legend = names(taxon_groups), fill = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1), bg = "white", xjust = 0, yjust = 1)
    if(group_legend_position == "top_right") graphics::legend(x = translate_to_stack_coordinates(x = x_range[2], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "x"] - legend_adjustment, y = translate_to_stack_coordinates(x = x_range[2], y = y_range[1], x_range = x_range, y_range = y_range, shear = shear, n_stacks = n_stacks, platform_size = platform_size)[, "y"] + (platform_addition * (n_stacks - 1)) - legend_adjustment, legend = names(taxon_groups), fill = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1), bg = "white", xjust = 1, yjust = 1)
  }

  # Plot axis labels:
  graphics::text(x = (100 * (1 - shear)) / 2, y = -0.5, pos = 1, srt = 0, labels = x_lab, xpd = NA) # Plots in middle front of bottom platform
  graphics::text(x = 103, y = 100 - ((100 - (100 / n_stacks) * platform_size) / 2), pos = 1, srt = 90, labels = y_lab, xpd = NA) # Plots at mid-height of back of stack
}
