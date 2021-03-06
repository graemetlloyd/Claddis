#' Assign taxa to time bins
#'
#' @description
#'
#' Given a set of first and last appearances assigns a set of taxa to a series of time bins.
#'
#' @param taxon_ages A matrix of taxon ages, with columns for first (\code{"fad"}) and last (\code{"lad"}) appearances and rownames correspodning to taxon names.
#' @param time_bins An object of class \code{timeBins}.
#'
#' @details
#'
#' The various disparity plotting functions (\link{plot_chronophylomorphospace}, \link{plot_morphospace_stack}, \link{plot_morphospace}, \link{plot_multi_morphospace}) are designed to allow assignment of taxa to named groups so that these groups may be assigned different colours when plotting. One way taxa may be grouped is temporally, by assignment to a series of time bins.
#'
#' There are many ways this may be automated and this function provides a very simple one: if the first and last appearance dates of a taxon overlap with a time bin then it can be assigned to that time bin. (In practice, taxa often have multiple occurrences with "ranges" that really represent uncertainty around their true age.)
#'
#' Note that it is recommended that time bins be named without special characters beyond letters and underscores.
#'
#' @return
#'
#' An object of class \code{taxonGroups}.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{plot_chronophylomorphospace}, \link{plot_morphospace_stack}, \link{plot_morphospace}, \link{plot_multi_morphospace}, \link{ordinate_cladistic_matrix}
#'
#' @examples
#'
#' # Build example time bins:
#' time_bins <- matrix(data = c(443.8, 358.9, 358.9, 298.9, 298.9, 251.9,
#'   251.9, 201.3, 201.3, 145.0, 145.0, 65.5, 65.5, 23.03), ncol = 2,
#'   byrow = TRUE, dimnames = list(c("Silurodevonian", "Carboniferous",
#'   "Permian", "Triassic", "Jurassic", "Cretaceous", "Paleogene"),
#'   c("fad", "lad")))
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Build example taxon ages:
#' taxon_ages <- matrix(data = c(385.3, 374.5, 407, 374.5, 251, 228, 385.3,
#'   251, 251, 251, 391.8, 251, 251, 228, 385.3, 391.8, 391.8, 385.3, 311.7,
#'   359.2, 359.2, 416, 407, 407, 407, 407, 385.3, 397.5, 385.3, 161.2, 385.3,
#'   345.3, 318.1, 385.3, 228, 385.3, 385.3, 385.3, 385.3, 385.3, 385.3, 385.3,
#'   385.3, 385.3, 391.8, 407, 391.8, 374.5, 407, 70.6, 311.7, 407, 145.5, 251,
#'   65.5, 251, 112, 374.5, 374.5, 374.5, 385.3, 311.7, 249.7, 359.2, 391.8,
#'   374.5, 385.3, 83.5, 418.7, 251, 385.3, 391.8, 374.5, 345.3, 385.3, 385.3,
#'   407, 411.2, 397.5, 345.3, 374.5, 407, 216.5, 326.4, 411.2, 411.2, 374.5,
#'   359.2, 391.8, 359.2, 245, 216.5, 374.5, 245, 245, 245, 385.3, 245, 245,
#'   199.6, 374.5, 385.3, 385.3, 374.5, 306.5, 345.3, 345.3, 411.2, 397.5,
#'   397.5, 397.5, 397.5, 374.5, 391.8, 374.5, 145.5, 374.5, 326.4, 311.7,
#'   374.5, 199.6, 374.5, 374.5, 374.5, 374.5, 374.5, 374.5, 374.5, 374.5,
#'   374.5, 385.3, 397.5, 385.3, 359.2, 397.5, 65.5, 306.5, 397.5, 99.6, 245,
#'   23.03, 245, 99.6, 359.2, 359.2, 359.2, 374.5, 306.5, 247.4, 318.1, 385.3,
#'   359.2, 374.5, 70.6, 416, 250.4, 374.5, 385.3, 359.2, 326.4, 374.5, 374.5,
#'   397.5, 407, 391.8, 326.4, 359.2, 397.5, 203.6, 318.1, 407, 407),
#'   ncol = 2, dimnames = list(c("Adololopas_moyasmithae", "Adelargo_schultzei",
#'   "Amadeodipterus_kencampbelli", "Andreyevichthys_epitomus",
#'   "Aphelodus_anapes", "Archaeoceratodus_avus", "Archaeonectes_pertusus",
#'   "Arganodus_atlantis", "Ariguna_formosa", "Asiatoceratodus_sharovi",
#'   "Barwickia_downunda", "Beltanodus_ambilobensis", "Ceratodus_formosa",
#'   "Ceratodus_latissimus", "Chirodipterus_australis",
#'   "Chirodipterus_onawwayensis", "Chirodipterus_rhenanus",
#'   "Chirodipterus_wildungensis", "Conchopoma_gadiforme", "Ctenodus_romeri",
#'   "Delatitia_breviceps", "Diabolepis_speratus", "Dipnorhynch_cathlesae",
#'   "Dipnorhynchus_sussmilchi", "Dipnorhynchus_kiandrensis",
#'   "Dipnorhynchus_kurikae", "Dipterus_cf_valenciennesi",
#'   "Dipterus_valenciennesi", "Eoctenodus_microsoma",
#'   "Ferganoceratodus_jurassicus", "Fleurantia_denticulata",
#'   "Ganopristodus_splendens", "Gnathorhiza_serrata", "Gogodipterus_paddyensis",
#'   "Gosfordia_truncata", "Griphognathus_minutidens", "Griphognathus_sculpta",
#'   "Griphognathus_whitei", "Grossipterus_crassus", "Holodipterus_elderae",
#'   "Holodipterus_gogoensis", "Robinsondipterus_longi",
#'   "Asthenorhynchus_meemannae", "Holodipterus_santacrucensis",
#'   "Howidipterus_donnae", "Ichnomylax_kurnai", "Iowadipterus_halli",
#'   "Jarvikia_arctica", "Jessenia_concentrica", "Lepidosiren_paradoxa",
#'   "Megapleuron_zangerli", "Melanognathus_canadensis",
#'   "Metaceratodus_wollastoni", "Microceratodus_angolensis",
#'   "Mioceratodus_gregoryi", "Namatozodia_pitikanta", "Neoceratodus_forsteri",
#'   "Nielsenia_nordica", "Oervigia_nordica", "Orlovichthys_limnatis",
#'   "Palaeodaphus_insignis", "Palaeophichthys_parvulus", "Paraceratodus_germaini",
#'   "Parasagenodus_sibiricus", "Pentlandia_macroptera",
#'   "Phaneropleuron_andersoni", "Pillararhynchus_longi", "Protopterus_annectens",
#'   "Psarolepis_romeri", "Ptychoceratodus_serratus", "Rhinodipterus_secans",
#'   "Rhinodipterus_ulrichi", "Rhynchodipterus_elginensis", "Sagenodus_inaequalis",
#'   "Scaumenacia_curta", "Soederberghia_groenlandica",
#'   "Sorbitorhynchus_deleaskitus", "Speonesydrion_iani", "Stomiahykus_thlaodus",
#'   "Straitonia_waterstoni", "Sunwapta_grandiceps", "Tarachomylax_oepiki",
#'   "Tellerodus_sturi", "Tranodis_castrensis", "Uranolophus_wyomingensis",
#'   "Westollrhynchus_lehmanni"), c("fad", "lad")))
#'
#' # Assign taxa to time bins:
#' assign_taxa_to_bins(taxon_ages = taxon_ages, time_bins = time_bins)
#'
#' @export assign_taxa_to_bins
assign_taxa_to_bins <- function(taxon_ages, time_bins) {
  
  # CHECK FOR TAXA NOT IN RANGE OF TIME BINS (SMALLEST TAXON FAD LESS THAN OLDEST TIME BIN LAD AND VICE VERSA
  # ENSURE TIME BINS DO NOT OVERLAP AND ARE CONTIGUOUS
  
  # Check time_bins is in a valid format and stop and warn user if not:
  if (!is.timeBins(x = time_bins)) stop(check_timeBins(time_bins = time_bins))
  
  # Ensure column names are lower case so will get called correctly later
  # (for consistency with paleotree and strap formats):
  colnames(x = taxon_ages) <- tolower(x = colnames(x = taxon_ages))
  colnames(x = time_bins) <- tolower(x = colnames(x = time_bins))
  
  # Assign taxa to bins:
  taxa_assigned_to_bins <- lapply(X = as.list(x = rownames(x = time_bins)), FUN = function(x) rownames(x = taxon_ages)[(c(taxon_ages[, "fad"] > time_bins[x, "lad"]) + c(taxon_ages[, "lad"] < time_bins[x, "fad"])) == 2])
  
  # Add group names from time bins:
  names(x = taxa_assigned_to_bins) <- rownames(x = time_bins)
  
  # Set class as taxonGroups:
  class(x = taxa_assigned_to_bins) <- "taxonGroups"
  
  # Return taxa assigned to bins:
  taxa_assigned_to_bins
}
