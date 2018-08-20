#==============================================================================
# Ratio Metrics
# These metrics should be avoided but are included due to their popularity
# in the literature.
#==============================================================================
#'Ratio of EPT to Chironomidae
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The ratio of EPT (Ephemeroptera, Plecoptera, and Trichoptera) to
#'chironomids (Family: Chironomidae). One is added to the count of EPT
#'and to the count of chironomids to avoid issues associate with the absence
#'of EPT or chironomids.  If no EPT individuals are observed the ratio
#'(0/# chironomids) will return 0 despite the number of chironomids observed.
#'If no chironomids were observed the ratio (# EPT/0) will return an inproper
#' fraction despite the number of EPT observed. I do not recommend using ratios to
#'develop indices of biotic integrity.  However, I included this metric becuase it
#'occurs frequently in the literature.  This metric is not included in the
#' all_metrics function.
#'@export

ratio_ept_chiro <- function(order.wide, family.wide) {
  EPT <- blank_col("EPHEMEROPTERA", order.wide) +
    blank_col("PLECOPTERA", order.wide) +
    blank_col("TRICHOPTERA", order.wide)
  chiro <- blank_col("CHIRONOMIDAE", family.wide[, 6:ncol(family.wide)])
  return((EPT + 1) / (shiro + 1))
}


