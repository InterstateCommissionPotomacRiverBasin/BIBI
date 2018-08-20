#==============================================================================
# Calculate All Metrics
#==============================================================================
#'Specific Metrics
#'
#'@param master.df = The master taxa list contains taxonomic ranks from Phylum
#'to species and known taxonomic attributes.
#'@param long.df = Taxonomic counts in a long data format.
#'@param taxa.rank = The lowest taxonomic rank ("ORDER", "FAMILY", or "GENUS")
#' used to calculate the metrics.  If the majority of your taxa are identified
#' to the family level, then it would be inappropriate to perform metric
#' calculations at the genus level.
#'@param rare = Sample counts can be rarefied or unaltered. The default,
#' "NON_RARE," does not alter sample counts. "ALL_RARE" and "DIV_RARE" will
#' rarefy sample counts to the pecified sample size (sample.size) using
#' expected rarefaction. "ALL_RARE" use the rarefied counts to calulate all
#' of the metrics. "DIV_RARE" will only use the rarefied sample counts to
#' calculate richness and diversity metrics, all other metrics will be calculated
#' with the original counts.  "VEG_RARE" will rarefy the data using the vegan
#' package function rrarefy.  rrarefy produces a random rarefied sample. I
#' argue that we can use hypergeometric formulas to provide the expected number
#' of each taxon that will be selected at specified sample size.  "ALL_RARE" and
#' "DIV_RARE" may produce sample counts abouve and below the specifed sample size.
#' This is a result of rounding the taxonomic counts to the nearest integer. The
#' expected rarefied value is rarely produces an integer and it is inappropriate
#' to report a portion of an individual (e.g., 1.52 individuals).  Therefore,
#' the final expected rarefied value is rounded to the nearest integer.
#'@param sample.size = The rarefaction subsample size.  The sample size should be
#'smaller than the majority of the sample counts.  If a sample count is less than
#'or equal to the specified sample size, the sample will not be altered.
#'@param bibi.standard = Indicate TRUE if you would like to standardize your
#'taxa following the practices used in the Chessie BIBI.  The default, FALSE,
#'will not attempt to standardize the taxa.
#'@param seed = if true set the seed for reproducible results.
#'@return Calculates all of applicable and available metrics in the package.
#'@export
#'

specific_metrics <- function(master.df, long.df, taxa.rank = "FAMILY", rare = "NON_RARE",
                        sample.size = 100, pct_un = NULL, bibi.standard = FALSE,
                        seed = FALSE, metrics.vec) {
  #Prep==========================================================================
  # Wide format data frames for necessary taxonomic levels.
  # Not rarefied
  if(rare %in% c("DIV_RARE", "NON_RARE")) class.wide <- wide(long.df, "CLASS")
  if(rare %in% c("DIV_RARE", "NON_RARE")) order.wide <- wide(long.df, "ORDER")
  if(rare %in% c("DIV_RARE", "NON_RARE")) family.wide <- wide(long.df, "FAMILY")
  if(rare %in% c("DIV_RARE", "NON_RARE") & taxa.rank %in% c("TRIBE", "GENUS", "SPECIES")){
    tribe.wide <- wide(long.df, "TRIBE")
  }else{
    tribe.wide <- NULL
  }
  
  if(rare %in% c("DIV_RARE", "NON_RARE") & taxa.rank %in% c("GENUS", "SPECIES")){
    genus.wide <- wide(long.df, "GENUS")
  }else{
    genus.wide <- NULL
  }
  # Specified Taxonomic level
  if(rare %in% c("DIV_RARE", "NON_RARE") & taxa.rank == "FAMILY"){
    taxa.rank.df <- family.wide
  }
  if(rare %in% c("DIV_RARE", "NON_RARE") & taxa.rank == "GENUS"){
    taxa.rank.df <- genus.wide
  }
  if(rare %in% c("DIV_RARE", "NON_RARE") & !(taxa.rank %in% c("FAMILY", "GENUS"))){
    taxa.rank.df <- wide(long.df, taxa.rank)
  }
  
  
  # Wide format rarefied
  if (seed == TRUE) set.seed(64)
  if(rare %in% "VEG_RARE") rare.long <- prep_rrarefy(long.df, taxa.rank, pct_un)
  if (seed == TRUE) set.seed(64)
  if(rare %in% c("DIV_RARE", "ALL_RARE")) rare.long <- prep_rare(long.df, master.df, taxa.rank, sample.size,
                                                                 pct_un, bibi.standard)
  if(rare %in% c("DIV_RARE", "ALL_RARE", "VEG_RARE")) rare.class.wide <- wide(rare.long, "CLASS")
  if(rare %in% c("DIV_RARE", "ALL_RARE", "VEG_RARE")) rare.order.wide <- wide(rare.long, "ORDER")
  if(rare %in% c("DIV_RARE", "ALL_RARE", "VEG_RARE")) rare.family.wide <- wide(rare.long, "FAMILY")
  if(rare %in% c("DIV_RARE", "ALL_RARE", "VEG_RARE") & taxa.rank %in% c("TRIBE", "GENUS", "SPECIES")){
    rare.tribe.wide <- wide(rare.long, "TRIBE")
  }else{
    rare.tribe.wide <- NULL
  }
  if(rare %in% c("DIV_RARE", "ALL_RARE", "VEG_RARE") & taxa.rank %in% c("GENUS", "SPECIES")){
    rare.genus.wide <- wide(rare.long, "GENUS")
  }else{
    rare.genus.wide <- NULL
  }
  
  # Specified Taxonomic level
  if(rare %in% c("DIV_RARE", "ALL_RARE", "VEG_RARE") & taxa.rank == "FAMILY"){
    rare.df <- rare.family.wide
  }
  if(rare %in% c("DIV_RARE", "ALL_RARE", "VEG_RARE") & taxa.rank == "GENUS"){
    rare.df <- rare.genus.wide
  }
  if(rare %in% c("DIV_RARE", "ALL_RARE", "VEG_RARE") & !(taxa.rank %in% c("FAMILY", "GENUS"))){
    rare.df <- wide(rare.long, taxa.rank)
  }
  
  
  #============================================================================
  print("Diversity Metrics")
  # Calculate diversity metrics
  if(rare %in% c("DIV_RARE", "ALL_RARE", "VEG_RARE")){
    all.div <- specific_diversity(master.df, rare.long,
                                  metrics.vec,  rare.class.wide,
                                  rare.order.wide, rare.family.wide, 
                                  rare.tribe.wide, rare.genus.wide,
                                  rare.df, taxa.rank)
  }else{
    if(rare %in% "NON_RARE"){
      all.div <- specific_diversity(master.df, long.df,
                                    metrics.vec,  class.wide,
                                    order.wide, family.wide, tribe.wide,
                                    genus.wide, taxa.rank.df, taxa.rank)
    }
  }
  
  #============================================================================
  print("Non-Diversity Metrics")
  # Calculate non-diversity metrics
  if(rare %in% c("ALL_RARE", "VEG_RARE")){
    all.non_div <- specific_non_diversity(master.df, rare.long, metrics.vec,
                                     rare.class.wide, rare.order.wide,
                                     rare.family.wide, rare.tribe.wide,
                                     rare.genus.wide, rare.df, taxa.rank)
  }else{
    if(rare %in% c("DIV_RARE", "NON_RARE")){
      all.non_div <- specific_non_diversity(master.df, long.df, metrics.vec, 
                                       class.wide, order.wide, family.wide,
                                       tribe.wide, genus.wide, taxa.rank.df,
                                       taxa.rank)
    }
  }
  
  #============================================================================
  # Merge the diversity metrics data frame with the non-diversity
  # metrics data frame
  merged <- merge(all.div, all.non_div, by = c("EVENT_ID", "STATION_ID",
                                               "DATE", "SAMPLE_NUMBER",
                                               "AGENCY_CODE"))
  #===========================================================================
  
  return(merged)
}

#==============================================================================
#'Run All Diversity Metrics
#'
#'@param master.df = The master taxa list containing all of the appropriate
#'taxonomic resolutions and the appropriate taxa attributes. A master taxa list
#'is contained with the BIBI package.
#'@param long.df = Taxonomic data in long format
#'@param metrics.vec = the vector of specified metrics to calculate.
#'@param class.wide = Class level wide taxonomic data frame
#'@param order.wide = Order level wide taxonomic data frame
#'@param family.wide = Family level wide taxonomic data frame
#'@param taxa.rank.df = Specified taxonomic level wide taxonomic data frame
#'@param taxa.rank = Taxonomic level ("FAMILY" or "GENUS") (Default = "FAMILY")
#'@return All diversity metric values
#'@export
#'

specific_diversity <- function(master.df, long.df, metrics.vec, 
                               class.wide, order.wide, family.wide, tribe.wide,
                               genus.wide, taxa.rank.df, taxa.rank){
  metrics <- data.frame(taxa.rank.df[, c("EVENT_ID", "STATION_ID", "DATE",
                                         "SAMPLE_NUMBER", "AGENCY_CODE")])
  colnames(metrics) <- c("EVENT_ID", "STATION_ID", "DATE",
                         "SAMPLE_NUMBER", "AGENCY_CODE")
  metrics <- with(metrics,  metrics[order(EVENT_ID), ])
  #============================================================================

  #============================================================================
  if (calc_metric("RICH", metrics.vec)) {
    metrics$RICH <- vegan::specnumber(taxa.rank.df[, 6:ncol(taxa.rank.df)]) 
  } 
  #----------------------------------------------------------------------------
  if (calc_metric("SHANNON", metrics.vec)) {
    metrics$SHANNON <- shannon(taxa.rank.df)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("EFFECTIVE_RICH_SHANNON", metrics.vec)) {
    metrics$EFFECTIVE_RICH_SHANNON <- effective_richness(taxa.rank.df,
                                                         "shannon")
  }
  #----------------------------------------------------------------------------
  if (calc_metric("SIMPSONS", metrics.vec)) {
    metrics$SIMPSONS <- simpsons(taxa.rank.df)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("EFFECTIVE_RICH_SIMPSON", metrics.vec)) {
    metrics$EFFECTIVE_RICH_SIMPSON <- effective_richness(taxa.rank.df,
                                                         "invsimpson")
  }
  #----------------------------------------------------------------------------
  if (calc_metric("HURLBERTS_PIE", metrics.vec)) {
    metrics$HURLBERTS_PIE <- hurlberts_pie(taxa.rank.df)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("MARGALEFS", metrics.vec)) {
    metrics$MARGALEFS <- margalefs(taxa.rank.df)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("MENHINICKS", metrics.vec)) {
    metrics$MENHINICKS <- menhinicks(taxa.rank.df)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PIELOU", metrics.vec)) {
    metrics$PIELOU <-  pielou(taxa.rank.df)
  }
  #============================================================================
  if(taxa.rank %in% c("FAMILY", "GENUS")){
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_EPHEMEROPTERA", metrics.vec)) {
      metrics$RICH_EPHEMEROPTERA <- rich_ephemeroptera(long.df, taxa.rank)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_PLECOPTERA", metrics.vec)) {
      metrics$RICH_PLECOPTERA <- rich_plecoptera(long.df, taxa.rank)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_TRICHOPTERA", metrics.vec)) {
      metrics$RICH_TRICHOPTERA <- rich_trichoptera(long.df, taxa.rank)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_EPT", metrics.vec)) {
      metrics$RICH_EPT <- rich_ept(long.df, taxa.rank)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("PCT_EPT_RICH", metrics.vec)) {
      metrics$PCT_EPT_RICH <- pct_ept_rich(long.df, taxa.rank)
    }
    #============================================================================
    # FFGs
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_COLLECT", metrics.vec)) {
      metrics$RICH_COLLECT <- rich_attribute(taxa.rank.df, master.df,
                                             attribute.column = "BIBI_FFG",
                                             attribute.interest = c("CG", "CF"),
                                             taxa.rank)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_GATHER", metrics.vec)) {
      metrics$RICH_GATHER <- rich_attribute(taxa.rank.df, master.df,
                                            attribute.column = "BIBI_FFG",
                                            attribute.interest = "CG", taxa.rank)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_FILTER", metrics.vec)) {
      metrics$RICH_FILTER <- rich_attribute(taxa.rank.df, master.df,
                                            attribute.column = "BIBI_FFG",
                                            attribute.interest = "CF", taxa.rank)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_SHRED", metrics.vec)) {
      metrics$RICH_SHRED <- rich_attribute(taxa.rank.df, master.df,
                                           attribute.column = "BIBI_FFG",
                                           attribute.interest = "SH", taxa.rank)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_SCRAPE", metrics.vec)) {
      metrics$RICH_SCRAPE <- rich_attribute(taxa.rank.df, master.df,
                                            attribute.column = "BIBI_FFG",
                                            attribute.interest = "SC", taxa.rank)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_PREDATOR", metrics.vec)) {
      metrics$RICH_PREDATOR <- rich_attribute(taxa.rank.df, master.df,
                                              attribute.column = "BIBI_FFG",
                                              attribute.interest = "PR", taxa.rank)
    }
    #==========================================================================
    # Habits
    #--------------------------------------------------------------------------
    if (calc_metric("RICH_CLIMB", metrics.vec)) {
      metrics$RICH_CLIMB <- rich_attribute(taxa.rank.df, master.df,
                                           attribute.column = "BIBI_HABIT",
                                           attribute.interest = "CB", taxa.rank)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("RICH_SWIM", metrics.vec)) {
      metrics$RICH_SWIM <- rich_attribute(taxa.rank.df, master.df,
                                          attribute.column = "BIBI_HABIT",
                                          attribute.interest = "SW", taxa.rank)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("RICH_CLING", metrics.vec)) {
      metrics$RICH_CLING <- rich_attribute(taxa.rank.df, master.df, 
                                           attribute.column = "BIBI_HABIT",
                                           attribute.interest = "CN", taxa.rank)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("RICH_BURROW", metrics.vec)) {
      metrics$RICH_BURROW <- rich_attribute(taxa.rank.df, master.df,
                                            attribute.column = "BIBI_HABIT",
                                            attribute.interest = "BU", taxa.rank)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("RICH_SPRAWL", metrics.vec)) {
      metrics$RICH_SPRAWL <- rich_attribute(taxa.rank.df, master.df,
                                            attribute.column = "BIBI_HABIT",
                                            attribute.interest = "SP", taxa.rank)
    }
    
    #--------------------------------------------------------------------------
    if (calc_metric("RICH_SKATE", metrics.vec)) {
      metrics$RICH_SKATE <- rich_attribute(taxa.rank.df, master.df,
                                           attribute.column = "BIBI_HABIT",
                                           attribute.interest = "SK", taxa.rank)
    }
  }
  #----------------------------------------------------------------------------
  if(taxa.rank %in% "GENUS"){
    if (calc_metric("RICH_NCO", metrics.vec)) {
      metrics$RICH_NCO <- rich_nco(long.df, taxa.rank)
    }
  }
  #============================================================================
  # Tolerance metrics need to be updated once more taxonomic levels are
  # added to the taxa attributes table
  
  if(taxa.rank %in% c("FAMILY", "GENUS")){
    if (calc_metric("RICH_INTOL", metrics.vec)) {
      metrics$RICH_INTOL <- rich_tolerance(taxa.rank.df, master.df, "BIBI_TV", 0, 3)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_MODTOL", metrics.vec)) {
      metrics$RICH_MODTOL <- rich_tolerance(taxa.rank.df, master.df, "BIBI_TV", 4, 6)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_TOL", metrics.vec)) {
      metrics$RICH_TOL <- rich_tolerance(taxa.rank.df, master.df, "BIBI_TV", 7, 10)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("RICH_EPT_NO_TOL", metrics.vec)) {
      metrics$RICH_EPT_NO_TOL <- ept_rich_no_tol(long.df, "FAMILY", master.df)
    }
    #----------------------------------------------------------------------------
    if (calc_metric("PCT_EPT_RICH_NO_TOL", metrics.vec)) {
      metrics$PCT_EPT_RICH_NO_TOL <- pct_ept_rich_no_tol(long.df, "FAMILY", master.df)
    }
  }
  #----------------------------------------------------------------------------
  if(taxa.rank %in% c("FAMILY")){
    if (calc_metric("BECKS_V1", metrics.vec)) {
      metrics$BECKS_V1 <- becks(taxa.rank.df,  taxa.rank, master.df, beck.version = 1)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("BECKS_V3", metrics.vec)) {
      becks(taxa.rank.df,  taxa.rank, master.df, beck.version = 3)
    }
  }
  #----------------------------------------------------------------------------
  if(taxa.rank %in% c("GENUS", "SPECIES")){
    if (calc_metric("RICH_EPHEM_EPEORUS", metrics.vec)) {
      metrics$RICH_EPHEM_EPEORUS <- rich_ephem_epeorus(long.df, genus.wide)
    }
  }
  #----------------------------------------------------------------------------
  Blank <- genus.wide
  # End select_diversity function.
  return(metrics)
  
}

#==============================================================================
#'Run All Non-Diversity Metrics
#'
#'@param master.df = Taxonomic Information
#'@param long.df = Taxonomic data in long format
#'@param metrics.vec = the vector of specified metrics to calculate.
#'@param class.wide = Class level wide taxonomic data frame
#'@param order.wide = Order level wide taxonomic data frame
#'@param family.wide = Family level wide taxonomic data frame
#'@param taxa.rank.df = Specified taxonomic level wide taxonomic data frame
#'@param taxa.rank = Taxonomic level ("FAMILY" or "GENUS") (Default = "FAMILY")
#'@return All non-diversity metric values
#'@export
#'

specific_non_diversity <- function(master.df, long.df, metrics.vec,
                              class.wide, order.wide,
                              family.wide, tribe.wide, genus.wide,
                              taxa.rank.df, taxa.rank){
  metrics <- data.frame(taxa.rank.df[, c("EVENT_ID", "STATION_ID", "DATE",
                                         "SAMPLE_NUMBER", "AGENCY_CODE")])
  colnames(metrics) <- c("EVENT_ID", "STATION_ID", "DATE",
                         "SAMPLE_NUMBER", "AGENCY_CODE")
  metrics <- with(metrics,  metrics[order(EVENT_ID), ])
  
  #============================================================================
  # These should be used for taxa attribute related metrics
  
  taxa <- c("PHYLUM", "SUBPHYLUM", "CLASS",
            "SUBCLASS", "ORDER", "SUBORDER",
            "FAMILY", "SUBFAMILY", "TRIBE",
            "GENUS", "SPECIES")
  
  master.fill <- fill_taxa(master.df)
  master.fill <- unique(master.fill[, c("TSN_R", taxa)])
  #test <- (master.fill[duplicated(master.fill$TSN_R), ])
  long.fill <- long.df
  long.fill <- long.fill[, !names(long.fill) %in% taxa]
  long.fill <- merge(long.fill, master.fill, by.x = "TSN",
                     by.y = "TSN_R", all.x = T)
  
  ord.fill <- wide(long.fill, "ORDER")
  if(taxa.rank %in% c("FAMILY", "GENUS")){
    fam.fill <- wide(long.fill, "FAMILY")
    if(taxa.rank %in% "GENUS") gen.fill <- wide(long.fill, "GENUS")
  }
  
  if(taxa.rank %in% "ORDER") taxa.rank.fill <- ord.fill
  if(taxa.rank %in% "FAMILY") taxa.rank.fill <- fam.fill
  if(taxa.rank %in% "GENUS") taxa.rank.fill <- gen.fill
  #============================================================================
  if (calc_metric("PCT_UNIDENTIFIED", metrics.vec)) {
    metrics$PCT_UNIDENTIFIED <- pct_unidentified(taxa.rank.df)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_DOM1", metrics.vec)) {
    metrics$PCT_DOM1 <- pct_dom(taxa.rank.df, 1)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_DOM2", metrics.vec)) {
    metrics$PCT_DOM2 <- pct_dom(taxa.rank.df, 2)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_DOM3", metrics.vec)) {
    metrics$PCT_DOM3 <- pct_dom(taxa.rank.df, 3)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_DOM4", metrics.vec)) {
    metrics$PCT_DOM4 <- pct_dom(taxa.rank.df, 4)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_DOM5", metrics.vec)) {
    metrics$PCT_DOM5 <- pct_dom(taxa.rank.df, 5)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("ABUNDANCE", metrics.vec)) {
    metrics$ABUNDANCE <- abundance(taxa.rank.df)
  }
  #============================================================================
  # Tolerance Metrics
  #============================================================================
  # Tolerance metrics need to be updated once more tolerance values are added
  # to the taxa attributes table
  print("Calculating Tolerance Metrics")
  if(taxa.rank %in% "FAMILY") {
    if (calc_metric("ASPT_MOD", metrics.vec)) {
      metrics$ASPT_MOD <- tol_index(long.fill, master.df, "ASPT", taxa.rank)
    }
  }
  #----------------------------------------------------------------------------
  if(taxa.rank %in% c("FAMILY", "GENUS")){
    if (calc_metric("HBI", metrics.vec)) {
      metrics$HBI <- tol_index(long.fill, master.df, "BIBI_TV", taxa.rank)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_URBAN_INTOL", metrics.vec)) {
      metrics$PCT_URBAN_INTOL <- pct_urban_intol(long.df, master.df)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_INTOL_0_3", metrics.vec)) {
      metrics$PCT_INTOL_0_3 <- pct_tol_val(taxa.rank.fill, master.df,
                                           "BIBI_TV", 0, 3)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_INTOL_0_4", metrics.vec)) {
      metrics$PCT_INTOL_0_4 <- pct_tol_val(taxa.rank.fill, master.df,
                                           "BIBI_TV", 0, 4)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_MOD_TOL_4_6", metrics.vec)) {
      metrics$PCT_MOD_TOL_4_6 <- pct_tol_val(taxa.rank.fill, master.df,
                                             "BIBI_TV", 4, 6)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_TOLERANT_7_10", metrics.vec)) {
      metrics$PCT_TOLERANT_7_10 <- pct_tol_val(taxa.rank.fill, master.df,
                                               "BIBI_TV", 7, 10)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_TOLERANT_5_10", metrics.vec)) {
      metrics$PCT_TOLERANT_5_10 <- pct_tol_val(taxa.rank.fill, master.df,
                                               "BIBI_TV", 5, 10)
    }
    }
  #============================================================================
  # Composition Metrics
  #============================================================================
  print("Calculating Composition Metrics")
  if(taxa.rank %in% c("FAMILY", "GENUS")){
    if (calc_metric("PCT_EPT_HYDRO_BAETID", metrics.vec)) {
      metrics$PCT_EPT_HYDRO_BAETID <- pct_ept_hydro_baetid(order.wide,
                                                           family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_EPT_NO_HYDRO", metrics.vec)) {
      metrics$PCT_EPT_NO_HYDRO <- pct_ept_no_hydro(order.wide, family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_NON_HYDROP_TRICHOPTERA", metrics.vec)) {
      metrics$PCT_NON_HYDROP_TRICHOPTERA <- pct_non_hydrop_trichoptera(order.wide,
                                                                       family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_HYDRO_TRICHOPTERA", metrics.vec)) {
      metrics$PCT_HYDRO_TRICHOPTERA <- pct_hydro_trichoptera(order.wide,
                                                             family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_EPHEMEROPTERA_NO_BAETID", metrics.vec)) {
      metrics$PCT_EPHEMEROPTERA_NO_BAETID <- pct_epmeroptera_no_baetid(order.wide,
                                                                       family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_HYDRO_EPT", metrics.vec)) {
      metrics$PCT_HYDRO_EPT <- pct_hydro_ept(order.wide, family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_RETREAT_CADDISFLY", metrics.vec)) {
      metrics$PCT_RETREAT_CADDISFLY <- pct_retreat_trichoptera(long.df)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_CORBICULIDAE", metrics.vec)) {
      metrics$PCT_CORBICULIDAE <- pct_corbiculidae(family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_SIMULIIDAE", metrics.vec)) {
      metrics$PCT_SIMULIIDAE <- pct_simuliidae(family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_CHIRONOMIDAE", metrics.vec)) {
      metrics$PCT_CHIRONOMIDAE <- pct_chironomidae(family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_OLIGO_CHIRO", metrics.vec)) {
      metrics$PCT_OLIGO_CHIRO <- pct_oligo_chiro(class.wide, family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_ANNELID_CHIRO", metrics.vec)) {
      phylum.wide <- wide(long.df, "PHYLUM")
      metrics$PCT_ANNELID_CHIRO <- pct_annelid_chiro(phylum.wide, family.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_LIMESTONE", metrics.vec)) {
      metrics$PCT_LIMESTONE <- pct_limestone(order.wide, family.wide)
    }
  }
  
  if(taxa.rank %in% c("GENUS")){
    if (calc_metric("PCT_CC_CHIRO", metrics.vec)) {
      metrics$PCT_CC_CHIRO <- pct_cc_chironomidae(family.wide, genus.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_EPT_CHEUMATOPSYCHE", metrics.vec)) {
      metrics$PCT_EPT_CHEUMATOPSYCHE <- pct_ept_cheumatopsyche(order.wide, genus.wide) 
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_EPT_HYDROPSYCHE", metrics.vec)) {
      metrics$PCT_EPT_HYDROPSYCHE <- pct_ept_hydropsyche(order.wide, genus.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_TANYTARSINI", metrics.vec)) {
      tribe.wide <- wide(long.df, "TRIBE")
      metrics$PCT_TANYTARSINI <- pct_tanytarsini(tribe.wide)
    }
    #--------------------------------------------------------------------------
    if (calc_metric("PCT_ORTHOCLADIINAE", metrics.vec)) {
      metrics$PCT_ORTHOCLADIINAE <- pct_orthocladiinae(long.df)
    }
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_EPT", metrics.vec)) {
    metrics$PCT_EPT <- pct_ept(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_EPHEMEROPTERA", metrics.vec)) {
    metrics$PCT_EPHEMEROPTERA <- pct_ephemeroptera(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_PLECOPTERA", metrics.vec)) {
    metrics$PCT_PLECOPTERA <- pct_plecoptera(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_TRICHOPTERA", metrics.vec)) {
    metrics$PCT_TRICHOPTERA <- pct_trichoptera(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_COLEOPTERA", metrics.vec)) {
    metrics$PCT_COLEOPTERA <- pct_coleoptera(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_ODONATA", metrics.vec)) {
    metrics$PCT_ODONATA <- pct_odonata(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_COTE", metrics.vec)) {
    metrics$PCT_COTE <- pct_cote(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_POTEC", metrics.vec)) {
    metrics$PCT_POTEC <- pct_potec(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_AMPHIPODA", metrics.vec)) {
    metrics$PCT_AMPHIPODA <- pct_amphipoda(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_BIVALVIA", metrics.vec)) {
    metrics$PCT_BIVALVIA <- pct_bivalvia(class.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_UNIONOIDA", metrics.vec)) {
    metrics$PCT_UNIONOIDA <- pct_unionoida(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_DIPTERA", metrics.vec)) {
    metrics$PCT_DIPTERA <- pct_diptera(order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("GOLD", metrics.vec)) {
    metrics$GOLD <- gold(class.wide, order.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_OLIGOCHAETA", metrics.vec)) {
    metrics$PCT_OLIGOCHAETA <- pct_oligochaeta(class.wide)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_NON_INSECT", metrics.vec)) {
    metrics$PCT_NON_INSECT <- pct_non_insect(class.wide)
  }
  #============================================================================
  # Functional Feeding Group (FFG) Metrics
  #============================================================================
  print("Calculating FFG Metrics")
  
  if (calc_metric(c("PCT_COLLECT", "PCT_FILTER", "PCT_GATHER", "PCT_PREDATOR",
                    "PCT_SCRAPE", "PCT_SHRED", "PCT_BURROW", "PCT_CLIMB",
                    "PCT_CLING", "PCT_SKATE", "SPCT_SPRAWL", "PCT_SWIM"),
                  metrics.vec)) {
    #if(!(taxa.rank %in% "FAMILY") | taxa.rank %in% "FAMILY") taxa.rank.att <- "FAMILY"
    if(taxa.rank %in% c("ORDER")) taxa.rank.att <- "ORDER"
    if(taxa.rank %in% c("FAMILY", "SUBFAMILY", "TRIBE")) taxa.rank.att <- "FAMILY"
    # Update once genus attributes added
    if(taxa.rank %in% c("GENUS", "SPECIES")) taxa.rank.att <- "GENUS"
    #metrics$PCT_FFG_UNASSIGNED <- pct_attribute(taxa.rank.df, master.df, "BIBI_FFG", "UA", taxa.rank.att)
    #metrics$PCT_DOM_FFG <- pct_dom1_group(long.df, master.df, "BIBI_FFG", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_COLLECT", metrics.vec)) {
    metrics$PCT_COLLECT <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG",
                                         c("CG", "CF"), taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_FILTER", metrics.vec)) {
    metrics$PCT_FILTER <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG",
                                        "CF", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_GATHER", metrics.vec)) {
    metrics$PCT_GATHER <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG",
                                        "CG", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_PREDATOR", metrics.vec)) {
    metrics$PCT_PREDATOR <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG",
                                          "PR", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_SCRAPE", metrics.vec)) {
    metrics$PCT_SCRAPE <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG",
                                        "SC", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_SHRED", metrics.vec)) {
    metrics$PCT_SHRED <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG",
                                       "SH", taxa.rank.att)
  }
  #============================================================================
  # Habit Metrics
  #============================================================================
  print("Calculating Habit Metrics")
  #metrics$PCT_DOM_Habit <- pct_dom1_group(long.df, master.df, "BIBI_HABIT", taxa.rank.att)
  #metrics$PCT_HABIT_UNASSIGNED <- pct_attribute(taxa.rank.df, master.df, "BIBI_HABIT", "UA", taxa.rank.att)
  
  if (calc_metric("PCT_BURROW", metrics.vec)) {
    metrics$PCT_BURROW <- pct_attribute(taxa.rank.fill, master.df,
                                        "BIBI_HABIT", "BU", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_CLIMB", metrics.vec)) {
    metrics$PCT_CLIMB <- pct_attribute(taxa.rank.fill, master.df,
                                       "BIBI_HABIT", "CB", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_CLING", metrics.vec)) {
    metrics$PCT_CLING <- pct_attribute(taxa.rank.fill, master.df,
                                       "BIBI_HABIT", "CN", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_SKATE", metrics.vec)) {
    metrics$PCT_SKATE <- pct_attribute(taxa.rank.fill, master.df,
                                       "BIBI_HABIT", "SK", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_SPRAWL", metrics.vec)) {
    metrics$PCT_SPRAWL <- pct_attribute(taxa.rank.fill, master.df,
                                        "BIBI_HABIT", "SP", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  if (calc_metric("PCT_SWIM", metrics.vec)) {
    metrics$PCT_SWIM <- pct_attribute(taxa.rank.fill, master.df,
                                      "BIBI_HABIT", "SW", taxa.rank.att)
  }
  #----------------------------------------------------------------------------
  return(metrics)
}

#==============================================================================
#'Should the Metric be Calculated?
#'
#'@param metric.name = the metric in question.
#'@param metrics.vec = the vector of specified metrics to calculate.
#'@return This function is used within the all_metrics function and is not 
#'necissarily useful outside of those functions. The purpose of the function is
#'to check that the metric in question (metric.name) is present within the
#'specifiec list of metrics to calculate (metrics.vec). If metrics.vec is 
#'set to "ALL" then all of the metrics will be calculated. This function was
#'created to speed up the calculation of metrics. Often calculating all of the
#'metrics is only useful during the exploratory stage of developing an index. 
#'Calculating all of the metrics, especially for large data sets, can be time 
#'consuming. This function helps to skip over metrics that the user is not 
#'interested in.
#'@export
#'

calc_metric <- function(metric.name, metrics.vec) {
  any(c("ALL", metric.name) %in% metrics.vec)
}

