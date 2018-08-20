#==============================================================================
# Composition Metrics
#==============================================================================
#'Proportion of Gastropoda, Oligochaeta, and Dipteran Individuals
#'
#'@param class.wide = Taxonomic counts aggregated at the class level
#'classification in a wide data format. Use the wide function to
#'prepare the data.
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return One minus the proportional number of gastropods (Class: Gastropoda),
#' oligochaetes (Class: Oligochaeta), and dipteran (Order: Diptera) individuals.
#' This metric typically decreases with degradation.
#'@export

gold <- function(class.wide, order.wide) {
  god <- blank_col("GASTROPODA", class.wide) +
         blank_col("OLIGOCHAETA", class.wide) +
         blank_col("DIPTERA", order.wide)
  return((1 - (god / rowSums(order.wide[, 6:ncol(order.wide)]))) * 100)
}

#==============================================================================
#'Percentage of Amphipod Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as amphipods
#' (Order: Amphipoda).  This metric will typically increase with degradation.
#'@export

pct_amphipoda <- function(order.wide) {
  return(blank_col("AMPHIPODA", order.wide) /
         rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of Chironomid Individuals
#'
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as chironomids
#'(Family: Chironomidae).
#'This metric will typically increase with degradation.
#'@export

pct_chironomidae <- function(family.wide) {
  return(blank_col("CHIRONOMIDAE", family.wide) /
         rowSums(family.wide[, 6:ncol(family.wide)]) * 100)
}

#==============================================================================
#'Percentage of Corbiculid Individuals
#'
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as corbiculids
#'(Family: Corbiculidae).
#'@export

pct_corbiculidae <- function(family.wide) {
  return(blank_col("CORBICULIDAE", family.wide) /
         rowSums(family.wide[, 6:ncol(family.wide)]) * 100)
}

#==============================================================================
#'Percentage of Bivalve Individuals
#'
#'@param class.wide = Taxonomic counts aggregated at the class level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as bivalves
#'(Class: Bivalvia).
#'@export

pct_bivalvia <- function(class.wide) {
  return(blank_col("BIVALVIA", class.wide) /
           rowSums(class.wide[, 6:ncol(class.wide)]) * 100)
}

#==============================================================================
#'Percentage of Unionoid Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as unionoids
#' (Order: Unionoida).
#'@export

pct_unionoida <- function(order.wide) {
  return(blank_col("UNIONOIDA", order.wide) /
           rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of Dipteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as dipterans
#' (Order: Diptera). This metric will typically increase with degradation.
#'@export

pct_diptera <- function(order.wide) {
  return(blank_col("DIPTERA", order.wide) /
         rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of Coleopteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as coleopterans
#'(Order: Coleoptera).
#'@export

pct_coleoptera <- function(order.wide) {
  return(blank_col("COLEOPTERA", order.wide) /
         rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of Odonate Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as odonates (Order: Odonata).
#'@export

pct_odonata <- function(order.wide) {
  return(blank_col("ODONATA", order.wide) /
         rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}


#==============================================================================
#'Percentage of Ephemeropteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as Ephemeropterans
#' (Order: Ephemeroptera).  This metric typically decreases with degradation.
#'@export
pct_ephemeroptera <- function(order.wide) {
  return(blank_col("EPHEMEROPTERA", order.wide) /
           rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of Ephemeropteran Individuals Minus Baetid
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as Ephemeropterans
#' (Order: Ephemeroptera) minus the percentage of baetid (Family: Baetidae)
#'  in the sample.  This metric typically decreases with degradation.
#'@export

pct_epmeroptera_no_baetid <- function(order.wide, family.wide) {
  baetid <- (blank_col("BAETIDAE", family.wide))
  ephem <- (blank_col("EPHEMEROPTERA", order.wide))

  final.vec <- ((ephem - baetid) / 
                  rowSums(order.wide[, 6:ncol(order.wide)])) * 100
  
  return(final.vec)
}

#==============================================================================
#'Percentage of EPT Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as ephemeropterans
#' (Order: Ephemeroptera), plecopterans (Order: Plecoptera),
#' and trichopterans (Order: Trichoptera) (EPT).  This metric will typically
#' decrease with degradation.
#'@export

pct_ept <- function(order.wide) {
  EPT <- blank_col("EPHEMEROPTERA", order.wide) +
         blank_col("PLECOPTERA", order.wide) +
         blank_col("TRICHOPTERA", order.wide)
  return(EPT / rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}


#==============================================================================
#'Percentage of EPT Individuals Minus Hydropsychid and Baetid Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as ephemeropterans
#' (Order: Ephemeroptera), plecopterans (Order: Plecoptera),
#' and trichopterans (Order: Trichoptera) (EPT) minus the percentage of
#' hydropsychids (Family: Hydropsychidae) and baetids (Family: Baetidae).
#' These two families are typically tolerant and can be hyperdominant,
#' resulting in elevated percentages of EPT that may not accurately
#' describe the degradation gradient of interest. This metric will
#' typically decrease with degradation.
#'@export

pct_ept_hydro_baetid <- function(order.wide, family.wide) {
  hydropsychidae <- blank_col("HYDROPSYCHIDAE", family.wide)
  baetidae <- blank_col("BAETIDAE", family.wide)
  ept <- (blank_col("EPHEMEROPTERA", order.wide) +
            blank_col("PLECOPTERA", order.wide) +
            blank_col("TRICHOPTERA", order.wide))
  final.vec <- ifelse(ept == 0, 0, (1 - ((hydropsychidae + baetidae) / ept)) * 100)
  return(final.vec)
}

#==============================================================================
#'Percentage of EPT Individuals Minus Hydropsychid Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as ephemeropterans
#' (Order: Ephemeroptera), plecopterans (Order: Plecoptera),
#' and trichopterans (Order: Trichoptera) (EPT) minus the percentage of
#' hydropsychids (Family: Hydropsychidae).
#' This family is typically tolerant and can be hyperdominant,
#' resulting in elevated percentages of EPT that may not accurately
#' describe the degradation gradient of interest. This metric will
#' typically decrease with degradation.
#'@export

pct_ept_no_hydro <- function(order.wide, family.wide) {
  hydropsychidae <- blank_col("HYDROPSYCHIDAE", family.wide)
  ept <- (blank_col("EPHEMEROPTERA", order.wide) +
            blank_col("PLECOPTERA", order.wide) +
            blank_col("TRICHOPTERA", order.wide))
  final.vec <- ifelse(ept == 0, 0, (1 - (hydropsychidae / ept)) * 100)
  return(final.vec)
}

#==============================================================================
#'Percentage EPT Individuals Identified as Hydropsychids
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return  The percentage of EPT (Orders: Ephemeroptera, Plecoptera,
#' and Trichoptera) individuals identified as hydropsychids
#' (Family: Hydropsychidae).  This metric typically increases with degradation.
#'@export

pct_hydro_ept <- function(order.wide, family.wide) {
  hydro <- (blank_col("HYDROPSYCHIDAE", family.wide))
  ept <- (blank_col("EPHEMEROPTERA", order.wide) +
            blank_col("PLECOPTERA", order.wide) +
            blank_col("TRICHOPTERA", order.wide))
  final.vec <- ifelse(ept == 0, 0, (hydro / ept) * 100)

  return(final.vec)
}

#==============================================================================
#'Percentage of Trichopterans Represented by Hydropsychids
#'
#'@param order.wide = Taxonomic counts aggregated at the class level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of trichopteran (Order:Trichoptera) individuals
#'that are hydropsychids (Family: Hydropsychidae). This metric typically
#'increases with degradation.
#'@export

pct_hydro_trichoptera <- function(order.wide, family.wide) {
  return(ifelse(order.wide[, "TRICHOPTERA"] == 0, 0,
                (blank_col("HYDROPSYCHIDAE", family.wide) /
                   blank_col("TRICHOPTERA", order.wide))  * 100))
}

#==============================================================================
#'Percentage of COTE Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as coleopterans
#' (Order: Coleoptera), odonates (Order: Odonata), trichopterans
#'  (Order: Trichoptera), and ephemeropterans (Order: Ephemeroptera) (COTE).
#'  These organisms were used by Kamman (2007) as a lentic benthic
#'  macroinverterbate metric.
#'
#'  Citation:
#'  KAMMAN, N. 2007. Development of Biocriteria for Vermont and
#'   New Hampshire Lakes Criteria Development for Macroinvertebrates
#'   for Three Lake Classes and Implementation Procedure for Biological
#'   Assessment of Vermont Lakes. Pages 1–21. Vermont Department
#'   of Environemental Conservation, Waterbury, VT.
#'@export

pct_cote <- function(order.wide) {
  COTE <- blank_col("COLEOPTERA", order.wide) +
          blank_col("ODONATA", order.wide) +
          blank_col("TRICHOPTERA", order.wide) +
          blank_col("EPHEMEROPTERA", order.wide)
  return(COTE / rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of POTEC Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as plecopterans
#' (Order: Plecoptera), odonates (Order: Odonata), trichopterans
#'  (Order: Trichoptera), ephemeropterans (Order: Ephemeroptera),
#'  and coleopterans (Order: Coleoptera) (POTEC).  This is an
#'  exploratory metric combining EPT taxa with COTE taxa.
#'  Coleopterans (Order: Coleoptera), odonates (Order: Odonata),
#'  trichopterans (Order: Trichoptera), and ephemeropterans
#'  (Order: Ephemeroptera) (COTE) were used by Kamman (2007) as
#'   a lentic benthic macroinverterbate metric.
#'
#'  Citation:
#'  KAMMAN, N. 2007. Development of Biocriteria for Vermont and
#'   New Hampshire Lakes Criteria Development for Macroinvertebrates
#'   for Three Lake Classes and Implementation Procedure for Biological
#'   Assessment of Vermont Lakes. Pages 1–21. Vermont Department
#'   of Environemental Conservation, Waterbury, VT.
#'@export

pct_potec <- function(order.wide) {
  POTEC <- blank_col("PLECOPTERA", order.wide) +
           blank_col("ODONATA", order.wide) +
           blank_col("TRICHOPTERA", order.wide) +
           blank_col("EPHEMEROPTERA", order.wide) +
           blank_col("COLEOPTERA", order.wide)
  return(POTEC / rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of Limestone Individuals
#'
#'@param order.wide = order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as isopods (Order: Isopoda),
#'amphipods (Order: Amphipoda), and ephemerellids (Family: Ephemerellidae).
#'@export

pct_limestone <- function(order.wide, family.wide) {
  IAE <- blank_col("ISOPODA", order.wide) +
         blank_col("AMPHIPODA", order.wide) +
         blank_col("EPHEMERELLIDAE", family.wide)
  return(IAE / rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of Retreat-Making Trichopteran Individuals
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@return The percentage of individuals identified as retreat-making
#' trichopterans (Suborder: Annulipalpia).
#'@export
pct_retreat_trichoptera <- function(long) {
  sub.ord <- wide(long, "SUBORDER")
  return(blank_col("ANNULIPALPIA", sub.ord) /
        rowSums(sub.ord[, 6:ncol(sub.ord)]) * 100)
}

#==============================================================================
#'Percentage of Hydropsychid Individuals
#'
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as hydropsychids
#' (Family: Hydropsychidae).
#'@export

pct_hydropsychidae <- function(family.wide) {
  return(blank_col("HYDROPSYCHIDAE", family.wide) /
           rowSums(family.wide[, 6:ncol(family.wide)]) * 100)
}

#==============================================================================
#'Percentage of Non-Insect Individuals
#'
#'@param class.wide = Taxonomic counts aggregated at the class level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals that were not identified as insects
#'(i.e., not identified as class Insecta).
#'@export

pct_non_insect <- function(class.wide){
  Non_Insect <- rowSums(class.wide[, 6:ncol(class.wide)]) -
    blank_col("INSECTA", class.wide)
  return(Non_Insect / rowSums(class.wide[, 6:ncol(class.wide)]) * 100)
}

#==============================================================================
#'Percentage of Oligochaet Individuals
#'
#'@param class.wide = Taxonomic counts aggregated at the class level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as oligochaets
#'(Class: Oligochaeta).
#'@export

pct_oligochaeta <- function(class.wide){
  return(blank_col("OLIGOCHAETA", class.wide) /
         rowSums(class.wide[, 6:ncol(class.wide)]) * 100)
}

#==============================================================================
#'Percentage of Plecopteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as plecopterans
#'(Order: Plecoptera).
#'@export

pct_plecoptera <- function(order.wide) {
  return(blank_col("PLECOPTERA", order.wide) /
         rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of Trichopteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as trichopterans
#'(Order: Trichoptera).
#'@export

pct_trichoptera <- function(order.wide) {
  return(blank_col("TRICHOPTERA", order.wide) /
         rowSums(order.wide[, 6:ncol(order.wide)]) * 100)
}

#==============================================================================
#'Percentage of Non-Hydropsychid Trichopteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as trichopterans
#'(Order:Trichoptera) excluding hydropsychids (Family: Hydropsychidae).
#'@export

pct_non_hydrop_trichoptera <- function(order.wide, family.wide) {
  trichop.all <- blank_col("TRICHOPTERA", order.wide)
  trichop.no.hydro <-  trichop.all - blank_col("HYDROPSYCHIDAE", family.wide)
  final.vec <- ifelse(trichop.all == 0, 0, (trichop.no.hydro / trichop.all) * 100)
  return(final.vec)
}

#==============================================================================
#'Percentage of Simuliid Individuals
#'
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as simuliids
#'(Family: Simuliidae).
#'@export

pct_simuliidae <- function(family.wide) {
  return(blank_col("SIMULIIDAE", family.wide) /
         rowSums(family.wide[, 6:ncol(family.wide)]) * 100)
}

#==============================================================================
#'Percentage of Oligochaet and Chironomid Individuals
#'
#'@param class.wide = Taxonomic counts aggregated at the class level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as oligochaets
#' (Class: Oligochaeta) and chironomids (Family: Chironomidae).
#'@export

pct_oligo_chiro <- function(class.wide, family.wide){
  chiro <- blank_col("CHIRONOMIDAE", family.wide)
  oligo <- blank_col("OLIGOCHAETA", class.wide)
  final.vec <- ((chiro + oligo) / rowSums(family.wide[, 6:ncol(family.wide)])) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Annelid and Chironomid Individuals
#'
#'@param phylum.wide = Taxonomic counts aggregated at the phylum level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as annelids (Phylum: Annelida)
#'and chironomids (Family: Chironomidae).  This metric typically increases with
#'degradation.
#'@export

pct_annelid_chiro <- function(phylum.wide, family.wide){
  chiro <- blank_col("CHIRONOMIDAE", family.wide)
  annelid <- blank_col("ANNELIDA", phylum.wide)
  final.vec <- ((chiro + annelid) /
                rowSums(family.wide[, 6:ncol(family.wide)])) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Unidentified Individuals
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@return The percentage of individuals that were not identified at the
#'specified taxonomic level.
#'@export

pct_unidentified <- function(taxa.wide) {
  return(blank_col("UNIDENTIFIED", taxa.wide) /
           rowSums(taxa.wide[, 6:ncol(taxa.wide)]) * 100)
}

#==============================================================================
#'Percentage of EPT Taxa Excluding Tolerant Taxa
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. You must
#'sepecify either 'FAMILY' or "GENUS.'
#'@param master = A master taxa list including taxonomic ranks Phylum through
#'the specified taxonomic rank (Family or Genus) and the an
#'associated list of tolerance values. The default is set to the master taxa
#'list included in the BIBI package.  The master taxa list can be viewed with
#'the following script: master.df <- BIBI::master
#'@param tolerance_value = The name of the column in the master taxon list
#'(specified using the master variable) that contains tolerance values on
#'a scale of 0-10.  Tolerant organisms are classified as organisms with a
#'tolerance value >= 7.  The defualt is set to the the BIBI tolerance values,
#'which are tolerance values summarized from multiple sources.
#'@return Percent of the assembalge represented by Trichoptera individuals,
#' excluding tolerant Trichoptera taxa.
#'@export

pct_trichoptera_no_tolerant <- function (long, rank, master = BIBI:master, tolerance_value = "BIBI_TV") {
  
  wide.df <- wide(long, rank)
  Order <- split(long[, rank], long$ORDER)
  trichop <- unique(Order$TRICHOPTERA)
  taxa.list <- trichop
  new.df <- wide.df[, names(wide.df) %in% taxa.list]
  new.df <- new.df[, !(names(new.df) %in% "UNIDENTIFIED")]
  
  # Find all of the tolerant taxa
  master$TOLERANCE <- ifelse(master[, tolerance_value] >= 7, "TOLERANT", NA)
  tol <- split(master[, rank], master$TOLERANCE)
  name.list <- as.character(unlist(unique(tol$TOLERANT)))
  # Remove any tolerant taxa from the trichop data frame
  no.tol.trichop <- data.frame(new.df[, !(names(new.df) %in% name.list)])
  
  # Caculate richness values using vegan
  sum.no.tol.trichop <- apply(no.tol.trichop[, 6:ncol(no.tol.trichop)], 1, sum)
  total.sum <- apply(new.df[, 6:ncol(new.df)], 1, sum)
  #sum.no.tol.trichop <- rowsum(no.tol.trichop[, 6:ncol(no.tol.trichop)])
  #total.sum <- rowsum(new.df[, 6:ncol(new.df)])
  final.vec <- ifelse(sum.no.tol.trichop == 0, 0, (sum.no.tol.trichop / total.sum) * 100)
  return(final.vec)
}
