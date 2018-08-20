#==============================================================================
#Prepare Data for Analysis
#==============================================================================
#==============================================================================
#'Clean Taxa
#'
#'@param taxa.long = a data frame, in a long data format, containing taxonomic
#'hierarchy columns.
#'@return Standardizes taxa according to the Chessie BIBI protocol.
#'@export

clean_taxa <- function(taxa.long){
  # Taxonomic hierarchy column names.
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
                 "ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                 "TRIBE", "GENUS", "SPECIES")
  #============================================================================
  # If any NAs exist in the taxonomic columns, replace the NAs with
  # "UNIDENTIFIED."
  taxa.long[, taxa.cols] <- apply(taxa.long[, taxa.cols], 2, function(x){
    ifelse(is.na(x), "UNIDENTIFIED", as.character(x))
  } )
  #============================================================================
  # Keep only taxa from the phyla Annelida, Arthropoda, Mollosca,
  # and Platyhelminthes.
  phy.keep <- c("ANNELIDA", "ARTHROPODA", "MOLLUSCA", "PLATYHELMINTHES")
  taxa.long <- taxa.long[taxa.long$PHYLUM %in% phy.keep, ]
  #============================================================================
  # Keep only taxa from the subphyla Clitellata, Crustacea, Hexapoda, 
  # Rhabditophora, and taxa unidentified at this level.
  subphy.keep <- c("CLITELLATA", "CRUSTACEA", "HEXAPODA", "RHABDITOPHORA",
                   "UNIDENTIFIED")
  taxa.long <- taxa.long[taxa.long$SUBPHYLUM %in% subphy.keep, ]
  #============================================================================
  # Remove any taxa from the class Branchiopoda, Maxillopoda, and Ostracoda.
  class.exc <- c("BRANCHIOPODA", "MAXILLOPODA", "OSTRACODA")
  taxa.long <- taxa.long[!(taxa.long$CLASS %in% class.exc), ]
  #============================================================================
  # Remove any taxa from the order Hymenoptera. Aquatic Hymenoptera are 
  # generally small, parasitic organisms. Therefore, they may go easily 
  # unnoticed during sorting. We decided it was best to remove these organisms
  # from the analysis.
  order.exc <- c("HYMENOPTERA")
  taxa.long <- taxa.long[!(taxa.long$ORDER  %in% order.exc), ]
  #============================================================================
  # Remove any taxa from the families Gerridae, Hebridae, Veliidae, 
  # Hydrometridae, and Saldidae.  These taxa are generally considerd 
  # semi-aquatic because they live on the surface of the water. Therefore,
  # it is not appropriate to include these organisims in a benthic IBI.
  family.exc <- c("GERRIDAE", "HEBRIDAE", "VELIIDAE", "HYDROMETRIDAE",
                  "SALDIDAE")
  taxa.long <- taxa.long[!(taxa.long$FAMILY  %in% family.exc), ]
  #============================================================================
  # Remove any taxa from the genus Stenus. These taxa are considered 
  # semi-aquatic, and thus, were removed from the analysis.
  genus.exc <- c("STENUS")
  taxa.long <- taxa.long[!(taxa.long$GENUS  %in% genus.exc), ]
  #============================================================================
  # These taxa were not consitently identified to the same taxonomic rank
  # by the agencies that contributed data. Therefore, the samples were rolled
  # up to the lowest common denominator.
  # These columns will be influenced by the common denominator taxa.
  class.spp <- c("CLASS", "SUBCLASS", "ORDER", "SUBORDER",
                 "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")
  taxa.long[taxa.long$CLASS %in% "BIVALVIA", class.spp] <- "BIVALVIA"
  taxa.long[taxa.long$CLASS %in% "GASTROPODA", class.spp] <- "GASTROPODA"
  taxa.long[taxa.long$CLASS %in% "OLIGOCHAETA", class.spp] <- "OLIGOCHAETA"
  taxa.long[taxa.long$CLASS %in% "TREPAXONEMATA", class.spp] <- "TREPAXONEMATA"
  #============================================================================
  # These taxa were not consitently identified to the same taxonomic rank
  # by the agencies that contributed data. Therefore, the samples were rolled
  # up to the lowest common denominator.
  # These columns will be influenced by the common denominator taxa.
  order.spp <- c("ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                 "TRIBE", "GENUS", "SPECIES")
  taxa.long[taxa.long$ORDER %in% "COLLEMBOLA", order.spp] <- "COLLEMBOLA"
  taxa.long[taxa.long$ORDER %in% "LEPIDOPTERA", order.spp] <- "LEPIDOPTERA"
  taxa.long[taxa.long$ORDER %in% "NEUROPTERA", order.spp] <- "NEUROPTERA"
  taxa.long[taxa.long$ORDER %in% "NEOOPHORA", order.spp] <- "NEOOPHORA"
  #============================================================================
  # End clean_taxa function.
  return(taxa.long)
}

#==============================================================================

#'Fill empty rows with lowest level of
#'
#'@param taxa.long = A long data frame containing multiple levels of the 
#'taxonomic hierarchy.
#'@return Fills in blank spaces with the previous lowest level of taxonomic
#' identification.
#'@export
#requires package zoo
#applies lowest identification to each row

fill_taxa <- function(taxa.long){
  # Taxonomic hierarchy columns.
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS",
                 "SUBCLASS", "ORDER", "SUBORDER",
                 "FAMILY", "SUBFAMILY", "TRIBE",
                 "GENUS", "SPECIES")
  # Keep only the taxonomic columns that exist in the 
  taxa.cols <- taxa.cols[taxa.cols %in% names(taxa.long)]
  
  # If there are any blanks ("") or values set to "UNIDENTIFIED", then
  # replace these values with NA. This will allow these values to be
  # replaced with the previous taxonomic rank.
  sub.taxa.cols <- taxa.cols[!taxa.cols %in% "PHYLUM"]
  taxa.long[, sub.taxa.cols] <- sapply(sub.taxa.cols, function(x){
    sub.df <- taxa.long[, x]
    ifelse(sub.df %in% c("", "UNIDENTIFIED"), NA, sub.df)
  })
  #taxa.long[taxa.long %in% ""]  <- NA
  #taxa.long[taxa.long == "UNIDENTIFIED"]  <- NA
  # Create a new final data.frame and remove any rows that contain NAs in 
  # lowest taxonomic resolution specified. For instance if there is an 
  # NA in the Phylum column, the zoo function na.locf with break.
  warning(paste("At least one NA was remove from column", taxa.cols[1],
                "which is assumed to be the lowest taxonomic resolution in
                your data frame. Please review/update these rows and rerun.
                If the lowest taxonomic resolution column contains NA,
                the zoo package function na.locf breaks."))
  final.df <- taxa.long[!is.na(taxa.long[, taxa.cols[1]]), ]
  # Apply the zoo package function na.locf (Last Observation Carried Forward)
  # to fill in the NAs.
  final.df[, taxa.cols] <- t(apply(final.df[, taxa.cols], 1, zoo::na.locf))
  #final.df[, taxa.cols] <- data.frame(t(apply(final.df[, taxa.cols], 1, zoo::na.locf)))
  # Make sure all column names are uppercase.
  names(final.df) <- toupper(colnames(final.df))
  # Make sure all taxonomic names are uppercase.
  final.df[, taxa.cols] <- lapply(final.df[, taxa.cols], toupper)
  # End fill_taxa function.
  return(final.df)
}

#==============================================================================

event.test <- exists("EVENT")
EVENT2 <- if (event.test == "TRUE"){
  t<-get0("EVENT")
} else{
  0
}


#==============================================================================
#'Long Data Frame
#'Transform taxa count data from wide to a long format.
#'@param x = Taxa count data
#'@param Level = Taxonomic Level (PHYLUM, CLASS, ORDER, FAMILY, GENUS)
#'@return Taxa counts in a long data format.
#'@export
long <- function (Wide.df, Merge_on = "FINAL_ID") {

  agg <- aggregate(Long$REPORTING_VALUE ~ Long$EVENT_ID + Long$STATION_ID +
                     Long[, colnames(Long) == Level], data = Long, FUN = "sum",
                   na.rm = TRUE)
  colnames(agg) <- c("EVENT_ID", "STATION_ID", Level, "REPORTING_VALUE")
  wide <- reshape(agg, v.names = "REPORTING_VALUE", idvar = "EVENT_ID",
                  timevar = Level, direction = "wide")

  colnames.removing.prefix <- function(df, prefix) {
    names <- colnames(df)
    indices <- (substr(names, 1, nchar(prefix)) == prefix)
    names[indices] <- substr(names[indices], nchar(prefix) + 1, nchar(names[indices]))
    return(names)
  }
  colnames(wide) <- colnames.removing.prefix(wide, "REPORTING_VALUE.")
  names(wide)[names(wide) == ''] <- "UNIDENTIFIED" # number of taxa not identified to family
  wide[is.na(wide)] <- 0 #NA = zero
  wide <- wide[rowSums(wide[, 3:ncol(wide)]) != 0, ]
  wide <- with(wide,  wide[order(EVENT_ID), ])
  names(wide) <- toupper(colnames(wide))
  return(wide)
}


#==============================================================================
#'Wide Data Frame
#'Transform taxa count data from long to wide format
#'@param long.df = Taxa count data
#'@param taxa.rank = Taxonomic Level (PHYLUM, CLASS, ORDER, FAMILY, GENUS)
#'@param pct.unidentified = A threshold can be established to exclude samples 
#'that do have too many unidentified taxa in a sample. Enter the percentage of
#' the sample that you are comfortable with being unidentified.  
#' If not specified the function will not excluded any samples.
#'@return Taxa counts in wide format
#'@import data.table
#'@export

wide <- function (long.df, taxa.rank, pct.unidentified = NULL) {
  
  #============================================================================
  #print("[1/2] Aggregating data for transformation")
  # Convert the data.frame to data.table to speed up aggregation below.
  long <- data.table::as.data.table(long.df)
  # Sum the reporting value by unique sample and the specified taxonomic rank.
  agg.df <- long[, sum(REPORTING_VALUE), by = list(EVENT_ID, STATION_ID, DATE,
                                                   SAMPLE_NUMBER, AGENCY_CODE,
                                                   long[[taxa.rank]])]
  # Update the column names.
  colnames(agg.df) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                       "AGENCY_CODE", taxa.rank, "REPORTING_VALUE")
  #============================================================================
  #print("[2/2] Transforming from long.df data format to wide data format.")
  # Reformat from a long to a wide data format.
  wide.df <- tidyr::spread_(agg.df, taxa.rank, "REPORTING_VALUE")

  #Fill all NA's with zeros
  wide.df[is.na(wide.df)] <- 0
  # Sort the data.table according to the order of the specified columns.
  wide.df <- wide.df[order(EVENT_ID, STATION_ID, DATE, SAMPLE_NUMBER,
                           AGENCY_CODE), ]
  # Make sure all column names are uppercase.
  names(wide.df) <- toupper(colnames(wide.df))
  # If the pct.unidentified is not NULL then identify how many and which taxa 
  # are not identified at the specified taxonomic rank and elminate samples with
  # percentages greater than the specified pct.unidentified threshold.
  if(!is.null(pct.unidentified) & "UNIDENTIFIED" %in% names(wide.df)){
    cat("Samples with >=", pct.unidentified, "% taxa unidentified at the specified taxonomic level were excluded from the data set (N = ",
        nrow(wide.df) - sum((wide.df$UNIDENTIFIED / 
                             rowSums(wide.df[, 6:ncol(wide.df)]) * 100 >= pct.unidentified)),
        "). \n The number of samples with >= ", pct.unidentified, "% unidentified taxa: ",
        sum((wide.df$UNIDENTIFIED / rowSums(wide.df[, 6:ncol(wide.df)]) *
             100 >= pct.unidentified)), " (N = ", nrow(wide.df), "; ",
        round((sum((wide.df$UNIDENTIFIED / rowSums(wide.df[, 6:ncol(wide.df)]) * 
               100 >= pct.unidentified)) / nrow(wide.df)) * 100, 2), "%)", sep ="")

    wide.df <- wide.df[!((wide.df$UNIDENTIFIED /
                            rowSums(wide.df[, 6:ncol(wide.df)])) * 100 >= pct.unidentified), ]
  }
  # Convert data.table to data.frame.
  final.df <- data.frame(wide.df)
  # End wide function.
  return(final.df)
}

#==============================================================================
