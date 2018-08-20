#==============================================================================
#Rarefaction prepration and related metrics

#'Rarify Taxonomic Counts
#'
#'@param Long = Taxa count data
#'@param Level = Taxonomic level to preform rarefaction on.
#'@param Sample = The rarefaction subsample size.
#'@param pct_un = A threshold can be established to exclude samples that do
#'have too many unidentified taxa in a sample. Enter the percentage of the
#'sample that you are comfortable with being unidentified.  If not specified
#'the function will not excluded any samples.
#'@param bibi.standard = Indicate TRUE if you would like to standardize your
#'taxa following the practices used in the Chessie BIBI.  The default, FALSE,
#'will not attempt to standardize the taxa.
#'@return Rarified counts
#'@export

prep_rare <- function(Long, master, Level = "FAMILY", Sample = 100, pct_un = NULL, bibi.standard = FALSE){
  #Transform the data to a wide format using the TSN as column names.
  TSN <- wide(Long, Level)

  rarefied.df <- if(min(rowSums(TSN[, 6:ncol(TSN)])) > Sample){
    # If all of the stations have more than 100 idividuals identified then
    # the vegan function rarefy is used to find richness if only
    # 100 individuals observed and the total number of individuals observed
    # is reported.
    large_comm <- function(tsn.df){
      #New data frame created with unique event ID's.
      test <- data.frame(unique(tsn.df[, c("EVENT_ID", "STATION_ID",
                                           "DATE","SAMPLE_NUMBER",
                                           "AGENCY_CODE")]))
      #The vegan function rarefy is used to find the richness
      #for the event subset to x (default = 100) number of individuals.
      #test$RICH <- round(vegan::rarefy(tsn.df[, 6:ncol(tsn.df)], Sample), 0)
      test$RICH <- vegan::rarefy(tsn.df[, 6:ncol(tsn.df)], Sample)
      #The total number of individuals identified per event.
      test$TOTAL <- rowSums(tsn.df[, 6:ncol(tsn.df)])
      return(test)
    }
    #Run Function
    final <-(large_comm(TSN))

  }else{
    #If at least one station has less than 100 individuals
    # identified then the data frame is divided into two data frames.
    # One data frame contains all of the stations with <= 100
    # individuals, while the other data frame contatins all of the
    # stations with > 100 individuals.
    small_comm <- function(tsn.df){
      #A new data frame is created for events with >= x (default = 100)
      # individuals identified during an event.
      TSN_101 <- tsn.df[(rowSums(tsn.df[, 6:ncol(tsn.df)]) > Sample), ]
      #rownames(TSN_101) <- TSN_101$EVENT_ID #NECESSARY???

      #A new data frame is created for events with < x (default = 100)
      # individuals identified during an event.
      TSN_0 <- tsn.df[(rowSums(tsn.df[, 6:ncol(tsn.df)]) <= Sample),]
      #rownames(TSN_101) <- TSN_101$EVENT_ID#NECESSARY???
      TS <- data.frame(unique(TSN_0[, c("EVENT_ID", "STATION_ID",
                                        "DATE","SAMPLE_NUMBER",
                                        "AGENCY_CODE")]))
      #The vegan function, specnumber, is used to find the richness in the
      # sample with <= x (defualt = 100) number of individuals.
      TS$RICH <- vegan::specnumber(TSN_0[, 6:ncol(TSN_0)])
      #The total number of individuals identified per event.
      TS$TOTAL <- rowSums(TSN_0[, 6:ncol(TSN_0)])

      #A new data frame is created for events with > x (default = 100)
      # individuals identified during an event.
      new <- data.frame(unique(TSN_101[, c("EVENT_ID", "STATION_ID",
                                           "DATE","SAMPLE_NUMBER",
                                           "AGENCY_CODE")]))
      #The vegan function rarefy is used to find the richness
      #for the event subset to x (default = 100) number of individuals.
      #new$RICH <-round(vegan::rarefy(TSN_101[, 6:ncol(TSN_101)], Sample), 0)
      new$RICH <- vegan::rarefy(TSN_101[, 6:ncol(TSN_101)], Sample)
      #The total number of individuals identified per event.
      new$TOTAL <- rowSums(TSN_101[, 6:ncol(TSN_101)])
      #Combine the data frame with <= x (default = 100) identified individuals
      # and the data frame with > x (default =100) identified individuals.
      return(rbind(new, TS))
    }
    #Run function
    final <- (small_comm(TSN))
  }
  #Reduce the number of columns in the original long data frame.
  new.df <- Long[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                     "AGENCY_CODE", Level, "REPORTING_VALUE")]
  new.df[, Level][is.na(new.df[, Level])] <- "UNIDENTIFIED"
  #Sum all counts at the specified taxonomic level (default = "FAMILY")
  # by unique event ID and station ID.
  #DT <- data.table::data.table(new.df)
  #names(DT) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
  #               "AGENCY_CODE", "FINAL_ID", "REPORTING_VALUE")
  #agg.rar <- DT[ , list(sum(REPORTING_VALUE)), by = list(EVENT_ID, STATION_ID,
  #                                                 DATE, SAMPLE_NUMBER,
  #                                                 AGENCY_CODE, FINAL_ID)]
  #names(agg.rar)[names(agg.rar) %in% "V1"] <- "REPORTING_VALUE"
  #names(agg.rar)[names(agg.rar) %in% "FINAL_ID"] <- Level
  agg.rar <- aggregate(REPORTING_VALUE ~ EVENT_ID + STATION_ID + DATE +
                         SAMPLE_NUMBER + AGENCY_CODE + new.df[, Level],
                      data = new.df, FUN = "sum")
  names(agg.rar)[names(agg.rar) %in% "new.df[, Level]"] <- Level
  #Order the taxa in descending order for each unique event ID + station ID.
  ord.rar <- agg.rar[order(agg.rar$EVENT_ID, agg.rar$STATION_ID, agg.rar$SAMPLE_NUMBER,
                           agg.rar$DATE, agg.rar$AGENCY_CODE,
                           - agg.rar$REPORTING_VALUE), ]
  #Merge the data frame organized in descending order with the corresponding
  # rarefied data frame for each event ID.
  merg.rar <- merge(ord.rar, rarefied.df, by = c("EVENT_ID", "STATION_ID",
                                                 "DATE", "SAMPLE_NUMBER",
                                                 "AGENCY_CODE"))
  merg.rar$RICH <- round(merg.rar$RICH, 0)

  #Create a new data frame containing the stations with <= 100 identified
  # individuals.
  less_x.rar <- data.frame(merg.rar[merg.rar$TOTAL <= Sample, ])
  colnames(less_x.rar) <- c("EVENT_ID", "STATION_ID",
                            "DATE", "SAMPLE_NUMBER",
                            "AGENCY_CODE", Level, "RV",
                            "RICH", "TOTAL")
  less_x.rar <- less_x.rar[, c("EVENT_ID", "STATION_ID",
                               "DATE","SAMPLE_NUMBER",
                               "AGENCY_CODE", Level, "RV",
                               "RICH", "TOTAL")]
  #remove the columns RICH and TOTAL.
  final_less_x.rar <- less_x.rar[, !(names(less_x.rar) %in% c("RICH", "TOTAL"))]

  #Create a new data frame containing the stations with > 100 identified
  # individuals.
  more_x.rar <- merg.rar[merg.rar$TOTAL > Sample, ]
 
  #This function samples the taxa with the highest probability of being
  # selected when only x (default = 100) number of individuals are subsampled.
  # The function creates a more stable rarefied community. Random sampling
  # without replacement is still necessary when last taxon selected in the
  # rareried sample has the same count as other taxa.
  rar_func <- function(x){
    #A new data frame using the richness column to subset each event.
    df.x <- head(x, n = round(mean(x$RICH)), 0)
    #Order the taxa in descending order for each unique event ID + station ID.
    df.x<- df.x[order(-df.x$REPORTING_VALUE), ]
    #The number of taxa with the same count value that will be included
    # in the final rarefied assemblage. Equation: # of taxa (richness)
    # with count values greater than the lowest count value in the
    # rarefied assemblage minus the richness value.
    num.extract <- as.numeric(as.character(abs(sum(df.x$REPORTING_VALUE >
                                        tail(df.x$REPORTING_VALUE, n = 1)) -
                                    round(mean(x$RICH)))))
    #Extract all of the taxa with the same count value as the last taxon
    # in the rarefied assemblage.
    same.num <- x[(tail(df.x$REPORTING_VALUE, n = 1) ==
                     x$REPORTING_VALUE), ]
    #Subsample (size is equal to the number needed to complete the
    # rarefied assemblage) the taxa with the same count value without
    # replacement.
    rand.num <- same.num[sample(nrow(same.num), size = num.extract,
                                replace = FALSE), ]
    #Remove all of the taxa with the same count value as the final
    # taxon in the rarefied assemblage.
    new.test <- df.x[df.x$REPORTING_VALUE !=
                       tail(df.x$REPORTING_VALUE, n = 1), ]
    #Bind the data frame the new data frame above with the random
    # subsample of taxa with the same count value.
    return(rbind(new.test, rand.num))

  }
  
  #A new data frame is created using the plyr package function "ddply"
  # to run the rar_func function by each unique event ID.
  new_rare.df <- plyr::ddply(more_x.rar, plyr::.(EVENT_ID, STATION_ID,
                                                 DATE, SAMPLE_NUMBER,
                                                 AGENCY_CODE), rar_func)

  #A new total count value is calculated for the rarefied assemblage.
  total_rare.df <- plyr::ddply(new_rare.df, plyr::.(EVENT_ID, STATION_ID,
                                                    DATE, SAMPLE_NUMBER,
                                                    AGENCY_CODE),
                               function(x) sum(x$REPORTING_VALUE))
  
  #dt.tr <- data.table::data.table(new_rare.df)
  #total_rare.df <- dt.tr[, sum(REPORTING_VALUE), by = EVENT_ID]
  colnames(total_rare.df) <- c("EVENT_ID", "STATION_ID",
                               "DATE", "SAMPLE_NUMBER",
                               "AGENCY_CODE", "SUM")
  #Merge the rarefied community with the new total count values.
  merg_rare.df <- merge(new_rare.df, total_rare.df, by = c("EVENT_ID", "STATION_ID",
                                                           "DATE", "SAMPLE_NUMBER",
                                                           "AGENCY_CODE"))
  
  #Create a new column of the proportion of the sample represented by each
  # taxon (default is out of 100).
  merg_rare.df$RV <- (merg_rare.df$REPORTING_VALUE / merg_rare.df$SUM) * Sample
  #Remove the unecessary columns
  merg_rare.df <- merg_rare.df[, -which(names(merg_rare.df) %in%
                                          c("TOTAL","RICH", "REPORTING_VALUE", "SUM"))]
  colnames(merg_rare.df) <- c("EVENT_ID","STATION_ID", "DATE", "SAMPLE_NUMBER",
                               "AGENCY_CODE", Level, "RV")

  #Combine the rarefied assemblage data frame with the data frame containing
  # events with fewer than x (default = 100) number of individuals identified.
  bound.df <- rbind(merg_rare.df, final_less_x.rar)
  names(bound.df)[names(bound.df) %in% "RV"] <- "REPORTING_VALUE"

  #The reshape2 package function melt and dcast enable 'Level' to be used
  # to specify a unique column. 'Level' provides flexability to the sricpt
  # allowing the user to specify which taxonomic level should be used
  # during the analysis. Could not use the tdyr package function spread because
  # the objects are called using unquoted names.

  #The reshape2 package function melt preps the data frame to be transformed
  # to a wide format
  melted_rare.df <- reshape2::melt(bound.df, id.vars = c("EVENT_ID", "STATION_ID", "DATE",
                                                         "SAMPLE_NUMBER", "AGENCY_CODE", Level))
  # The reshape2 package function dcast transforms the data frame to a wide
  # format.
  final_rare.df<- reshape2::dcast(melted_rare.df, EVENT_ID + STATION_ID + DATE +
                                    SAMPLE_NUMBER + AGENCY_CODE ~ melted_rare.df[, Level])
  #Fill in all NA's with zero.
  final_rare.df[is.na(final_rare.df)] <- 0
  #Round all numbers to the nearest integer
  final_rare.df[, 6:ncol(final_rare.df)] <- ceiling(final_rare.df[, 6:ncol(final_rare.df)])
 
#  if(!is.null(pct_un) & "UNIDENTIFIED" %in% names(final_rare.df)){
#    cat("Samples with >=", pct_un, "% taxa unidentified at the specified taxonomic level were excluded from the data set (N = ",
#        nrow(final_rare.df) - sum((final_rare.df$UNIDENTIFIED / 
#                               rowSums(final_rare.df[, 6:ncol(final_rare.df)]) * 100 >= pct_un)),
#        "). \n The number of samples with >= ", pct_un, "% unidentified taxa:",
#        sum((final_rare.df$UNIDENTIFIED / rowSums(final_rare.df[, 6:ncol(final_rare.df)]) *
#               100 >= pct_un)), "(N = ", nrow(final_rare.df), "; ",
#        round((sum((final_rare.df$UNIDENTIFIED / rowSums(final_rare.df[, 6:ncol(final_rare.df)]) * 
#                      100 >= pct_un)) / nrow(final_rare.df)) * 100, 2), "%)", sep ="")
    
#    final_rare.df <- final_rare.df[!((final_rare.df$UNIDENTIFIED /
#                            rowSums(final_rare.df[, 6:ncol(final_rare.df)])) * 100 >= pct_un), ]
#  }

 gather.df <- tidyr::gather(final_rare.df, FINAL_ID, REPORTING_VALUE, 6:ncol(final_rare.df))
 final.df <- gather.df[gather.df$REPORTING_VALUE > 0, ]
 
 master.u <- unique(master[, !(names(master) %in% c("TAXON_LEVEL", "TSN_FINAL", "TSN_DB",
                                                    "TSN_R",
                                             "TSN_AGENCY_ID", "AGENCY_ID",
                                              "SYN", "VALIDATION", "NOTES"))])
 master.u$FINAL_ID <- as.character(master.u$FINAL_ID)
 #test3 <- master[master$FINAL_ID %in% "OLIGOCHAETA",]

 final.df <- dplyr::left_join(final.df, master.u, by = "FINAL_ID")
 
 if(bibi.standard == TRUE){
   final.df <- BIBI::clean_taxa(final.df)
 }
  
  return(final.df)
}

#==============================================================================
#'rrarefied counts in a long data frame
#'
#'@param Long = Taxa count data
#'@param Level = The taxonomic level to preform the analysis.
#'@return The vegan package function rrarefy is used to subsample the
#' orginal assemblage. Level is used to specify the lowest taxonomic resolution.
#' Rarefaction is performed at the rank specified by Level. The new wide data
#'  frame is merged with the orginal Long data frame to represent that
#'  rarified counts and all of the applicable taxonomic ranks.
#'@export


prep_rrarefy <- function(Long.df, Level){
  wide.df <- wide(Long.df, Level)
  rare.df <- data.frame(vegan::rrarefy(wide.df[, 6:ncol(wide.df)], sample = 100))
  bound <- cbind(wide.df[, 1:5], rare.df)
  long.rare <- reshape2::melt(bound, id.vars = c("EVENT_ID", "STATION_ID",
                                                 "DATE", "AGENCY_CODE",
                                                 "SAMPLE_NUMBER"))
  long.rare <- long.rare[long.rare$value > 0, ]
  names(long.rare) <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE",
                        "SAMPLE_NUMBER", Level, "REPORTING_VALUE")
  master.taxa <- unique(master[, c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
                                   "ORDER",
                                   if(Level %in% c("FAMILY", "GENUS", "SPECIES")) "SUBORDER",
                                   if(Level %in% c("FAMILY", "GENUS", "SPECIES")) "FAMILY",
                                   if(Level %in% c("GENUS", "SPECIES")) "SUBFAMILY",
                                   if(Level %in% c("GENUS", "SPECIES")) "TRIBE",
                                   if(Level %in% c("GENUS", "SPECIES")) "GENUS",
                                   if(Level %in% c("GENUS", "SPECIES"))"SPECIES")])
  final.df <- merge(long.rare, master.taxa, by = Level)
  return(final.df)
}




#==============================================================================
#'Most Probable Rarefied Assemblage
#'
#'@param Long = Taxa count data
#'@param Sample = The rarefied sample size
#'@return Selects the assemblage with the highest probability of being selected
#'during a random rarefaction.
#'@export

prarefy <- function(Long, Sample = 100){
  ID <- c("EVENT_ID", "STATION_ID", "DATE",
          "AGENCY_CODE", "SAMPLE_NUMBER")

  agg.df <- data.frame(aggregate(REPORTING_VALUE ~ EVENT_ID + STATION_ID +
                                   DATE + AGENCY_CODE + SAMPLE_NUMBER,
                                 data = Long, FUN = sum))
  names(agg.df) <- c(ID, "CSUM")

  frarefy <- function(Long, agg.df, Sample){
    merged <- merge(Long, agg.df, by = ID)
    merged$REPORTING_VALUE <- round((merged$REPORTING_VALUE /
                                       merged$CSUM) * Sample, digits = 0)
    merged <- merged[!(merged$REPORTING_VALUE == 0),]
    merged <- merged[,!(names(merged) %in% "CSUM")]
    return(merged)
  }

  if(min(agg.df$CSUM) <= Sample){
    small.list <- agg.df[agg.df$CSUM <= Sample, ID]
    small.df <- merge(small.list, Long, by = ID)
    large.list <- agg.df[agg.df$CSUM > Sample, ID]
    large.df <- merge(large.list, Long, by = ID)
    rare.large <- frarefy(large.df, agg.df, Sample)
    final.df <- rbind(rare.large, small.df)
  }

  if(min(agg.df$CSUM) > Sample){
    final.df <- frarefy(Long, agg.df, Sample)
  }


  return(final.df)
}
