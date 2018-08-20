#==============================================================================
# Discrimination Efficiency
#==============================================================================
#'Optimal Discrimination Efficiency
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@param quant.df = Data frame containing upper.class quantile values.
#'@param upper.class = the site class that represents the better condition.
#'@param lower.class = the site class that represents the poorer condition.
#'@param ref.df = a data frame of only the upper.class values.
#'@param quant.ref = a data frame of reference quantile values.
#'@return Determines the threshold at which a metric best categorizes
#'reference and degraded stations.
#'@export
#'

ode <- function(metrics.df, quant.df, upper.class, lower.class, ref.df, quant.ref){
  #Transform the metrics data frame from a wide format to a long data format.
  melted <- reshape2::melt(metrics.df, id.vars = c("EVENT_ID", "STATION_ID", "CATEGORY",
                                                   "DATE", "SAMPLE_NUMBER", "AGENCY_CODE"))
  colnames(melted) <- c("EVENT_ID", "STATION_ID", "CATEGORY", "DATE",
                        "SAMPLE_NUMBER", "AGENCY_CODE",  "METRICS", "VALUE")
  #Merge the new long format metrics data frame with the quantile (percentile) values.
  if(any(melted$METRICS %in% quant.df$METRICS)){
    long.df <- merge(melted, quant.df, by = "METRICS")
    
    #Create a new data frame of just reference values
    long.ref <- long.df[long.df$CATEGORY %in% upper.class, ]
    #Column numbers can change easily. The script below specifies column "0%" to column "100%."
    # These columns represent the percentile values.
    ref.columns <-  which(colnames(long.ref) == "0%") : which(colnames(long.ref) == "100%")
    #Looking for the percentage of sites correctly identified as reference or degraded based
    # on each threshold.  The thresholds are defined by the reference percentiles. If a site
    # is correctly identified as reference, then a 1 is returned. If the site is incorrectly
    # identified as degraded, then a 0 is returned.  Essentially, 1 is equivalent to "yes"
    # and 0 is equivalent to "no."  This ifelse statement is specific to reference sites.
    # Below the ifelse statement is specific to degraded sites. If the metric decreases
    # with disturbance and the raw metric score for a sampling event
    # is greater than the percentile value, then the site was correctly identified as
    # a reference site and a 1 is returned.
    
    
    long.ref[, ref.columns] <- ifelse((long.ref$VALUE >= long.ref[, ref.columns] &
                                         long.ref$CATEGORY == upper.class &
                                         long.ref$DISTURBANCE == "DECREASE") |
                                        (long.ref$VALUE <= long.ref[, ref.columns] &
                                           long.ref$CATEGORY == upper.class &
                                           long.ref$DISTURBANCE == "INCREASE"), 1, 0)
   
    #Transform the long reference data frame to a wide format.
    melted.ref <- reshape2::melt(long.ref, id.vars = c("METRICS", "EVENT_ID",
                                                       "STATION_ID", "CATEGORY", "DATE",
                                                       "SAMPLE_NUMBER", "AGENCY_CODE",
                                                       "VALUE", "DISTURBANCE"))
    
    #Aggregate the values (1's and 0's) by distrubance, metric, and variable (percentile).
    # The aggregation function finds the mean of the values (1's and 0's) and multiplies
    # the mean by 100. This value represents the percentage of reference sites correctly
    # identified as reference sites for a particular metric at a particular
    # threshold (percentile).
    pct.ref <- aggregate(value ~ DISTURBANCE + METRICS + variable, data = melted.ref,
                         function(x) mean(x) * 100)
    colnames(pct.ref) <- c("DISTURBANCE", "METRICS", "PERCENTILE", "PCT_REF")
    
    #Create a new data frame of just degraded values.
    long.deg <- long.df[long.df$CATEGORY %in% lower.class, ]
    #Column numbers can change easily. The script below specifies column "0%" to column "100%."
    # These columns represent the percentile values.
    deg.columns <-  which(colnames(long.deg) == "0%") : which(colnames(long.deg) == "100%")
    #See above for further description. Same process for identifing the number
    # of reference sites correctly identified at each threshold but the script
    # below is for correctly identified degraded sites.
    long.deg[, deg.columns] <- ifelse((long.deg$VALUE < long.deg[, deg.columns] &
                                         long.deg$DISTURBANCE == "DECREASE") |
                                        (long.deg$VALUE > long.deg[, deg.columns] &
                                           long.deg$DISTURBANCE == "INCREASE"), 1, 0)

    #Transform the long degraded data frame to a wide format.
    melted.deg <- reshape2::melt(long.deg, id.vars = c("METRICS", "EVENT_ID",
                                                       "STATION_ID", "CATEGORY", "DATE",
                                                       "SAMPLE_NUMBER", "AGENCY_CODE",
                                                       "VALUE", "DISTURBANCE"))
    #Aggregate the values (1's and 0's) by distrubance, metric, and variable (percentile).
    # The aggregation function finds the mean of the values (1's and 0's) and multiplies
    # the mean by 100. This value represents the percentage of degraded sites correctly
    # identified as degraded sites for a particular metric at a particular
    # threshold (percentile).
    pct.deg <- aggregate(value ~ DISTURBANCE + METRICS + variable , data = melted.deg,
                         function(x) mean(x) * 100)
    #dt <- data.table::data.table(melted.deg)
    #pct.deg <- dt[, mean(value) * 100, by = list(DISTURBANCE, METRICS, variable)]
    colnames(pct.deg) <- c("DISTURBANCE", "METRICS", "PERCENTILE", "PCT_DEG")
    
    #Merge the two tables containing the percent of reference and the percent of
    # degraded correctly identified.
    merge.pct <- merge(pct.ref, pct.deg, by = c("DISTURBANCE", "METRICS", "PERCENTILE"))
    #Calculate the discrimination efficiency of each threshold.
    # DE = (Correctly identified Reference + Correctly identified Degraded) / 2
    merge.pct$SENSITIVITY <- merge.pct$PCT_REF + merge.pct$PCT_DEG
    
    tp <- (merge.pct$PCT_REF * nrow(unique(long.df[long.df$CATEGORY %in% upper.class, c("EVENT_ID", "CATEGORY")]))) / 100
    fn <- nrow(unique(long.df[long.df$CATEGORY %in% upper.class, c("EVENT_ID", "CATEGORY")])) - tp
    tn <- (merge.pct$PCT_DEG * nrow(unique(long.df[long.df$CATEGORY %in% lower.class, c("EVENT_ID", "CATEGORY")]))) / 100
    fp <- nrow(unique(long.df[long.df$CATEGORY %in% lower.class, c("EVENT_ID", "CATEGORY")])) - tn
    
    merge.pct$ACCURACY <- ((tp + tn) / (tp + tn + fp + fn)) * 100
    merge.pct$TPR <- tp / (tp + fn)
    merge.pct$FPR <- fp / (fp + tn)
    merge.pct$FNR <- fn / (tp + fn)
    merge.pct$New_SENSITIVITY <- (merge.pct$TPR + (merge.pct$ACCURACY) / 100) - (merge.pct$FPR + merge.pct$FNR + (abs(merge.pct$PCT_REF - merge.pct$PCT_DEG) / 100))
   
    #Aggregate the table to select the best DE score.  That is the max DE score for
    # each metric.
    agg.df <- aggregate(SENSITIVITY ~ METRICS , FUN = max, data= merge.pct)
    
    
    
    # Merge the aggregated data frame to the merged.pct data frame.  Represents the
    # percentile with the best DE score. However, several percentiles may have the same
    # DE score. Hence, "almost_best."
    almost_best.df <- unique(merge(agg.df, merge.pct, by = c("METRICS", "SENSITIVITY")))
    
    #Create a table to count the number of metrics with multiple thresholds.
    metrics.table <- data.frame(table(almost_best.df$METRICS))
    colnames(metrics.table) <- c("Metric", "Count")
    
    #A data frame containing only metrics with multiple thresholds.
    metric.repeats <- metrics.table[metrics.table$Count > 1, ]
    
    #A data frame containing only metrics with a single threshold.
    fine <- almost_best.df[!almost_best.df$METRICS %in% metric.repeats$Metric, ]
    
    #Use the list of metrics from metric.repeats to further inspect
    # thresholds in almost_best.df.
    mult.metrics <- almost_best.df[almost_best.df$METRICS %in% metric.repeats$Metric, ]
    #The balance issue occurs infrequently but must be accounted for.
    # Sometimes the same DE score can be attained from different pct_ref and pct_deg values.
    # For example, 90% pct_ref and 10% pct_deg is equal to 50% pct_ref and 50% pct_deg.
    # In these cases we prefer the more balanced solution (i.e., 50% and 50%) because both
    # groups are better represented.  If the DE score is a product of 90% and 10%, then the
    # majority of the data is being binned into a single group.  Therefore, there is poor
    # distinction between reference and degraded conditions.
    # Balance = | pct_ref - pct_deg |
    # The smaller the value the better the balance between the two conditions.
    mult.metrics$BALANCE <- abs(mult.metrics$PCT_REF - mult.metrics$PCT_DEG)
    balanced.df <- plyr::ddply(mult.metrics, plyr::.(METRICS),
                               function(x) x[which(x$BALANCE == min(x$BALANCE)), ])
    balanced.df <- balanced.df[ , -which(names(balanced.df) %in% c("BALANCE"))]
    
    #Create new data frames for each disturbance category. Disturbance indicates
    # how the metric responds to disturbance.
    dec.df <- balanced.df[balanced.df$DISTURBANCE == "DECREASE", ]
    inc.df <- balanced.df[balanced.df$DISTURBANCE == "INCREASE", ]
    equ.df <-balanced.df[balanced.df$DISTURBANCE == "EQUAL", ]
    
    #Most of the time the balance check does not uliminate multiple thresholds.
    # This final step selects the single best threshold for each metric.
    # If the metric decreases with disturbance, the lowest percentile is
    # selected as the threshold.  If the metric increases with disturbance
    # the largest percentile is selected as the threshold. If the metric
    # cannot distinquish between reference and degraded then no threshold
    # is returned and all values equal zero.
    final.dec <- plyr::ddply(dec.df, plyr::.(METRICS), function(x) x[which.min(x$PERCENTILE), ])
    final.inc <- plyr::ddply(inc.df, plyr::.(METRICS), function(x) x[which.max(x$PERCENTILE), ])
    final.equ <- plyr::ddply(equ.df, plyr::.(METRICS), function(x) x[which.max(x$PERCENTILE), ])
    #Join all of the final threshold values together.
    bound.df <- rbind(fine, final.dec, final.inc, final.equ)
    
    #Use the ref.df data frame created in the beginning of the function
    # to report the reference median value.
    med.df <- data.frame(sapply(ref.df[, 7:ncol(ref.df)], quantile, 0.50))
    #Remove .50% from row names
    med.df$METRICS <- gsub(".50%", "", row.names(med.df))
    colnames(med.df) <- c("MEDIAN", "METRICS") #Rename the columns
    #Join the reference median values for each metric to the bound
    # data frame with the SENSITIVITYs and percentile thresholds
    # for each metric. af = almost final
    af.df <- merge(bound.df, med.df, by = "METRICS")
    
    #Use the quant.ref data frame created in the beginning of the function
    # to represent the actual threshold value that the chosen percentile
    # represents.
    #Create a new column from row names.
    quant.ref$PERCENTILE <- row.names(quant.ref)
    
    #Transform the data frame from a wide to long format.
    melted.pct <- reshape2::melt(quant.ref, id.vars = "PERCENTILE")
    names(melted.pct) <- c("PERCENTILE", "METRICS", "THRESHOLD")#change column names
    #Merge the threshold values with the almost final data frame (af.df).
    final.df <- merge(af.df, melted.pct, by = c("METRICS", "PERCENTILE"), all.x = TRUE)
    #Round the values to the hundreths place.
    final.df[, c("THRESHOLD", "SENSITIVITY", "PCT_REF", "PCT_DEG",
                 "MEDIAN")] <- round(final.df[, c("THRESHOLD", "SENSITIVITY", "PCT_REF",
                                                  "PCT_DEG", "MEDIAN")], digits = 2)
    final.df <- final.df[, c("METRICS", "DISTURBANCE", "SENSITIVITY", "PERCENTILE",
                             "PCT_REF", "PCT_DEG", "MEDIAN", "THRESHOLD")]
    return(final.df)
  }else{
    return(NULL)
  }
  
}

#==============================================================================
#'CMA
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@param quant.df = Data frame containing upper.class quantile values.
#'@param upper.class = the site class that represents the better condition.
#'@param lower.class = the site class that represents the poorer condition.
#'@param ref.df = a data frame of only the upper.class values.
#'@param quant.ref = a data frame of reference quantile values.
#'@return Determines the threshold at which a metric best categorizes
#'reference and degraded stations.
#'@export
#'

cma <- function(metrics.df, quant.df, upper.class, lower.class, ref.df, quant.ref){
  #Transform the metrics data frame from a wide format to a long data format.
  melted <- reshape2::melt(metrics.df, id.vars = c("EVENT_ID", "STATION_ID", "CATEGORY",
                                                   "DATE", "SAMPLE_NUMBER", "AGENCY_CODE"))
  colnames(melted) <- c("EVENT_ID", "STATION_ID", "CATEGORY", "DATE",
                        "SAMPLE_NUMBER", "AGENCY_CODE",  "METRICS", "VALUE")
  #Merge the new long format metrics data frame with the quantile (percentile) values.
  if(any(melted$METRICS %in% quant.df$METRICS)){
    long.df <- merge(melted, quant.df, by = "METRICS")
    
    #Create a new data frame of just reference values
    long.ref <- long.df[long.df$CATEGORY %in% upper.class, ]
    #Column numbers can change easily. The script below specifies column "0%" to column "100%."
    # These columns represent the percentile values.
    ref.columns <-  which(colnames(long.ref) == "0%") : which(colnames(long.ref) == "100%")
    #Looking for the percentage of sites correctly identified as reference or degraded based
    # on each threshold.  The thresholds are defined by the reference percentiles. If a site
    # is correctly identified as reference, then a 1 is returned. If the site is incorrectly
    # identified as degraded, then a 0 is returned.  Essentially, 1 is equivalent to "yes"
    # and 0 is equivalent to "no."  This ifelse statement is specific to reference sites.
    # Below the ifelse statement is specific to degraded sites. If the metric decreases
    # with disturbance and the raw metric score for a sampling event
    # is greater than the percentile value, then the site was correctly identified as
    # a reference site and a 1 is returned.
    
    
    long.ref[, ref.columns] <- ifelse((long.ref$VALUE >= long.ref[, ref.columns] &
                                         long.ref$CATEGORY == upper.class &
                                         long.ref$DISTURBANCE == "DECREASE") |
                                        (long.ref$VALUE <= long.ref[, ref.columns] &
                                           long.ref$CATEGORY == upper.class &
                                           long.ref$DISTURBANCE == "INCREASE"), 1, 0)
    
    #Transform the long reference data frame to a wide format.
    melted.ref <- reshape2::melt(long.ref, id.vars = c("METRICS", "EVENT_ID",
                                                       "STATION_ID", "CATEGORY", "DATE",
                                                       "SAMPLE_NUMBER", "AGENCY_CODE",
                                                       "VALUE", "DISTURBANCE"))
    
    #Aggregate the values (1's and 0's) by distrubance, metric, and variable (percentile).
    # The aggregation function finds the mean of the values (1's and 0's) and multiplies
    # the mean by 100. This value represents the percentage of reference sites correctly
    # identified as reference sites for a particular metric at a particular
    # threshold (percentile).
    pct.ref <- aggregate(value ~ DISTURBANCE + METRICS + variable, data = melted.ref,
                         function(x) mean(x) * 100)
    colnames(pct.ref) <- c("DISTURBANCE", "METRICS", "PERCENTILE", "PCT_REF")
    
    #Create a new data frame of just degraded values.
    long.deg <- long.df[long.df$CATEGORY %in% lower.class, ]
    #Column numbers can change easily. The script below specifies column "0%" to column "100%."
    # These columns represent the percentile values.
    deg.columns <-  which(colnames(long.deg) == "0%") : which(colnames(long.deg) == "100%")
    #See above for further description. Same process for identifing the number
    # of reference sites correctly identified at each threshold but the script
    # below is for correctly identified degraded sites.
    long.deg[, deg.columns] <- ifelse((long.deg$VALUE < long.deg[, deg.columns] &
                                         long.deg$DISTURBANCE == "DECREASE") |
                                        (long.deg$VALUE > long.deg[, deg.columns] &
                                           long.deg$DISTURBANCE == "INCREASE"), 1, 0)
    
    #Transform the long degraded data frame to a wide format.
    melted.deg <- reshape2::melt(long.deg, id.vars = c("METRICS", "EVENT_ID",
                                                       "STATION_ID", "CATEGORY", "DATE",
                                                       "SAMPLE_NUMBER", "AGENCY_CODE",
                                                       "VALUE", "DISTURBANCE"))
    #Aggregate the values (1's and 0's) by distrubance, metric, and variable (percentile).
    # The aggregation function finds the mean of the values (1's and 0's) and multiplies
    # the mean by 100. This value represents the percentage of degraded sites correctly
    # identified as degraded sites for a particular metric at a particular
    # threshold (percentile).
    pct.deg <- aggregate(value ~ DISTURBANCE + METRICS + variable , data = melted.deg,
                         function(x) mean(x) * 100)
    #dt <- data.table::data.table(melted.deg)
    #pct.deg <- dt[, mean(value) * 100, by = list(DISTURBANCE, METRICS, variable)]
    colnames(pct.deg) <- c("DISTURBANCE", "METRICS", "PERCENTILE", "PCT_DEG")
    
    #Merge the two tables containing the percent of reference and the percent of
    # degraded correctly identified.
    merge.pct <- merge(pct.ref, pct.deg, by = c("DISTURBANCE", "METRICS", "PERCENTILE"))
    #Calculate the discrimination efficiency of each threshold.
    # DE = (Correctly identified Reference + Correctly identified Degraded) / 2
    merge.pct$SENSITIVITY <- (merge.pct$PCT_REF + merge.pct$PCT_DEG) / 2
    
    tp <- (merge.pct$PCT_REF * nrow(unique(long.df[long.df$CATEGORY %in% upper.class, c("EVENT_ID", "CATEGORY")]))) / 100
    fn <- nrow(unique(long.df[long.df$CATEGORY %in% upper.class, c("EVENT_ID", "CATEGORY")])) - tp
    tn <- (merge.pct$PCT_DEG * nrow(unique(long.df[long.df$CATEGORY %in% lower.class, c("EVENT_ID", "CATEGORY")]))) / 100
    fp <- nrow(unique(long.df[long.df$CATEGORY %in% lower.class, c("EVENT_ID", "CATEGORY")])) - tn
    
    merge.pct$ACCURACY <- ((tp + tn) / (tp + tn + fp + fn)) * 100
    merge.pct$TPR <- tp / (tp + fn)
    merge.pct$FPR <- fp / (fp + tn)
    merge.pct$FNR <- fn / (tp + fn)
    #merge.pct$NEW_SENSITIVITY <- (merge.pct$TPR + (merge.pct$ACCURACY) / 100) - (merge.pct$FPR + merge.pct$FNR + (abs(merge.pct$PCT_REF - merge.pct$PCT_DEG) / 100))
    merge.pct$NEW_SENSITIVITY <- ((merge.pct$PCT_REF + merge.pct$PCT_DEG) / 2) - abs(merge.pct$PCT_REF - merge.pct$PCT_DEG)
    #Aggregate the table to select the best DE score.  That is the max DE score for
    # each metric.
    agg.df <- aggregate(NEW_SENSITIVITY ~ METRICS , FUN = max, data= merge.pct)
    
    
    
    # Merge the aggregated data frame to the merged.pct data frame.  Represents the
    # percentile with the best DE score. However, several percentiles may have the same
    # DE score. Hence, "almost_best."
    almost_best.df <- unique(merge(agg.df, merge.pct, by = c("METRICS", "NEW_SENSITIVITY")))
    
    #Create a table to count the number of metrics with multiple thresholds.
    metrics.table <- data.frame(table(almost_best.df$METRICS))
    colnames(metrics.table) <- c("Metric", "Count")
    
    #A data frame containing only metrics with multiple thresholds.
    metric.repeats <- metrics.table[metrics.table$Count > 1, ]
    
    #A data frame containing only metrics with a single threshold.
    fine <- almost_best.df[!almost_best.df$METRICS %in% metric.repeats$Metric, ]
    
    #Use the list of metrics from metric.repeats to further inspect
    # thresholds in almost_best.df.
    mult.metrics <- almost_best.df[almost_best.df$METRICS %in% metric.repeats$Metric, ]
    #The balance issue occurs infrequently but must be accounted for.
    # Sometimes the same DE score can be attained from different pct_ref and pct_deg values.
    # For example, 90% pct_ref and 10% pct_deg is equal to 50% pct_ref and 50% pct_deg.
    # In these cases we prefer the more balanced solution (i.e., 50% and 50%) because both
    # groups are better represented.  If the DE score is a product of 90% and 10%, then the
    # majority of the data is being binned into a single group.  Therefore, there is poor
    # distinction between reference and degraded conditions.
    # Balance = | pct_ref - pct_deg |
    # The smaller the value the better the balance between the two conditions.
    mult.metrics$BALANCE <- abs(mult.metrics$PCT_REF - mult.metrics$PCT_DEG)
    balanced.df <- plyr::ddply(mult.metrics, plyr::.(METRICS),
                               function(x) x[which(x$BALANCE == min(x$BALANCE)), ])
    balanced.df <- balanced.df[ , -which(names(balanced.df) %in% c("BALANCE"))]
    
    #Create new data frames for each disturbance category. Disturbance indicates
    # how the metric responds to disturbance.
    dec.df <- balanced.df[balanced.df$DISTURBANCE == "DECREASE", ]
    inc.df <- balanced.df[balanced.df$DISTURBANCE == "INCREASE", ]
    equ.df <-balanced.df[balanced.df$DISTURBANCE == "EQUAL", ]
    
    #Most of the time the balance check does not uliminate multiple thresholds.
    # This final step selects the single best threshold for each metric.
    # If the metric decreases with disturbance, the lowest percentile is
    # selected as the threshold.  If the metric increases with disturbance
    # the largest percentile is selected as the threshold. If the metric
    # cannot distinquish between reference and degraded then no threshold
    # is returned and all values equal zero.
    final.dec <- plyr::ddply(dec.df, plyr::.(METRICS), function(x) x[which.min(x$PERCENTILE), ])
    final.inc <- plyr::ddply(inc.df, plyr::.(METRICS), function(x) x[which.max(x$PERCENTILE), ])
    final.equ <- plyr::ddply(equ.df, plyr::.(METRICS), function(x) x[which.max(x$PERCENTILE), ])
    #Join all of the final threshold values together.
    bound.df <- rbind(fine, final.dec, final.inc, final.equ)
    
    #Use the ref.df data frame created in the beginning of the function
    # to report the reference median value.
    if(ncol(ref.df) > 7){
      med.df <- data.frame(sapply(ref.df[, 7:ncol(ref.df)], quantile, 0.50, na.rm = TRUE))
      #Remove .50% from row names
      med.df$METRICS <- gsub(".50%", "", row.names(med.df))
      colnames(med.df) <- c("MEDIAN", "METRICS") #Rename the columns
    } 
    if(ncol(ref.df) <= 7){
      med.df <- data.frame(METRICS = colnames(ref.df)[7])
      med.df$MEDIAN <- quantile(ref.df[, 7], 0.50, na.rm = TRUE)
    } 
    
    #Join the reference median values for each metric to the bound
    # data frame with the SENSITIVITYs and percentile thresholds
    # for each metric. af = almost final
    af.df <- merge(bound.df, med.df, by = "METRICS")
    
    #Use the quant.ref data frame created in the beginning of the function
    # to represent the actual threshold value that the chosen percentile
    # represents.
    #Create a new column from row names.
    quant.ref$PERCENTILE <- row.names(quant.ref)
    
    #Transform the data frame from a wide to long format.
    melted.pct <- reshape2::melt(quant.ref, id.vars = "PERCENTILE")
    names(melted.pct) <- c("PERCENTILE", "METRICS", "THRESHOLD")#change column names
    #Merge the threshold values with the almost final data frame (af.df).
    final.df <- merge(af.df, melted.pct, by = c("METRICS", "PERCENTILE"), all.x = TRUE)
    #Round the values to the hundreths place.
    final.df[, c("THRESHOLD", "NEW_SENSITIVITY", "SENSITIVITY", "PCT_REF", "PCT_DEG",
                 "MEDIAN")] <- round(final.df[, c("THRESHOLD", "NEW_SENSITIVITY",
                                                  "SENSITIVITY","PCT_REF",
                                                  "PCT_DEG", "MEDIAN")], digits = 2)
    final.df <- final.df[, c("METRICS", "DISTURBANCE", "SENSITIVITY",
                             "NEW_SENSITIVITY", "PERCENTILE",
                             "PCT_REF", "PCT_DEG", "MEDIAN", "THRESHOLD")]
    final.df$BOUND <- ifelse(final.df$DISTURBANCE %in% "DECREASE",
                             final.df$THRESHOLD - abs(final.df$MEDIAN - final.df$THRESHOLD),
                             ifelse(final.df$DISTURBANCE %in% "INCREASE",
                                    final.df$THRESHOLD + abs(final.df$MEDIAN - final.df$THRESHOLD), "ERROR"))
    
    return(final.df)
  }else{
    return(NULL)
  }
  
}


#==============================================================================
#'Discrimination Efficiency
#'
#'@param deg.df = a data frame of only the lower.class values.
#'@param quant.df = Data frame containing upper.class quantile values.
#'@return Determines the threshold at which a metric best categorizes
#'reference and degraded stations.
#'@export

d_e <- function(deg.df, quant.df){
  ID <- c("EVENT_ID", "CATEGORY",
          "STATION_ID", "DATE",
          "AGENCY_CODE", "SAMPLE_NUMBER")

  deg.long <- data.frame(t(deg.df[, 7:ncol(deg.df)]))
  names(deg.long) <- deg.df$EVENT_ID
  deg.long$METRICS <- row.names(deg.long)
  merged <- merge(quant.df, deg.long, by = "METRICS")

  melted <- reshape2::melt(deg.df, id.vars= ID)
  names(melted) <- c(ID, "METRICS", "VALUES")
  merged <- merge(quant.df, melted, by = "METRICS")

  merged$A <- ifelse(merged$DISTURBANCE == "DECREASE" &
                       merged$VALUES < merged$`25%`, 1,
                     ifelse(merged$DISTURBANCE == "DECREASE" &
                              merged$VALUES >= merged$`25%`, 0,
                            ifelse(merged$DISTURBANCE == "INCREASE" &
                                     merged$VALUES > merged$`75%`, 1,
                                   ifelse(merged$DISTURBANCE == "INCREASE" &
                                            merged$VALUES <= merged$`75%`, 0,
                                          ifelse(merged$DISTURBANCE == "EQUAL", 0, 1000)))))
  if(!is.na(sum(merged$A))){
    agg <- aggregate(A ~ METRICS + CATEGORY + DISTURBANCE, data = merged, FUN = sum)
    agg$SENSITIVITY <- (agg$A / nrow(deg.df)) * 100
    
    final.df <- agg[,!(names(agg) %in% c("CATEGORY", "A"))]
  }else{
    final.df <- data.frame(METRICS = names(deg.df)[!names(deg.df) %in% c("EVENT_ID", "STATION_ID",
                                                    "SAMPLE_NUMBER", "AGENCY_CODE",
                                                    "DATE", "CATEGORY")])
    final.df$DISTURBANCE <- NA
    final.df$SENSITIVITY <- NA
  }

  return(final.df)
}

#==============================================================================
# Discrimination Efficiency
#==============================================================================
#'Optimal Discrimination Efficiency
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@param quant.df = Data frame containing upper.class quantile values.
#'@param upper.class = the site class that represents the better condition.
#'@param lower.class = the site class that represents the poorer condition.
#'@param ref.df = a data frame of only the upper.class values.
#'@param quant.ref = a data frame of reference quantile values.
#'@return Determines the threshold at which a metric best categorizes
#'reference and degraded stations.
#'@export
#'


sse <- function(metrics.df, quant.df, upper.class, lower.class, ref.df, quant.ref){
  #Transform the metrics data frame from a wide format to a long data format.
  melted <- reshape2::melt(metrics.df, id.vars = c("EVENT_ID", "STATION_ID", "CATEGORY",
                                                   "DATE", "SAMPLE_NUMBER", "AGENCY_CODE"))
  colnames(melted) <- c("EVENT_ID", "STATION_ID", "CATEGORY", "DATE",
                        "SAMPLE_NUMBER", "AGENCY_CODE",  "METRICS", "VALUE")
  #Merge the new long format metrics data frame with the quantile (percentile) values.
  long.df <- merge(melted, quant.df, by = "METRICS")

  #Create a new data frame of just reference values
  long.ref <- long.df[long.df$CATEGORY %in% upper.class, ]
  #Column numbers can change easily. The script below specifies column "0%" to column "100%."
  # These columns represent the percentile values.
  ref.columns <-  which(colnames(long.ref) == "0%") : which(colnames(long.ref) == "100%")
  #Looking for the percentage of sites correctly identified as reference or degraded based
  # on each threshold.  The thresholds are defined by the reference percentiles. If a site
  # is correctly identified as reference, then a 1 is returned. If the site is incorrectly
  # identified as degraded, then a 0 is returned.  Essentially, 1 is equivalent to "yes"
  # and 0 is equivalent to "no."  This ifelse statement is specific to reference sites.
  # Below the ifelse statement is specific to degraded sites. If the metric decreases
  # with disturbance and the raw metric score for a sampling event
  # is greater than the percentile value, then the site was correctly identified as
  # a reference site and a 1 is returned.

  long.ref[,ref.columns] <- ifelse((long.ref$VALUE > long.ref[, ref.columns] &
                                      long.ref$CATEGORY == upper.class &
                                      long.ref$DISTURBANCE == "DECREASE") |
                                     (long.ref$VALUE < long.ref[, ref.columns] &
                                        long.ref$CATEGORY == upper.class &
                                        long.ref$DISTURBANCE == "INCREASE"), 1, 0)

  #Transform the long reference data frame to a wide format.
  melted.ref <- reshape2::melt(long.ref, id.vars = c("METRICS", "EVENT_ID",
                                                     "STATION_ID", "CATEGORY", "DATE",
                                                     "SAMPLE_NUMBER", "AGENCY_CODE",
                                                     "VALUE", "DISTURBANCE"))
  #Aggregate the values (1's and 0's) by distrubance, metric, and variable (percentile).
  # The aggregation function finds the mean of the values (1's and 0's) and multiplies
  # the mean by 100. This value represents the percentage of reference sites correctly
  # identified as reference sites for a particular metric at a particular
  # threshold (percentile).
  pct.ref <- aggregate(value ~ DISTURBANCE + METRICS + variable, data = melted.ref,
                       function(x) mean(x) * 100)
  colnames(pct.ref) <- c("DISTURBANCE", "METRICS", "PERCENTILE", "PCT_REF")

  #Create a new data frame of just degraded values.
  long.deg <- long.df[long.df$CATEGORY %in% lower.class, ]
  #Column numbers can change easily. The script below specifies column "0%" to column "100%."
  # These columns represent the percentile values.
  deg.columns <-  which(colnames(long.deg) == "0%") : which(colnames(long.deg) == "100%")
  #See above for further description. Same process for identifing the number
  # of reference sites correctly identified at each threshold but the script
  # below is for correctly identified degraded sites.
  long.deg[, deg.columns] <- ifelse((long.deg$VALUE < long.deg[, deg.columns] &
                                       long.deg$DISTURBANCE == "DECREASE") |
                                      (long.deg$VALUE > long.deg[, deg.columns] &
                                         long.deg$DISTURBANCE == "INCREASE"), 1, 0)
  #Transform the long degraded data frame to a wide format.
  melted.deg <- reshape2::melt(long.deg, id.vars = c("METRICS", "EVENT_ID",
                                                     "STATION_ID", "CATEGORY", "DATE",
                                                     "SAMPLE_NUMBER", "AGENCY_CODE",
                                                     "VALUE", "DISTURBANCE"))
  #Aggregate the values (1's and 0's) by distrubance, metric, and variable (percentile).
  # The aggregation function finds the mean of the values (1's and 0's) and multiplies
  # the mean by 100. This value represents the percentage of degraded sites correctly
  # identified as degraded sites for a particular metric at a particular
  # threshold (percentile).
  pct.deg <- aggregate(value ~ DISTURBANCE + METRICS + variable , data = melted.deg, FUN = sum)
  colnames(pct.deg) <- c("DISTURBANCE", "METRICS", "PERCENTILE", "PCT_DEG")

  #Merge the two tables containing the percent of reference and the percent of
  # degraded correctly identified.
  merge.pct <- merge(pct.ref, pct.deg, by = c("DISTURBANCE", "METRICS", "PERCENTILE"))
  #Calculate the discrimination efficiency of each threshold.
  # DE = (Correctly identified Reference + Correctly identified Degraded) / 2




  #  #merge.pct$SENSITIVITY <- (merge.pct$PCT_REF + merge.pct$PCT_DEG)/2

  #=================================================================================
  merge.pct$SENSITIVITY <- abs((merge.pct$PCT_REF(merge.pct$PCT_DEG +
                           (nrow(metrics.df[metrics.df$CATEGORY %in% lower.class, ])))) -
                           (merge.pct$PCT_DEG(merge.pct$PCT_REF +
                           (nrow(metrics.df[metrics.df$CATEGORY %in% upper.class, ])))))
  #=================================================================================
  #Aggregate the table to select the best DE score.  That is the max DE score for
  # each metric.
  agg.df <- aggregate(SENSITIVITY ~ METRICS , FUN = min, data= merge.pct)
  # Merge the aggregated data frame to the merged.pct data frame.  Represents the
  # percentile with the best DE score. However, several percentiles may have the same
  # DE score. Hence, "almost_best."
  almost_best.df <- unique(merge(agg.df, merge.pct, by = c("METRICS", "SENSITIVITY")))

  #Create a table to count the number of metrics with multiple thresholds.
  metrics.table <- data.frame(table(almost_best.df$METRICS))
  colnames(metrics.table) <- c("Metric", "Count")

  #A data frame containing only metrics with multiple thresholds.
  metric.repeats <- metrics.table[metrics.table$Count > 1, ]

  #A data frame containing only metrics with a single threshold.
  fine <- almost_best.df[!almost_best.df$METRICS %in% metric.repeats$Metric, ]

  #Use the list of metrics from metric.repeats to further inspect
  # thresholds in almost_best.df.
  mult.metrics <- almost_best.df[almost_best.df$METRICS %in% metric.repeats$Metric, ]
  #The balance issue occurs infrequently but must be accounted for.
  # Sometimes the same DE score can be attained from different pct_ref and pct_deg values.
  # For example, 90% pct_ref and 10% pct_deg is equal to 50% pct_ref and 50% pct_deg.
  # In these cases we prefer the more balanced solution (i.e., 50% and 50%) because both
  # groups are better represented.  If the DE score is a product of 90% and 10%, then the
  # majority of the data is being binned into a single group.  Therefore, there is poor
  # distinction between reference and degraded conditions.
  # Balance = | pct_ref - pct_deg |
  # The smaller the value the better the balance between the two conditions.
  mult.metrics$BALANCE <- abs(mult.metrics$PCT_REF - mult.metrics$PCT_DEG)
  balanced.df <- plyr::ddply(mult.metrics, plyr::.(METRICS),
                             function(x) x[which(x$BALANCE == min(x$BALANCE)), ])
  balanced.df <- balanced.df[ , -which(names(balanced.df) %in% c("BALANCE"))]

  #Create new data frames for each disturbance category. Disturbance indicates
  # how the metric responds to disturbance.
  dec.df <- balanced.df[balanced.df$DISTURBANCE == "DECREASE", ]
  inc.df <- balanced.df[balanced.df$DISTURBANCE == "INCREASE", ]
  equ.df <-balanced.df[balanced.df$DISTURBANCE == "EQUAL", ]

  #Most of the time the balance check does not uliminate multiple thresholds.
  # This final step selects the single best threshold for each metric.
  # If the metric decreases with disturbance, the lowest percentile is
  # selected as the threshold.  If the metric increases with disturbance
  # the largest percentile is selected as the threshold. If the metric
  # cannot distinquish between reference and degraded then no threshold
  # is returned and all values equal zero.
  final.dec <- plyr::ddply(dec.df, plyr::.(METRICS), function(x) x[which.min(x$PERCENTILE), ])
  final.inc <- plyr::ddply(inc.df, plyr::.(METRICS), function(x) x[which.max(x$PERCENTILE), ])
  final.equ <- plyr::ddply(equ.df, plyr::.(METRICS), function(x) x[which.max(x$PERCENTILE), ])
  #Join all of the final threshold values together.
  bound.df <- rbind(fine, final.dec, final.inc, final.equ)

  #Use the ref.df data frame created in the beginning of the function
  # to report the reference median value.
  med.df <- data.frame(sapply(ref.df[, 7:ncol(ref.df)], quantile, 0.50))
  #Remove .50% from row names
  med.df$METRICS <- gsub(".50%", "", row.names(med.df))
  colnames(med.df) <- c("MEDIAN", "METRICS") #Rename the columns
  #Join the reference median values for each metric to the bound
  # data frame with the SENSITIVITYs and percentile thresholds
  # for each metric. af = almost final
  af.df <- merge(bound.df, med.df, by = "METRICS")

  #Use the quant.ref data frame created in the beginning of the function
  # to represent the actual threshold value that the chosen percentile
  # represents.
  #Create a new column from row names.
  quant.ref$PERCENTILE <- row.names(quant.ref)

  #Transform the data frame from a wide to long format.
  melted.pct <- reshape2::melt(quant.ref, id.vars = "PERCENTILE")
  names(melted.pct) <- c("PERCENTILE", "METRICS", "THRESHOLD")#change column names
  #Merge the threshold values with the almost final data frame (af.df).
  final.df <- merge(af.df, melted.pct, by = c("METRICS", "PERCENTILE"), all.x = TRUE)
  #Round the values to the hundreths place.
  final.df[, c("THRESHOLD", "SENSITIVITY", "PCT_REF", "PCT_DEG",
               "MEDIAN")] <- round(final.df[, c("THRESHOLD", "SENSITIVITY", "PCT_REF",
                                                "PCT_DEG", "MEDIAN")], digits = 2)
  final.df <- final.df[, c("METRICS", "DISTURBANCE", "SENSITIVITY", "PERCENTILE",
                           "PCT_REF", "PCT_DEG", "MEDIAN", "THRESHOLD")]
  return(final.df)
}

