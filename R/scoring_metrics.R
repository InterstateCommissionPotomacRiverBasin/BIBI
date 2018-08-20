#==============================================================================
#Scoring metrics
#==============================================================================
#'Score Samples
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param index.metrics = a data frame containing the metrics selected for
#'the index.
#'@param method = Designates the the method(s) used to score the samples.
#'@return Scores the samples.
#'@export

score_samples <- function(metrics.df, index.metrics, method = "all"){
  ID <- c("EVENT_ID", "STATION_ID","DATE", "AGENCY_CODE",
          "SAMPLE_NUMBER", "CATEGORY")
  best.df <- metrics.df[, names(metrics.df) %in% index.metrics$METRICS]
  best.df <- cbind(metrics.df[, ID], best.df)

  ref.df <- best.df[best.df$CATEGORY %in% "REF", ]
  ref.quant <- data.frame(t(sapply(ref.df[, 7:ncol(ref.df)], quantile,
                                   c(0.05, 0.25, 0.50, 0.75, 0.95))))
  names(ref.quant) <- c("5%", "25%", "50%", "75%", "95%")
  ref.quant$METRICS <- row.names(ref.quant)

  score.info <- merge(index.metrics, ref.quant, by = "METRICS")
  long.df <- reshape2::melt(metrics.df,
                            id.vars = c("EVENT_ID", "STATION_ID",
                                        "DATE", "AGENCY_CODE",
                                        "SAMPLE_NUMBER", "CATEGORY"))
  names(long.df) <- c(ID, "METRICS", "REPORTING_VALUE")
  long.df <- merge(long.df, score.info, by = "METRICS")

  if(method %in% "ref_categorical") long.df <- ref_categorical(long.df)
  if(method %in% "ref_threshold" & !(is.null(long.df$THRESHOLD))) long.df <- ref_threshold(long.df)
  if(method %in% "ref_threshold_mod" & !(is.null(long.df$THRESHOLD))) long.df <- ref_threshold_mod(long.df)
  if(method %in% "ref_gradient") long.df <- ref_gradient(long.df)
  if(method %in% "all_gradient") long.df <- all_gradient(metrics.df, index.metrics)
  if(method %in% "all"){
    ag.df <- all_gradient(metrics.df, index.metrics)
    long.df$`ALL_5%` <- ag.df$`ALL_5%`
    long.df$`ALL_95%` <- ag.df$`ALL_95%`
    long.df$ALL_GRADIENT <- ag.df$ALL_GRADIENT
    long.df <- ref_gradient(long.df)
    if(!is.null(long.df$THRESHOLD)) long.df <- ref_threshold(long.df)
    if(!is.null(long.df$THRESHOLD)) long.df <- ref_threshold_mod(long.df)
    long.df <- ref_categorical(long.df)
    long.df
  }

  return(long.df)
}

#==============================================================================
#'Scoring Method: Reference Categorical
#'
#'@param long.df <- a data frame of samples to be scored.
#'@return Scores the samples.
#'@export
ref_categorical <- function(long.df){
  long.df$REF_CATEGORICAL <- ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                      long.df$REPORTING_VALUE <= long.df$`25%`, 1,
                                    ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                             long.df$REPORTING_VALUE > long.df$`25%` &
                                             long.df$REPORTING_VALUE < long.df$`50%` , 3,
                                           ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                                    long.df$REPORTING_VALUE >= long.df$`50%`, 5,
                                                  ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                           long.df$REPORTING_VALUE >= long.df$`75%`, 1,
                                                         ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                                  long.df$REPORTING_VALUE < long.df$`75%` &
                                                                  long.df$REPORTING_VALUE > long.df$`50%` , 3,
                                                                ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                                         long.df$REPORTING_VALUE <= long.df$`50%`, 5, 100000))))))
  long.df$REF_CATEGORICAL <- round(long.df$REF_CATEGORICAL, digits = 2)
  return(long.df)
}

#==============================================================================
#'Scoring Method: Reference Threshold
#'
#'@param long.df <- a data frame of samples to be scored.
#'@return Scores the samples.
#'@export

ref_threshold <- function(long.df){
  long.df$REF_THRESHOLD <- ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                    long.df$REPORTING_VALUE <= long.df$THRESHOLD, 0,
                                  ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                           long.df$REPORTING_VALUE > long.df$THRESHOLD &
                                           long.df$REPORTING_VALUE < long.df$`50%` ,
                                         ((long.df$REPORTING_VALUE - long.df$THRESHOLD) /
                                            (long.df$`50%` - long.df$THRESHOLD)) * 100,
                                         ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                                  long.df$REPORTING_VALUE >= long.df$`50%`, 100,
                                                ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                         long.df$REPORTING_VALUE >= long.df$THRESHOLD, 100,
                                                       ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                                long.df$REPORTING_VALUE < long.df$THRESHOLD &
                                                                long.df$REPORTING_VALUE > long.df$`50%` ,
                                                              ((long.df$THRESHOLD - long.df$REPORTING_VALUE) /
                                                                 (long.df$THRESHOLD - long.df$`50%`)) * 100,
                                                              ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                                       long.df$REPORTING_VALUE <= long.df$`50%`, 0, 100000))))))
  long.df$REF_THRESHOLD <- round(long.df$REF_THRESHOLD, digits = 2)
  return(long.df)
}

#==============================================================================
#'Scoring Method: Reference Threshold Modified
#'
#'@param long.df <- a data frame of samples to be scored.
#'@return Scores the samples.
#'@export

ref_threshold_mod <- function(long.df){
  long.df$LOW_THRESH <- long.df$THRESHOLD - abs(long.df$'50%' - long.df$THRESHOLD)
  long.df$REF_THRESHOLD_MOD <- ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                    long.df$REPORTING_VALUE <= long.df$LOW_THRESH, 0,
                                  ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                           long.df$REPORTING_VALUE > long.df$LOW_THRESH &
                                           long.df$REPORTING_VALUE < long.df$`50%` ,
                                         ((long.df$REPORTING_VALUE - long.df$LOW_THRESH) /
                                            (long.df$`50%` - long.df$LOW_THRESH)) * 100,
                                         ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                                  long.df$REPORTING_VALUE >= long.df$`50%`, 100,
                                                ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                         long.df$REPORTING_VALUE >= long.df$LOW_THRESH, 100,
                                                       ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                                long.df$REPORTING_VALUE < long.df$LOW_THRESH &
                                                                long.df$REPORTING_VALUE > long.df$`50%` ,
                                                              ((long.df$LOW_THRESH - long.df$REPORTING_VALUE) /
                                                                 (long.df$LOW_THRESH - long.df$`50%`)) * 100,
                                                              ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                                       long.df$REPORTING_VALUE <= long.df$`50%`, 0, 100000))))))
  long.df$REF_THRESHOLD <- round(long.df$REF_THRESHOLD, digits = 2)
  return(long.df)
}
#==============================================================================
#'Scoring Method: Reference Gradient
#'
#'@param long.df <- a data frame of samples to be scored.
#'@return Scores the samples.
#'@export

ref_gradient <- function(long.df){
  long.df$REF_GRADIENT <- ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                   long.df$REPORTING_VALUE <= long.df$`5%`, 0,
                                 ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                          long.df$REPORTING_VALUE > long.df$`5%` &
                                          long.df$REPORTING_VALUE < long.df$`95%` ,
                                        ((long.df$REPORTING_VALUE - long.df$`5%`) /
                                           (long.df$`95%` - long.df$`5%`)) * 100,
                                        ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                                 long.df$REPORTING_VALUE >= long.df$`95%`, 100,
                                               ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                        long.df$REPORTING_VALUE <= long.df$`5%`, 100,
                                                      ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                               long.df$REPORTING_VALUE > long.df$`5%` &
                                                               long.df$REPORTING_VALUE < long.df$`95%` ,
                                                             ((long.df$`95%` - long.df$REPORTING_VALUE) /
                                                                (long.df$`95%` - long.df$`5%`)) * 100,
                                                             ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                                      long.df$REPORTING_VALUE >= long.df$`95%`, 0, 100000))))))
  long.df$REF_GRADIENT <- round(long.df$REF_GRADIENT, digits = 2)
  return(long.df)
}

#==============================================================================
#'Scoring Method: All Site Classes Gradient
#'
#'@param metrics.df = data frame of metric values for each station
#'@return Scores the samples.
#'@export

all_gradient <- function(metrics.df, index.metrics){
  ID <- c("EVENT_ID", "STATION_ID","DATE", "AGENCY_CODE",
          "SAMPLE_NUMBER", "CATEGORY")
  all.quant <- data.frame(t(sapply(metrics.df[, 7:ncol(metrics.df)], quantile,
                                   c(0.05, 0.95))))
  names(all.quant) <- c("ALL_5%", "ALL_95%")
  all.quant$METRICS <- row.names(all.quant)
  score.info <- merge(index.metrics, all.quant, by = "METRICS")
  long.df <- reshape2::melt(metrics.df,
                            id.vars = c("EVENT_ID", "STATION_ID",
                                        "DATE", "AGENCY_CODE",
                                        "SAMPLE_NUMBER", "CATEGORY"))
  names(long.df) <- c(ID, "METRICS", "REPORTING_VALUE")
  long.df <- merge(long.df, score.info, by = "METRICS")

  long.df$ALL_GRADIENT <- ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                   long.df$REPORTING_VALUE <= long.df$`ALL_5%`, 0,
                                 ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                          long.df$REPORTING_VALUE > long.df$`ALL_5%` &
                                          long.df$REPORTING_VALUE < long.df$`ALL_95%` ,
                                        ((long.df$REPORTING_VALUE - long.df$`ALL_5%`) /
                                           (long.df$`ALL_95%` - long.df$`ALL_5%`)) * 100,
                                        ifelse(long.df$DISTURBANCE %in% "DECREASE" &
                                                 long.df$REPORTING_VALUE >= long.df$`ALL_95%`, 100,
                                               ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                        long.df$REPORTING_VALUE <= long.df$`ALL_5%`, 100,
                                                      ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                               long.df$REPORTING_VALUE > long.df$`ALL_5%` &
                                                               long.df$REPORTING_VALUE < long.df$`ALL_95%` ,
                                                             ((long.df$`ALL_95%` - long.df$REPORTING_VALUE) /
                                                                (long.df$`ALL_95%` - long.df$`ALL_5%`)) * 100,
                                                             ifelse(long.df$DISTURBANCE %in% "INCREASE" &
                                                                      long.df$REPORTING_VALUE >= long.df$`ALL_95%`, 0, 100000))))))
  long.df$ALL_GRADIENT <- round(long.df$ALL_GRADIENT, digits = 2)
  return(long.df)
}

#==============================================================================
#'Select the best metrics available
#'
#'@param Minus_Redund.df = data frame after redundancy analysis
#'@param best_type =
#'@param more_metrics =
#'@return Determines the best metric to represent each metric type
#'(i.e. Diversity, Tolerance, Composition, FFG, Habit) and includes any additional
#'metrics if the metrics discrimination efficiency value is >= 75.
#'@export


best_metrics <- function(Minus_Redund.df, best_type = 55, more_metrics = 60){
  Diversity <- c("RICH", "SHANNON", "SIMPSONS", "HURLBERTS_PIE", "MARGALEFS", "MENHINICKS",
                 "PCT_DOM1", "PCT_DOM2", "PCT_DOM3", "ABUNDANCE", "PCT_EPT_TAXA_RICH",
                 "PCT_EPT_TAXA_RICH_100", "EPT_CHIRO_100", "SCRAPE_FILTER_100",
                 "RICH_100", "SHANNON_100", "SIMPSONS_100", "HURLBERTS_PIE_100",
                 "MARGALEFS_100", "MENHINICKS_100",
                 "PCT_DOM1_100", "PCT_DOM2_100", "PCT_DOM3_100", "ABUNDANCE_100")
  Tolerance <- c("ASPT_MOD", "FBI", "PCT_URBAN_INTOL", "PCT_SENSITIVE",
                 "PCT_MOD_TOL", "PCT_TOLERANT", "BECKS")
  Composition <- c("PCT_EPT", "PCT_EPHEMEROPTERA", "PCT_PLECOPTERA",
                   "PCT_TRICHOPTERA", "PCT_NET_CADDISFLY", "PCT_TRICHOPTERA_NO_TOL",
                   "PCT_COLEOPTERA", "PCT_ODONATA", "PCT_COTE", "PCT_POTEC",
                   "PCT_HYDRO_TRICHOPTERA", "PCT_HYDRO_EPT", "PCT_CORBICULIDAE",
                   "PCT_DIPTERA", "PCT_SIMULIIDAE", "PCT_CHIRONOMIDAE", "PCT_NCO",
                   "GOLD", "PCT_NON_INSECT", "PCT_OLIGOCHAETA", "PCT_LIMESTONE",
                   "EPT_CHIRO", "PCT_RETREAT_CADDISFLY")
  FFG <- c("PCT_COLLECT", "PCT_GATHER", "PCT_FILTER", "PCT_PREDATOR",
           "PCT_SCRAPE", "PCT_SHRED", "SCRAPE_FILTER")
  Habit <- c("PCT_BURROW", "PCT_CLIMB", "PCT_CLING", "PCT_SWIM",
             "PCT_SKATE", "PCT_SPRAWL")


  Minus_Redund.df$TYPE <- ifelse(Minus_Redund.df$METRICS %in% Diversity, "DIVERSITY",
                          ifelse(Minus_Redund.df$METRICS %in% Tolerance, "TOLERANCE",
                          ifelse(Minus_Redund.df$METRICS %in% Composition, "COMPOSTION",
                          ifelse(Minus_Redund.df$METRICS %in% FFG, "FFG",
                          ifelse(Minus_Redund.df$METRICS %in% Habit, "HABIT", "NO MATCH")))))

  best.type <- plyr::ddply(Minus_Redund.df, c("TYPE"), function(x) x[which.max(x$SENSITIVITY),])
  bt <- best.type[best.type$SENSITIVITY >=  best_type, ]
  bt <- best.type[best.type$SENSITIVITY >= best_type, ]
  best_metrics <- Minus_Redund.df[Minus_Redund.df$SENSITIVITY >= more_metrics, ]
  final_metrics <- unique(rbind(bt, best_metrics))
  return(final_metrics)
}

#==============================================================================
#'Assigns Index Rating
#'
#'@param Scores = Scores for each metric in the index.
#'@param Metrics = Raw metric scores.  Includes only metrics within the index.
#'@return Prepares a the raw metric values, the metric score, and the sampling
#'event's rating (i.e. "Very Poor", "Poor", "Fair", "Good", or "Excellent).
#'@export
prep_score <- function(Scores, Metrics){
  Scores[, 6:ncol(Scores)] <- sapply(Scores[, 6:ncol(Scores)], function(x) as.numeric(as.character(x)))
  Scores$Mean <- rowMeans(Scores[, 6:ncol(Scores)])
  
  

  Scores$Rating <- ifelse (Scores$Mean < 17, "VeryPoor",
                           ifelse (Scores$Mean < 30, "Poor",
                                   ifelse (Scores$Mean < 50, "Fair",
                                           ifelse (Scores$Mean < 67, "Good",
                                                   ifelse (Scores$Mean >= 67 & Scores$Mean <= 100,
                                                           "Excellent", "ERROR")))))
  Scores <- with(Scores, Scores[order(EVENT_ID), ])

  melted.score <- (reshape2::melt(Scores,id.vars = c("EVENT_ID", "STATION_ID",
                                                     "DATE", "SAMPLE_NUMBER",
                                                     "AGENCY_CODE",
                                                     "Mean", "Rating")))
  colnames(melted.score) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                              "AGENCY_CODE", "MEAN", "RATING",
                              "METRIC", "METRIC_SCORE")
  melted.score <- with(melted.score, melted.score[order(EVENT_ID), ])

  melted.metrics <- (reshape2::melt(Metrics, id.vars = c("EVENT_ID","STATION_ID",
                                                         "DATE", "SAMPLE_NUMBER",
                                                         "AGENCY_CODE")))
  colnames(melted.metrics) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                "AGENCY_CODE", "METRIC","METRIC_VALUE")
  melted.metrics <- with(melted.metrics, melted.metrics[order(EVENT_ID), ])

  new.data <- merge(melted.metrics, melted.score,
                    by = c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                           "AGENCY_CODE", "METRIC"))
  return(new.data)
}

#==============================================================================
#'BIBI Rating
#'
#'@param metrics.df = data frame with each metric scored based on
#'@param index = 2011 or 2016 rating procedure.
#'@return The index rating for each site. The rating is based on
#'the average score of the metrics in the index.
#'@export

ibi_rating <- function(Scored.df, index = 2011){
  Mean_Score <- aggregate(Scored.df$METRIC_SCORE ~ Scored.df$EVENT_ID +
                            Scored.df$STATION_ID, FUN = mean)

  colnames(Mean_Score) <- c("EVENT_ID", "STATION_ID", "MEAN")
  if(index %in% 2011){
    Mean_Score$RATING <- ifelse (Mean_Score$MEAN < 17, "VeryPoor",
                                 ifelse (Mean_Score$MEAN < 30, "Poor",
                                         ifelse (Mean_Score$MEAN < 50, "Fair",
                                                 ifelse (Mean_Score$MEAN < 67, "Good",
                                                         ifelse (Mean_Score$MEAN >= 67 & Mean_Score$MEAN <= 100,
                                                                 "Excellent", "ERROR")))))
  }else{
    if(index %in% 2016){
      Mean_Score$RATING <- ifelse (Mean_Score$MEAN < 20, "VeryPoor",
                                   ifelse (Mean_Score$MEAN < 40, "Poor",
                                           ifelse (Mean_Score$MEAN < 61, "Fair",
                                                   ifelse (Mean_Score$MEAN < 80, "Good",
                                                           ifelse (Mean_Score$MEAN >= 80 & Mean_Score$MEAN <= 100,
                                                                   "Excellent", "ERROR")))))
    }
  }

  
  mean.score <- with(Mean_Score, Mean_Score[order(EVENT_ID), ])

  metric.scores <- Scored.df[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                 "AGENCY_CODE", "CATEGORY", "VALUE", "METRICS",
                                 "METRIC_SCORE", "THRESHOLD", "MEDIAN")]
  melted.score <- (reshape2::melt(metric.scores,
                                  id.vars = c("EVENT_ID", "STATION_ID", "DATE",
                                              "SAMPLE_NUMBER", "AGENCY_CODE",
                                              "CATEGORY", "VALUE", "METRICS",
                                              "THRESHOLD", "MEDIAN")))
  melted.score <- melted.score[ , -which(names(melted.score) %in% c("variable"))]
  names(melted.score)[names(melted.score) == "value"] <-  "METRIC_SCORE"
  melted.score <- with(melted.score, melted.score[order(EVENT_ID), ])

  new.df <- merge(melted.score, mean.score, by = c("EVENT_ID", "STATION_ID"))
  final.df <- subset(new.df, select = c(EVENT_ID:CATEGORY, METRICS, THRESHOLD,
                                        MEDIAN, VALUE, METRIC_SCORE, MEAN,
                                        RATING))
  return(final.df)
}

#==============================================================================
#'Rate Samples
#'
#'@param Metrics.df = data frame with all of the raw metric scores.
#'@param Final_Metrics.df = a list of metrics selected/identified as
#'sensitive to a disturbance gradient.
#'@return The index rating for each sampling event. The metrics are scored on
#'a scale of 0-100 based on descrimination efficiency values between each of
#'the designated site classes.  The mean score of all of the metrics for each
#'sampling event is used to rate sampling station.
#'@export

rate_samples <- function(Metrics.df, Final_Metrics.df){
  melted <- reshape2::melt(Metrics.df, id.vars = c("EVENT_ID", "STATION_ID", "DATE",
                                                   "SAMPLE_NUMBER", "AGENCY_CODE", "CATEGORY"))
  colnames(melted) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                        "AGENCY_CODE","CATEGORY", "METRICS", "VALUE")
  merged.df <- merge(melted, Final_Metrics.df, by = "METRICS", all.y = TRUE)

  merged.df$METRIC_SCORE <- ifelse(merged.df$VALUE <= merged.df$THRESH_REF_NEAR &
                                     merged.df$DISTURBANCE == "INCREASE", 100,
                                   ifelse(merged.df$VALUE > merged.df$THRESH_REF_NEAR &
                                            merged.df$VALUE <= merged.df$THRESH_NEAR_MOD &
                                            merged.df$DISTURBANCE == "INCREASE", 66,
                                          ifelse(merged.df$VALUE > merged.df$THRESH_NEAR_MOD &
                                                   merged.df$VALUE <= merged.df$THRESH_MOD_SEV &
                                                   merged.df$DISTURBANCE == "INCREASE", 33,
                                                 ifelse(merged.df$VALUE > merged.df$THRESH_MOD_SEV &
                                                          merged.df$DISTURBANCE == "INCREASE", 0,
                                                        ifelse(merged.df$VALUE >= merged.df$THRESH_REF_NEAR &
                                                                 merged.df$DISTURBANCE == "DECREASE", 100,
                                                               ifelse(merged.df$VALUE < merged.df$THRESH_REF_NEAR &
                                                                        merged.df$VALUE >= merged.df$THRESH_NEAR_MOD &
                                                                        merged.df$DISTURBANCE == "DECREASE",66,
                                                                      ifelse(merged.df$VALUE < merged.df$THRESH_NEAR_MOD &
                                                                               merged.df$VALUE >= merged.df$THRESH_MOD_SEV &
                                                                               merged.df$DISTURBANCE == "DECREASE", 33,
                                                                             ifelse(merged.df$VALUE < merged.df$THRESH_MOD_SEV &
                                                                                      merged.df$DISTURBANCE == "DECREASE", 0, 10000))))))))



  merged.df$METRIC_SCORE <- round(as.numeric(as.character(merged.df$METRIC_SCORE)), digits = 2)


  Scored.df <- merged.df

  Mean_Score <- aggregate(Scored.df$METRIC_SCORE ~ Scored.df$EVENT_ID +
                            Scored.df$STATION_ID, FUN = mean)

  colnames(Mean_Score) <- c("EVENT_ID", "STATION_ID", "MEAN")

  Mean_Score$RATING <- ifelse (Mean_Score$MEAN <= 25, "VeryPoor",
                               ifelse (Mean_Score$MEAN < 25, "Poor",
                                       ifelse (Mean_Score$MEAN < 50, "Fair",
                                               ifelse (Mean_Score$MEAN < 75, "Good",
                                                       ifelse (Mean_Score$MEAN >= 75 & Mean_Score$MEAN <= 100,
                                                               "Excellent", "ERROR")))))
  mean.score <- with(Mean_Score, Mean_Score[order(EVENT_ID), ])

  metric.scores <- Scored.df[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                 "AGENCY_CODE", "CATEGORY", "VALUE", "METRICS",
                                 "METRIC_SCORE")]
  melted.score <- (reshape2::melt(metric.scores,
                                  id.vars = c("EVENT_ID", "STATION_ID", "DATE",
                                              "SAMPLE_NUMBER", "AGENCY_CODE",
                                              "CATEGORY", "VALUE", "METRICS")))
  melted.score <- melted.score[ , -which(names(melted.score) %in% c("variable"))]
  names(melted.score)[names(melted.score) == "value"] <-  "METRIC_SCORE"
  melted.score <- with(melted.score, melted.score[order(EVENT_ID), ])

  new.df <- merge(melted.score, mean.score, by = c("EVENT_ID", "STATION_ID"))
  final.df <- subset(new.df, select = c(EVENT_ID:CATEGORY, METRICS,
                                        VALUE, METRIC_SCORE, MEAN,
                                        RATING))
  return(final.df)
}

#==============================================================================
#'Rate Samples 2
#'
#'@param Metrics.df = data frame with all of the raw metric scores.
#'@param Final_Metrics.df = a list of metrics selected/identified as
#'sensitive to a disturbance gradient.
#'@return The index rating for each sampling event. The metrics are scored on
#'a scale of 0-100 based on descrimination efficiency values between each of
#'the designated site classes.  The mean score of all of the metrics for each
#'sampling event is used to rate sampling station.
#'@export

rate_samples2 <- function(Metrics.df, Final_Metrics.df){
  melted <- reshape2::melt(Metrics.df, id.vars = c("EVENT_ID", "STATION_ID", "DATE",
                                                   "SAMPLE_NUMBER", "AGENCY_CODE", "CATEGORY"))
  colnames(melted) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                        "AGENCY_CODE","CATEGORY", "METRICS", "VALUE")
  merged.df <- merge(melted, Final_Metrics.df, by = "METRICS", all.y = TRUE)

  merged.df$METRIC_SCORE <- ifelse(merged.df$VALUE <= merged.df$THRESH_REF_NEAR &
                                     merged.df$DISTURBANCE == "INCREASE", 100,
                                   ifelse(merged.df$VALUE > merged.df$THRESH_REF_NEAR &
                                            merged.df$VALUE <= merged.df$THRESH_NEAR_MIN &
                                            merged.df$DISTURBANCE == "INCREASE", 100,
                                          ifelse(merged.df$VALUE > merged.df$THRESH_NEAR_MIN &
                                                   merged.df$VALUE <= merged.df$THRESH_MIN_MOD &
                                                   merged.df$DISTURBANCE == "INCREASE", 50,
                                                 ifelse(merged.df$VALUE > merged.df$THRESH_MIN_MOD &
                                                          merged.df$VALUE <= merged.df$THRESH_MOD_SEV &
                                                          merged.df$DISTURBANCE == "INCREASE", 0,
                                                 ifelse(merged.df$VALUE > merged.df$THRESH_MOD_SEV &
                                                          merged.df$DISTURBANCE == "INCREASE", 0,
                                                        ifelse(merged.df$VALUE >= merged.df$THRESH_REF_NEAR &
                                                                 merged.df$DISTURBANCE == "DECREASE", 100,
                                                               ifelse(merged.df$VALUE < merged.df$THRESH_REF_NEAR &
                                                                        merged.df$VALUE >= merged.df$THRESH_NEAR_MIN &
                                                                        merged.df$DISTURBANCE == "DECREASE",100,
                                                                      ifelse(merged.df$VALUE < merged.df$THRESH_NEAR_MIN &
                                                                               merged.df$VALUE >= merged.df$THRESH_MIN_MOD &
                                                                               merged.df$DISTURBANCE == "DECREASE", 50,
                                                                             ifelse(merged.df$VALUE < merged.df$THRESH_MIN_MOD &
                                                                                      merged.df$VALUE >= merged.df$THRESH_MOD_SEV &
                                                                                      merged.df$DISTURBANCE == "DECREASE", 0,
                                                                             ifelse(merged.df$VALUE < merged.df$THRESH_MOD_SEV &
                                                                                      merged.df$DISTURBANCE == "DECREASE", 0, 10000))))))))))



  merged.df$METRIC_SCORE <- round(as.numeric(as.character(merged.df$METRIC_SCORE)), digits = 2)


  Scored.df <- merged.df

  Mean_Score <- aggregate(Scored.df$METRIC_SCORE ~ Scored.df$EVENT_ID +
                            Scored.df$STATION_ID, FUN = mean)

  colnames(Mean_Score) <- c("EVENT_ID", "STATION_ID", "MEAN")

  Mean_Score$RATING <- ifelse (Mean_Score$MEAN <= 25, "VeryPoor",
                               ifelse (Mean_Score$MEAN < 25, "Poor",
                                       ifelse (Mean_Score$MEAN < 50, "Fair",
                                               ifelse (Mean_Score$MEAN < 75, "Good",
                                                       ifelse (Mean_Score$MEAN >= 75 & Mean_Score$MEAN <= 100,
                                                               "Excellent", "ERROR")))))
  mean.score <- with(Mean_Score, Mean_Score[order(EVENT_ID), ])

  metric.scores <- Scored.df[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                 "AGENCY_CODE", "CATEGORY", "VALUE", "METRICS",
                                 "METRIC_SCORE")]
  melted.score <- (reshape2::melt(metric.scores,
                                  id.vars = c("EVENT_ID", "STATION_ID", "DATE",
                                              "SAMPLE_NUMBER", "AGENCY_CODE",
                                              "CATEGORY", "VALUE", "METRICS")))
  melted.score <- melted.score[ , -which(names(melted.score) %in% c("variable"))]
  names(melted.score)[names(melted.score) == "value"] <-  "METRIC_SCORE"
  melted.score <- with(melted.score, melted.score[order(EVENT_ID), ])

  new.df <- merge(melted.score, mean.score, by = c("EVENT_ID", "STATION_ID"))
  final.df <- subset(new.df, select = c(EVENT_ID:CATEGORY, METRICS,
                                        VALUE, METRIC_SCORE, MEAN,
                                        RATING))
  return(final.df)
}

#==============================================================================
#'Rate Samples 3
#'
#'@param Metrics.df = data frame with all of the raw metric scores.
#'@param Final_Metrics.df = a list of metrics selected/identified as
#'sensitive to a disturbance gradient.
#'@return The index rating for each sampling event. The metrics are scored on
#'a scale of 0-100 based on descrimination efficiency values between each of
#'the designated site classes.  The mean score of all of the metrics for each
#'sampling event is used to rate sampling station.
#'@export

rate_samples3 <- function(Metrics.df, Final_Metrics.df){
  melted <- reshape2::melt(Metrics.df, id.vars = c("EVENT_ID", "STATION_ID", "DATE",
                                                   "SAMPLE_NUMBER", "AGENCY_CODE", "CATEGORY"))
  colnames(melted) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                        "AGENCY_CODE","CATEGORY", "METRICS", "VALUE")
  merged.df <- merge(melted, Final_Metrics.df, by = "METRICS", all.y = TRUE)
  
  merged.df$METRIC_SCORE <- ifelse(merged.df$VALUE <= merged.df$THRESH_REF_MIN &
                                     merged.df$DISTURBANCE == "INCREASE", 100,
                                   ifelse(merged.df$VALUE > merged.df$THRESH_REF_MIN &
                                            merged.df$VALUE <= merged.df$THRESH_MIN_MOD &
                                            merged.df$DISTURBANCE == "INCREASE", 50,
                                          ifelse(merged.df$VALUE > merged.df$THRESH_MIN_MOD &
                                                   merged.df$VALUE <= merged.df$THRESH_MOD_SEV &
                                                   merged.df$DISTURBANCE == "INCREASE", 25,
                                                        ifelse(merged.df$VALUE > merged.df$THRESH_MOD_SEV &
                                                                 merged.df$DISTURBANCE == "INCREASE", 0,
                                                               ifelse(merged.df$VALUE >= merged.df$THRESH_REF_MIN &
                                                                        merged.df$DISTURBANCE == "DECREASE", 100,
                                                                      ifelse(merged.df$VALUE < merged.df$THRESH_REF_MIN &
                                                                               merged.df$VALUE >= merged.df$THRESH_MIN_MOD &
                                                                               merged.df$DISTURBANCE == "DECREASE", 75,
                                                                             ifelse(merged.df$VALUE < merged.df$THRESH_MIN_MOD &
                                                                                      merged.df$VALUE >= merged.df$THRESH_MOD_SEV &
                                                                                      merged.df$DISTURBANCE == "DECREASE", 25,
                                                                                           ifelse(merged.df$VALUE < merged.df$THRESH_MOD_SEV &
                                                                                                    merged.df$DISTURBANCE == "DECREASE", 0, 10000))))))))
  
  
  
  merged.df$METRIC_SCORE <- round(as.numeric(as.character(merged.df$METRIC_SCORE)), digits = 2)
  
  
  Scored.df <- merged.df
  
  Mean_Score <- aggregate(Scored.df$METRIC_SCORE ~ Scored.df$EVENT_ID +
                            Scored.df$STATION_ID, FUN = mean)
  
  colnames(Mean_Score) <- c("EVENT_ID", "STATION_ID", "MEAN")
  
  Mean_Score$RATING <- ifelse (Mean_Score$MEAN <= 60, "VeryPoor",
                               ifelse (Mean_Score$MEAN < 60, "Poor",
                                       ifelse (Mean_Score$MEAN < 70, "Fair",
                                               ifelse (Mean_Score$MEAN < 80, "Good",
                                                       ifelse (Mean_Score$MEAN >= 80 & Mean_Score$MEAN <= 100,
                                                               "Excellent", "ERROR")))))

                                                               
  mean.score <- with(Mean_Score, Mean_Score[order(EVENT_ID), ])
  
  metric.scores <- Scored.df[, c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                 "AGENCY_CODE", "CATEGORY", "VALUE", "METRICS",
                                 "METRIC_SCORE")]
  melted.score <- (reshape2::melt(metric.scores,
                                  id.vars = c("EVENT_ID", "STATION_ID", "DATE",
                                              "SAMPLE_NUMBER", "AGENCY_CODE",
                                              "CATEGORY", "VALUE", "METRICS")))
  melted.score <- melted.score[ , -which(names(melted.score) %in% c("variable"))]
  names(melted.score)[names(melted.score) == "value"] <-  "METRIC_SCORE"
  melted.score <- with(melted.score, melted.score[order(EVENT_ID), ])
  
  new.df <- merge(melted.score, mean.score, by = c("EVENT_ID", "STATION_ID"))
  final.df <- subset(new.df, select = c(EVENT_ID:CATEGORY, METRICS,
                                        VALUE, METRIC_SCORE, MEAN,
                                        RATING))
  return(final.df)
}

#==============================================================================
#'Aggregate Scores
#'
#'@param Data frame containing multiple scoring procedures.
#'@return Summarize Reference Gradient, Reference Categorical, and 
#'All Gradient scoring procedures.
#'@export

agg_score <- function(scores.df){
  Am6rg <- aggregate(REF_GRADIENT ~  EVENT_ID + STATION_ID + DATE + 
                       AGENCY_CODE + SAMPLE_NUMBER + CATEGORY, data = scores.df, FUN = mean)
  Am6rc <- aggregate(REF_CATEGORICAL ~  EVENT_ID + STATION_ID + DATE + 
                       AGENCY_CODE + SAMPLE_NUMBER + CATEGORY, data = scores.df, FUN = mean)
  Am6rtm <- aggregate(REF_THRESHOLD_MOD ~  EVENT_ID + STATION_ID + DATE + 
                       AGENCY_CODE + SAMPLE_NUMBER + CATEGORY, data = scores.df, FUN = mean)
  Am6 <- aggregate(ALL_GRADIENT ~  EVENT_ID + STATION_ID + DATE + 
                     AGENCY_CODE + SAMPLE_NUMBER + CATEGORY, data = scores.df, FUN = mean)
  Am6$REF_GRADIENT <- Am6rg$REF_GRADIENT
  Am6$REF_CATEGORICAL <- Am6rc$REF_CATEGORICAL
  Am6$REF_THRESHOLD_MOD <- Am6rtm$REF_THRESHOLD_MOD
  return(Am6)
}


#==============================================================================
#'Build a new IBI
#'
#'@param Info = taxonomic information
#'@param Long = taxonomic data in a long data format
#'@param Level = the lowest taxonomic level that metrics should be calculated
#'@param Category = each sampling event assigned a site category (i.e. Reference,
#'Mixed, or Degraded)
#'@param all.metrics = the function that should be used to calculate
#'all of the metrics
#'@return Build an IBI from scratch.
#'@export

build_ibi <- function(Info, Long, Level = "FAMILY", Category, all.metrics = all_metrics_div_rare){
  calc.metrics <- all.metrics(Info, Long, Level)
  metrics.category <- merge(Category, calc.metrics, by = "EVENT_ID")
  de.thresh <- de_thresh2(metrics.category)
  redund.metrics <- redundancy3(metrics.category, de.thresh)
  best.metrics <- best_metrics(redund.metrics)
  score.metrics <- ibi_score(metrics.category, best.metrics)
  rate.events <- ibi_rating(score.metrics)
  return(rate.events)
}
