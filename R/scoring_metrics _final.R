#==============================================================================
#
# Scoring Metrics
#
#==============================================================================
#'Score Metrics
#'
#'@param metrics.df = a data frame containing raw metric scores.
#'@param metrics.summary = a data frame of the output from the function 
#'metric_summary.
#'@param bioregion = the bioregion of interest.
#'@param scoring.method = the scoring method of interest.
#'@return Update
#'@export

score_metrics <- function(metrics.df, metrics.summary, bioregion, scoring.method, bound.limits = TRUE){
  if(!is.data.frame((metrics.summary))){
    metrics.summary <- unique(metrics_summary2(metrics.df, bioregion))
  }
  
  if(bound.limits == TRUE){
    metrics.summary$BOUND <- ifelse(metrics.summary$BOUND_BI_CMA < 0, 0,
                          ifelse(metrics.summary$BOUND_BI_CMA > 100, 100,
                                 as.numeric(as.character(metrics.summary$BOUND_BI_CMA))))
  }
  
  if(scoring.method %in% c("ALL", "ALL_GRADIENT", "TIGHT_GRADIENT", "UNEVEN_GRADIENT")){

    #============================================================================
    ref_deg.df <- metrics.df[metrics.df$CATEGORY %in% c("REF", "SEV"), ]
    #============================================================================
    quantile.df <- data.frame(METRICS = names(metrics.df[, 8:ncol(metrics.df)]))
    quantile.df$"5th" <- apply(ref_deg.df[, 8:ncol(ref_deg.df)], 2, function(x) quantile(x, 0.05))
    quantile.df$"95th" <- apply(ref_deg.df[, 8:ncol(ref_deg.df)], 2, function(x) quantile(x, 0.95))
    #============================================================================
    merge.df <- merge(metrics.summary, quantile.df, by = "METRICS")
    #merge.df <- merge.df[, c("METRICS", "DISTURBANCE", "MIDPOINT_BI_CMA", "5th", "95th")]
    #names(merge.df) <- c("METRICS", "DISTURBANCE", "BSP", "5th", "95th")
    #============================================================================
    long.df <- tidyr::gather(metrics.df, "METRICS", "REPORTING_VALUE", 8:ncol(metrics.df))
    all.long.df <- merge(long.df, merge.df, by = "METRICS")
    #all.long.df <- tidyr::gather(metrics.df, METRICS, REPORTING_VALUE, 8:ncol(metrics.df))
    #all.long.df <- merge(all.long.df, score.info, by = "METRICS")
  }

  #============================================================================
  if(scoring.method %in% c("ALL", "REF_CATEGORICAL", "REF_GRADIENT")){
    ref.df <- metrics.df[metrics.df$CATEGORY %in% "REF", ]
    ref.quant <- data.frame(t(sapply(ref.df[, 8:ncol(ref.df)], quantile,
                                     c(0.05, 0.25, 0.50, 0.75, 0.95))))
    names(ref.quant) <- c("5%", "25%", "50%", "75%", "95%")
    ref.quant$METRICS <- row.names(ref.quant)
    
    #score.info <- merge(index.metrics, ref.quant, by = "METRICS")
    ref.long.df <- reshape2::melt(metrics.df,
                                  id.vars = c("EVENT_ID", "STATION_ID",
                                              "DATE", "AGENCY_CODE",
                                              "SAMPLE_NUMBER", "CATEGORY"))
    ref.long.df <- tidyr::gather(metrics.df, METRICS, REPORTING_VALUE, 8:ncol(metrics.df))
    #names(ref.long.df) <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE",
    #                    "SAMPLE_NUMBER", "CATEGORY", "METRICS", "REPORTING_VALUE")
    #ref.long.df <- merge(ref.long.df, score.info, by = "METRICS")
  }
  
  #============================================================================
  # UNEVEN GRADIENT SCORING
  if(scoring.method  %in% c("ALL", "UNEVEN_GRADIENT")){
    final.df <- uneven_gradient(all.long.df)
  }
  #============================================================================
  # TIGHT GRADIENT SCORING
  if(scoring.method  %in% c("ALL", "TIGHT_GRADIENT")){
    final.df <- tight_gradient(all.long.df)
  }
  #============================================================================
  # ALL GRADIENT SCORING
  if(scoring.method  %in% c("ALL", "ALL_GRADIENT")){
    final.df <- all_gradient2(all.long.df)
  }
  #============================================================================
  #if(scoring.method %in% "ref_categorical") long.df <- ref_categorical(long.df)
  #if(scoring.method %in% "ref_threshold" & !(is.null(long.df$THRESHOLD))) long.df <- ref_threshold(long.df)
  #if(scoring.method %in% "ref_threshold_mod" & !(is.null(long.df$THRESHOLD))) long.df <- ref_threshold_mod(long.df)
  #if(scoring.method %in% "ref_gradient") long.df <- ref_gradient(long.df)
  #============================================================================
  # REF CATEGORICAL
  if(scoring.method  %in% c("ALL", "REF_CATEGORICAL")){
    final.df <- ref_categorical(all.long.df)
  }
  #============================================================================
  return(final.df)
  
}

#==============================================================================
#'Uneven Gradient Scoring
#'
#'@param metrics.df = a data frame containing raw metric scores.
#'@param metrics.summary = a data frame of the output from the function 
#'metric_summary.
#'@param bioregion = the bioregion of interest.
#'@return Update
#'@export

uneven_gradient <- function(metrics.df){
 
  metrics.df$SCORE <- ifelse(metrics.df$REPORTING_VALUE == metrics.df$MIDPOINT_BI_CMA, 50,
                             ifelse(metrics.df$DISTURBANCE %in% "EQUAL", 50,
                            ifelse(metrics.df$REPORTING_VALUE <= metrics.df$`5th` & metrics.df$DISTURBANCE %in% "DECREASE", 0,
                                   ifelse(metrics.df$REPORTING_VALUE <= metrics.df$`5th` & metrics.df$DISTURBANCE %in% "INCREASE", 100,       
                                          ifelse(metrics.df$REPORTING_VALUE >= metrics.df$`95th` & metrics.df$DISTURBANCE %in% "DECREASE", 100,
                                                 ifelse(metrics.df$REPORTING_VALUE >= metrics.df$`95th` & metrics.df$DISTURBANCE %in% "INCREASE", 0,
                                                        ifelse(metrics.df$REPORTING_VALUE < metrics.df$MIDPOINT_BI_CMA & metrics.df$REPORTING_VALUE > metrics.df$`5th` & metrics.df$DISTURBANCE %in% "DECREASE",
                                                               ((metrics.df$REPORTING_VALUE - metrics.df$`5th`) / (metrics.df$MIDPOINT_BI_CMA - metrics.df$`5th`)) * 50,
                                                               ifelse(metrics.df$REPORTING_VALUE < metrics.df$MIDPOINT_BI_CMA & metrics.df$REPORTING_VALUE > metrics.df$`5th` & metrics.df$DISTURBANCE %in% "INCREASE",
                                                                      (((metrics.df$MIDPOINT_BI_CMA - metrics.df$REPORTING_VALUE) / (metrics.df$MIDPOINT_BI_CMA - metrics.df$`5th`)) * 50) + 50,
                                                                      ifelse(metrics.df$REPORTING_VALUE > metrics.df$MIDPOINT_BI_CMA & metrics.df$REPORTING_VALUE < metrics.df$`95th` & metrics.df$DISTURBANCE %in% "DECREASE",
                                                                             (((metrics.df$REPORTING_VALUE - metrics.df$MIDPOINT_BI_CMA) / (metrics.df$`95th` - metrics.df$MIDPOINT_BI_CMA)) * 50) + 50,
                                                                             ifelse(metrics.df$REPORTING_VALUE > metrics.df$MIDPOINT_BI_CMA & metrics.df$REPORTING_VALUE < metrics.df$`95th` & metrics.df$DISTURBANCE %in% "INCREASE",
                                                                                    ((metrics.df$`95th` - metrics.df$REPORTING_VALUE) / (metrics.df$`95th` - metrics.df$MIDPOINT_BI_CMA)) * 50, 10000))))))))))
  metrics.df$SCORE <- round(metrics.df$SCORE, 2)
  
  metrics.df <- metrics.df[, c("EVENT_ID", "CATEGORY", "STATION_ID", "SAMPLE_NUMBER",
                             "AGENCY_CODE", "DATE", "METRICS", "SCORE")]
  
  final.df <- tidyr::spread(metrics.df, METRICS, SCORE)
  return(final.df)
  
}


#==============================================================================
#'Tight Gradient Scoring
#'
#'@param metrics.df = a data frame containing raw metric values.
#'@return Update
#'@export

tight_gradient <- function(metrics.df){
  
  metrics.df$SCORE <- ifelse(metrics.df$DISTURBANCE %in% "EQUAL", 50,
                                    ifelse(metrics.df$REPORTING_VALUE < metrics.df$BOUND_BI_CMA & metrics.df$DISTURBANCE %in% "DECREASE", 0,
                                           ifelse(metrics.df$REPORTING_VALUE < metrics.df$REF_MEDIAN & metrics.df$DISTURBANCE %in% "INCREASE", 0,       
                                                  ifelse(metrics.df$REPORTING_VALUE >= metrics.df$REF_MEDIAN & metrics.df$DISTURBANCE %in% "DECREASE", 100,
                                                         ifelse(metrics.df$REPORTING_VALUE >= metrics.df$BOUND_BI_CMA & metrics.df$DISTURBANCE %in% "INCREASE", 100,
                                                                ifelse(metrics.df$REPORTING_VALUE < metrics.df$REF_MEDIAN & metrics.df$REPORTING_VALUE >= metrics.df$BOUND_BI_CMA & metrics.df$DISTURBANCE %in% "DECREASE",
                                                                       ((metrics.df$REPORTING_VALUE - metrics.df$BOUND_BI_CMA) / (metrics.df$REF_MEDIAN - metrics.df$BOUND_BI_CMA)) * 100,
                                                                       ifelse(metrics.df$REPORTING_VALUE < metrics.df$BOUND_BI_CMA & metrics.df$REPORTING_VALUE >= metrics.df$REF_MEDIAN & metrics.df$DISTURBANCE %in% "INCREASE",
                                                                              ((metrics.df$BOUND_BI_CMA - metrics.df$REPORTING_VALUE) / (metrics.df$BOUND_BI_CMA - metrics.df$REF_MEDIAN)) * 100,
                                                                              1000)))))))
  metrics.df$SCORE <- round(metrics.df$SCORE, 2)
  
  
  metrics.df <- metrics.df[, c("EVENT_ID", "CATEGORY", "STATION_ID", "SAMPLE_NUMBER",
                        "AGENCY_CODE", "DATE", "METRICS", "SCORE")]
  
  final.df <- tidyr::spread(metrics.df, METRICS, SCORE)
  return(final.df)
}





#==============================================================================
#'Scoring Method: All Gradient2
#'
#'@param metrics.df = data frame of metric values for each station
#'@return Scores the samples.
#'@export

all_gradient2 <- function(metrics.df){

  metrics.df$SCORE <- ifelse(metrics.df$DISTURBANCE %in% "EQUAL", 50,
                             ifelse(metrics.df$DISTURBANCE %in% "DECREASE" &
                                   metrics.df$REPORTING_VALUE <= metrics.df$`5th`, 0,
                                 ifelse(metrics.df$DISTURBANCE %in% "DECREASE" &
                                          metrics.df$REPORTING_VALUE > metrics.df$`5th` &
                                          metrics.df$REPORTING_VALUE < metrics.df$`95th` ,
                                        ((metrics.df$REPORTING_VALUE - metrics.df$`5th`) /
                                           (metrics.df$`95th` - metrics.df$`5th`)) * 100,
                                        ifelse(metrics.df$DISTURBANCE %in% "DECREASE" &
                                                 metrics.df$REPORTING_VALUE >= metrics.df$`95th`, 100,
                                               ifelse(metrics.df$DISTURBANCE %in% "INCREASE" &
                                                        metrics.df$REPORTING_VALUE <= metrics.df$`5th`, 100,
                                                      ifelse(metrics.df$DISTURBANCE %in% "INCREASE" &
                                                               metrics.df$REPORTING_VALUE > metrics.df$`5th` &
                                                               metrics.df$REPORTING_VALUE < metrics.df$`95th` ,
                                                             ((metrics.df$`95th` - metrics.df$REPORTING_VALUE) /
                                                                (metrics.df$`95th` - metrics.df$`5th`)) * 100,
                                                             ifelse(metrics.df$DISTURBANCE %in% "INCREASE" &
                                                                      metrics.df$REPORTING_VALUE >= metrics.df$`95th`, 0, 100000)))))))
  metrics.df$SCORE <- round(metrics.df$SCORE, digits = 2)
  metrics.df <- metrics.df[, c("EVENT_ID", "CATEGORY", "STATION_ID", "SAMPLE_NUMBER",
                               "AGENCY_CODE", "DATE", "METRICS", "SCORE")]
  
  final.df <- tidyr::spread(metrics.df, METRICS, SCORE)
  return(metrics.df)
}

#==============================================================================

#==============================================================================
#'Old Scoring 4
#'
#'@param metrics.df
#'@param bioregion
#'@return Scores the samples.
#'@export


old_scoring4 <- function(metrics.df, bioregion,
                         zero_null = TRUE, bound.lim = TRUE, metric.summary = NULL,
                         scoring.method = "UNEVEN_GRADIENT"){
  new.df <- metrics.df[metrics.df$BIOREGION %in% bioregion, ]
  #new.df <- new.df[new.df$ABUNDANCE >= 70, ]
  
  bi_de <- unique(chunk_sensitivity(new.df[, !names(new.df) %in% "BIOREGION"], "REF", "SEV", "DE"))
  names(bi_de) <- c("METRICS", "DISTURBANCE", "BINARY_de")
  
  ref.df <- new.df[new.df$CATEGORY %in% "REF", ]
  
  if(zero_null == TRUE){
    ref.zero <- ref.df
    ref.zero[, 8:ncol(ref.zero)] <- ifelse(ref.zero[, 8:ncol(ref.zero)] == 0, 1, 0)
    
    zero.cols <- names(ref.zero[, 8:ncol(ref.zero)])[(colSums(ref.zero[, 8:ncol(ref.zero)]) / nrow(ref.zero[, 8:ncol(ref.zero)])) >= 0.05]
    if(any(bi_de$DISTURBANCE %in% "EQUAL")){
      equal.cols <- unlist(list(bi_de[bi_de$DISTURBANCE %in% "EQUAL", "METRICS"]))
      zero.cols <- unlist(list(zero.cols, equal.cols))
    }
    zero.keep <- unlist(list(names(new.df[, 1:7]), zero.cols, "ABUNDANCE"))
    
    new.with.zero <- new.df[, !names(new.df) %in% zero.cols]
    new.without.zero <- new.df[, names(new.df) %in% zero.keep]
    new.without.zero[new.without.zero == 0] <- NA
  }
  
  
  
  
  if(zero_null == TRUE){
    if(ncol(new.with.zero) > 7 & ncol(new.without.zero) > 7){
      final.with <- score_metrics(new.with.zero, metric.summary, bioregion, scoring.method)
      final.without <- score_metrics(new.without.zero, metric.summary, bioregion, scoring.method)
      
      #test22 <- colSums(new.without.zero[, 8:ncol(new.without.zero)], na.rm = T) == 0
      final.df <- merge(final.with, final.without, by = c("EVENT_ID", "CATEGORY",
                                                          "STATION_ID", "SAMPLE_NUMBER",
                                                          "AGENCY_CODE", "DATE"))
    }
    
    if(ncol(new.with.zero) > 7 & ncol(new.without.zero) <= 7){
      final.df <- score_metrics(new.with.zero, metric.summary, bioregion, scoring.method)
    }
    
    if(ncol(new.with.zero) <= 7 & ncol(new.without.zero) > 7){
      final.df <- score_metrics(new.without.zero, metric.summary, bioregion, scoring.method)
    }
  }else{
    final.df <- score_metrics(new.df, metric.summary, bioregion, scoring.method)
  }
  
  
  
  return(final.df)
  
}
