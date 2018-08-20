#==============================================================================
# Redundancy Analysis
#==============================================================================
#'Redundancy Analysis
#'
#'@param metrics.df = data frame of metric values for each station
#'@param analysis = the method for descriminating between two redundant methods.
#'Two choices are available. "wilcox" will run a Wilcoxin Rank Sum Test on the
#' designated upper.class and lower.class for each metric. If two metrics are
#' redundant the metric with the smaller p-value is kept for further analysis.
#' "bibi" requires a data from the "ODE" metric sensitivity analysis. The
#' Descrimination Efficiency threshold from the "ODE" analysis is used to
#'  discriminate between two redundant metrics.  The metric with the larger
#'  threshold is kept for further analysis.
#'@param sensitivity.df = data frame of metrics scored using the Barbour et al. 1996 method
#'for testing metric sensitivity
#'@param method.corr = correlation method to be used.
#'@param upper.coef = the upper correlation coeficient to be used in the
#' correlation analysis.
#'@param lower.coef = the lower correlation coeficient to be used in the
#' correlation analysis.
#'@param upper.class = Only enter this class if using analysis = "wilcox".
#'The reference class to be used in the the Wilcoxin Rank Sum Test.
#'@param lower.class = Only enter this class if using analysis = "wilcox".
#'The degraded class to be used in the the Wilcoxin Rank Sum Test.
#'@return A data frame of metrics scored using the Barbour et al. 1996 method
#'for testing metric sensitivity excluding the "weaker" redundant metric(s)
#'using the interquartiles.
#'@export

old_redundancy <- function(metrics.df, analysis = "bibi",  sensitivity.df,
                       method.corr = "spearman",
                       upper.coef = 0.7, lower.coef = -0.7,
                       upper.class = NULL, lower.class = NULL){

  #Remove the highest scoring metrics if correlated (r < -0.7 or r > 0.7)
  Corr.df <- cor(metrics.df[, 7:ncol(metrics.df)], method = method.corr)
  Corr.df[upper.tri(Corr.df, diag = "TRUE")] <- NA
  Corr.df <- data.frame(Corr.df)
  Corr.df$Metrics <- rownames(Corr.df)
  melted <- reshape2::melt(Corr.df, id.vars = c("Metrics"))
  melted <- melted[!is.na(melted$value), ]
  sub.melted <- na.omit(melted[c(melted$value >= upper.coef | melted$value <= lower.coef),])


  if(analysis == "wilcox"){
    #Greg Pond suggested using the test statistic instead.
    #We both agreed that a combination of the test statistic 
    #and the p-value might be a better option for selecting the better metric.
    
    M2 <- metrics.df[metrics.df$CATEGORY %in% c(upper.class, lower.class), ]
    nt <- lapply(M2[, 7:ncol(M2)], function(x) wilcox.test(x ~ M2[, "CATEGORY"]))
    df <- data.frame(names(nt), matrix(unlist(nt), ncol = 6, byrow = T), stringsAsFactors = FALSE)
    df <- df[, c(1, 3)]
    colnames(df) <- c("METRICS", "P_VALUE")

    wilcox.x <- merge(sub.melted, df, by.x = "Metrics", by.y = "METRICS")
    wilcox.xy <- merge(wilcox.x, df, by.x = "variable", by.y = "METRICS")
    wilcox.xy <- wilcox.xy[, c("Metrics", "variable", "value", "P_VALUE.x", "P_VALUE.y")]
    names(wilcox.xy) <- c("METRICS_X", "METRICS_Y", "CORRELATION", "P_VALUE.X", "P_VALUE.Y")

    wilcox.xy$Remove <- ifelse(wilcox.xy$P_VALUE.X < wilcox.xy$P_VALUE.Y,
                               as.character(wilcox.xy$METRICS_Y),
                               ifelse(wilcox.xy$P_VALUE.X > wilcox.xy$P_VALUE.Y,
                                      as.character(wilcox.xy$METRICS_X),
                                      ifelse(wilcox.xy$P_VALUE.X == wilcox.xy$P_VALUE.Y,
                                             as.character(wilcox.xy$METRICS_X), "ERROR")))
    remove.metric <- unique(wilcox.xy$Remove)
    
    de.df <- sensitivity.df[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    score.df <- de.df[!de.df$METRICS %in% remove.metric, ]
    final.df <- score.df[order(-score.df$SENSITIVITY),]
    final.df

  }else{
    if(analysis == "bibi"){

      de.df <- sensitivity.df[, c("METRICS", "SENSITIVITY")]
      long.y <- merge(sub.melted, de.df, by.x = "Metrics", by.y = "METRICS")
      long.xy <- merge(long.y, de.df , by.x = "variable", by.y = "METRICS")
      colnames(long.xy) <- c("Metrics_X", "Metrics_Y",
                             "Correlation",
                             "DE_SCORE_Y", "DE_SCORE_X")
      long.xy <- long.xy[, c("Metrics_X", "Metrics_Y",
                             "Correlation",
                             "DE_SCORE_X", "DE_SCORE_Y")]
      long.xy$Remove <- ifelse(long.xy$DE_SCORE_X > long.xy$DE_SCORE_Y,
                               as.character(long.xy$Metrics_Y),
                               ifelse(long.xy$DE_SCORE_X < long.xy$DE_SCORE_Y,
                                      as.character(long.xy$Metrics_X),
                                      ifelse(long.xy$DE_SCORE_X == long.xy$DE_SCORE_Y,
                                             as.character(long.xy$Metrics_X), "ERROR")))

      remove.metric <- unique(long.xy$Remove)
      if(is.null(sensitivity.df$THRESHOLD)){
        med_thresh <- sensitivity.df[, c("METRICS", "DISTURBANCE", "THRESH_REF_SEV",
                                         "THRESH_REF_SEV",
                                         "THRESH_REF_MIN",
                                         #"THRESH_REF_NEAR", "THRESH_NEAR_MIN",
                                         "THRESH_MIN_MOD",
                                         "THRESH_MOD_SEV")]
        names(med_thresh) <- c("METRICS", "DISTURBANCE", "THRESHOLD",
                               "THRESH_REF_SEV",
                               "THRESH_REF_MIN",
                               #"THRESH_REF_NEAR", "THRESH_NEAR_MIN",
                               "THRESH_MIN_MOD", 
                               "THRESH_MOD_SEV")
        med_thresh
      }
      if(!is.null(sensitivity.df$THRESHOLD)){
        med_thresh <- sensitivity.df[, c("METRICS", "THRESHOLD", "DISTURBANCE")]
      }

      score.df <- de.df[!de.df$METRICS %in% remove.metric, ]
      final.df <- merge(score.df, med_thresh, by = "METRICS")
      final.df <- final.df[order(-final.df$SENSITIVITY),]
      final.df
    }
  }

  return(final.df)
}

#==============================================================================
# Redundancy Analysis
#==============================================================================
#'Redundancy Analysis
#'
#'@param metrics.df = data frame of metric values for each station
#'@param analysis = the method for descriminating between two redundant methods.
#'Two choices are available. "wilcox" will run a Wilcoxin Rank Sum Test on the
#' designated upper.class and lower.class for each metric. If two metrics are
#' redundant the metric with the smaller p-value is kept for further analysis.
#' "bibi" requires a data from the "ODE" metric sensitivity analysis. The
#' Descrimination Efficiency threshold from the "ODE" analysis is used to
#'  discriminate between two redundant metrics.  The metric with the larger
#'  threshold is kept for further analysis.
#'@param sensitivity.df = data frame of metrics scored using the Barbour et al. 1996 method
#'for testing metric sensitivity
#'@param method.corr = correlation method to be used.
#'@param upper.coef = the upper correlation coeficient to be used in the
#' correlation analysis.
#'@param lower.coef = the lower correlation coeficient to be used in the
#' correlation analysis.
#'@param upper.class = Only enter this class if using analysis = "wilcox".
#'The reference class to be used in the the Wilcoxin Rank Sum Test.
#'@param lower.class = Only enter this class if using analysis = "wilcox".
#'The degraded class to be used in the the Wilcoxin Rank Sum Test.
#'@return A data frame of metrics scored using the Barbour et al. 1996 method
#'for testing metric sensitivity excluding the "weaker" redundant metric(s)
#'using the interquartiles.
#'@export

redundancy_work <- function(metrics.df, analysis = "bibi",  sensitivity.df,
                       method.corr = "spearman",
                       upper.coef = 0.7, lower.coef = -0.7,
                       upper.class = NULL, lower.class = NULL){
  
  corr_long <- function(corr.df, sensitivity.df) {
    de.df <- sensitivity.df[, c("METRICS", "SENSITIVITY")]
    de.df$METRICS <- as.character(de.df$METRICS)
    #de.df <- sensitivity.df[, c("METRICS", "BDE_SENSITIVITY")]
    long.y <- merge(corr.df, de.df, by.x = "Metrics", by.y = "METRICS", all = TRUE)
    final.df <- merge(long.y, de.df, by = "METRICS")
    colnames(final.df) <- c("Metrics_X", "Metrics_Y",
                            "Correlation",
                            "DE_SCORE_Y", "DE_SCORE_X")
    final.df <- final.df[, c("Metrics_X", "Metrics_Y",
                             "Correlation",
                             "DE_SCORE_X", "DE_SCORE_Y")]
    final.df$Metrics_X <- as.character(final.df$Metrics_X)
    final.df$Metrics_Y <- as.character(final.df$Metrics_Y)
    #========================================================================
    final.df$KEEP <- ifelse(final.df$DE_SCORE_X > final.df$DE_SCORE_Y,
                            as.character(final.df$Metrics_X),
                            ifelse(final.df$DE_SCORE_X < final.df$DE_SCORE_Y,
                                   as.character(final.df$Metrics_Y),
                                   ifelse(final.df$DE_SCORE_X == final.df$DE_SCORE_Y &
                                            final.df$Metrics_X < final.df$Metrics_Y,
                                          as.character(final.df$Metrics_Y),
                                          ifelse(final.df$DE_SCORE_X == final.df$DE_SCORE_Y &
                                                   final.df$Metrics_X > final.df$Metrics_Y,
                                                 as.character(final.df$Metrics_X), "ERROR"))))
    
    keep.metric <- c(unique(final.df$KEEP))
    #========================================================================
    final.df$Remove <- ifelse(final.df$DE_SCORE_X > final.df$DE_SCORE_Y,
                              as.character(final.df$Metrics_Y),
                              ifelse(final.df$DE_SCORE_X < final.df$DE_SCORE_Y,
                                     as.character(final.df$Metrics_X),
                                     ifelse(final.df$DE_SCORE_X == final.df$DE_SCORE_Y &
                                              final.df$Metrics_X < final.df$Metrics_Y,
                                            as.character(final.df$Metrics_X),
                                            ifelse(final.df$DE_SCORE_X == final.df$DE_SCORE_Y &
                                                     final.df$Metrics_X > final.df$Metrics_Y,
                                                   as.character(final.df$Metrics_Y), "ERROR"))))
    
    remove.metric <- unique(final.df$Remove)
    return(final.df)
  }
  #============================================================================
  corr_df <- function(metrics.df, method.corr) {
    complete.metrics <- metrics.df[complete.cases(metrics.df[, 7:ncol(metrics.df)]), ]
    corr.df <- cor(complete.metrics[, 7:ncol(complete.metrics)], method = method.corr)
    corr.df[upper.tri(corr.df, diag = "TRUE")] <- NA
    corr.df <- data.frame(corr.df)
    corr.df$Metrics <- rownames(corr.df)
    final.df <- tidyr::gather(corr.df, METRICS, COEF, -Metrics)
    final.df <- final.df[!is.na(final.df$COEF), ]
    return(final.df)
  }
  #============================================================================
  rm.equal.metrics <- unique(as.character(sensitivity.df[sensitivity.df$DISTURBANCE %in% "EQUAL", "METRICS"]))
  metrics.df <- metrics.df[, !names(metrics.df) %in% rm.equal.metrics]
  sensitivity.df <- sensitivity.df[!sensitivity.df$DISTURBANCE %in% "EQUAL", ]
  sensitivity.df <- sensitivity.df[sensitivity.df$METRICS %in% names(metrics.df), ]
  #KEEP the highest scoring metrics if correlated (r < -0.7 or r > 0.7)
  melted <- corr_df(metrics.df, method.corr)
  corr.melted <- na.omit(melted[c(melted$COEF >= upper.coef | melted$COEF <= lower.coef),])
  corr.metrics <- c(unique(corr.melted$Metrics), unique(corr.melted$METRICS))
  uncorr.melted <- na.omit(melted[c(melted$COEF < upper.coef & melted$COEF > lower.coef),])
  uncorr.metrics <- c(unique(uncorr.melted$Metrics), unique(uncorr.melted$METRICS))
  uncorr.metrics <- unique(uncorr.metrics[!uncorr.metrics %in% corr.metrics])
  if(analysis %in% "wilcox"){
      #Greg Pond suggested using the test statistic instead.
      #We both agreed that a combination of the test statistic 
      #and the p-value might be a better option for selecting the better metric.
      
      M2 <- metrics.df[metrics.df$CATEGORY %in% c(upper.class, lower.class), ]
      nt <- lapply(M2[, 7:ncol(M2)], function(x) wilcox.test(x ~ M2[, "CATEGORY"]))
      df <- data.frame(names(nt), matrix(unlist(nt), ncol = 6, byrow = T), stringsAsFactors = FALSE)
      df <- df[, c(1, 3)]
      #colnames(df) <- c("METRICS", "P_VALUE")
      colnames(df) <- c("METRICS", "SENSITIVITY")
      long.xy <- corr_long(corr.melted, df)
      remove.metrics <- unique(long.xy$Remove)
      keep.metrics <- unique(long.xy$KEEP)
      remove.cols <- remove.metrics[!remove.metrics %in% keep.metrics]
      #========================================================================
      de.df <- sensitivity.df[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      score.df <- de.df[!de.df$METRICS %in% remove.cols, ]
      final.df <- score.df[order(-score.df$SENSITIVITY),]
      #final.df
      
      #wilcox.x <- merge(corr.melted, df, by.x = "Metrics", by.y = "METRICS")
      #wilcox.xy <- merge(wilcox.x, df, by = "METRICS")
      #wilcox.xy <- wilcox.xy[, c("Metrics", "METRICS", "value", "P_VALUE.x", "P_VALUE.y")]
      #names(wilcox.xy) <- c("METRICS_X", "METRICS_Y", "CORRELATION", "P_VALUE.X", "P_VALUE.Y")
      
      #wilcox.xy$KEEP <- ifelse(wilcox.xy$P_VALUE.X < wilcox.xy$P_VALUE.Y,
      #                         as.character(wilcox.xy$METRICS_X),
      #                         ifelse(wilcox.xy$P_VALUE.X > wilcox.xy$P_VALUE.Y,
      #                                as.character(wilcox.xy$METRICS_Y),
      #                                ifelse(wilcox.xy$P_VALUE.X == wilcox.xy$P_VALUE.Y,
      #                                       as.character(wilcox.xy$METRICS_X), "ERROR")))
      #keep.metric <- unique(wilcox.xy$KEEP)
      
      #de.df <- sensitivity.df[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      #score.df <- de.df[!de.df$METRICS %in% keep.metric, ]
      #final.df <- score.df[order(-score.df$SENSITIVITY),]
      #final.df
      
    
  }else{
    if(analysis %in% "bibi"){
      long.xy <- corr_long(corr.melted, sensitivity.df)
      remove.metrics <- unique(long.xy$Remove)
      keep.metrics <- unique(long.xy$KEEP)
      remove.cols <- remove.metrics[!remove.metrics %in% keep.metrics]
      #========================================================================
      
      #if(is.null(sensitivity.df$THRESHOLD)){
      if(is.null(sensitivity.df$BOUND_BI_CMA)){
        med_thresh <- sensitivity.df[, c("METRICS", "DISTURBANCE", "THRESH_REF_SEV",
                                         "THRESH_REF_SEV",
                                         "THRESH_REF_MIN",
                                         #"THRESH_REF_NEAR", "THRESH_NEAR_MIN",
                                         "THRESH_MIN_MOD",
                                         "THRESH_MOD_SEV")]
        names(med_thresh) <- c("METRICS", "DISTURBANCE", "THRESHOLD",
                               "THRESH_REF_SEV",
                               "THRESH_REF_MIN",
                               #"THRESH_REF_NEAR", "THRESH_NEAR_MIN",
                               "THRESH_MIN_MOD", 
                               "THRESH_MOD_SEV")
        med_thresh
      }
      #if(!is.null(sensitivity.df$THRESHOLD)){
      if(!is.null(sensitivity.df$BOUND_BI_CMA)){
        #med_thresh <- sensitivity.df[, c("METRICS", "THRESHOLD", "DISTURBANCE")]
        med_thresh <- sensitivity.df[, c("METRICS", "BOUND_BI_CMA", "DISTURBANCE")]
      }
      
      de.df <- sensitivity.df[, c("METRICS", "SENSITIVITY")]
      score.df <- de.df[!de.df$METRICS %in% remove.cols, ]
      final.df <- merge(score.df, med_thresh, by = "METRICS")
      final.df <- final.df[order(-final.df$SENSITIVITY),]
      final.df$METRICS <- as.character(final.df$METRICS)
      #final.df
    }
  }
  
  return(final.df)
}

#==============================================================================
# Redundancy Analysis
#==============================================================================
#'Redundancy Analysis
#'
#'@param metrics.df = data frame of metric values for each station
#'@param analysis = the method for descriminating between two redundant methods.
#'Two choices are available. "wilcox" will run a Wilcoxin Rank Sum Test on the
#' designated upper.class and lower.class for each metric. If two metrics are
#' redundant the metric with the smaller p-value is kept for further analysis.
#' "bibi" requires a data from the "ODE" metric sensitivity analysis. The
#' Descrimination Efficiency threshold from the "ODE" analysis is used to
#'  discriminate between two redundant metrics.  The metric with the larger
#'  threshold is kept for further analysis.
#'@param sensitivity.df = data frame of metrics scored using the Barbour et al. 1996 method
#'for testing metric sensitivity
#'@param method.corr = correlation method to be used.
#'@param upper.coef = the upper correlation coeficient to be used in the
#' correlation analysis.
#'@param lower.coef = the lower correlation coeficient to be used in the
#' correlation analysis.
#'@param upper.class = Only enter this class if using analysis = "wilcox".
#'The reference class to be used in the the Wilcoxin Rank Sum Test.
#'@param lower.class = Only enter this class if using analysis = "wilcox".
#'The degraded class to be used in the the Wilcoxin Rank Sum Test.
#'@return A data frame of metrics scored using the Barbour et al. 1996 method
#'for testing metric sensitivity excluding the "weaker" redundant metric(s)
#'using the interquartiles.
#'@export

redundancy <- function(metrics.df, analysis = "bibi",  sensitivity.df,
                            method.corr = "spearman",
                            upper.coef = 0.7, lower.coef = -0.7,
                            upper.class = NULL, lower.class = NULL){
  for (i in 1:25){
    #print(i)
    redund.df <- redundancy_work(metrics.df, analysis, sensitivity.df,
                            method.corr,
                            upper.coef, lower.coef,
                            upper.class, lower.class)
    keep.metrics <- unique(as.character(redund.df$METRICS))
    keep.metrics
    keep.cols <- names(metrics.df[, 1:6])
    keep.cols <- c(keep.cols, keep.metrics)
    metrics.df <- metrics.df[, keep.cols]
    #keep.cols[!keep.cols %in% names(metrics.df)]
  }
  return(redund.df)
}

