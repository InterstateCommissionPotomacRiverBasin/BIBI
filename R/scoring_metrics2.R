#==============================================================================
#New Scoring metrics
#==============================================================================
#'Metrics Summary
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param bioregion = the bioregion of interest.
#'@return Scores the samples.
#'@export

metrics_summary <- function(metrics.df, bioregion, zero = TRUE){
  #metrics.df <- metrics.df[metrics.df$ABUNDANCE >= 70, ]
  metrics.df <- metrics.df[metrics.df$BIOREGION %in% bioregion, ]
  metrics.df <- metrics.df[, !names(metrics.df) %in% "BIOREGION"]
  
  
  
  ref.df <- metrics.df[metrics.df$CATEGORY %in% "REF", ]
  if(zero == TRUE){
    ref.zero <- ref.df
    ref.zero[, 7:ncol(ref.zero)] <- ifelse(ref.zero[, 7:ncol(ref.zero)] == 0, 1, 0)
    
  }
  
  bi_de <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "DE"))
  names(bi_de) <- c("METRICS", "DISTURBANCE", "BINARY_DE")
  
  if(zero == TRUE){
    if(ncol(ref.zero) > 7){
      zero.cols <- names(ref.zero[, 7:ncol(ref.zero)])[(colSums(ref.zero[, 7:ncol(ref.zero)]) /
                                                          nrow(ref.zero[, 7:ncol(ref.zero)])) >= 0.05]
    }
    
    if(ncol(ref.zero) == 7){
      zero.cols <- names(ref.zero[, 7])[(sum(ref.zero[, 7], na.rm = TRUE) / nrow(ref.zero[, 7])) >= 0.05]
    } 
    
    if(any(bi_de$DISTURBANCE %in% "EQUAL")){
      equal.cols <- unlist(list(bi_de[bi_de$DISTURBANCE %in% "EQUAL", "METRICS"]))
      zero.cols <- unlist(list(zero.cols, equal.cols))
    }
    zero.keep <- unlist(list(names(metrics.df[, 1:6]), zero.cols, "ABUNDANCE"))
    
    new.without.zero <- metrics.df[, names(metrics.df) %in% zero.keep]
    if(ncol(new.without.zero) > 6) new.without.zero[new.without.zero == 0] <- NA
    
    new.with.zero <- metrics.df[, !names(metrics.df) %in% zero.cols]
  }else{
    new.with.zero <- metrics.df
  }
  
  
  
  sub_ms <- function(metrics.df){
    #pair.cma <- unique(pairwise_sensitivity(metrics.df, "ODE"))
    #names(pair.cma)[names(pair.cma) %in% "SENSITIVITY"] <- "PAIRWISE_CMA"
    bi.cma <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "CMA"))
    bi.cma$BOUND <- ifelse(bi.cma$DISTURBANCE %in% c("DECREASE", "EQUAL"),
                           bi.cma$THRESHOLD - abs(bi.cma$MEDIAN - bi.cma$THRESHOLD),
                           ifelse(bi.cma$DISTURBANCE %in% "INCREASE",
                                  bi.cma$THRESHOLD + abs(bi.cma$MEDIAN - bi.cma$THRESHOLD),
                                  100000))
    names(bi.cma) <- c("METRICS", "DISTURBANCE", "BINARY_CMA",
                       "BALANCED_BINARY_CMA", "PRECENTILE_BINARY_CMA",
                       "PCT_REF_BI_CMA", "PCT_DEG_BI_CMA",
                       "REF_MEDIAN", "MIDPOINT_BI_CMA",
                       "BOUND_BI_CMA")
    bi.cma$SENSITIVITY <- (as.numeric(as.character(bi.cma$PCT_REF)) + as.numeric(as.character(bi.cma$PCT_DEG))) / 2
    
    bi.de <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "DE"))
    names(bi.de) <- c("METRICS", "DISTURBANCE", "BINARY_DE")
    bi_barbour <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "BARBOUR"))
    names(bi_barbour) <- c("METRICS", "DISTURBANCE", "BINARY_BARBOUR")
    range.var <- unique(range_variability(metrics.df))
    names(range.var) <- c("METRICS", "REF_MIN", "REF_MAX", "REF_RANGE_VALUE",
                          "REF_RANGE_CLASS", "REF_Q25", "REF_Q75",
                          "REF_VARIABILITY_VALUE", "REF_VARIABILITY_CLASS")
    #zero.inflate <- zero_inflate(metrics.df, bi_barbour)
    final.df <- plyr::join_all(list(#pair.cma,
      bi.cma, 
      bi.de[, c(1, 3)],
      bi_barbour[, c(1, 3)], 
      range.var),
      #zero.inflate[,c(1, 3:5)]),
      "METRICS")
    return(final.df)
  }
  
  new.with.zero <- new.with.zero[, !names(new.with.zero) %in% "ABUNDANCE"]
  if(zero == TRUE){
    new.without.zero <- new.without.zero[, !names(new.without.zero) %in% "ABUNDANCE"]
    
    if(ncol(new.with.zero) >= 7 & ncol(new.without.zero) >= 7){
      final_with_zero <- sub_ms(new.with.zero)
      final_without_zero <- sub_ms(new.without.zero)
      final.df <- rbind(final_with_zero, final_without_zero)
    }
    
    if(ncol(new.with.zero) >= 7 & ncol(new.without.zero) < 7){
      final.df <- sub_ms(new.with.zero)
    }
    
    if(ncol(new.with.zero) < 7 & ncol(new.without.zero) >= 7){
      final.df <- sub_ms(new.without.zero)
    }
  }
  if(zero != TRUE){
    
    if(ncol(new.with.zero) >= 7){
      final.df <- sub_ms(new.with.zero)
    }
    
  }
 
  
  return(final.df)
}

#==============================================================================
#New Scoring metrics
#==============================================================================
#'Metrics Summary
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param bioregion = the bioregion of interest.
#'@return Scores the samples.
#'@export

metrics_summary2 <- function(metrics.df, bioregion){
  #metrics.df <- metrics.df[metrics.df$ABUNDANCE >= 70, ]
  metrics.df <- metrics.df[metrics.df$BIOREGION %in% bioregion, ]
  metrics.df <- metrics.df[, !names(metrics.df) %in% "BIOREGION"]
  
  
  ref.df <- metrics.df[metrics.df$CATEGORY %in% "REF", ]
  #ref.zero <- ref.df
  #ref.zero[, 7:ncol(ref.zero)] <- ifelse(ref.zero[, 7:ncol(ref.zero)] == 0, 1, 0)
  
  bi_de <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "DE"))
  names(bi_de) <- c("METRICS", "DISTURBANCE", "BINARY_DE")
  
 
  
  if(any(bi_de$DISTURBANCE %in% "EQUAL")){
    equal.cols <- unlist(list(bi_de[bi_de$DISTURBANCE %in% "EQUAL", "METRICS"]))
    #zero.cols <- unlist(list(zero.cols, equal.cols))
  }
  #zero.keep <- unlist(list(names(metrics.df[, 1:6]), zero.cols, "ABUNDANCE"))
  
  new.with.zero <- metrics.df
  #new.without.zero <- metrics.df[, names(metrics.df) %in% zero.keep]
  #if(ncol(new.without.zero) > 6) new.without.zero[new.without.zero == 0] <- NA
  
  sub_ms <- function(metrics.df){
    #pair.cma <- unique(pairwise_sensitivity(metrics.df, "ODE"))
    #names(pair.cma)[names(pair.cma) %in% "SENSITIVITY"] <- "PAIRWISE_CMA"
    bi.cma <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "CMA"))
    bi.cma$BOUND <- ifelse(bi.cma$DISTURBANCE %in% c("DECREASE", "EQUAL"),
                           bi.cma$THRESHOLD - abs(bi.cma$MEDIAN - bi.cma$THRESHOLD),
                           ifelse(bi.cma$DISTURBANCE %in% "INCREASE",
                                  bi.cma$THRESHOLD + abs(bi.cma$MEDIAN - bi.cma$THRESHOLD),
                                  100000))
    names(bi.cma) <- c("METRICS", "DISTURBANCE", "BINARY_CMA",
                       "BALANCED_BINARY_CMA", "PRECENTILE_BINARY_CMA",
                       "PCT_REF_BI_CMA", "PCT_DEG_BI_CMA",
                       "REF_MEDIAN", "MIDPOINT_BI_CMA",
                       "BOUND_BI_CMA")
    bi.cma$SENSITIVITY <- (as.numeric(as.character(bi.cma$PCT_REF)) + as.numeric(as.character(bi.cma$PCT_DEG))) / 2
    
    bi.de <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "DE"))
    names(bi.de) <- c("METRICS", "DISTURBANCE", "BINARY_DE")
    bi_barbour <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "BARBOUR"))
    names(bi_barbour) <- c("METRICS", "DISTURBANCE", "BINARY_BARBOUR")
    range.var <- unique(range_variability(metrics.df))
    names(range.var) <- c("METRICS", "REF_MIN", "REF_MAX", "REF_RANGE_VALUE",
                          "REF_RANGE_CLASS", "REF_Q25", "REF_Q75",
                          "REF_VARIABILITY_VALUE", "REF_VARIABILITY_CLASS")
    #zero.inflate <- zero_inflate(metrics.df, bi_barbour)
    final.df <- plyr::join_all(list(#pair.cma,
      bi.cma, 
      bi.de[, c(1, 3)],
      bi_barbour[, c(1, 3)], 
      range.var),
      #zero.inflate[,c(1, 3:5)]),
      "METRICS")
    return(final.df)
  }

  
    final.df <- sub_ms(new.with.zero)
    
    #rand.num <- sample(nrow(new.with.zero), 20)
    #test <- new.with.zero[rand.num, ]
    #sub_ms(test)

  return(final.df)
}

#==============================================================================
#'Zero Inflation
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param bi.barbour = binary barbour output.
#'@return Scores the samples.
#'@export

zero_inflate <- function(metrics.df, bi.barbour){
 # metrics.orig <- metrics.df
  #metrics.df <- metrics.orig[1:20,]
  
  ref.df <- metrics.df[metrics.df$CATEGORY %in% "REF", ]
  if(ncol(ref.df) > 7){
    ref.df <- ref.df[, c(names(ref.df[, 1:6]), sort(names(ref.df[, 7:ncol(ref.df)])))]
  }
  
  deg.df <- metrics.df[metrics.df$CATEGORY %in% "SEV", ]
  if(ncol(deg.df) > 7){
    deg.df <- deg.df[, c(names(deg.df[, 1:6]), sort(names(deg.df[, 7:ncol(deg.df)])))]
  }
  
  
  barb <- bi.barbour[, c("METRICS", "DISTURBANCE")]
  new.df <- data.frame(METRICS = names(metrics.df)[7])
  if(ncol(metrics.df) > 7) new.df <- data.frame(METRICS = names(metrics.df[, 7:ncol(metrics.df)]))
  # all.x = TRUE added 4/13/17
  new.df <- merge(new.df , barb, by = "METRICS", all.x = TRUE)
  
  if(ncol(ref.df) > 7){
    new.df$PCT_0_REF <- apply(ref.df[, 7:ncol(ref.df)], 2, function(x){
      round((sum(x == 0) / length(x)) * 100, 0)
    })
    
    new.df$PCT_0_DEG <- apply(deg.df[, 7:ncol(deg.df)], 2, function(x){
      round((sum(x == 0) / length(x)) * 100, 0)
    })
  }else{
    new.df$PCT_0_REF <- round((sum(ref.df[, 7] == 0) / length(ref.df[, 7])) * 100, 0)
    new.df$PCT_0_DEG <- round((sum(deg.df[, 7] == 0) / length(deg.df[, 7])) * 100, 0)
  }
  
  
  
  
  new.df$ZERO_INFLATE <- ifelse(new.df$PCT_0_REF > 10 & new.df$PCT_0_REF <= 50 &
                                  new.df$PCT_0_DEG > 10 & new.df$PCT_0_DEG <= 50, "REVIEW",
                                ifelse(new.df$PCT_0_REF > 10 & new.df$PCT_0_REF <= 50 & new.df$PCT_0_DEG > 50, "REVIEW",
                                       ifelse(new.df$PCT_0_REF > 50 & new.df$PCT_0_DEG > 10 & new.df$PCT_0_DEG <= 50, "REVIEW",
                                              ifelse(new.df$PCT_0_REF > 50 & new.df$PCT_0_DEG > 50, "POOR",
                                                     ifelse(new.df$PCT_0_REF <= 10 & new.df$PCT_0_DEG <= 10, "GOOD",
                                                            ifelse(new.df$DISTURBANCE %in% "DECREASE" &
                                                                     new.df$PCT_0_REF > 10 & new.df$PCT_0_DEG <= 10, "POOR",
                                                                   ifelse(new.df$DISTURBANCE %in% "DECREASE" &
                                                                            new.df$PCT_0_REF <= 10 & new.df$PCT_0_DEG > 10, "GOOD",
                                                                          ifelse(new.df$DISTURBANCE %in% "INCREASE" &
                                                                                   new.df$PCT_0_REF > 10 & new.df$PCT_0_DEG <= 10, "GOOD",
                                                                                 ifelse(new.df$DISTURBANCE %in% "INCREASE" &
                                                                                          new.df$PCT_0_REF <= 10 & new.df$PCT_0_DEG > 10, "POOR", "ERROR")))))))))
  
  return(new.df)
}

#==============================================================================
#'Metric Class Break
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param class = metric class.
#'@return Scores the samples.
#'@export

metric_class_break <- function(metric_class, class){
  
  diversity <- metric_class[metric_class$METRIC_CLASS %in% "DIVERSITY", ]
  FFG <- metric_class[metric_class$METRIC_CLASS %in% "FFG", ]
  HABIT <- metric_class[metric_class$METRIC_CLASS %in% "HABIT", ]
  TOL <- metric_class[metric_class$METRIC_CLASS %in% "TOLERANCE", ]
  scores3 <- unlist(list(as.character(diversity$METRICS),
                         as.character(FFG$METRICS),
                         as.character(HABIT$METRICS),
                         as.character(TOL$METRICS),
                         "EVENT_ID", "CATEGORY", "STATION_ID",
                         "SAMPLE_NUMBER", "AGENCY_CODE", "DATE"))
  scores4 <- unlist(scores3)
  
  div <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(diversity$METRICS)))
  ffg <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(FFG$METRICS)))
  habit <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(HABIT$METRICS)))
  tol <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(TOL$METRICS)))
  comp.exc <- unlist(list(as.character(diversity$METRICS),
                          as.character(FFG$METRICS),
                          as.character(HABIT$METRICS),
                          as.character(TOL$METRICS)))
  comp <- names(metrics.df[, !names(metrics.df) %in% comp.exc])
  comp <- comp[!comp %in% c("BIOREGION", "ABUNDANCE", "EFFECTIVE_RICH_SIMPSON", "NO_MATCH")]
  
  if(class %in% "DIVERSITY") final.vec <- div
  if(class %in% "FFG") final.vec <- ffg
  if(class %in% "HABIT") final.vec <- habit
  if(class %in% "TOLERANCE") final.vec <- tol
  if(class %in% "COMPOSITION") final.vec <- comp
  
  return(final.vec)
}

#==============================================================================
#'Best Distribution H
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param metric.summary
#'@param scores.df
#'@param m.c
#'@return Scores the samples.
#'@export

best.distribution.H <- function(metrics.df, metric.summary, scores.df, m.c,
                                range_var = TRUE, redund = FALSE, best.metrics){
  if(range_var == TRUE){
    metric.summary <- metric.summary[(metric.summary$REF_RANGE_CLASS %in% "HIGH" &
                                        metric.summary$REF_VARIABILITY_CLASS %in% "LOW"), ]
  }
  
  if(redund == TRUE & nrow(metric.summary) >= 2){
    list.metrics <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                as.character(metric.summary $METRICS)))
    R_test <- redundancy(metrics.df[, names(metrics.df) %in% list.metrics],
                         analysis = "wilcox", sensitivity.df = metric.summary,
                         lower.class = "SEV", upper.class = "REF",
                         lower.coef = -0.85, upper.coef = 0.85)
    metric.summary  <- metric.summary[metric.summary$METRICS %in% R_test$METRICS, ]
  }
  
  datalist = list()
  for (n in 1:100) {
    # ... make some data
    ms.n <- metric.summary[metric.summary$SENSITIVITY >= n, ]
    ms.bm <- metric.summary[metric.summary$METRICS %in% best.metrics$METRICS, ]
    ms.df <- unique(rbind(ms.n, ms.bm))
    #if(range_var == TRUE){
    #  ms.df <- ms.df[(ms.df$REF_RANGE_CLASS %in% "HIGH" & ms.df$REF_VARIABILITY_CLASS %in% "LOW"), ]
    #}
    
    if(nrow(ms.df) >= 5){
      metrics.vec <- ms.df$METRICS
      metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                 "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                 "BIOREGION", as.character(metrics.vec)))
      #sub.metrics <- metrics.df[, metrics.vec]
      #os.df <- old_scoring3(sub.metrics, i, FALSE)
      os.df <- scores.df[, names(scores.df) %in% metrics.vec]
      if(nrow(ms.df) > 1){
        os.df[, 7:ncol(os.df)] <- apply(os.df[, 7:ncol(os.df)], 2, function(x) as.numeric(as.character(x)))
        os.df$FINAL_SCORE <- apply(os.df[, 7:ncol(os.df)], 1, mean)
      }else{
        os.df[, 7] <- as.numeric(as.character(os.df[, 7]))
        os.df$FINAL_SCORE <- os.df[, 7]
      }
      test <- scree(os.df)
      
    }else{
      test <- data.frame(NA, NA, NA, NA, NA, NA)
      names(test) <- c("PERCENTILE", "PCT_REF", "PCT_DEG",
                       "SENSITIVITY", "NEW_SENSITIVITY", "PERCENTILE_VALUE")
    }
    test$THRESH <- n
    datalist[[n]] <- test # add it to your list
  }
  
  big_data <- unique(do.call(rbind, datalist))
  big_data <- big_data[order(-big_data$SENSITIVITY, -big_data$THRESH), ]
  #============================================================================
  thresh <- big_data[1, 7]
  ms.n <- metric.summary[metric.summary$SENSITIVITY >= thresh, ]
  ms.bm <- metric.summary[metric.summary$METRICS %in% best.metrics$METRICS, ]
  final.df <- unique(rbind(ms.n, ms.bm))
  return(final.df)
}
#==============================================================================
#'Select the Best Metrics
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param metric.summary
#'@param scores.df
#'@param m.c
#'@return Select the best metrics from 4 of the 5 metric categories
#'@export

select_best <- function(metrics.df, metric.summary, m.c, range_var, redund){
  if(range_var == TRUE){
    metric.summary <- metric.summary[(metric.summary$REF_RANGE_CLASS %in% "HIGH" &
                                        metric.summary$REF_VARIABILITY_CLASS %in% "LOW"), ]
  }
  
  if(redund == TRUE & nrow(metric.summary) >= 2){
    list.metrics <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                as.character(metric.summary $METRICS)))
    R_test <- redundancy(metrics.df[, names(metrics.df) %in% list.metrics],
                         analysis = "wilcox", sensitivity.df = metric.summary,
                         lower.class = "SEV", upper.class = "REF",
                         lower.coef = -0.85, upper.coef = 0.85)
    metric.summary  <- metric.summary[metric.summary$METRICS %in% R_test$METRICS, ]
  }
  
  ms.t <- merge(m.c, metric.summary, by = "METRICS", all.y = T)
  ms.t$METRIC_CLASS <- ifelse(is.na(ms.t$METRIC_CLASS), "COMPOSITION", as.character(ms.t$METRIC_CLASS))
  
  ms.div <- ms.t[ms.t$METRIC_CLASS %in% "DIVERSITY" &
                   ms.t$BINARY_CMA == max(ms.t[ms.t$METRIC_CLASS %in% "DIVERSITY", "BINARY_CMA"]), ]
  if(nrow(ms.div) > 1){
    ms.div$CHECK <- ms.div$REF_RANGE_VALUE - ms.div$REF_VARIABILITY_VALUE
    ms.div <- ms.div[ms.div$CHECK == max(ms.div$CHECK), ]
    ms.div <- ms.div[, !names(ms.div) %in% "CHECK"]
  } 
  #============================================================================
  ms.ffg <- ms.t[ms.t$METRIC_CLASS %in% "FFG" &
                   ms.t$BINARY_CMA == max(ms.t[ms.t$METRIC_CLASS %in% "FFG", "BINARY_CMA"]), ]
  if(nrow(ms.ffg) > 1){
    ms.ffg$CHECK <- ms.ffg$REF_RANGE_VALUE - ms.ffg$REF_VARIABILITY_VALUE
    ms.ffg <- ms.ffg[ms.ffg$CHECK == max(ms.ffg$CHECK), ]
    ms.ffg <- ms.ffg[, !names(ms.ffg) %in% "CHECK"]
  } 
  #============================================================================
  ms.hab <- ms.t[ms.t$METRIC_CLASS %in% "HABIT" &
                   ms.t$BINARY_CMA == max(ms.t[ms.t$METRIC_CLASS %in% "HABIT", "BINARY_CMA"]), ]
  if(nrow(ms.hab) > 1){
    ms.hab$CHECK <- ms.hab$REF_RANGE_VALUE - ms.hab$REF_VARIABILITY_VALUE
    ms.hab <- ms.hab[ms.hab$CHECK == max(ms.hab$CHECK), ]
    ms.hab <- ms.hab[, !names(ms.hab) %in% "CHECK"]
  } 
  #============================================================================
  ms.tol <- ms.t[ms.t$METRIC_CLASS %in% "TOLERANCE" &
                   ms.t$BINARY_CMA == max(ms.t[ms.t$METRIC_CLASS %in% "TOLERANCE", "BINARY_CMA"]), ]
  if(nrow(ms.tol) > 1){
    ms.tol$CHECK <- ms.tol$REF_RANGE_VALUE - ms.tol$REF_VARIABILITY_VALUE
    ms.tol <- ms.tol[ms.tol$CHECK == max(ms.tol$CHECK), ]
    ms.tol <- ms.tol[, !names(ms.tol) %in% "CHECK"]
  } 
  #============================================================================
  ms.comp <- ms.t[ms.t$METRIC_CLASS %in% "COMPOSITION" &
                    ms.t$BINARY_CMA == max(ms.t[ms.t$METRIC_CLASS %in% "COMPOSITION", "BINARY_CMA"]), ]
  if(nrow(ms.comp) > 1){
    ms.comp$CHECK <- ms.comp$REF_RANGE_VALUE - ms.comp$REF_VARIABILITY_VALUE
    ms.comp <- ms.comp[ms.comp$CHECK ==max(ms.comp$CHECK), ]
    ms.comp <- ms.comp[, !names(ms.comp) %in% "CHECK"]
  } 
  #============================================================================
  best.metrics <- rbind(ms.div, ms.ffg, ms.hab, ms.tol, ms.comp)
  best.metrics <- best.metrics[order(best.metrics$BINARY_CMA, decreasing =  T), ]
  final.df <- best.metrics[1:4, ]
  
  return(final.df)
}

#==============================================================================
#'Best Distribution 
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param metric.summary
#'@param scores.df
#'@param m.c
#'@return Scores the samples.
#'@export


best.distribution0 <- function(metrics.df, metric.summary, scores.df, m.c, range_var = TRUE, redund = FALSE){
  datalist = list()
  
  for (n in 1:100) {
    # ... make some data
    ms.df <- metric.summary[metric.summary$SENSITIVITY >= n, ]
    if(range_var == TRUE){
      ms.df <- ms.df[(ms.df$REF_RANGE_CLASS %in% "HIGH" & ms.df$REF_VARIABILITY_CLASS %in% "LOW"), ]
    }
    
    if(redund == TRUE & nrow(ms.df) >= 2){
      list.metrics <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                  "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                  as.character(ms.df$METRICS)))
      R_test <- redundancy(metrics.df[, names(metrics.df) %in% list.metrics],
                           analysis = "wilcox", sensitivity.df = ms.df,
                           lower.class = "SEV", upper.class = "REF",
                           lower.coef = -0.85, upper.coef = 0.85)
      ms.df <- ms.df[ms.df$METRICS %in% R_test$METRICS, ]
    }
    
    if(nrow(ms.df) >= 5){
      metrics.vec <- ms.df$METRICS
      metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                 "BIOREGION", as.character(metrics.vec)))
      #sub.metrics <- metrics.df[, metrics.vec]
      #os.df <- old_scoring3(sub.metrics, i, FALSE)
      os.df <- scores.df[, names(scores.df) %in% metrics.vec]
      if(nrow(ms.df) > 1){
        os.df[, 7:ncol(os.df)] <- apply(os.df[, 7:ncol(os.df)], 2, function(x) as.numeric(as.character(x)))
        os.df$FINAL_SCORE <- apply(os.df[, 7:ncol(os.df)], 1, mean)
      }else{
        os.df[, 7] <- as.numeric(as.character(os.df[, 7]))
        os.df$FINAL_SCORE <- os.df[, 7]
      }
      test <- scree(os.df)
      
    }else{
      test <- data.frame(NA, NA, NA, NA, NA, NA)
      names(test) <- c("PERCENTILE", "PCT_REF", "PCT_DEG",
                       "SENSITIVITY", "NEW_SENSITIVITY", "PERCENTILE_VALUE")
    }
    test$THRESH <- n
    datalist[[n]] <- test # add it to your list
  }
  
  big_data <- unique(do.call(rbind, datalist))
  big_data <- big_data[order(-big_data$SENSITIVITY, -big_data$THRESH), ]
  return(big_data)
}
#==============================================================================
#'Best Distribution 
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param metric.summary
#'@param scores.df
#'@param m.c
#'@return Scores the samples.
#'@export


best.distribution <- function(metrics.df, metric.summary, scores.df, m.c, range_var = TRUE){
  datalist = list()
  
  for (n in 1:100) {
    # ... make some data
    ms.df <- metric.summary[metric.summary$SENSITIVITY >= n, ]
    if(range_var == TRUE){
      ms.df <- ms.df[(ms.df$REF_RANGE_CLASS %in% "HIGH" & ms.df$REF_VARIABILITY_CLASS %in% "LOW"), ]
    }
    
    if(nrow(ms.df) >= 5){
      scored <- scoring(metrics.df, scores.df, ms.df, n, m.c)
      test <- scree(scored)
    }else{
      test <- data.frame(NA, NA, NA, NA, NA, NA)
      names(test) <- c("PERCENTILE", "PCT_REF", "PCT_DEG",
                       "SENSITIVITY", "NEW_SENSITIVITY", "PERCENTILE_VALUE")
    }
    test$THRESH <- n
    datalist[[n]] <- test # add it to your list
  }
  
  big_data <- do.call(rbind, datalist)
  big_data <- big_data[order(-big_data$SENSITIVITY, -big_data$THRESH), ]
  return(big_data)
}

#==============================================================================
#'Best Distribution 2
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param metric.summary
#'@param scores.df
#'@param m.c
#'@return Scores the samples.
#'@export


best.distribution2 <- function(metrics.df, metric.summary, scores.df, m.c,
                               Family = FALSE, master, range_var = TRUE){
  
  ms_div <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "DIVERSITY", "METRICS"]), ]
 
  ms_ffg <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "FFG", "METRICS"]), ]
  
  ms_hab <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "HABIT", "METRICS"]), ]
  
  ms_tol <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "TOLERANCE", "METRICS"]), ]
  
  ms_comp <- metric.summary[metric.summary$METRICS %in% 
                              (m.c[m.c$METRIC_CLASS %in% "COMPOSITION", "METRICS"]), ]
  
  if(Family == TRUE){
    fam.list <- paste("PCT", na.omit(unique(master$FAMILY)), sep = "_")
    ms_fam <- metric.summary[metric.summary$METRICS %in% fam.list, ]
    ms_comp <- ms_comp[!ms_comp$METRICS %in% fam.list, ]
  }
  
  if(range_var == TRUE){
    ms_div <- ms_div[(ms_div$REF_RANGE_CLASS %in% "HIGH" &
                        ms_div$REF_VARIABILITY_CLASS %in% "LOW"), ]
    ms_ffg <- ms_ffg[(ms_ffg$REF_RANGE_CLASS %in% "HIGH" &
                        ms_ffg$REF_VARIABILITY_CLASS %in% "LOW"), ]
    ms_hab <- ms_hab[(ms_hab$REF_RANGE_CLASS %in% "HIGH" &
                        ms_hab$REF_VARIABILITY_CLASS %in% "LOW"), ]
    ms_tol <- ms_tol[(ms_tol$REF_RANGE_CLASS %in% "HIGH" &
                        ms_tol$REF_VARIABILITY_CLASS %in% "LOW"), ]
    ms_comp <- ms_comp[(ms_comp$REF_RANGE_CLASS %in% "HIGH" &
                          ms_comp$REF_VARIABILITY_CLASS %in% "LOW"), ]
    if(Family == TRUE){
      ms_fam <- ms_fam[(ms_fam$REF_RANGE_CLASS %in% "HIGH" &
                          ms_fam$REF_VARIABILITY_CLASS %in% "LOW"), ]
    } 
  }
  
  diversity <- m.c[m.c$METRIC_CLASS %in% "DIVERSITY", ]
  FFG <- m.c[m.c$METRIC_CLASS %in% "FFG", ]
  HABIT <- m.c[m.c$METRIC_CLASS %in% "HABIT", ]
  TOL <- m.c[m.c$METRIC_CLASS %in% "TOLERANCE", ]
  scores3 <- unlist(list(as.character(diversity$METRICS),
                         as.character(FFG$METRICS),
                         as.character(HABIT$METRICS),
                         as.character(TOL$METRICS),
                         "EVENT_ID", "CATEGORY", "STATION_ID",
                         "SAMPLE_NUMBER", "AGENCY_CODE", "DATE"))
  scores4 <- unlist(scores3)
  
  div <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(diversity$METRICS)))
  ffg <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(FFG$METRICS)))
  habit <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(HABIT$METRICS)))
  tol <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(TOL$METRICS)))
  comp.exc <- unlist(list(as.character(diversity$METRICS),
                          as.character(FFG$METRICS),
                          as.character(HABIT$METRICS),
                          as.character(TOL$METRICS)))
  comp <- names(metrics.df[, !names(metrics.df) %in% comp.exc])
  comp <- comp[!comp %in% c("BIOREGION", "ABUNDANCE", "EFFECTIVE_RICH_SIMPSON", "NO_MATCH")]
  

  scree_classes <- function(metrics.df, scores.df, metric.summary, class){
    datalist = list()
    for (n in 50:100) {
      
      ms.df <- metric.summary[metric.summary$SENSITIVITY >= n, ]
      if(range_var == TRUE){
        ms.df <- ms.df[(ms.df$REF_RANGE_CLASS %in% "HIGH" & ms.df$REF_VARIABILITY_CLASS %in% "LOW"), ]
      }
      
      if(length(metrics.df[, names(metrics.df) %in% class]) > 7){
        R_test <- redundancy(metrics.df[, names(metrics.df) %in% class],
                            analysis = "wilcox", 
                            sensitivity.df = metric.summary[metric.summary$METRICS %in% class, ],
                            lower.class = "SEV", upper.class = "REF",
                            lower.coef = -0.85, upper.coef = 0.85)
      }else{
        R_test <- metric.summary[metric.summary$METRICS %in% class, ]
        R_test <- R_test[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      }
      
      ms.df <- ms.df[ms.df$METRICS %in% R_test$METRICS, ]
      
      
      if(nrow(ms.df) >= 1){
        scored <- scoring(metrics.df, scores.df, ms.df, n, m.c)
        scree.df  <- scree(scored)
        
      }else{
  
          scree.df <- data.frame(NA, NA, NA, NA, NA, NA)
          names(scree.df) <- c("PERCENTILE", "PCT_REF", "PCT_DEG",
                               "SENSITIVITY", "NEW_SENSITIVITY",
                               "PERCENTILE_VALUE")
 
      }
      scree.df$THRESH <- n
      datalist[[n]] <- scree.df # add it to your list
      
      big_data <- do.call(rbind, datalist)
      big_data <-big_data[order(-big_data$SENSITIVITY,
                                -big_data$NEW_SENSITIVITY,
                                -big_data$THRESH), ]
    }
    return(big_data)
  }
  
  scree.div <- scree_classes(metrics.df, scores.df, ms_div, div)
  div.metrics <- ms_div[ms_div$SENSITIVITY >= scree.div$THRESH[1], ]
  scree.ffg <- scree_classes(metrics.df, scores.df, ms_ffg, ffg)
  ffg.metrics <- ms_ffg[ms_ffg$SENSITIVITY >= scree.ffg$THRESH[1], ]
  scree.hab <- scree_classes(metrics.df, scores.df, ms_hab, habit)
  hab.metrics <- ms_hab[ms_hab$SENSITIVITY >= scree.hab$THRESH[1], ]
  scree.tol <- scree_classes(metrics.df, scores.df, ms_tol, tol)
  tol.metrics <- ms_tol[ms_tol$SENSITIVITY >= scree.tol$THRESH[1], ]
  scree.comp <- scree_classes(metrics.df, scores.df, ms_comp, comp)
  comp.metrics <- ms_comp[ms_comp$SENSITIVITY >= scree.comp$THRESH[1], ]
  
  if(length(metrics.df[, names(metrics.df) %in% div]) > 7){
    R_div <- redundancy(metrics.df[, names(metrics.df) %in% div], analysis = "wilcox", sensitivity.df = ms_div[ms_div$METRICS %in% div, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_div <- ms_div[ms_div$METRICS %in% div, ]
    R_div <- R_div[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% ffg]) > 7){
    R_ffg <- redundancy(metrics.df[, names(metrics.df) %in% ffg], analysis = "wilcox", sensitivity.df = ms_ffg[ms_ffg$METRICS %in% ffg, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_ffg <- ms_ffg[ms_ffg$METRICS %in% ffg, ]
    R_ffg <- R_ffg[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% habit]) > 7){
    R_habit <- redundancy(metrics.df[, names(metrics.df) %in% habit], analysis = "wilcox", sensitivity.df = ms_hab[ms_hab$METRICS %in% habit, ],
                          lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_habit <- ms_hab[ms_hab$METRICS %in% habit, ]
    R_habit <- R_habit[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% tol]) > 7){
    R_tol <- redundancy(metrics.df[, names(metrics.df) %in% tol], analysis = "wilcox", sensitivity.df = ms_tol[ms_tol$METRICS %in% tol, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_tol <- ms_tol[ms_tol$METRICS %in% tol, ]
    R_tol <- R_tol[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% comp]) > 7){
    R_comp <- redundancy(metrics.df[, names(metrics.df) %in% comp], analysis = "wilcox", sensitivity.df = ms_comp[ms_comp$METRICS %in% comp, ],
                         lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_comp <- ms_comp[ms_comp$METRICS %in% comp, ]
    R_comp <- R_comp[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  R_test <- rbind(R_div, R_ffg, R_habit, R_tol, R_comp)
  
  
  final.metrics.summary <- unique(rbind(div.metrics, ffg.metrics, hab.metrics, tol.metrics, comp.metrics))
  
  final.df <- final.metrics.summary[final.metrics.summary$METRICS %in% R_test$METRICS, ]
  return(final.df)
}

#==============================================================================
#'Best Distribution 3
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param metric.summary
#'@param scores.df
#'@param m.c
#'@return Scores the samples.
#'@export


best.distribution3 <- function(metrics.df, metric.summary, scores.df, m.c,
                               Family = FALSE, master, range_var = TRUE){
  
  ms_div <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "DIVERSITY", "METRICS"]), ]
  
  ms_ffg <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "FFG", "METRICS"]), ]
  
  ms_hab <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "HABIT", "METRICS"]), ]
  
  ms_tol <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "TOLERANCE", "METRICS"]), ]
  
  ms_comp <- metric.summary[metric.summary$METRICS %in% 
                              (m.c[m.c$METRIC_CLASS %in% "COMPOSITION", "METRICS"]), ]
  
  if(Family == TRUE){
    fam.list <- paste("PCT", na.omit(unique(master$FAMILY)), sep = "_")
    ms_fam <- metric.summary[metric.summary$METRICS %in% fam.list, ]
    ms_comp <- ms_comp[!ms_comp$METRICS %in% fam.list, ]
  }
  
  if(range_var == TRUE){
    ms_div <- ms_div[(ms_div$REF_RANGE_CLASS %in% "HIGH" &
                        ms_div$REF_VARIABILITY_CLASS %in% "LOW"), ]
    #ms_ffg <- ms_ffg[(ms_ffg$REF_RANGE_CLASS %in% "HIGH" &
    #                    ms_ffg$REF_VARIABILITY_CLASS %in% "LOW"), ]
    #ms_hab <- ms_hab[(ms_hab$REF_RANGE_CLASS %in% "HIGH" &
    #                    ms_hab$REF_VARIABILITY_CLASS %in% "LOW"), ]
    #ms_tol <- ms_tol[(ms_tol$REF_RANGE_CLASS %in% "HIGH" &
    #                    ms_tol$REF_VARIABILITY_CLASS %in% "LOW"), ]
    ms_comp <- ms_comp[(ms_comp$REF_RANGE_CLASS %in% "HIGH" &
                          ms_comp$REF_VARIABILITY_CLASS %in% "LOW"), ]
    if(Family == TRUE){
      ms_fam <- ms_fam[(ms_fam$REF_RANGE_CLASS %in% "HIGH" &
                          ms_fam$REF_VARIABILITY_CLASS %in% "LOW"), ]
    } 
  }
  
  diversity <- m.c[m.c$METRIC_CLASS %in% "DIVERSITY", ]
  FFG <- m.c[m.c$METRIC_CLASS %in% "FFG", ]
  HABIT <- m.c[m.c$METRIC_CLASS %in% "HABIT", ]
  TOL <- m.c[m.c$METRIC_CLASS %in% "TOLERANCE", ]
  scores3 <- unlist(list(as.character(diversity$METRICS),
                         as.character(FFG$METRICS),
                         as.character(HABIT$METRICS),
                         as.character(TOL$METRICS),
                         "EVENT_ID", "CATEGORY", "STATION_ID",
                         "SAMPLE_NUMBER", "AGENCY_CODE", "DATE"))
  scores4 <- unlist(scores3)
  
  div <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(diversity$METRICS)))
  ffg <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(FFG$METRICS)))
  ffg <- ffg[!ffg %in% "PCT_COLLECT"]
  habit <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(HABIT$METRICS)))
  tol <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(TOL$METRICS)))
  tol <- tol[tol %in% c("PCT_INTOL_0_3", "PCT_MOD_TOL_4_6", "PCT_TOLERANT_7_10")]
  comp.exc <- unlist(list(as.character(diversity$METRICS),
                          as.character(FFG$METRICS),
                          as.character(HABIT$METRICS),
                          as.character(TOL$METRICS)))
  comp <- names(metrics.df[, !names(metrics.df) %in% comp.exc])
  comp <- comp[!comp %in% c("BIOREGION", "ABUNDANCE", "EFFECTIVE_RICH_SIMPSON", "NO_MATCH")]
  
  
  scree_classes <- function(metrics.df, scores.df, metric.summary, class){
    datalist = list()
    for (n in 50:100) {
      
      ms.df <- metric.summary[metric.summary$SENSITIVITY >= n, ]
      if(range_var == TRUE){
        ms.df <- ms.df[(ms.df$REF_RANGE_CLASS %in% "HIGH" & ms.df$REF_VARIABILITY_CLASS %in% "LOW"), ]
      }
      
      if(length(metrics.df[, names(metrics.df) %in% class]) > 7){
        R_test <- redundancy(metrics.df[, names(metrics.df) %in% class],
                             analysis = "wilcox", sensitivity.df = ms.df[ms.df$METRICS %in% class, ],
                             lower.class = "SEV", upper.class = "REF",
                             lower.coef = -0.85, upper.coef = 0.85)
      }else{
        R_test <- ms.df[ms.df$METRICS %in% class, ]
        R_test <- R_test[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      }
      
      ms.df <- ms.df[ms.df$METRICS %in% R_test$METRICS, ]
      
      
      if(nrow(ms.df) >= 1){
        scored <- scoring(metrics.df, scores.df, ms.df, n, m.c)
        scree.df  <- scree(scored)
        
      }else{
        
        scree.df <- data.frame(NA, NA, NA, NA, NA, NA)
        names(scree.df) <- c("PERCENTILE", "PCT_REF", "PCT_DEG",
                             "SENSITIVITY", "NEW_SENSITIVITY",
                             "PERCENTILE_VALUE")
        
      }
      scree.df$THRESH <- n
      datalist[[n]] <- scree.df # add it to your list
      
      big_data <- do.call(rbind, datalist)
      big_data <-big_data[order(-big_data$SENSITIVITY,
                                -big_data$NEW_SENSITIVITY,
                                -big_data$THRESH), ]
    }
    return(big_data)
  }
  
  scree.div <- scree_classes(metrics.df, scores.df, ms_div, div)
  div.metrics <- ms_div[ms_div$SENSITIVITY >= scree.div$THRESH[1], ]
  #scree.ffg <- scree_classes(metrics.df, scores.df, ms_ffg, ffg)
  ffg.metrics <- ms_ffg
  #scree.hab <- scree_classes(metrics.df, scores.df, ms_hab, habit)
  hab.metrics <- ms_hab
  #scree.tol <- scree_classes(metrics.df, scores.df, ms_tol, tol)
  tol.metrics <- ms_tol
  scree.comp <- scree_classes(metrics.df, scores.df, ms_comp, comp)
  comp.metrics <- ms_comp[ms_comp$SENSITIVITY >= scree.comp$THRESH[1], ]
  
  if(length(metrics.df[, names(metrics.df) %in% div]) > 7){
    R_div <- redundancy(metrics.df[, names(metrics.df) %in% div],
                        analysis = "wilcox", sensitivity.df = ms_div[ms_div$METRICS %in% div, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_div <- ms_div[ms_div$METRICS %in% div, ]
    R_div <- R_div[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  

    R_ffg <- ms_ffg[ms_ffg$METRICS %in% ffg, ]
    R_ffg <- R_ffg[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
 
  

    R_habit <- ms_hab[ms_hab$METRICS %in% habit, ]
    R_habit <- R_habit[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  
  
 
    R_tol <- ms_tol[ms_tol$METRICS %in% tol, ]
    R_tol <- R_tol[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
 
  
  if(length(metrics.df[, names(metrics.df) %in% comp]) > 7){
    R_comp <- redundancy(metrics.df[, names(metrics.df) %in% comp],
                         analysis = "wilcox", sensitivity.df = ms_comp[ms_comp$METRICS %in% comp, ],
                         lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_comp <- ms_comp[ms_comp$METRICS %in% comp, ]
    R_comp <- R_comp[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  R_test <- rbind(R_div, R_ffg, R_habit, R_tol, R_comp)
  
  
  final.metrics.summary <- unique(rbind(div.metrics, ffg.metrics, hab.metrics, tol.metrics, comp.metrics))
  
  final.df <- final.metrics.summary[final.metrics.summary$METRICS %in% R_test$METRICS, ]
  return(final.df)
}
#==============================================================================
#'Final Score
#'
#'@param metrics.df = data frame of metric values for each station.
#'@param final.metric.summary
#'@param scores.df
#'@param m.c
#'@return Scores the samples.
#'@export


final.score <- function(final.metrics.summary, metrics.df, scores.df, m.c){
  final.metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                   "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                   as.character(final.metrics.summary$METRICS)))
  datalist = list()
  for (i in 50:100) {
    
    if(nrow(final.metrics.summary[final.metrics.summary$SENSITIVITY >= i, ]) > 1){
      scored.df <- scoring2(metrics.df[, names(metrics.df) %in% final.metrics.vec],
                            scores.df[, names(scores.df) %in% final.metrics.vec],
                            final.metrics.summary, i, m.c)
      scree.df <- scree(scored.df)
      
    }else{
      scree.df <- data.frame(NA, NA, NA, NA, NA)
      names(scree.df) <- c("PERCENTILE", "PCT_REF", "PCT_DEG",
                           "SENSITIVITY", "NEW_SENSITIVITY")
    }
    scree.df$THRESH <- i
    datalist[[i]] <- scree.df # add it to your list
    
    final.big.data <- do.call(rbind, datalist)
    final.big.data <- final.big.data[order(-final.big.data$SENSITIVITY, -final.big.data$THRESH), ]
  }
  
  return(final.big.data)
  
}

#==============================================================================
#'Scree
#'
#'@param scored
#'@param upper.class
#'@param lower.class
#'@return Scores the samples.
#'@export

scree <- function(scored, upper.class = "REF", lower.class = "SEV"){  
  #Provide the each reference percentile value for each metric.
  quant.ref <- data.frame(PERCENTILE = quantile(scored$FINAL_SCORE, probs = seq(0, 1, by = 0.01), na.rm = TRUE))
  
  scored2 <- scored[, c("EVENT_ID", "STATION_ID", "CATEGORY",
                        "DATE", "SAMPLE_NUMBER", "AGENCY_CODE", "FINAL_SCORE")]
  quant.ref.t <- t(quant.ref)
  row.names(quant.ref.t) <- NULL 
  quant.df <- cbind(data.frame(scored2$EVENT_ID), quant.ref.t)
  names(quant.df)[1] <- "EVENT_ID" #Rename column 1
  
  long.df <- merge(scored2, quant.df, by = "EVENT_ID")
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
  
  
  long.ref[, ref.columns] <- ifelse(long.ref$FINAL_SCORE >= long.ref[, ref.columns], 1, 0)
  
  #Transform the long reference data frame to a wide format.
  melted.ref <- reshape2::melt(long.ref, id.vars = c("EVENT_ID",
                                                     "STATION_ID", "CATEGORY", "DATE",
                                                     "SAMPLE_NUMBER", "AGENCY_CODE",
                                                     "FINAL_SCORE"))
  
  #Aggregate the values (1's and 0's) by distrubance, metric, and variable (percentile).
  # The aggregation function finds the mean of the values (1's and 0's) and multiplies
  # the mean by 100. This value represents the percentage of reference sites correctly
  # identified as reference sites for a particular metric at a particular
  # threshold (percentile).
  pct.ref <- aggregate(value ~ variable, data = melted.ref,
                       function(x) mean(x) * 100)
  colnames(pct.ref) <- c("PERCENTILE", "PCT_REF")
  
  #Create a new data frame of just degraded values.
  long.deg <- long.df[long.df$CATEGORY %in% lower.class, ]
  #Column numbers can change easily. The script below specifies column "0%" to column "100%."
  # These columns represent the percentile values.
  deg.columns <-  which(colnames(long.deg) == "0%") : which(colnames(long.deg) == "100%")
  #See above for further description. Same process for identifing the number
  # of reference sites correctly identified at each threshold but the script
  # below is for correctly identified degraded sites.
  long.deg[, deg.columns] <- ifelse(long.deg$FINAL_SCORE < long.deg[, deg.columns], 1, 0)
  
  #Transform the long degraded data frame to a wide format.
  melted.deg <- reshape2::melt(long.deg, id.vars = c("EVENT_ID",
                                                     "STATION_ID", "CATEGORY", "DATE",
                                                     "SAMPLE_NUMBER", "AGENCY_CODE",
                                                     "FINAL_SCORE"))
  #Aggregate the values (1's and 0's) by distrubance, metric, and variable (percentile).
  # The aggregation function finds the mean of the values (1's and 0's) and multiplies
  # the mean by 100. This value represents the percentage of degraded sites correctly
  # identified as degraded sites for a particular metric at a particular
  # threshold (percentile).
  pct.deg <- aggregate(value ~ variable , data = melted.deg,
                       function(x) mean(x) * 100)
  #dt <- data.table::data.table(melted.deg)
  #pct.deg <- dt[, mean(value) * 100, by = list(DISTURBANCE, METRICS, variable)]
  colnames(pct.deg) <- c("PERCENTILE", "PCT_DEG")
  
  #Merge the two tables containing the percent of reference and the percent of
  # degraded correctly identified.
  merge.pct <- merge(pct.ref, pct.deg, by = c("PERCENTILE"))
  #Calculate the discrimination efficiency of each threshold.
  # DE = (Correctly identified Reference + Correctly identified Degraded) / 2
  merge.pct$SENSITIVITY <- (merge.pct$PCT_REF + merge.pct$PCT_DEG) / 2
  
  
  merge.pct$NEW_SENSITIVITY <- ((merge.pct$PCT_REF + merge.pct$PCT_DEG) / 2) - abs(merge.pct$PCT_REF - merge.pct$PCT_DEG)
  #Aggregate the table to select the best DE score.  That is the max DE score for
  # each metric.
  final.df <- merge.pct[merge.pct$NEW_SENSITIVITY %in% max(merge.pct$NEW_SENSITIVITY), ]
  if(nrow(final.df) > 1){
    final.df <- final.df[order(final.df$PERCENTILE, decreasing = TRUE), ]
    final.df <- final.df[1, ]
  } 
  final.df$PERCENTILE_VALUE <- quant.ref[row.names(quant.ref) %in% final.df$PERCENTILE, ]
  
  
  
  return(final.df)
  
}

#==============================================================================
#'Solid Metrics
#'
#'@param metrics.df
#'@param scores
#'@param sensitivity
#'@param m.c
#'@return Scores the samples.
#'@export

solid_metrics <- function(metrics.df, scores, ms, sensitivity, m.c){
  #names(ms)[names(ms) %in% "BINARY_CMA" ] <- "SENSITIVITY"
  ms2 <- ms[ms$SENSITIVITY < sensitivity, ]
  ms3 <- ms[!(ms$REF_RANGE_CLASS %in% "HIGH" & ms$REF_VARIABILITY_CLASS %in% "LOW"),]
  ms.remove <- unique(unlist(list(ms2$METRICS, ms3$METRICS)))
  
  test <- dplyr::anti_join(ms2, ms3, by = "METRICS")
  keep <- c("EVENT_ID", "CATEGORY")
  scores2 <- scores[, !names(scores) %in%  ms2$METRICS]
  scores2 <- scores2[, !names(scores2) %in%  ms3$METRICS]

  scores2[, 7:ncol(scores2)] <- apply(scores2[,7:ncol(scores2)], 2, function(x) as.numeric(as.character(x)))
  #scores2$SCORE <- apply(scores2[,7:ncol(scores2)], 1, function(x) mean(x, na.rm = TRUE))
  scores2$CATEGORY <- factor(scores2$CATEGORY, levels= c("REF", "MIN", "MOD", "SEV", "MIX"))
  
  #write.csv(scores, "scores_unp.csv", row.names = FALSE)
  metrics.df <- metrics.df[ , !names(metrics.df) %in%  ms.remove]
  m.c <- m.c[!m.c$METRICS %in% ms.remove, ]
  
  diversity <- m.c[m.c$METRIC_CLASS %in% "DIVERSITY", ]
  FFG <- m.c[m.c$METRIC_CLASS %in% "FFG", ]
  HABIT <- m.c[m.c$METRIC_CLASS %in% "HABIT", ]
  TOL <- m.c[m.c$METRIC_CLASS %in% "TOLERANCE", ]
  scores3 <- unlist(list(as.character(diversity$METRICS),
                         as.character(FFG$METRICS),
                         as.character(HABIT$METRICS),
                         as.character(TOL$METRICS),
                         "EVENT_ID", "CATEGORY", "STATION_ID",
                         "SAMPLE_NUMBER", "AGENCY_CODE", "DATE"))
  scores4 <- unlist(scores3)
  
  div <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(diversity$METRICS)))
  ffg <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(FFG$METRICS)))
  habit <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(HABIT$METRICS)))
  tol <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(TOL$METRICS)))
  comp.exc <- unlist(list(as.character(diversity$METRICS),
                          as.character(FFG$METRICS),
                          as.character(HABIT$METRICS),
                          as.character(TOL$METRICS)))
  comp <- names(metrics.df[, !names(metrics.df) %in% comp.exc])
  comp <- comp[!comp %in% c("BIOREGION", "ABUNDANCE", "EFFECTIVE_RICH_SIMPSON", "NO_MATCH")]
  
  
  #rr.ffg <- metrics.df[, names(metrics.df) %in% test.this]
  if(length(metrics.df[, names(metrics.df) %in% div]) > 7){
    R_div <- redundancy(metrics.df[, names(metrics.df) %in% div], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% div, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_div <- ms[ms$METRICS %in% div, ]
    R_div <- R_div[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% ffg]) > 7){
    R_ffg <- redundancy(metrics.df[, names(metrics.df) %in% ffg], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% ffg, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_ffg <- ms[ms$METRICS %in% ffg, ]
    R_ffg <- R_ffg[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% habit]) > 7){
    R_habit <- redundancy(metrics.df[, names(metrics.df) %in% habit], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% habit, ],
                          lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_habit <- ms[ms$METRICS %in% habit, ]
    R_habit <- R_habit[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% tol]) > 7){
    R_tol <- redundancy(metrics.df[, names(metrics.df) %in% tol], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% tol, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_tol <- ms[ms$METRICS %in% tol, ]
    R_tol <- R_tol[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% comp]) > 7){
    R_comp <- redundancy(metrics.df[, names(metrics.df) %in% comp], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% comp, ],
                         lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_comp <- ms[ms$METRICS %in% comp, ]
    R_comp <- R_comp[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  final.df <- rbind(R_div, R_ffg, R_habit, R_tol, R_comp)
  return(final.df)
}


#==============================================================================
#'Scoring
#'
#'@param metrics.df
#'@param scores
#'@param sensitivity
#'@param m.c
#'@return Scores the samples.
#'@export

scoring <- function(metrics.df, scores, ms, sensitivity, m.c){
  #names(ms)[names(ms) %in% "BINARY_CMA" ] <- "SENSITIVITY"
  ms2 <- ms[ms$SENSITIVITY >= sensitivity, ]
  ms3 <- ms2[(ms2$REF_RANGE_CLASS %in% "HIGH" & ms2$REF_VARIABILITY_CLASS %in% "LOW"),]
  ms.keep <- unique(unlist(list("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE",
                                  "SAMPLE_NUMBER", "CATEGORY",
                                  as.character(ms2$METRICS),
                                as.character(ms3$METRICS))))
  
  test <- dplyr::anti_join(ms2, ms3, by = "METRICS")
  keep <- c("EVENT_ID", "CATEGORY")
  scores2 <- scores[, names(scores) %in%  ms.keep]
  #scores2 <- scores2[, names(scores2) %in%  ms.keep]
  
  if(ncol(scores2) > 7){
    scores2[, 7:ncol(scores2)] <- apply(scores2[, 7:ncol(scores2)], 2, function(x) as.numeric(as.character(x)))
  }
  if(ncol(scores2) == 7){
    scores2[, 7] <- as.numeric(as.character(scores2[, 7]))
  }
  
  #scores2$SCORE <- apply(scores2[,7:ncol(scores2)], 1, function(x) mean(x, na.rm = TRUE))
  scores2$CATEGORY <- factor(scores2$CATEGORY, levels= c("REF", "MIN", "MOD", "SEV", "MIX"))
  
  #write.csv(scores, "scores_unp.csv", row.names = FALSE)
  metrics.df <- metrics.df[ , names(metrics.df) %in%  ms.keep]
  m.c <- m.c[m.c$METRICS %in% ms.keep, ]
  
  diversity <- m.c[m.c$METRIC_CLASS %in% "DIVERSITY", ]
  FFG <- m.c[m.c$METRIC_CLASS %in% "FFG", ]
  HABIT <- m.c[m.c$METRIC_CLASS %in% "HABIT", ]
  TOL <- m.c[m.c$METRIC_CLASS %in% "TOLERANCE", ]
  scores3 <- unlist(list(as.character(diversity$METRICS),
                         as.character(FFG$METRICS),
                         as.character(HABIT$METRICS),
                         as.character(TOL$METRICS),
                         "EVENT_ID", "CATEGORY", "STATION_ID",
                         "SAMPLE_NUMBER", "AGENCY_CODE", "DATE"))
  scores4 <- unlist(scores3)
  
  div <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(diversity$METRICS)))
  ffg <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(FFG$METRICS)))
  habit <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(HABIT$METRICS)))
  tol <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(TOL$METRICS)))
  comp.exc <- unlist(list(as.character(diversity$METRICS),
                          as.character(FFG$METRICS),
                          as.character(HABIT$METRICS),
                          as.character(TOL$METRICS)))
  comp <- names(metrics.df[, !names(metrics.df) %in% comp.exc])
  comp <- comp[!comp %in% c("BIOREGION", "ABUNDANCE", "EFFECTIVE_RICH_SIMPSON", "NO_MATCH")]
  
  
  #rr.ffg <- metrics.df[, names(metrics.df) %in% test.this]
  if(length(metrics.df[, names(metrics.df) %in% div]) > 7){
    R_div <- redundancy(metrics.df[, names(metrics.df) %in% div], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% div, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_div <- ms[ms$METRICS %in% div, ]
    R_div <- R_div[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% ffg]) > 7){
    R_ffg <- redundancy(metrics.df[, names(metrics.df) %in% ffg], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% ffg, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_ffg <- ms[ms$METRICS %in% ffg, ]
    R_ffg <- R_ffg[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% habit]) > 7){
    R_habit <- redundancy(metrics.df[, names(metrics.df) %in% habit], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% habit, ],
                          lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_habit <- ms[ms$METRICS %in% habit, ]
    R_habit <- R_habit[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% tol]) > 7){
    R_tol <- redundancy(metrics.df[, names(metrics.df) %in% tol], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% tol, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_tol <- ms[ms$METRICS %in% tol, ]
    R_tol <- R_tol[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  if(length(metrics.df[, names(metrics.df) %in% comp]) > 7){
    R_comp <- redundancy(metrics.df[, names(metrics.df) %in% comp], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% comp, ],
                         lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_comp <- ms[ms$METRICS %in% comp, ]
    R_comp <- R_comp[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  
  solid.metrics <- rbind(R_div, R_ffg, R_habit, R_tol, R_comp)
  metrics.vec <- unlist(list(names(scores2[, 1:6]), as.character(solid.metrics$METRICS)))
  # scores2 <- scores2[, !names(scores2) %in% "PCT_SCRAPE" ]
  
  
  #scores2$COMP <- apply(scores2[,names(scores2) %in% R_comp$METRICS], 1, function(x) mean(x, na.rm = TRUE))
  # scores2$diversity <- apply(scores2[, names(scores2) %in% R_div$METRICS], 1, function(x) mean(x, na.rm = TRUE))

  scores2 <- scores2[, names(scores2) %in% metrics.vec]
  
  if(is.na(table(names(scores2) %in% R_div$METRICS)[2])){
    scores2$div <- NA
  }else{
    if(table(names(scores2) %in% R_div$METRICS)[2] > 1){
      scores2$div <- apply(scores2[, names(scores2) %in% R_div$METRICS], 1, function(x) mean(as.numeric(as.character(x)), na.rm = TRUE))
    }else{
      scores2$div <- scores2[, names(scores2) %in% R_div$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_comp$METRICS)[2])){
    scores2$comp <- NA
  }else{
    if(table(names(scores2) %in% R_comp$METRICS)[2] > 1){
      scores2$comp <- apply(scores2[, names(scores2) %in% R_comp$METRICS], 1, function(x) mean(as.numeric(as.character(x)), na.rm = TRUE))
    }else{
      scores2$comp <- scores2[, names(scores2) %in% R_comp$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_ffg$METRICS)[2])){
    scores2$FFG <- NA
  }else{
    if(table(names(scores2) %in% R_ffg$METRICS)[2] > 1){
      scores2$FFG <- apply(scores2[, names(scores2) %in% R_ffg$METRICS], 1, function(x) mean(as.numeric(as.character(x)), na.rm = TRUE))
    }else{
      scores2$FFG <- scores2[, names(scores2) %in% R_ffg$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_habit$METRICS)[2])){
    scores2$HABIT <- NA
  }else{
    if(table(names(scores2) %in% R_habit$METRICS)[2] > 1){
      scores2$HABIT <- apply(scores2[, names(scores2) %in% R_habit$METRICS], 1, function(x) mean(as.numeric(as.character(x)), na.rm = TRUE))
    }else{
      scores2$HABIT <- scores2[, names(scores2) %in% R_habit$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_tol$METRICS)[2])){
    scores2$tol <- NA
  }else{
    if(table(names(scores2) %in% R_tol$METRICS)[2] > 1){
      scores2$tol <- apply(scores2[, names(scores2) %in% R_tol$METRICS], 1, function(x) mean(as.numeric(as.character(x)), na.rm = TRUE))
    }else{
      scores2$tol <- scores2[, names(scores2) %in% R_tol$METRICS]
    }
  }
  
  #scores2$HABIT <- apply(scores2[, names(scores2) %in% R_habit$METRICS], 1, function(x) mean(x, na.rm = TRUE))
  #scores2$TOL <- apply(scores2[, names(scores2) %in% R_tol$METRICS], 1, function(x) mean(x, na.rm = TRUE))
  scores2$FINAL_SCORE <- apply(scores2[, names(scores2) %in% c("comp", "div",
                                                               "FFG", "HABIT", "tol")], 1,
                               function(x) mean(x, na.rm = TRUE))
  return(scores2)
}

#==============================================================================
#'Old Scoring
#'
#'@param metrics.df
#'@param bioregion
#'@return Scores the samples.
#'@export


old_scoring <- function(metrics.df, bioregion, zero = TRUE){
  new.df <- metrics.df[metrics.df$BIOREGION %in% bioregion, ]
  new.df <- new.df[new.df$ABUNDANCE >= 70, ]
  
  bi_de <- unique(chunk_sensitivity(new.df[, !names(new.df) %in% "BIOREGION"], "REF", "SEV", "DE"))
  names(bi_de) <- c("METRICS", "DISTURBANCE", "BINARY_de")
  
  ref.df <- new.df[new.df$CATEGORY %in% "REF", ]
  if(zero == TRUE){
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
  }else{
    new.with.zero <- new.df
  }
 
  
  #new.df <- new.without.zero
  scoring_process <- function(new.df, bioregion, zero){
    ms <- metrics_summary(new.df, bioregion, zero)
    m.sum <- ms[,c("METRICS", "DISTURBANCE", "PCT_0_REF", "PCT_0_DEG",
                   "REF_MEDIAN", "BOUND_BI_CMA")]
    #test <- new.without.zero[,!names(new.without.zero) %in% ms$METRICS]
    names(m.sum)[names(m.sum) %in% "BOUND_BI_CMA"] <- "BOUND"
    m.sum$BOUND <- as.numeric(as.chracter(m.sum$BOUND))
    #m.sum$BOUND <- ifelse(m.sum$REF_MEDIAN == 0, 0, m.sum$BOUND)
    #m.sum$BOUND <- ifelse(m.sum$BOUND < 0, 0, m.sum$BOUND)
    
    long.df <- tidyr::gather(new.df, METRICS, VALUE, 8:ncol(new.df))
    
    merge.df <- merge(m.sum, long.df, by = "METRICS")
    
    dec.case.1 <- ifelse(merge.df$VALUE < merge.df$BOUND, 0,
                         ifelse(merge.df$VALUE >= merge.df$REF_MEDIAN, 1,
                                ifelse(merge.df$VALUE < merge.df$REF_MEDIAN &
                                         merge.df$VALUE >= merge.df$BOUND,
                                       (merge.df$VALUE - merge.df$BOUND) / 
                                         (merge.df$REF_MEDIAN - merge.df$BOUND), "ERROR")))
    
    inc.case.1 <- ifelse(merge.df$VALUE > merge.df$BOUND, 0,
                         ifelse(merge.df$VALUE <= merge.df$REF_MEDIAN, 1,
                                ifelse(merge.df$VALUE > merge.df$REF_MEDIAN &
                                         merge.df$VALUE <= merge.df$BOUND,
                                       (merge.df$BOUND - merge.df$VALUE) / 
                                         (merge.df$BOUND - merge.df$REF_MEDIAN), "ERROR")))
    merge.df$TEST_SCORE <- ifelse(merge.df$DISTURBANCE %in% c("EQUAL", "DECREASE"), dec.case.1,
                                  ifelse(merge.df$DISTURBANCE %in% "INCREASE", inc.case.1,
                                         "ERROR"))
    
    newer <- merge.df[, c("EVENT_ID", "CATEGORY", "STATION_ID", "SAMPLE_NUMBER",
                          "AGENCY_CODE", "DATE", "METRICS", "TEST_SCORE")]
    final.df <- tidyr::spread(newer, METRICS, TEST_SCORE)
    return(final.df)
  }
  
  final.with <- scoring_process(new.with.zero, bioregion)
  if(zero == TRUE){
    final.without <- scoring_process(new.without.zero, bioregion)
    test22 <- colSums(new.without.zero[, 8:ncol(new.without.zero)], na.rm = T) ==0
    final.df <- merge(final.with, final.without, by = c("EVENT_ID", "CATEGORY", "STATION_ID", "SAMPLE_NUMBER",
                                                        "AGENCY_CODE", "DATE"))
  }else{
    final.df <- scoring_process(new.without.zero, bioregion)
  }
 
  return(final.df)
  
}

#==============================================================================
#'Old Scoring 2
#'
#'@param metrics.df
#'@param bioregion
#'@return Scores the samples.
#'@export


old_scoring2 <- function(metrics.df, bioregion, bound.lim = TRUE){
  new.df <- metrics.df[metrics.df$BIOREGION %in% bioregion, ]
  new.df <- new.df[new.df$ABUNDANCE >= 70, ]
  
  #bi_de <- unique(chunk_sensitivity(new.df[, !names(new.df) %in% "BIOREGION"], "REF", "SEV", "DE"))
  #names(bi_de) <- c("METRICS", "DISTURBANCE", "BINARY_de")
  
  ref.df <- new.df[new.df$CATEGORY %in% "REF", ]
  #ref.zero <- ref.df
  #ref.zero[, 8:ncol(ref.zero)] <- ifelse(ref.zero[, 8:ncol(ref.zero)] == 0, 1, 0)
  
  #zero.cols <- names(ref.zero[, 8:ncol(ref.zero)])[(colSums(ref.zero[, 8:ncol(ref.zero)]) / nrow(ref.zero[, 8:ncol(ref.zero)])) >= 0.05]
  #if(any(bi_de$DISTURBANCE %in% "EQUAL")){
  #  equal.cols <- unlist(list(bi_de[bi_de$DISTURBANCE %in% "EQUAL", "METRICS"]))
  #  zero.cols <- unlist(list(zero.cols, equal.cols))
  #}
  #zero.keep <- unlist(list(names(new.df[, 1:7]), zero.cols, "ABUNDANCE"))
  
  #new.with.zero <- new.df[, !names(new.df) %in% zero.cols]
  #new.without.zero <- new.df[, names(new.df) %in% zero.keep]
  #new.without.zero[new.without.zero == 0] <- NA
  
  #new.df <- new.without.zero
  scoring_process <- function(new.df, bioregion, bound.limits){
    ms <- metrics_summary(new.df, bioregion)
    m.sum <- ms[,c("METRICS", "DISTURBANCE", "PCT_0_REF", "PCT_0_DEG",
                   "REF_MEDIAN", "BOUND_BI_CMA")]
    #test <- new.without.zero[,!names(new.without.zero) %in% ms$METRICS]
    names(m.sum)[names(m.sum) %in% "BOUND_BI_CMA"] <- "BOUND"
    m.sum$BOUND <- as.numeric(as.character(m.sum$BOUND))
    if(bound.limits == TRUE){
      m.sum$BOUND <- ifelse(m.sum$BOUND < 0, 0,
                            ifelse(m.sum$BOUND > 100, 100, as.numeric(as.character(m.sum$BOUND))))
    }
    #m.sum$BOUND <- ifelse(m.sum$REF_MEDIAN == 0, 0, m.sum$BOUND)
    #m.sum$BOUND <- ifelse(m.sum$BOUND < 0, 0, m.sum$BOUND)
    
    long.df <- tidyr::gather(new.df, METRICS, VALUE, 8:ncol(new.df))
    
    merge.df <- merge(m.sum, long.df, by = "METRICS")
    
    dec.case.1 <- ifelse(merge.df$VALUE < merge.df$BOUND, 0,
                         ifelse(merge.df$VALUE >= merge.df$REF_MEDIAN, 1,
                                ifelse(merge.df$VALUE < merge.df$REF_MEDIAN &
                                         merge.df$VALUE >= merge.df$BOUND,
                                       (merge.df$VALUE - merge.df$BOUND) / 
                                         (merge.df$REF_MEDIAN - merge.df$BOUND), "ERROR")))
    
    inc.case.1 <- ifelse(merge.df$VALUE > merge.df$BOUND, 0,
                         ifelse(merge.df$VALUE <= merge.df$REF_MEDIAN, 1,
                                ifelse(merge.df$VALUE > merge.df$REF_MEDIAN &
                                         merge.df$VALUE <= merge.df$BOUND,
                                       (merge.df$BOUND - merge.df$VALUE) / 
                                         (merge.df$BOUND - merge.df$REF_MEDIAN), "ERROR")))
    merge.df$TEST_SCORE <- ifelse(merge.df$DISTURBANCE %in% c("EQUAL", "DECREASE"), dec.case.1,
                                  ifelse(merge.df$DISTURBANCE %in% "INCREASE", inc.case.1,
                                         "ERROR"))
    
    newer <- merge.df[, c("EVENT_ID", "CATEGORY", "STATION_ID", "SAMPLE_NUMBER",
                          "AGENCY_CODE", "DATE", "METRICS", "TEST_SCORE")]
    final.df <- tidyr::spread(newer, METRICS, TEST_SCORE)
    return(final.df)
  }
  
  final.df <- scoring_process(new.df, bioregion, bound.lim)
  #final.without <- scoring_process(new.without.zero, bioregion)
  #test22 <- colSums(new.without.zero[, 8:ncol(new.without.zero)], na.rm = T) ==0
  #final.df <- merge(final.with, final.without, by = c("EVENT_ID", "CATEGORY", "STATION_ID", "SAMPLE_NUMBER",
  #                                                    "AGENCY_CODE", "DATE"))
  return(final.df)
  
}

#==============================================================================
#'Old Scoring 3
#'
#'@param metrics.df
#'@param bioregion
#'@return Scores the samples.
#'@export


old_scoring3 <- function(metrics.df, bioregion, scoring.method,
                         zero_null = TRUE, bound.lim = TRUE, metric.summary = NULL){
  new.df <- metrics.df[metrics.df$BIOREGION %in% bioregion, ]
  #new.df <- new.df[new.df$ABUNDANCE >= 70, ]
  
  ref.df <- new.df[new.df$CATEGORY %in% "REF", ]
  
  if(zero_null == TRUE){
    ref.zero <- ref.df
    ref.zero[, 8:ncol(ref.zero)] <- ifelse(ref.zero[, 8:ncol(ref.zero)] == 0, 1, 0)
    bi_de <- unique(chunk_sensitivity(new.df[, !names(new.df) %in% "BIOREGION"], "REF", "SEV", "DE"))
    names(bi_de) <- c("METRICS", "DISTURBANCE", "BINARY_de")
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
 
  #new.df <- new.without.zero
  scoring_process <- function(new.df, bioregion, bound.limits, metric.sum){
    if(!is.data.frame((metric.sum))){
      ms <- unique(metrics_summary2(new.df, bioregion))
    }else{
      ms <- metric.sum
    }
    
    #m.sum <- ms[,c("METRICS", "DISTURBANCE", "PCT_0_REF", "PCT_0_DEG",
    #               "REF_MEDIAN", "BOUND_BI_CMA")]
    m.sum <- ms[,c("METRICS", "DISTURBANCE", "REF_MEDIAN", "BOUND_BI_CMA")]
    #test <- new.without.zero[,!names(new.without.zero) %in% ms$METRICS]
    names(m.sum)[names(m.sum) %in% "BOUND_BI_CMA"] <- "BOUND"
    m.sum$BOUND <- as.numeric(as.character(m.sum$BOUND))
    if(bound.limits == TRUE){
      m.sum$BOUND <- ifelse(m.sum$BOUND < 0, 0,
                            ifelse(m.sum$BOUND > 100, 100, as.numeric(as.character(m.sum$BOUND))))
    }
    #m.sum$BOUND <- ifelse(m.sum$REF_MEDIAN == 0, 0, m.sum$BOUND)
    #m.sum$BOUND <- ifelse(m.sum$BOUND < 0, 0, m.sum$BOUND)
    
    long.df <- tidyr::gather(new.df, METRICS, VALUE, 8:ncol(new.df))
    
    merge.df <- merge(m.sum, long.df, by = "METRICS")
    
    dec.case.1 <- ifelse(merge.df$VALUE < merge.df$BOUND, 0,
                         ifelse(merge.df$VALUE >= merge.df$REF_MEDIAN, 1,
                                ifelse(merge.df$VALUE < merge.df$REF_MEDIAN &
                                         merge.df$VALUE >= merge.df$BOUND,
                                       (merge.df$VALUE - merge.df$BOUND) / 
                                         (merge.df$REF_MEDIAN - merge.df$BOUND), "ERROR")))
    
    inc.case.1 <- ifelse(merge.df$VALUE > merge.df$BOUND, 0,
                         ifelse(merge.df$VALUE <= merge.df$REF_MEDIAN, 1,
                                ifelse(merge.df$VALUE > merge.df$REF_MEDIAN &
                                         merge.df$VALUE <= merge.df$BOUND,
                                       (merge.df$BOUND - merge.df$VALUE) / 
                                         (merge.df$BOUND - merge.df$REF_MEDIAN), "ERROR")))
    merge.df$TEST_SCORE <- ifelse(merge.df$DISTURBANCE %in% c("EQUAL", "DECREASE"), dec.case.1,
                                  ifelse(merge.df$DISTURBANCE %in% "INCREASE", inc.case.1,
                                         "ERROR"))
    
    newer <- merge.df[, c("EVENT_ID", "CATEGORY", "STATION_ID", "SAMPLE_NUMBER",
                          "AGENCY_CODE", "DATE", "METRICS", "TEST_SCORE")]
    
    final.df <- tidyr::spread(newer, METRICS, TEST_SCORE)
    return(final.df)
  }
  
  if(zero_null == TRUE){
    if(ncol(new.with.zero) > 7 & ncol(new.without.zero) > 7){
      final.with <- scoring_process(new.with.zero, bioregion, bound.lim, metric.summary)
      final.without <- scoring_process(new.without.zero, bioregion, bound.lim, metric.summary)
      
      #test22 <- colSums(new.without.zero[, 8:ncol(new.without.zero)], na.rm = T) == 0
      final.df <- merge(final.with, final.without, by = c("EVENT_ID", "CATEGORY",
                                                          "STATION_ID", "SAMPLE_NUMBER",
                                                          "AGENCY_CODE", "DATE"))
    }
    
    if(ncol(new.with.zero) > 7 & ncol(new.without.zero) <= 7){
      final.df <- scoring_process(new.with.zero, bioregion, bound.lim, metric.summary)
    }
    
    if(ncol(new.with.zero) <= 7 & ncol(new.without.zero) > 7){
      final.df <- scoring_process(new.without.zero, bioregion, bound.lim, metric.summary)
    }
  }else{
    final.df <- scoring_process(new.df, bioregion, bound.lim, metric.summary)
  }
  
  return(final.df)
  
}
#==============================================================================
#'Scoring 2
#'
#'@param metrics.df
#'@param scores
#'@param sensitivity
#'@param m.c
#'@return Scores the samples.
#'@export

scoring2 <- function(metrics.df, scores, ms, m.c, Family = FALSE, master = master, redund, method = ibi.method){
  #names(ms)[names(ms) %in% "BINARY_CMA" ] <- "SENSITIVITY"
  #ms2 <- ms[ms$SENSITIVITY < sensitivity, ]
  ms2 <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(ms$METRICS)))
  keep <- c("EVENT_ID", "CATEGORY")
  #scores2 <- scores[, !names(scores) %in%  ms2$METRICS]
  scores2 <- scores
  scores2[, 7:ncol(scores2)] <- apply(scores2[,7:ncol(scores2)], 2, function(x) as.numeric(as.character(x)))
  #scores2$SCORE <- apply(scores2[,7:ncol(scores2)], 1, function(x) mean(x, na.rm = TRUE))
  scores2$CATEGORY <- factor(scores2$CATEGORY, levels= c("REF", "MIN", "MOD", "SEV", "MIX"))
  
  #write.csv(scores, "scores_unp.csv", row.names = FALSE)
  metrics.df <- metrics.df[ , names(metrics.df) %in%  ms2]
  m.c <- m.c[m.c$METRICS %in% ms2, ]
  
  diversity <- m.c[m.c$METRIC_CLASS %in% "DIVERSITY", ]
  FFG <- m.c[m.c$METRIC_CLASS %in% "FFG", ]
  HABIT <- m.c[m.c$METRIC_CLASS %in% "HABIT", ]
  TOL <- m.c[m.c$METRIC_CLASS %in% "TOLERANCE", ]
  scores3 <- unlist(list(as.character(diversity$METRICS),
                         as.character(FFG$METRICS),
                         as.character(HABIT$METRICS),
                         as.character(TOL$METRICS),
                         "EVENT_ID", "CATEGORY", "STATION_ID",
                         "SAMPLE_NUMBER", "AGENCY_CODE", "DATE"))
  scores4 <- unlist(scores3)
  
  div <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(diversity$METRICS)))
  ffg <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(FFG$METRICS)))
  habit <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(HABIT$METRICS)))
  tol <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(TOL$METRICS)))
  comp.exc <- unlist(list(as.character(diversity$METRICS),
                          as.character(FFG$METRICS),
                          as.character(HABIT$METRICS),
                          as.character(TOL$METRICS)))
  comp <- names(metrics.df[, !names(metrics.df) %in% comp.exc])
  comp <- comp[!comp %in% c("BIOREGION", "ABUNDANCE", "EFFECTIVE_RICH_SIMPSON", "NO_MATCH")]
  
  if(Family == TRUE){
    fam.metrics <- na.omit(unique(master$FAMILY))
    fam.list <- paste("PCT_", fam.metrics, sep = "")
    fam <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                              "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", fam.list))
    comp <- comp[!comp %in% fam.list]
  }
  
  if(redund == TRUE){
    #rr.ffg <- metrics.df[, names(metrics.df) %in% test.this]
    if(length(metrics.df[, names(metrics.df) %in% div]) > 7){
      R_div <- redundancy(metrics.df[, names(metrics.df) %in% div], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% div, ],
                          lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
    }else{
      R_div <- ms[ms$METRICS %in% div, ]
      R_div <- R_div[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    }
    
    if(length(metrics.df[, names(metrics.df) %in% ffg]) > 7 & !method %in% "B"){
      R_ffg <- redundancy(metrics.df[, names(metrics.df) %in% ffg], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% ffg, ],
                          lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
    }else{
      R_ffg <- ms[ms$METRICS %in% ffg, ]
      R_ffg <- R_ffg[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    }
    
    if(length(metrics.df[, names(metrics.df) %in% habit]) > 7 & !method %in% "B"){
      R_habit <- redundancy(metrics.df[, names(metrics.df) %in% habit], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% habit, ],
                            lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
    }else{
      R_habit <- ms[ms$METRICS %in% habit, ]
      R_habit <- R_habit[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    }
    
    if(length(metrics.df[, names(metrics.df) %in% tol]) > 7 & !method %in% "B"){
      R_tol <- redundancy(metrics.df[, names(metrics.df) %in% tol], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% tol, ],
                          lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
    }else{
      R_tol <- ms[ms$METRICS %in% tol, ]
      R_tol <- R_tol[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    }
    
    if(length(metrics.df[, names(metrics.df) %in% comp]) > 7){
      R_comp <- redundancy(metrics.df[, names(metrics.df) %in% comp], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% comp, ],
                           lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
    }else{
      R_comp <- ms[ms$METRICS %in% comp, ]
      R_comp <- R_comp[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    }
    
    if(Family == TRUE){
      if(length(metrics.df[, names(metrics.df) %in% fam]) > 7){
        R_fam <- redundancy(metrics.df[, names(metrics.df) %in% fam], analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% fam, ],
                            lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
      }else{
        R_fam <- ms[ms$METRICS %in% fam, ]
        R_fam <- R_fam[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      }
    }
    }else{
      R_div <- ms[ms$METRICS %in% div, ]
      R_div <- R_div[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      
      R_ffg <- ms[ms$METRICS %in% ffg, ]
      R_ffg <- R_ffg[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      
      R_habit <- ms[ms$METRICS %in% habit, ]
      R_habit <- R_habit[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      
      R_tol <- ms[ms$METRICS %in% tol, ]
      R_tol <- R_tol[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      
      R_comp <- ms[ms$METRICS %in% comp, ]
      R_comp <- R_comp[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
      
      if(Family == TRUE){
        R_fam <- ms[ms$METRICS %in% fam, ]
        R_fam <- R_fam[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    }
  }
  
 
  
  # scores2 <- scores2[, !names(scores2) %in% "PCT_SCRAPE" ]
  
  
  #scores2$COMP <- apply(scores2[,names(scores2) %in% R_comp$METRICS], 1, function(x) mean(x, na.rm = TRUE))
  # scores2$diversity <- apply(scores2[, names(scores2) %in% R_div$METRICS], 1, function(x) mean(x, na.rm = TRUE))
  
  
  if(is.na(table(names(scores2) %in% R_div$METRICS)[2])){
    scores2$div <- NA
  }else{
    if(table(names(scores2) %in% R_div$METRICS)[2] > 1){
      scores2$div <- apply(scores2[, names(scores2) %in% R_div$METRICS], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$div <- scores2[, names(scores2) %in% R_div$METRICS]
    }
  }
  
 
  if(is.na(table(names(scores2) %in% R_ffg$METRICS)[2])){
    scores2$FFG <- NA
  }else{
    if(table(names(scores2) %in% R_ffg$METRICS)[2] > 1){
      scores2$FFG <- apply(scores2[, names(scores2) %in% R_ffg$METRICS], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$FFG <- scores2[, names(scores2) %in% R_ffg$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_habit$METRICS)[2])){
    scores2$HABIT <- NA
  }else{
    if(table(names(scores2) %in% R_habit$METRICS)[2] > 1){
      scores2$HABIT <- apply(scores2[, names(scores2) %in% R_habit$METRICS], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$HABIT <- scores2[, names(scores2) %in% R_habit$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_tol$METRICS)[2])){
    scores2$tol <- NA
  }else{
    if(table(names(scores2) %in% R_tol$METRICS)[2] > 1){
      scores2$tol <- apply(scores2[, names(scores2) %in% R_tol$METRICS], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$tol <- scores2[, names(scores2) %in% R_tol$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_comp$METRICS)[2])){
    scores2$comp <- NA
  }else{
    if(table(names(scores2) %in% R_comp$METRICS)[2] > 1){
      scores2$comp <- apply(scores2[, names(scores2) %in% R_comp$METRICS], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$comp <- scores2[, names(scores2) %in% R_comp$METRICS]
    }
  }
  
  if(Family == TRUE){
    if(is.na(table(names(scores2) %in% R_fam$METRICS)[2])){
      scores2$fam <- NA
    }else{
      if(table(names(scores2) %in% R_fam$METRICS)[2] > 1){
        scores2$fam <- apply(scores2[, names(scores2) %in% R_fam$METRICS], 1, function(x) mean(x, na.rm = TRUE))
      }else{
        scores2$fam <- scores2[, names(scores2) %in% R_fam$METRICS]
      }
    }
  }
  
  
  
  #scores2$HABIT <- apply(scores2[, names(scores2) %in% R_habit$METRICS], 1, function(x) mean(x, na.rm = TRUE))
  #scores2$TOL <- apply(scores2[, names(scores2) %in% R_tol$METRICS], 1, function(x) mean(x, na.rm = TRUE))
  if(Family == FALSE){
    scores2$FINAL_SCORE <- apply(scores2[, names(scores2) %in% c("comp", "div",
                                                                 "FFG", "HABIT", "tol")], 1,
                                 function(x) mean(x, na.rm = TRUE))
  }else{
    scores2$FINAL_SCORE <- apply(scores2[, names(scores2) %in% c("comp", "div",
                                                                 "FFG", "HABIT", "tol", "fam")], 1,
                                 function(x) mean(x, na.rm = TRUE))
  }
  
  return(scores2)
}

#==============================================================================
#'Scoring 3
#'
#'@param metrics.df
#'@param scores
#'@param sensitivity
#'@param m.c
#'@return Scores the samples.
#'@export

scoring3 <- function(metrics.df, scores, ms, sensitivity, m.c){
  #names(ms)[names(ms) %in% "BINARY_CMA" ] <- "SENSITIVITY"
  #ms2 <- ms[ms$SENSITIVITY < sensitivity, ]
  ms2 <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(ms$METRICS)))
  keep <- c("EVENT_ID", "CATEGORY")
  #scores2 <- scores[, !names(scores) %in%  ms2$METRICS]
  scores2 <- scores
  scores2[, 7:ncol(scores2)] <- apply(scores2[,7:ncol(scores2)], 2, function(x) as.numeric(as.character(x)))
  #scores2$SCORE <- apply(scores2[,7:ncol(scores2)], 1, function(x) mean(x, na.rm = TRUE))
  scores2$CATEGORY <- factor(scores2$CATEGORY, levels= c("REF", "MIN", "MOD", "SEV", "MIX"))
  
  #write.csv(scores, "scores_unp.csv", row.names = FALSE)
  metrics.df <- metrics.df[ , names(metrics.df) %in%  ms2]
  m.c <- m.c[m.c$METRICS %in% ms2, ]
  
  diversity <- m.c[m.c$METRIC_CLASS %in% "DIVERSITY", ]
  FFG <- m.c[m.c$METRIC_CLASS %in% "FFG", ]
  HABIT <- m.c[m.c$METRIC_CLASS %in% "HABIT", ]
  TOL <- m.c[m.c$METRIC_CLASS %in% "TOLERANCE", ]
  scores3 <- unlist(list(as.character(diversity$METRICS),
                         as.character(FFG$METRICS),
                         as.character(HABIT$METRICS),
                         as.character(TOL$METRICS),
                         "EVENT_ID", "CATEGORY", "STATION_ID",
                         "SAMPLE_NUMBER", "AGENCY_CODE", "DATE"))
  scores4 <- unlist(scores3)
  
  div <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(diversity$METRICS)))
  ffg <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(FFG$METRICS)))
  habit <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(HABIT$METRICS)))
  tol <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(TOL$METRICS)))
  comp.exc <- unlist(list(as.character(diversity$METRICS),
                          as.character(FFG$METRICS),
                          as.character(HABIT$METRICS),
                          as.character(TOL$METRICS)))
  comp <- names(metrics.df[, !names(metrics.df) %in% comp.exc])
  comp <- comp[!comp %in% c("BIOREGION", "ABUNDANCE", "EFFECTIVE_RICH_SIMPSON", "NO_MATCH")]
  
  
  
  
  if(is.na(table(names(scores2) %in% R_div$METRICS)[2])){
    scores2$div <- NA
  }else{
    if(table(names(scores2) %in% R_div$METRICS)[2] > 1){
      scores2$div <- apply(scores2[, names(scores2) %in% div], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$div <- scores2[, names(scores2) %in% R_div$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_comp$METRICS)[2])){
    scores2$comp <- NA
  }else{
    if(table(names(scores2) %in% R_comp$METRICS)[2] > 1){
      scores2$comp <- apply(scores2[, names(scores2) %in% R_comp$METRICS], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$comp <- scores2[, names(scores2) %in% R_comp$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_ffg$METRICS)[2])){
    scores2$FFG <- NA
  }else{
    if(table(names(scores2) %in% R_ffg$METRICS)[2] > 1){
      scores2$FFG <- apply(scores2[, names(scores2) %in% R_ffg$METRICS], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$FFG <- scores2[, names(scores2) %in% R_ffg$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_habit$METRICS)[2])){
    scores2$HABIT <- NA
  }else{
    if(table(names(scores2) %in% R_habit$METRICS)[2] > 1){
      scores2$HABIT <- apply(scores2[, names(scores2) %in% R_habit$METRICS], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$HABIT <- scores2[, names(scores2) %in% R_habit$METRICS]
    }
  }
  
  if(is.na(table(names(scores2) %in% R_tol$METRICS)[2])){
    scores2$tol <- NA
  }else{
    if(table(names(scores2) %in% R_tol$METRICS)[2] > 1){
      scores2$tol <- apply(scores2[, names(scores2) %in% R_tol$METRICS], 1, function(x) mean(x, na.rm = TRUE))
    }else{
      scores2$tol <- scores2[, names(scores2) %in% R_tol$METRICS]
    }
  }
  
  #scores2$HABIT <- apply(scores2[, names(scores2) %in% R_habit$METRICS], 1, function(x) mean(x, na.rm = TRUE))
  #scores2$TOL <- apply(scores2[, names(scores2) %in% R_tol$METRICS], 1, function(x) mean(x, na.rm = TRUE))
  scores2$FINAL_SCORE <- apply(scores2[, names(scores2) %in% c("comp", "div",
                                                               "FFG", "HABIT", "tol")], 1,
                               function(x) mean(x, na.rm = TRUE))
  return(scores2)
}

#==============================================================================
#'Balanced Discrimination Efficiency
#'
#'@param scored.df = final scored index.
#'@param bioregion = the bioregion of interest.
#'@return Update
#'@export

bde <- function(scored.df, bioregion){
  datalist = list()
  for(i in 1:100){
    ref <- scored.df[scored.df$CATEGORY %in% "REF", ]
    deg <- scored.df[scored.df$CATEGORY %in% "SEV", ]
    
    pct.ref <- ifelse(ref$FINAL_SCORE * 100 >= i, 1, 0)
    pct.deg <- ifelse(deg$FINAL_SCORE * 100 < i, 1, 0)
    
    final.df <- data.frame(BIOREGION = bioregion)
    final.df$PCT_REF <- (sum(pct.ref) / length(pct.ref)) * 100
    final.df$PCT_DEG <- (sum(pct.deg) / length(pct.deg)) * 100
    final.df$BDE <- ((final.df$PCT_REF + final.df$PCT_DEG) / 2) - abs(final.df$PCT_REF - final.df$PCT_DEG)
    final.df$CE <- ((final.df$PCT_REF + final.df$PCT_DEG) / 2)
    final.df$THRESHOLD <- i
    datalist[[i]] <- final.df # add it to your list
  }
  
  big_data <- do.call(rbind, datalist)
  big_data <- big_data[order(-big_data$BDE, -big_data$THRESHOLD), ]
}

#==============================================================================
#'Redundancy analysis by metric class
#'
#'@param metrics.df = a data frame with raw metric scores
#'@param metric.summary = metric summary data frame
#'@param m.c = data frame containing metric classes
#'@param method = the ibi.method specified
#'@return Update
#'@export


metric_class_redund <- function(metrics.df, metric.summary, m.c, method = ibi.method){
  
  ms_div <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "DIVERSITY", "METRICS"]), ]
  
  ms_ffg <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "FFG", "METRICS"]), ]
  
  ms_hab <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "HABIT", "METRICS"]), ]
  
  ms_tol <- metric.summary[metric.summary$METRICS %in% 
                             (m.c[m.c$METRIC_CLASS %in% "TOLERANCE", "METRICS"]), ]
  
  ms_comp <- metric.summary[metric.summary$METRICS %in% 
                              (m.c[m.c$METRIC_CLASS %in% "COMPOSITION", "METRICS"]), ]
  #============================================================================
  diversity <- m.c[m.c$METRIC_CLASS %in% "DIVERSITY", ]
  FFG <- m.c[m.c$METRIC_CLASS %in% "FFG", ]
  HABIT <- m.c[m.c$METRIC_CLASS %in% "HABIT", ]
  TOL <- m.c[m.c$METRIC_CLASS %in% "TOLERANCE", ]
  #============================================================================
  div <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(diversity$METRICS)))
  ffg <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(FFG$METRICS)))
  habit <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(HABIT$METRICS)))
  tol <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                     "SAMPLE_NUMBER", "AGENCY_CODE", "DATE", as.character(TOL$METRICS)))
  comp.exc <- unlist(list(as.character(diversity$METRICS),
                          as.character(FFG$METRICS),
                          as.character(HABIT$METRICS),
                          as.character(TOL$METRICS)))
  comp <- names(metrics.df[, !names(metrics.df) %in% comp.exc])
  comp <- comp[!comp %in% c("BIOREGION", "ABUNDANCE", "EFFECTIVE_RICH_SIMPSON", "NO_MATCH")]
  #============================================================================
  
  if(length(metrics.df[, names(metrics.df) %in% div]) > 7){
    R_div <- redundancy(metrics.df[, names(metrics.df) %in% div], analysis = "wilcox",
                        sensitivity.df = ms_div[ms_div$METRICS %in% div, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_div <- ms_div[ms_div$METRICS %in% div, ]
    R_div <- R_div[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  #============================================================================
  if(length(metrics.df[, names(metrics.df) %in% ffg]) > 7 & method %in% "A"){
    R_ffg <- redundancy(metrics.df[, names(metrics.df) %in% ffg], analysis = "wilcox", sensitivity.df = ms_ffg[ms_ffg$METRICS %in% ffg, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_ffg <- ms_ffg[ms_ffg$METRICS %in% ffg, ]
    R_ffg <- R_ffg[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  #============================================================================
  if(length(metrics.df[, names(metrics.df) %in% habit]) > 7 & method %in% "A"){
    R_habit <- redundancy(metrics.df[, names(metrics.df) %in% habit], analysis = "wilcox", sensitivity.df = ms_hab[ms_hab$METRICS %in% habit, ],
                          lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_habit <- ms_hab[ms_hab$METRICS %in% habit, ]
    R_habit <- R_habit[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  #============================================================================
  if(length(metrics.df[, names(metrics.df) %in% tol]) > 7 & method %in% "A"){
    R_tol <- redundancy(metrics.df[, names(metrics.df) %in% tol], analysis = "wilcox", sensitivity.df = ms_tol[ms_tol$METRICS %in% tol, ],
                        lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_tol <- ms_tol[ms_tol$METRICS %in% tol, ]
    R_tol <- R_tol[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  #============================================================================
  if(length(metrics.df[, names(metrics.df) %in% comp]) > 7){
    R_comp <- redundancy(metrics.df[, names(metrics.df) %in% comp], analysis = "wilcox", sensitivity.df = ms_comp[ms_comp$METRICS %in% comp, ],
                         lower.class = "SEV", upper.class = "REF", lower.coef = -0.85, upper.coef = 0.85)
  }else{
    R_comp <- ms_comp[ms_comp$METRICS %in% comp, ]
    R_comp <- R_comp[, c("METRICS", "SENSITIVITY", "DISTURBANCE")]
  }
  #============================================================================
  final.df <- rbind(R_div, R_ffg, R_habit, R_tol, R_comp)
  
  return(final.df)
}













