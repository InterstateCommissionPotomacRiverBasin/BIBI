#==============================================================================
#Metric Sensitivity
#==============================================================================
#'Metric Sensitivity
#'
#'@param metrics.df = data frame of metric values for each station
#'@param upper.class = The site classification that represents better
#'environmental conditions.
#'@param lower.class = The site classification that represents the degraded
#'environmental conditions.
#'@param method = the sensitivity function to be used during the assessment.
#'@return Determines the threshold at which a metric best categorizes
#'reference and degraded stations.
#'@export
#'
sensitivity <- function(metrics.df, upper.class, lower.class, method = "ODE"){
  if("PCT_UNIDENTIFIED"  %in% names(metrics.df)){
    metrics.df <- metrics.df[, !grepl("PCT_UNIDENTIFIED", names(metrics.df))]
  }
  
  if(any(grepl(".y", names(metrics.df)))){
    metrics.df <- metrics.df[, !grepl(".y", names(metrics.df))]
  }

  if(any(grepl(".x", names(metrics.df)))){
    names(metrics.df) <- gsub(".x", "", names(metrics.df))
  }
  #Create new data frames specific for Degraded and Reference sites
  deg.df <- metrics.df[metrics.df$CATEGORY %in% lower.class, ]
  ref.df <- metrics.df[metrics.df$CATEGORY %in% upper.class, ]

  if(method == "BARBOUR"){
    final.df <- barbour(metrics.df, ref.df, deg.df)
  }

  #Calculate the median values for the reference and degraded distributions.
  if(length(ref.df) < 8){
    ref_50 <- quantile(ref.df[, 7], 0.50, na.rm = TRUE)
    deg_50 <- quantile(deg.df[, 7], 0.50, na.rm = TRUE)
    
    #Provide the each reference percentile value for each metric.
    quant.ref <- data.frame(quantile(ref.df[, 7], probs = seq(0, 1, by = 0.01), na.rm = TRUE))
    colnames(quant.ref) <- colnames(ref.df)[7]
    #Create a column listing all of the metrics and join the reference percentile values
    quant.df <- cbind(data.frame(colnames(metrics.df[7])), t(quant.ref))
    names(quant.df)[1] <- "METRICS" #Rename column 1
    
    
    
    #quant.df$DISTURBANCE <- ifelse(ref_50 > deg_50, "DECREASE",
    #                              ifelse(ref_50 < deg_50, "INCREASE", "EQUAL"))
    #quant.df <- quant.df[!(quant.df$DISTURBANCE %in% "EQUAL"), ]
    #quant.df <- quant.df[rowSums(quant.df[, 2:102]) > 0, ]
    
    #Create new data frames specific for Degraded and Reference sites
    severe.df <- metrics.df[metrics.df$CATEGORY %in% lower.class, ]
    reference.df <- metrics.df[metrics.df$CATEGORY %in% upper.class, ]
    #Calculate the median values for the reference and degraded distributions.
    reference_50 <- quantile(reference.df[, 7],  0.50, na.rm = TRUE)
    severe_50 <- quantile(severe.df[, 7],  0.50, na.rm = TRUE)
    #Insert a column to suggest how the metric reacts to disturbance. If the reference median
    # is greater than the degraded median, the metric decreases with distrubance. If the reference
    # median is less than the degraded median, the metric increases with disturbance.  If the
    # medians are equal, equal is return to indicate that this metric shows no distinction between
    # reference and degraded contions.
    quant.df$DISTURBANCE <- ifelse(reference_50 > severe_50, "DECREASE",
                                   ifelse(reference_50 < severe_50, "INCREASE", "EQUAL"))
    
  }
  
  if(length(ref.df) >= 8){
    ref_50 <- sapply(ref.df[, 7:ncol(ref.df)], quantile, 0.50, na.rm = TRUE)
    deg_50 <- sapply(deg.df[, 7:ncol(deg.df)], quantile, 0.50, na.rm = TRUE)
    
    #Provide the each reference percentile value for each metric.
    quant.ref <- data.frame(apply(ref.df[, 7:ncol(ref.df)], 2, function(x){
      quantile(x, probs = seq(0, 1, by = 0.01), na.rm = TRUE)
    } ))
    #Create a column listing all of the metrics and join the reference percentile values
    quant.df <- cbind(data.frame(colnames(metrics.df[7:ncol(metrics.df)])), t(quant.ref))
    names(quant.df)[1] <- "METRICS" #Rename column 1
    
    
    
    #quant.df$DISTURBANCE <- ifelse(ref_50 > deg_50, "DECREASE",
    #                              ifelse(ref_50 < deg_50, "INCREASE", "EQUAL"))
    #quant.df <- quant.df[!(quant.df$DISTURBANCE %in% "EQUAL"), ]
    #quant.df <- quant.df[rowSums(quant.df[, 2:102]) > 0, ]
    
    #Create new data frames specific for Degraded and Reference sites
    severe.df <- metrics.df[metrics.df$CATEGORY == lower.class, ]
    reference.df <- metrics.df[metrics.df$CATEGORY == upper.class, ]
    #Calculate the median values for the reference and degraded distributions.
    reference_50 <- sapply(reference.df[, 7:ncol(reference.df)], quantile, 0.50, na.rm = TRUE)
    severe_50 <- sapply(severe.df[, 7:ncol(severe.df)], quantile, 0.50, na.rm = TRUE)
    #Insert a column to suggest how the metric reacts to disturbance. If the reference median
    # is greater than the degraded median, the metric decreases with distrubance. If the reference
    # median is less than the degraded median, the metric increases with disturbance.  If the
    # medians are equal, equal is return to indicate that this metric shows no distinction between
    # reference and degraded contions.
    quant.df$DISTURBANCE <- ifelse(reference_50 > severe_50, "DECREASE",
                                   ifelse(reference_50 < severe_50, "INCREASE", "EQUAL"))
    
  }


  if(method == "DE"){
    final.df <- d_e(deg.df, quant.df)
  }

  if(method == "ODE"){
    final.df <- ode(metrics.df, quant.df, upper.class, lower.class,
                    ref.df, quant.ref)
  }
  
  if(method == "CMA"){
    final.df <- cma(metrics.df, quant.df, upper.class, lower.class,
                    ref.df, quant.ref)
  }
  
  if(method == "SSE"){
    final.df <- sse(metrics.df, quant.df, upper.class, lower.class,
                    ref.df, quant.ref)
  }

  return(final.df)
}

#==============================================================================
#'
#'Chunk Sensitivity
#'
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@param upper.class = the site class that represents the better condition.
#'@param lower.class = the site class that represents the poorer condition.
#'@param method = the sensitivity function to be used during the assessment.
#'@return Determines the threshold at which a metric best categorizes
#'two defined environmental conditions.
#'@export

chunk_sensitivity <- function(metrics.df, upper.class = "REF", lower.class = "SEV", method){
  
  metrics.list <- break.me(metrics.df, 100, 6)
  #============================================================================
  datalist = list()
  for(j in 1:length(metrics.list)){
    sub.metrics <- metrics.list[[j]]
    de.thresh <- sensitivity(sub.metrics, upper.class, lower.class, method)
    datalist[[j]] <- de.thresh
  }
  #============================================================================
  final.df <- do.call(rbind, datalist)
  
  return(final.df)

}

#==============================================================================
#'Pairwise Sensitivity
#'
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@param method = the sensitivity function to be used during the assessment.
#'@return Determines the threshold at which a metric best categorizes
#'two defined environmental conditions.
#'@export

pairwise_sensitivity <- function(metrics.df, method){
  #rn.df <- chunk_sensitivity(metrics.df, "REF", "NEAR", method)
  #nmin.df <- chunk_sensitivity(metrics.df, "NEAR", "MIN", method)
  rm.df <- chunk_sensitivity(metrics.df, "REF", "MIN", method)
  mm.df <- chunk_sensitivity(metrics.df, "MIN", "MOD", method)
  ms.df <- chunk_sensitivity(metrics.df, "MOD", "SEV", method)
  rs.df <- chunk_sensitivity(metrics.df, "REF", "SEV", method)


  #rs.df <- chunk_sensitivity(metrics.df, "NEAR", "SEV", method)
  #rn.df <- chunk_sensitivity(metrics.df, "REF", "NEAR", method)
  #nmin.df <- chunk_sensitivity(metrics.df, "REF", "MIN", method)
  #mm.df <- chunk_sensitivity(metrics.df, "REF", "MOD", method)
  #ms.df <- chunk_sensitivity(metrics.df, "REF", "SEV", method)


  if(method %in% c("ODE", "SSE")){
    #rn.df <- rn.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    #names(rn.df) <- c("METRICS", "SENSITIVITY_REF_NEAR", "THRESH_REF_NEAR")

    #nmin.df <- nmin.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    #names(nmin.df) <- c("METRICS", "SENSITIVITY_NEAR_MIN", "THRESH_NEAR_MIN")
    
    rm.df <- rm.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    names(rm.df) <- c("METRICS", "SENSITIVITY_REF_MIN", "THRESH_REF_MIN")
    
    mm.df <- mm.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    names(mm.df) <- c("METRICS", "SENSITIVITY_MIN_MOD", "THRESH_MIN_MOD")

    ms.df <- ms.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    names(ms.df) <- c("METRICS", "SENSITIVITY_MOD_SEV", "THRESH_MOD_SEV")

    rs.df <- rs.df[,c("METRICS", "SENSITIVITY", "THRESHOLD",
                      "DISTURBANCE")]
    names(rs.df) <- c("METRICS", "SENSITIVITY_REF_SEV",
                      "THRESH_REF_SEV", "DISTURBANCE")

  }

  if(method %in% c("DE", "BARBOUR")){
    rn.df <- rn.df[, c("METRICS", "SENSITIVITY")]
    names(rn.df) <- c("METRICS", "SENSITIVITY_REF_NEAR")

    nmin.df <- nmin.df[, c("METRICS", "SENSITIVITY")]
    names(nmin.df) <- c("METRICS", "SENSITIVITY_NEAR_MIN")

    mm.df <- mm.df[, c("METRICS", "SENSITIVITY")]
    names(mm.df) <- c("METRICS", "SENSITIVITY_MIN_MOD")

    ms.df <- ms.df[, c("METRICS", "SENSITIVITY")]
    names(ms.df) <- c("METRICS", "SENSITIVITY_MOD_SEV")

    rs.df <- rs.df[,c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    names(rs.df) <- c("METRICS", "SENSITIVITY_REF_SEV", "DISTURBANCE")
  }


  if(method %in% c("ODE", "SSE", "CMA")){
    m3 <- cbind(rm.df, 
                #rn.df, nmin.df[, c(2,3)], 
                mm.df[, c(2,3)],
                ms.df[, c(2,3)], rs.df[, 2:4])
  }

  if(method %in% c("DE", "BARBOUR")){
    m3 <- cbind(rn.df, nmin.df[, c(2)], mm.df[, c(2)], ms.df[, c(2)], rs.df[, 2:3])
    names(m3) <- c("METRICS",
                   #"SENSITIVITY_REF_NEAR", "SENSITIVITY_NEAR_MIN",
                   "SENSITIVITY_REF_MIN",
                   "SENSITIVITY_MIN_MOD","SENSITIVITY_MOD_SEV",
                   "SENSITIVITY_REF_SEV", "DISTURBANCE")
  }

  m3$SENSITIVITY <- (rowSums(m3[, c("SENSITIVITY_REF_MIN",
                                    #"SENSITIVITY_REF_NEAR", "SENSITIVITY_NEAR_MIN",
                                   "SENSITIVITY_MIN_MOD","SENSITIVITY_MOD_SEV",
                                   "SENSITIVITY_REF_SEV")])) / 4

  return(m3)
}

#==============================================================================
#'Pairwise Sensitivity 2
#'
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@param method = the sensitivity function to be used during the assessment.
#'@return Determines the threshold at which a metric best categorizes
#'two defined environmental conditions.
#'@export

pairwise_sensitivity2 <- function(metrics.df, method){
  ns.df <- chunk_sensitivity(metrics.df, "NEAR", "SEV", method)
  rm.df <- chunk_sensitivity(metrics.df, "REF", "MOD", method)
  rs.df <- chunk_sensitivity(metrics.df, "REF", "SEV", method)
  
  
  #rs.df <- chunk_sensitivity(metrics.df, "NEAR", "SEV", method)
  #rn.df <- chunk_sensitivity(metrics.df, "REF", "NEAR", method)
  #nmin.df <- chunk_sensitivity(metrics.df, "REF", "MIN", method)
  #mm.df <- chunk_sensitivity(metrics.df, "REF", "MOD", method)
  #ms.df <- chunk_sensitivity(metrics.df, "REF", "SEV", method)
  
  
  if(method %in% c("ODE", "SSE")){
    ns.df <- ns.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    names(ns.df) <- c("METRICS", "SENSITIVITY_NEAR_SEV", "THRESH_NEAR_SEV")
    
    rm.df <- rm.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    names(rm.df) <- c("METRICS", "SENSITIVITY_REF_MOD", "THRESH_REF_MOD")
    
    rs.df <- rs.df[,c("METRICS", "SENSITIVITY", "THRESHOLD",
                      "DISTURBANCE")]
    names(rs.df) <- c("METRICS", "SENSITIVITY_REF_SEV",
                      "THRESH_REF_SEV", "DISTURBANCE")
    
  }
  
  if(method %in% c("DE", "BARBOUR")){
    ns.df <- ns.df[, c("METRICS", "SENSITIVITY")]
    names(ns.df) <- c("METRICS", "SENSITIVITY_NEAR_SEV")
    
    rm.df <- rm.df[, c("METRICS", "SENSITIVITY")]
    names(rm.df) <- c("METRICS", "SENSITIVITY_REF_MOD")
    
    rs.df <- rs.df[,c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    names(rs.df) <- c("METRICS", "SENSITIVITY_REF_SEV", "DISTURBANCE")
  }
  
  
  if(method %in% c("ODE", "SSE")){
    m3 <- cbind(ns.df, rm.df[, c(2,3)], rs.df[, 2:4])
  }
  
  if(method %in% c("DE", "BARBOUR")){
    m3 <- cbind(rs.df, rm.df[, c(2)], rs.df[, 2:3])
    names(m3) <- c("METRICS", "SENSITIVITY_NEAR_SEV", "SENSITIVITY_REF_MOD",
                   "SENSITIVITY_REF_SEV", "DISTURBANCE")
  }
  
  m3$SENSITIVITY <- rowSums(m3[, c("SENSITIVITY_NEAR_SEV", "SENSITIVITY_REF_MOD",
                                    "SENSITIVITY_REF_SEV")]) / 3

  return(m3)
}

#==============================================================================
#'Pairwise Sensitivity 3
#'
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@param method = the sensitivity function to be used during the assessment.
#'@return Determines the threshold at which a metric best categorizes
#'two defined environmental conditions.
#'@export

pairwise_sensitivity3 <- function(metrics.df, method){
  #rn.df <- chunk_sensitivity(metrics.df, "REF", "NEAR", method)
  #nmin.df <- chunk_sensitivity(metrics.df, "NEAR", "MIN", method)
  rm.df <- chunk_sensitivity(metrics.df, "REF", "MIN", method)
  mm.df <- chunk_sensitivity(metrics.df, "MIN", "MOD", method)
  ms.df <- chunk_sensitivity(metrics.df, "MOD", "SEV", method)
  rs.df <- chunk_sensitivity(metrics.df, "REF", "SEV", method)
  
  
  #rs.df <- chunk_sensitivity(metrics.df, "NEAR", "SEV", method)
  #rn.df <- chunk_sensitivity(metrics.df, "REF", "NEAR", method)
  #nmin.df <- chunk_sensitivity(metrics.df, "REF", "MIN", method)
  #mm.df <- chunk_sensitivity(metrics.df, "REF", "MOD", method)
  #ms.df <- chunk_sensitivity(metrics.df, "REF", "SEV", method)
  
  
  if(method %in% c("ODE", "SSE")){
    #rn.df <- rn.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    #names(rn.df) <- c("METRICS", "SENSITIVITY_REF_NEAR", "THRESH_REF_NEAR")
    
    #nmin.df <- nmin.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    #names(nmin.df) <- c("METRICS", "SENSITIVITY_NEAR_MIN", "THRESH_NEAR_MIN")
    
    rm.df <- rm.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    names(rm.df) <- c("METRICS", "SENSITIVITY_REF_MIN", "THRESH_REF_MIN")
    
    mm.df <- mm.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    names(mm.df) <- c("METRICS", "SENSITIVITY_MIN_MOD", "THRESH_MIN_MOD")
    
    ms.df <- ms.df[, c("METRICS", "SENSITIVITY", "THRESHOLD")]
    names(ms.df) <- c("METRICS", "SENSITIVITY_MOD_SEV", "THRESH_MOD_SEV")
    
    rs.df <- rs.df[,c("METRICS", "SENSITIVITY", "THRESHOLD",
                      "DISTURBANCE")]
    names(rs.df) <- c("METRICS", "SENSITIVITY_REF_SEV",
                      "THRESH_REF_SEV", "DISTURBANCE")
    
  }
  
  if(method %in% c("DE", "BARBOUR")){
    rn.df <- rn.df[, c("METRICS", "SENSITIVITY")]
    names(rn.df) <- c("METRICS", "SENSITIVITY_REF_NEAR")
    
    nmin.df <- nmin.df[, c("METRICS", "SENSITIVITY")]
    names(nmin.df) <- c("METRICS", "SENSITIVITY_NEAR_MIN")
    
    mm.df <- mm.df[, c("METRICS", "SENSITIVITY")]
    names(mm.df) <- c("METRICS", "SENSITIVITY_MIN_MOD")
    
    ms.df <- ms.df[, c("METRICS", "SENSITIVITY")]
    names(ms.df) <- c("METRICS", "SENSITIVITY_MOD_SEV")
    
    rs.df <- rs.df[,c("METRICS", "SENSITIVITY", "DISTURBANCE")]
    names(rs.df) <- c("METRICS", "SENSITIVITY_REF_SEV", "DISTURBANCE")
  }
  
  
  if(method %in% c("ODE", "SSE")){
    m3 <- cbind(rm.df, 
                #rn.df, nmin.df[, c(2,3)], 
                mm.df[, c(2,3)],
                ms.df[, c(2,3)], rs.df[, 2:4])
  }
  
  if(method %in% c("DE", "BARBOUR")){
    m3 <- cbind(rn.df, nmin.df[, c(2)], mm.df[, c(2)], ms.df[, c(2)], rs.df[, 2:3])
    names(m3) <- c("METRICS",
                   #"SENSITIVITY_REF_NEAR", "SENSITIVITY_NEAR_MIN",
                   "SENSITIVITY_REF_MIN",
                   "SENSITIVITY_MIN_MOD","SENSITIVITY_MOD_SEV",
                   "SENSITIVITY_REF_SEV", "DISTURBANCE")
  }
  
  m3$SENSITIVITY <- (rowSums(m3[, c("SENSITIVITY_REF_MIN",
                                    #"SENSITIVITY_REF_NEAR", "SENSITIVITY_NEAR_MIN",
                                    "SENSITIVITY_MIN_MOD","SENSITIVITY_MOD_SEV",
                                    "SENSITIVITY_REF_SEV")]) + m3$SENSITIVITY_REF_SEV) / 4
  
  return(m3)
}



#==============================================================================
#'Range and Variability Test
#'
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@return Tests that the range of the reference condition is not too low and
#'that variability is not too high.
#'@export

range_variability <- function(metrics.df){
  if("NO_MATCH" %in% names(metrics.df)){
    metrics.df <- metrics.df[, !(names(metrics.df) %in% "NO_MATCH")]
  } 
  if("EFFECTIVE_RICH_SIMPSON" %in% names(metrics.df)){
    metrics.df <- metrics.df[, !(names(metrics.df) %in% "EFFECTIVE_RICH_SIMPSON")]
  } 
  
  ref <- metrics.df[metrics.df$CATEGORY %in% "REF", ]
  
  if("PIELOU" %in% names(ref)){
    ref$PIELOU <- ref$PIELOU * 100
  } 
  if("HURLBERTS_PIE" %in% names(ref)){
    ref$HURLBERTS_PIE <- ref$HURLBERTS_PIE * 100
  } 
  if("SIMPSONS" %in% names(ref)){
    ref$SIMPSONS <- ref$SIMPSONS * 100
  } 
  
  if(ncol(ref) > 7){
    df <- data.frame(METRICS = names(ref[, 7:ncol(ref)]))
    #df <- merge(df, sensitivity.df[, c("METRICS", "DISTURBANCE")], by = "METRICS", all.x = TRUE, sort = FALSE)
    df$MIN <- apply(ref[, 7:ncol(ref)], 2, function(x) quantile(x, probs = 0.05, na.rm = TRUE))
    df$MAX <- apply(ref[, 7:ncol(ref)], 2, function(x) quantile(x, probs = 0.95, na.rm = TRUE))
  }
  
  if(ncol(ref) == 7){
    df <- data.frame(METRICS = names(ref)[7])
    df$MIN <- quantile(ref[, 7], probs = 0.05, na.rm = TRUE)
    df$MAX <- quantile(ref[, 7], probs = 0.95, na.rm = TRUE)
  }

  df$DIFF <- abs(df$MIN - df$MAX)
  
  pct.m <- paste(c("PCT", "PIELOU", "GOLD", "SIMPSON", "HURLBERT"), collapse = "|")
  rich.m <- paste(c("RICH", "BECK"), collapse = "|")
  div.m <- paste(c("SHANNON", "MENHINICKS", "MARGALEFS"), collapse = "|")
  tol.m <- paste(c("HBI", "ASPT"), collapse = "|")
  
  df$RANGE <-  ifelse(grepl(pct.m, df$METRIC) & df$DIFF <= 10, "LOW",
                      ifelse(grepl(pct.m, df$METRIC) & df$DIFF > 10, "HIGH",
                             ifelse(grepl(div.m, df$METRIC) & df$DIFF < 1, "LOW",
                                    ifelse(grepl(div.m, df$METRIC) & df$DIFF >= 1, "HIGH",
                                           ifelse(grepl(tol.m, df$METRIC) & df$DIFF < 2, "LOW",
                                                  ifelse(grepl(tol.m, df$METRIC) & df$DIFF >= 2, "HIGH",
                                                         ifelse(grepl(rich.m, df$METRIC) & df$DIFF < 3, "LOW",
                                                                ifelse(grepl(rich.m, df$METRIC) & df$DIFF >= 3, "HIGH",
                                                                       ifelse(!grepl(pct.m, df$METRICS) & !grepl(rich.m, df$METRICS) &
                                                                                !grepl(tol.m, df$METRICS) & !grepl(div.m, df$METRICS), "Not Measured", "ERROR")))))))))
  if(ncol(ref) > 7){
    df$Q25 <- round(apply(ref[, 7:ncol(ref)], 2, function(x) quantile(x, probs = 0.25, na.rm = TRUE)), 0)
    df$Q75 <- round(apply(ref[, 7:ncol(ref)], 2, function(x) quantile(x, probs = 0.75, na.rm = TRUE)), 0)
  }
  
  if(ncol(ref) == 7){
    df$Q25 <- round(quantile(ref[, 7], probs = 0.25, na.rm = TRUE), 0)
    df$Q75 <- round(quantile(ref[, 7], probs = 0.75, na.rm = TRUE), 0)
  }
  
  df$Q_DIFF <- df$Q75 - df$Q25
  #df$VARIABILITY <- ifelse((df$Q_DIFF) == 0, "LOW",
  #                         ifelse((df$Q_DIFF / df$Q25) > 1, "HIGH",
  #                                ifelse((df$Q_DIFF / df$Q25) <= 1, "LOW", "ERROR")))
  df$VARIABILITY <- ifelse((df$Q_DIFF) == 0, "LOW",
                            ifelse((df$Q_DIFF / df$Q25) > 3, "HIGH",
                                   ifelse((df$Q_DIFF / df$Q25) <= 3, "LOW", "ERROR")))
  
  return(df)
  
}

#==============================================================================
#'Summary of Metric Tests
#'
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@param bioregion = the bioregion to perform the analysis.
#'@return Summarizes multiple metric tests into a single table.
#'@export

metrics_summary <- function(metrics.df, bioregion, de.method = "CMA"){
  metrics.df <- metrics.df[metrics.df$ABUNDANCE >= 70, ]
  metrics.df <- metrics.df[metrics.df$BIOREGION %in% bioregion, ]
  metrics.df <- metrics.df[, !names(metrics.df) %in% "BIOREGION"]
  #pair.cma <- unique(pairwise_sensitivity(metrics.df, de.method))
  #names(pair.cma)[names(pair.cma) %in% "SENSITIVITY"] <- "PAIRWISE_CMA"
  bi.cma <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", de.method))
  names(bi.cma) <- c("METRICS", "DISTURBANCE", "BINARY_CMA",
                     "PRECENTILE_BINARY_CMA",
                     "PCT_REF_BI_CMA", "PCT_DEG_BI_CMA",
                     "REF_MEDIAN", "THRESHOLD_BI_CMA", "BOUND_BI_CMA")
  bi.de <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "DE"))
  names(bi.de) <- c("METRICS", "DISTURBANCE", "BINARY_DE")
  bi_barbour <- unique(chunk_sensitivity(metrics.df, "REF", "SEV", "BARBOUR"))
  names(bi_barbour) <- c("METRICS", "DISTURBANCE", "BINARY_BARBOUR")
  range.var <- unique(range_variability(metrics.df))
  names(range.var) <- c("METRICS", "REF_MIN", "REF_MAX", "REF_RANGE_VALUE",
                        "REF_RANGE_CLASS", "REF_Q25", "REF_Q75",
                        "REF_VARIABILITY_VALUE", "REF_VARIABILITY_CLASS")
  zero.inflate <- zero_inflate(metrics.df, bi_barbour)
  final.df <- plyr::join_all(list(#pair.cma,
                                  bi.cma[, c(1, 3:9)], 
                                  bi.de[, c(1, 3)],
                                  bi_barbour[, c(1, 3)], 
                                  range.var, zero.inflate), "METRICS")
  
  final.df$QUALITY <- ifelse(#final.df$SENSITIVITY >= 70 &
                             final.df$BINARY_CMA >= 70 &
                             final.df$BINARY_BARBOUR >= 2 &
                             final.df$REF_RANGE_CLASS %in% "HIGH" &
                             final.df$REF_VARIABILITY_CLASS %in% "LOW" &
                             final.df$ZERO_INFLATE %in% "GOOD", "HIGH",
                             ifelse(#final.df$SENSITIVITY >= 70 &
                                       final.df$BINARY_CMA >= 70 &
                                       final.df$BINARY_BARBOUR >= 2 &
                                       final.df$REF_RANGE_CLASS %in% "HIGH" &
                                       final.df$REF_VARIABILITY_CLASS %in% "LOW" &
                                       final.df$ZERO_INFLATE %in% "REVIEW", "REVIEW", "POOR"))
  final.df <- final.df[!final.df$METRICS %in% "EFFECTIVE_RICH_SIMPSON", ]
  
  return(final.df)
}

#==============================================================================
#'Zero Inflation Test
#'
#'@param metrics.df = data frame of metric values for each station with site
#'a column of site classes defined by environmental variables.
#'@param bi.barbour = a data frame created within another function and used for
#'the calculated disturbance value.
#'@return Tests the influence of zeros on the results.
#'@export

zero_inflate <- function(metrics.df, bi.barbour){
  ref.df <- metrics.df[metrics.df$CATEGORY %in% "REF", ]
  if(ncol(ref.df) > 7){
    ref.df <- ref.df[, c(names(ref.df[, 1:6]), sort(names(ref.df[, 7:ncol(ref.df)])))]
  }
  deg.df <- metrics.df[metrics.df$CATEGORY %in% "SEV", ]
  if(ncol(deg.df) > 7){
    deg.df <- deg.df[, c(names(deg.df[, 1:6]), sort(names(deg.df[, 7:ncol(deg.df)])))]
  }
  
  barb <- bi.barbour[, c("METRICS", "DISTURBANCE")]
  new.df <- data.frame(METRICS = names(metrics.df[, 7:ncol(metrics.df)]))
  new.df <- merge(new.df , barb, by = "METRICS")
  
  if(ncol(ref.df) > 7){
    new.df$PCT_0_REF <- apply(ref.df[, 7:ncol(ref.df)], 2, function(x){
      round((sum(x == 0) / length(x)) * 100, 0)
    })
    
    new.df$PCT_0_DEG <- apply(deg.df[, 7:ncol(deg.df)], 2, function(x){
      round((sum(x == 0) / length(x)) * 100, 0)
    })
  }
 
  if(ncol(ref.df) == 7){
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