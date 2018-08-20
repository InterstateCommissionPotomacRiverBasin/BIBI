#==============================================================================
#
# Jackknife Functions
#
#==============================================================================
#==============================================================================
#'Delete-d Jackknife Select Training Data Set
#'
#'@param ref.df = a data frame containing only Reference samples.
#'@param deg.df = a data frame containing only Degraded samples.
#'@param ref_deg = a data frame containing only Reference and Degraded samples.
#'@param runs = number of jackknife resampling events.
#'@param keep = the proportion of the reference population that should be kept.
#'@return Performes a Delete-d Jackknife to create a unique subset of the orginal
#'data set for the number of specified runs.
#'@export

jack_train_samp <- function(ref.df, deg.df, ref_deg, runs, keep){
  #ref.df <- metrics.df[metrics.df$CATEGORY %in% "REF", ]
  #not_ref.df <- metrics.df[!metrics.df$CATEGORY %in% "REF", ]
  #deg.df <- metrics.df[metrics.df$CATEGORY %in% "SEV", ]
  #ref_deg <- rbind(ref.df, deg.df)
  
  
  choose.int <- function(x, n, k) {
    if(n <= k) return(rep(TRUE, k))
    u <- choose(n-1, k-1)
    pick <- x < u
    if (pick) y <- choose.int(x, n-1, k-1) else y <- choose.int(x-u, n-1, k)
    return(c(pick, y))
  }
  
  #REFERENCE subsample
  n <- nrow(ref.df)
  k <- floor(nrow(ref.df) * keep)
  if (n > 1000){
    n <- 1000
    k <- 1000 * keep
  }
  bf <- gmp::chooseZ(n, k)
  rand.samp <- gmp::urand.bigz(nb = runs, size = (log(bf) / log(2)), seed = 0)
  ref.sub.samp <- data.frame(sapply(as.numeric(as.character(rand.samp)), choose.int, n = n, k = k))
  #DEGRADED subsample
  n <- nrow(deg.df)
  k <- floor(nrow(deg.df) * keep)
  if (n > 1000){
    n <- 1000
    k <- 1000 * keep
  }
  bf <- gmp::chooseZ(n, k)
  
  rand.samp <- gmp::urand.bigz(nb = runs, size = (log(bf) / log(2)), seed = 0)
  deg.sub.samp <- data.frame(sapply(as.numeric(as.character(rand.samp)), choose.int, n = n, k = k))
  
  train.samp <- rbind(ref.sub.samp, deg.sub.samp)
  return(train.samp)
}

#==============================================================================
#'Delete-d Jackknife Calculations
#'
#'@param metrics.df = a data frame containing the raw metric values. Only the
#' metrics selected for the final index are included.
#'@param runs = number of jackknife resampling events.
#'@param keep = the proportion of the reference population that should be kept.
#'@param bioregion = the bioregion of interest
#'@param m.c = a data frame containing all of the metrics and thier associated
#'metric category.
#'@param Fam = If TRUE, all family metrics are grouped and scored together 
#'prior to the final scoring procedure.  If FALSE, family metrics are grouped
#'with composition metrics.
#'@param Master = a data frame containing taxonomic traits and attributes.
#'@param zero.null = If TRUE, all zeros are converted to NA for all metrics
#' with > 5% of Reference values equal to zero. If FALSE, metric values are
#' unaltered.
#'@param redund.df = IF TRUE, metrics are evaluated for redundancy. If FALSE,
#'a redundancy analysis is not preformed.
#'@param method = specify the method being preformed. FIGURE OUT A BETTER WAY!!!!
#'@return Recomputes the final IBI threshold after removing a specified
#' percentage of the data.
#'@export

jack_sim <- function(train.samp, ref_deg, bioregion, zero.null, Fam, Master,
                     redund.df, m.c = metric.class,
                     metric_types = FALSE, method){
  data.list <- list()
  for(i in 1:ncol(train.samp)){
    print(paste("Start Jackknife Iteration:", i))
    train_sites.df <- ref_deg[unlist(train.samp[, i]), ]
    
    
    #all_sites.df <- rbind(ref, not_ref.df)
    
    ms.train <- metrics_summary2(train_sites.df, bioregion)
    os.train <- old_scoring3(train_sites.df, bioregion, zero_null = zero.null)
    #bd.test <- best.distribution2(all_sites.df, ms.test, os.test, m.c)
    #solid.m <- solid_metrics(new.df, os2, ms, bd2[1, "THRESH"], m.c)
    if(metric_types == TRUE){
      final.metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                       as.character(ms.train$METRICS)))
      scored.train <- scoring2(metrics.df = train_sites.df[, names(train_sites.df) %in% final.metrics.vec],
                               scores = os.train[, names(os.train) %in% final.metrics.vec],
                               ms = ms.train, m.c,
                               Family = Fam,
                               master = Master,
                               redund = redund.df,
                               method)
      
      bde.train <- bde(scored.train, bioregion)
    }else{
      if(ncol(os.train) > 7){
        os.train[, 7:ncol(os.train)] <- apply(os.train[, 7:ncol(os.train)], 2, function(x) as.numeric(as.character(x)))
        os.train$FINAL_SCORE <- apply(os.train[, 7:ncol(os.train)], 1, mean, na.rm = TRUE)
      }else{
        os.train[, 7] <- as.numeric(as.character(os.train[, 7]))
        os.train$FINAL_SCORE <- os.train[, 7]
      }
      bde.train <- bde(os.train, bioregion)
    }
    
    
    val.samp <- ifelse(train.samp[, i] == TRUE, FALSE, TRUE)
    val_sites.df <- ref_deg[unlist(val.samp), ]
    os.val <- old_scoring3(val_sites.df, bioregion, zero_null = zero.null, metric.summary = ms.train)
    #tt <- ms[ms$METRICS %in% final.metrics.vec, ]
    if(metric_types == TRUE){
      scored.val <- scoring2(metrics.df =val_sites.df[, names(val_sites.df) %in% final.metrics.vec],
                             scores = os.val[, names(os.val) %in% final.metrics.vec],
                             ms = ms.train,
                             m.c,
                             Family = Fam,
                             master = Master,
                             redund = redund.df, 
                             method)
    }else{
      if(ncol(os.val) > 7){
        os.val[, 7:ncol(os.val)] <- apply(os.val[, 7:ncol(os.val)], 2, function(x) as.numeric(as.character(x)))
        os.val$FINAL_SCORE <- apply(os.val[, 7:ncol(os.val)], 1, mean, na.rm = TRUE)
      }else{
        os.val[, 7] <- as.numeric(as.character(os.val[, 7]))
        os.val$FINAL_SCORE <- os.val[, 7]
      }
      
      scored.val <- os.val
    }
    
    score_ref_deg <- function(scored.df, bde.train){
      scored.ref <- scored.df[scored.df$CATEGORY %in% "REF", ]
      ref.vec <- ifelse(scored.ref$FINAL_SCORE * 100 >= bde.train$THRESHOLD[1], 1, 0) 
      pct.ref <- (sum(ref.vec) / length(ref.vec)) * 100
      scored.deg <- scored.df[scored.df$CATEGORY %in% "SEV", ]
      deg.vec <- ifelse(scored.deg$FINAL_SCORE * 100 < bde.train$THRESHOLD[1], 1, 0)
      pct.deg <- (sum(deg.vec) / length(deg.vec)) * 100
      
      final.value <- (pct.ref + pct.deg) / 2
      return(final.value)
    }
    
    comb.df <- data.frame(RUN = i)
    comb.df$VAL_CE <- score_ref_deg(scored.val, bde.train)
    comb.df$TRAIN_CE <- score_ref_deg(os.train, bde.train)
    comb.df$BSP <- bde.train$THRESHOLD[1]
    
    gathered.df <- tidyr::gather(comb.df[, 2:ncol(comb.df)], MEASURE, VALUE)
    gathered.df$RUN <- i
    gathered.df <- gathered.df[, c(ncol(gathered.df), 1:(ncol(gathered.df) - 1))]
    
    data.list[[i]] <- gathered.df
    
  }
  
  final.df <- do.call(rbind, data.list)
  
  
  return(final.df)
}          

#==============================================================================
#'Delete-d Jackknife Select Training Data Set
#'
#'@param ref.df = a data frame containing only Reference samples.
#'@param deg.df = a data frame containing only Degraded samples.
#'@param bioregion = the bioregion or region of interest.
#'@param jack.sim = output from the jack_sim function.
#'@return Summarizes the Validation-set CE, Training-set CE, and BSP.
#'@export

jack_summarize <- function(jack.sim, bioregion, ref.df, deg.df){
  
  n.ref <- nrow(ref.df)
  n.deg <- nrow(deg.df)
  n <- n.ref + n.deg
  d.ref <- round(n.ref * 0.25, 0)
  d.deg <- round(n.deg * 0.25, 0)
  d <- d.ref + d.deg
  #N <- length(jack.train)
  
  data.list <- list()
  for(i in unique(jack.sim$MEASURE)){
    sub.jack <- jack.sim[jack.sim$MEASURE %in% i, ]
    jack <- data.frame(BIOREGION = bioregion)
    jack$MEASURE <- i
    jack$MEAN_SIM_VALUE <- mean(sub.jack$VALUE)
    
    resid.jack <- (sub.jack$VALUE) - jack$MEAN_SIM_VALUE
    jack$MSE <- (sum(resid.jack ^ 2)) / nrow(sub.jack)
    jack$RMSE <- sqrt(jack$MSE)
    jack$BIAS <- sum(jack$MEAN_SIM_VALUE - sub.jack$VALUE) / nrow(sub.jack)
    
    data.list[[i]] <- jack
  }
  final.df <- do.call(rbind, data.list)
  
  return(final.df)
}

#==============================================================================
#'Delete-d Jackknife
#'
#'@param metrics.df = a data frame containing the raw metric values. Only the
#' metrics selected for the final index are included.
#'@param runs = number of jackknife resampling events.
#'@param keep = the proportion of the reference population that should be kept.
#'@param bioregion = the bioregion of interest
#'@param m.c = a data frame containing all of the metrics and thier associated
#'metric category.
#'@param Fam = If TRUE, all family metrics are grouped and scored together 
#'prior to the final scoring procedure.  If FALSE, family metrics are grouped
#'with composition metrics.
#'@param Master = a data frame containing taxonomic traits and attributes.
#'@param zero.null = If TRUE, all zeros are converted to NA for all metrics
#' with > 5% of Reference values equal to zero. If FALSE, metric values are
#' unaltered.
#'@param redund.df = IF TRUE, metrics are evaluated for redundancy. If FALSE,
#'a redundancy analysis is not preformed.
#'@param method = specify the method being preformed. FIGURE OUT A BETTER WAY!!!!
#'@param seed = set.seed for repeatable results if TRUE.
#'@return Recomputes the final IBI threshold after removing a specified
#' percentage of the data.
#'@export
#'
bibi_jackknife <- function(metrics.df, runs, keep, bioregion, m.c, Fam = FALSE, Master = master,
                           zero.null = TRUE, metric_types = TRUE, redund.df, method = ibi.method,
                           seed = FALSE){
  ref.df <- metrics.df[metrics.df$CATEGORY %in% "REF", ]
  not_ref.df <- metrics.df[!metrics.df$CATEGORY %in% "REF", ]
  deg.df <- metrics.df[metrics.df$CATEGORY %in% "SEV", ]
  ref_deg <- rbind(ref.df, deg.df)
  if(seed == TRUE) set.seed(64)
  jack.train <- jack_train_samp(ref.df, deg.df, ref_deg, runs, keep)
  jack.sim <- jack_sim(jack.train, ref_deg, bioregion, zero.null, Family, Master,
                       redund.df, m.c, metric_types, method)
  final.df <- jack_summarize(jack.sim, bioregion, ref.df, deg.df)
  return(final.df)
}




