#==============================================================================
# IBI Methodologies
#==============================================================================
#==============================================================================
#'Run Multiple IBI Methodologies
#'
#'@param raw.data
#'@param master.df
#'@param ibi.method
#'@param todays_date
#'@param redundant
#'@param range.var
#'@param zero.inflate
#'@param metric.groups
#'@param fam.group
#'@param jack.iterate
#'@param metric.class
#'@return 
#'@export


test_ibi_method <- function(raw.data, master.df, ibi.method, todays_date, redundant,
                            range.var, zero.inflate, metric.groups, fam.group,
                            jack.iterate, metric.class = m.c){
  
  for(i in unique(raw.data$BIOREGION)){
    #============================================================================
    # Select only data from the bioregion of interest
    new.df <- raw.data[raw.data$BIOREGION %in% i, ]
    # Samples must contain more than 70 observed individuals.
    # Smaller sample sizes may result in odd percentage metric values.
    new.df <- new.df[new.df$ABUNDANCE > 70, ]
    
    remove.cols.1 <- c("EFFECTIVE_RICH_SHANNON", "EFFECTIVE_RICH_SIMPSON",
                     "PCT_UNIDENTIFIED", "NO_MATCH")
    rc <- names(new.df)[colSums(new.df[, 8:ncol(new.df)]) ==0]
    remove.cols.2 <- 
    new.df <- new.df[, !names(new.df) %in% remove.cols.1]
    
    if(ibi.method %in% c("B", "C", "D", "E")){
      if(ibi.method %in% c("C", "D", "E")){
        reg.metrics <- c("PCT_EPT_RICH_NO_TOL", "RICH", "RICH_CLING",
                         "RICH_INTOL", "SIMPSONS", "BECKS_V3", "HBI",
                         "PCT_DOM3", "PCT_INTOL_0_3", "PCT_COLLECT",
                         "PCT_PREDATOR", "PCT_SCRAPE", "PCT_SHRED",
                         "PCT_BURROW", "PCT_CLIMB", "PCT_CLING",
                         "PCT_SPRAWL", "PCT_SWIM", "PCT_NON_INSECT",
                         "PCT_ANNELID_CHIRO", "PCT_COTE", "PCT_DIPTERA",
                         "PCT_EPHEMEROPTERA_NO_BAETID", "PCT_EPT")
        if(ibi.method %in% c("D", "E")){
          fam.metrics <- na.omit(unique(master.df$FAMILY))
          fam.metrics <- paste("PCT_", fam.metrics, sep = "")
          keep.metrics <- c("EVENT_ID", "BIOREGION", "CATEGORY",
                            "STATION_ID", "SAMPLE_NUMBER",
                            "AGENCY_CODE", "DATE",
                            reg.metrics, fam.metrics)
        }else{
          keep.metrics <- c("EVENT_ID", "BIOREGION", "CATEGORY",
                            "STATION_ID", "SAMPLE_NUMBER",
                            "AGENCY_CODE", "DATE",
                            reg.metrics)
        }
        
      }else{
        standard.metrics <- c("EVENT_ID", "BIOREGION", "CATEGORY","STATION_ID", "SAMPLE_NUMBER",
                         "AGENCY_CODE", "DATE","PCT_INTOL_0_3", "PCT_MOD_TOL_4_6", "PCT_TOLERANT_7_10",
                         "PCT_GATHER", "PCT_FILTER", "PCT_PREDATOR", "PCT_SCRAPE", "PCT_SHRED",
                         "PCT_BURROW", "PCT_CLIMB", "PCT_CLING", "PCT_SPRAWL", "PCT_SWIM")
        div_comp.metrics <- as.character(m.c[m.c$METRIC_CLASS %in% c("DIVERSITY", "COMPOSITION"), "METRICS"])
        keep.metrics <- c(standard.metrics, div_comp.metrics)
        
      }
      
      new.df <- new.df[, names(new.df) %in% keep.metrics]
    }
    #============================================================================
    # Metric Assessment
    #keep.cols <- new.df[, 8:ncol(new.df)]
    ms <- metrics_summary(new.df, i, zero = zero.inflate)
    os <- old_scoring3(new.df, i, zero_null = zero.inflate, bound.lim = TRUE) #, metric.summary = ms)
    remove.cols <- c("ABUNDANCE", "EFFECTIVE_RICH_SHANNON", "EFFECTIVE_RICH_SIMPSON",
                     "PCT_UNIDENTIFIED", "NO_MATCH")
    os <- os[, !grepl(paste0(remove.cols, collapse = "|"), names(os))]
    if(ibi.method %in% c("A", "B", "H")){
      if(ibi.method %in% "H"){
        #Remove metrics that cannot identify disturbance better than a random flip of a coin
        new.df <- new.df[, !names(new.df) %in% ms[ms$SENSITIVITY < 50, "METRICS"]]
        ms <- ms[ms$SENSITIVITY >= 50, ]
        os <- os[, !names(os) %in% ms[ms$SENSITIVITY < 50, "METRICS"]]
        #===================================================================================
        b.metrics <- select_best(metrics.df = new.df, ms, m.c = metric.class,
                                    range_var = range.var,
                                    redund = redundant)
        bd <- best.distribution.H(new.df, ms, os, m.c = metric.class,
                                 range_var = range.var,
                                 redund = redundant,
                                 best.metrics = b.metrics)
        
        metrics.vec2 <- bd$METRICS
      }
      if(ibi.method %in% "A"){
        bd <- best.distribution2(metrics.df = new.df,
                                 metric.summary = ms,
                                 scores.df = os,
                                 m.c = metric.class,
                                 Family = fam.group,
                                 master = master.df,
                                 range_var = range.var)
        metrics.vec2 <- bd$METRICS
        
      }
      if(ibi.method %in% "B"){

        bd <- best.distribution3(new.df, ms, os, m.c= metric.class,
                                 Family = fam.group, 
                                 master <- master.df,
                                 range_var = range.var)
        metrics.vec2 <- bd$METRICS
      }
      
      
      write.csv(bd, paste(ibi.method, "_", i, "_ms_", todays_date, ".csv", sep = ""))
      #============================================================================
      # Creat a vectorof  metrics with sensitivity values (DE) >= the best 
      # threshold provided by the best.distribution function.
      #metrics.vec <- ms[ms$SENSITIVITY >= bd[1, 7], "METRICS"]
      
      
      
      final.metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                       "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                       as.character(metrics.vec2)))
      final.metrics.vec2 <- final.metrics.vec[!final.metrics.vec %in% remove.cols]
      if(redundant == TRUE & length(final.metrics.vec) > 7){
        if(ibi.method %in% c("A", "B")){
        R_test <- metric_class_redund(new.df[, names(new.df) %in% final.metrics.vec],
                            ms[ms$METRICS %in% final.metrics.vec,],
                            m.c, method = ibi.method)
        final.metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                         "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                         as.character(R_test$METRICS)))
        
      }else{
        R_test <- redundancy(new.df[, names(new.df) %in% final.metrics.vec],
                             analysis = "wilcox", sensitivity.df = ms[ms$METRICS %in% final.metrics.vec, ],
                             lower.class = "SEV", upper.class = "REF",
                             lower.coef = -0.85, upper.coef = 0.85)
        final.metrics.vec <- unlist(list("EVENT_ID", "CATEGORY", "STATION_ID",
                                         "SAMPLE_NUMBER", "AGENCY_CODE", "DATE",
                                         as.character(R_test$METRICS)))
      }
      }
      
      
      if(metric.groups != TRUE){
        #============================================================================
        # Select only the metrics selected for the final index from the os (contains scored values)
        scored <- os[, final.metrics.vec]
        
        if(length(final.metrics.vec) > 7){
          scored[, 7:ncol(scored)] <- apply(scored[, 7:ncol(scored)], 2, function(x) as.numeric(as.character(x)))
          scored$FINAL_SCORE <- apply(scored[, 7:ncol(scored)], 1, mean, na.rm = TRUE)
        }else{
          scored[, 7] <- as.numeric(as.character(scored[, 7]))
          scored$FINAL_SCORE <- scored[, 7]
        }

      }else{
        scored <- scoring2(metrics.df = new.df[, names(new.df) %in% final.metrics.vec],
                           scores = os[, names(os) %in% final.metrics.vec],
                           ms,
                           m.c = metric.class,
                           Family = fam.group,
                           master = master.df,
                           redund = redundant,
                           method = ibi.method)
      }
      
    }else{
      if(metric.groups == TRUE){
        scored <- scoring2(metrics.df = new.df,
                           scores = os,
                           ms,
                           m.c = metric.class,
                           Family = fam.group,
                           master = master.df,
                           redund = redundant,
                           method = ibi.method)
        
      }else{
        scored <- os
        scored[, 7:ncol(scored)] <- apply(scored[, 7:ncol(scored)], 2, function(x) as.numeric(as.character(x)))
        scored$FINAL_SCORE <- apply(scored[, 7:ncol(scored)], 1, mean, na.rm = TRUE)
      }
    }
    
    
    write.csv(scored, paste(ibi.method, "_", i, "_metrics_", todays_date, ".csv", sep = ""))
    #============================================================================
    # Balanced Descrimination Efficiency
    t.bde <- bde(scored, i)
    write.csv(scored, paste(ibi.method, "_", i, "_bde_", todays_date, ".csv", sep = ""))
    #============================================================================
    # Create a pdf of the final index distributions
    pdf(paste(ibi.method, "_", i, "_Index_", todays_date, ".pdf", sep = ""))
    scored$CATEGORY <- factor(scored$CATEGORY,
                              levels = c("REF", "MIN", "MOD", "SEV", "MIX"))
    boxplot(scored[, "FINAL_SCORE"] ~ scored$CATEGORY,
            a <- i,
            data = raw.data,
            col="white", las = 0,
            main = substitute(paste(a)))
    
    dev.off()
    
    #============================================================================
    # Create a pdf of all of the individual metric used in the final index
    pdf(paste(ibi.method, "_", i, "_raw_metrics_", todays_date, ".pdf", sep = ""))
    par(mfrow=c(2,2), mar=c(4,2,2,2), oma=c(2,2,2,2))
    plot.me <- new.df[,names(new.df) %in% names(scored)]
    #***
    col.names <- names(plot.me)
    col.names <- col.names[!col.names %in% c("div", "FFG", "HABIT",
                                             "tol", "comp", "fam", "FINAL_SCORE")]
    for (j in col.names[!col.names %in% c("EVENT_ID", "CATEGORY", "STATION_ID", 
                                          "SAMPLE_NUMBER", "AGENCY_CODE", "DATE")]){
      a <- j
      #c_list <- list('REF', 'NEAR', 'MIN', 'MOD', 'SEV', 'MIX')
      plot.me$CATEGORY <- factor(plot.me$CATEGORY,
                                 levels = c("REF", "MIN", "MOD", "SEV", "MIX"))
      
      boxplot(plot.me[, j] ~ plot.me$CATEGORY,
              #data = raw.data,
              #names = c_list,
              col="white", las = 0,
              main = substitute(paste(a)))
    }
    dev.off()
    #============================================================================
    # Create a pdf of all of the individual metric used in the final index
    pdf(paste(ibi.method, "_", i, "_scored_metrics_", todays_date, ".pdf", sep = ""))
    par(mfrow=c(2,2), mar=c(4,2,2,2), oma=c(2,2,2,2))
    plot.me <- scored
    #***
    col.names <- names(plot.me)
    col.names <- col.names[!col.names %in% c("div", "FFG", "HABIT",
                                             "tol", "comp", "fam", "FINAL_SCORE")]
    for (j in col.names[!col.names %in% c("EVENT_ID", "CATEGORY", "STATION_ID", 
                                          "SAMPLE_NUMBER", "AGENCY_CODE", "DATE")]){
      a <- j
      #c_list <- list('REF', 'NEAR', 'MIN', 'MOD', 'SEV', 'MIX')
      plot.me$CATEGORY <- factor(plot.me$CATEGORY,
                                 levels = c("REF", "MIN", "MOD", "SEV", "MIX"))
      
      boxplot(plot.me[, j] ~ plot.me$CATEGORY,
              #data = raw.data,
              #names = c_list,
              col="white", las = 0,
              main = substitute(paste(a)))
    }
    dev.off()
    #============================================================================
    #============================================================================
    # Create a data frame for jackknife percision and valdiation
    fmv <- unlist(list("BIOREGION", as.character(col.names)))
    
    final_metrics.df <- new.df[, names(new.df) %in% fmv]
    
    #============================================================================
    # Jackknife Percision
    # should always be false because we don't want the index to remove the
    # metrics selected for by the program
    system.time(test1 <- jack.calc(final_metrics.df, jack.iterate, 0.75, i,
                                   m.c = metric.class,
                                   Fam = fam.group,
                                   Master = master.df,
                                   zero.null = zero.inflate,
                                   metric_types = metric.groups,
                                   redund = FALSE,
                                   method = ibi.method))
    
    
    jack.summary <- function(scored.df, jack.output, percentile.value, bioregion, measure){
      jack <- data.frame(BIOREGION = bioregion)
      jack$MEASURE <- measure
      n.ref <- nrow(scored.df[scored.df$CATEGORY %in% "REF", ])
      n.deg <- nrow(scored.df[scored.df$CATEGORY %in% "SEV", ])
      n <- n.ref + n.deg
      d.ref <- round(n.ref * 0.25, 0)
      d.deg <- round(n.deg * 0.25, 0)
      d <- d.ref + d.deg
      N <- ncol(jack.output)
      
      jack$ORGINAL_VALUE <- percentile.value[1]
      jack$MEAN_SIM_VALUE <- mean(jack.output)
      resid.jack <- (jack.output) - jack$MEAN_SIM_VALUE
      jack$MSE <- ((n - d) / (d * length(jack.output))) * (sum(resid.jack ^ 2))
      jack$RMSE <- sqrt(jack$MSE)
      
      jack$BIAS <- sum(jack$MEAN_SIM_VALUE - jack.output) / length(jack.output)
      
      return(jack)
    }
    
    j1 <- jack.summary(scored, unlist(lapply(test1, "[[", "THRESHOLD")),
                       t.bde$THRESHOLD, i, "THRESHOLD")
    #j1[, 3:ncol(j1)] <- j1[, 3:ncol(j1)] * 100
    j2 <- jack.summary(scored, unlist(lapply(test1, "[[", "CE")), t.bde$CE, i, "CE")
    j3 <- jack.summary(scored, unlist(lapply(test1, "[[", "BDE")), t.bde$BDE, i, "BDE")
    j4 <- jack.summary(scored, unlist(lapply(test1, "[[", "PCT_REF")), t.bde$PCT_REF, i, "PCT_REF")
    j5 <- jack.summary(scored, unlist(lapply(test1, "[[", "PCT_DEG")), t.bde$PCT_DEG, i, "PCT_DEG")
    
    jack.final <- rbind(j1, j2, j3, j4, j5)
    write.csv(jack.final, paste(ibi.method, "_", i, "_jack_precision_", todays_date, ".csv", sep = ""))
    
    #============================================================================
    # Jackknife Validation
    system.time(jack_val <- jack.val(final_metrics.df, jack.iterate, 0.75, i,
                                     m.c = metric.class,
                                     Fam = fam.group,
                                     Master = master.df,
                                     zero.null = zero.inflate,
                                     metric_types = metric.groups,
                                     redund = FALSE,
                                     method = ibi.method))
    jack_val$ORIGINAL_CE <- t.bde$CE[1]
    write.csv(jack_val, paste(ibi.method, "_", i, "_jack_validation_", todays_date, ".csv", sep = ""))
    
  }
}
