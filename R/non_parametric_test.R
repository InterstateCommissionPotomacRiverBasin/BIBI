#==============================================================================
#Nonparametric Test Functions
#==============================================================================
#'Kruskal Wallis p-value and test Statistic
#'
#'@param Prep_Taxa = Taxanomic counts in a long data format.
#'@param Info = Taxanomic attributes table.
#'@param Tab_Project = Table containing Project_ID and agency name.
#'@param Tab_Event = Table containing Project_ID and Event_ID
#'@return Creates a table of the Kruskal Wallis output. The table includes
#'the EVENT_ID, test statistic, and the p-values comparing agencies.
#'@export
kruskal_table <- function(Prep_Taxa, Info, Tab_Project, Tab_Event){
  sub_pro <- Tab_Project[, c("PROJECT_ID", "AGENCY_CODE")]
  sub_event <- Tab_Event[, c("EVENT_ID", "PROJECT_ID")]

  Bioregion_Metrics <- all_metrics(Info, Prep_Taxa)
  merged1 <- merge(sub_event, Bioregion_Metrics, by = "EVENT_ID", all.y = TRUE)

  merged2 <- merge(sub_pro, merged1, by = "PROJECT_ID", all.y = TRUE)

  nt <- lapply(merged2[, 5:ncol(merged2)], function(x) kruskal.test(x ~ merged2$AGENCY_CODE))
  df <- data.frame(matrix(unlist(nt), nrow =  length(nt), byrow=T),stringsAsFactors=FALSE)
  colnames(df) <- c("kw_Statistic", "df", "kw_p_value", "Test", "Description")

  new.df <- df[, 1:3]
  new.df$kw_Statistic <- round(as.numeric(as.character(new.df$kw_Statistic)), digits = 4)
  new.df$kw_p_value <- round(as.numeric(as.character(new.df$kw_p_value)), digits = 4)

  metrics.df <- data.frame(colnames(merged2[, 5:ncol(merged2)]))
  colnames(metrics.df) <- "Metric"
  cbound <- cbind(metrics.df, new.df)
  return(cbound)
}

#==============================================================================
#'Post hoc Dunn's test p-value and test Statistic
#'
#'@param Prep_Taxa = Taxanomic counts in a long data format.
#'@param Info = Taxanomic attributes table.
#'@param Tab_Project = Table containing Project_ID and agency name.
#'@param Tab_Event = Table containing Project_ID and Event_ID
#'@return Creates a table of the post hoc Dunn's testoutput.
#'The table includes the EVENT_ID, test statistic, and the p-values
#'comparing agencies.
#'@export
dunn_table <- function(Prep_Taxa, Info, Tab_Project, Tab_Event){
  sub_pro <- Tab_Project[, c("PROJECT_ID", "AGENCY_CODE")]
  sub_event <- Tab_Event[, c("EVENT_ID", "PROJECT_ID")]

  Bioregion_Metrics <- all_metrics(Info, Prep_Taxa)
  merged1 <- merge(sub_event, Bioregion_Metrics, by = "EVENT_ID", all.y = TRUE)

  merged2 <- merge(sub_pro, merged1, by = "PROJECT_ID", all.y = TRUE)
  merged2$AGENCY_CODE <- as.character(merged2$AGENCY_CODE)

  merged3 <- merged2[5:ncol(merged2)]
  merged4 <- merged3[, !colSums(merged3) %in% 0]
  merged5 <- cbind(merged2[, 1:5], merged4)

  nt <- lapply(merged5[, 5:ncol(merged5)], function(x) dunn.test(x, merged5$AGENCY_CODE, kw = TRUE, label = TRUE))
  df <- data.frame(names(nt), matrix(unlist(nt),
                                     ncol = 1 + 4 * length(unique(nt$RICH$comparisons)),
                                     byrow = T), stringsAsFactors = FALSE)
  colnames(df) <- c("Metric", "Dunn_Chi2", paste("Z", nt$RICH$comparisons, sep="_"), paste("Dunn_P_Value", nt$RICH$comparisons, sep="_") , paste("P_Adjust", nt$RICH$comparisons, sep="_"))
  grepped <- grep("Dunn_P_Value", colnames(df))
  dunn_pvalue <- df[,c(1, grepped)]
  if(length(grepped) > 1){
    dunn_pvalue[, 2:ncol(dunn_pvalue)] <- apply(dunn_pvalue[, 2:ncol(dunn_pvalue)], 2,
                                                function(x) as.numeric(as.character(x)))
  }else{
    dunn_pvalue[, 2:ncol(dunn_pvalue)] <- as.numeric(as.character(dunn_pvalue[, 2:ncol(dunn_pvalue)]))
  }


  dunn_pvalue[,2:ncol(dunn_pvalue)] <- round(dunn_pvalue[,2:ncol(dunn_pvalue)], digits = 4)
  return(dunn_pvalue)
}
#==============================================================================
#'Kruskal Wallis p-value and Dunn's test
#'
#'@param Prep_Taxa = Taxanomic counts in a long data format.
#'@param Info = Taxanomic attributes table.
#'@param Tab_Project = Table containing Project_ID and agency name.
#'@param Tab_Event = Table containing Project_ID and Event_ID
#'@return Creates a table of the Kruskal Wallis and Dunn's test output.
#'The table includes the EVENT_ID, test statistic, the
#'Kruskal Wallis p-values, and the Dunn's test p-values comparing agencies.
#'@export
kd_table <- function(Prep_Taxa, Info, Tab_Project, Tab_Event){
  dunn.df <- dunn_table(Prep_Taxa, Info, Tab_Project, Tab_Event)
  kruskal.df <- kruskal_table(Prep_Taxa, Info, Tab_Project, Tab_Event)
  both.df <- merge(kruskal.df, dunn.df, by = "Metric")
  return(both.df)
}

