kruskal_bias <- function(Prep_Taxa, Info, Tab_Project, Tab_Event){
  sub_pro <- Tab_Project[, c("PROJECT_ID", "AGENCY_CODE")]
  sub_event <- Tab_Event[, c("EVENT_ID", "PROJECT_ID")]

  Bioregion_Metrics <- all_metrics(Info, Prep_Taxa)
  merged1 <- merge(sub_event, Bioregion_Metrics, by = "EVENT_ID", all.y = TRUE)

  merged2 <- merge(sub_pro, merged1, by = "PROJECT_ID", all.y = TRUE)

  nt <- lapply(merged2[, 5:ncol(merged2)], function(x) kruskal.test(x ~ merged2$AGENCY_CODE))
  df <- data.frame(matrix(unlist(nt), nrow =  length(nt), byrow=T),stringsAsFactors=FALSE)
  colnames(df) <- c("kw_Statistic", "df", "kw_p_value", "Test", "Description")

  new.df <- df[, 1:3]
  new.df$kw_Statistic <- round(as.numeric(as.character(new.df$kw_statistic)), digits = 4)
  new.df$kw_p_value <- round(as.numeric(as.character(new.df$kw_p_value)), digits = 4)

  metrics.df <- data.frame(colnames(merged2[, 5:ncol(merged2)]))
  colnames(metrics.df) <- "Metric"
  cbound <- cbind(metrics.df, new.df)
  return(cbound)
}
