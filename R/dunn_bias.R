dunn_bias <- function(Prep_Taxa, Info, Tab_Project, Tab_Event){
  sub_pro <- Tab_Project[, c("PROJECT_ID", "AGENCY_CODE")]
  sub_event <- Tab_Event[, c("EVENT_ID", "PROJECT_ID")]

  Bioregion_Metrics <- all_metrics(Info, Prep_Taxa)
  merged1 <- merge(sub_event, Bioregion_Metrics, by = "EVENT_ID", all.y = TRUE)

  merged2 <- merge(sub_pro, merged1, by = "PROJECT_ID", all.y = TRUE)
  merged2$AGENCY_CODE <- as.character(merged2$AGENCY_CODE)

  merged3 <- merged2[5:ncol(merged2)]
  merged4 <- merged3[, colSums(merged3) > 0]
  merged5 <- cbind(merged2[, 1:5], merged4)

  nt <- lapply(merged5[, 5:ncol(merged5)], function(x) dunn.test(x, merged5$AGENCY_CODE, kw = TRUE, label = TRUE))
  df <- data.frame(names(nt), matrix(unlist(nt),
                                     ncol = 1 + 4 * length(unique(nt$RICH$comparisons)),
                                     byrow = T), stringsAsFactors = FALSE)
  colnames(df) <- c("Metric", "Dunn_Chi2", paste("Z", nt$RICH$comparisons, sep="_"), paste("P_Value", nt$RICH$comparisons, sep="_") , paste("P_Adjust", nt$RICH$comparisons, sep="_"))
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
