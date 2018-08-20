#==============================================================================
# Functions Testing for Correlations
#==============================================================================
#'Taxa Correlations
#'
#'@param Long = Long data formate containing taxonomic counts.
#'@param Level = The taxonomic level used in the analysis.
#'@param Coef = The correlation coeficient to be used. The value will be used
#'as the positive and negative threshold.
#'@return Taxa positively or negatively correlated with each other.  These
#'taxa may work well togther as a percent metric.
#'@export

corr_taxa <- function(Long, Level, value = 0.50){
  wide.df <- prep_rare2(Long, Level)
  new.df <- data.frame(names(wide.df[, 6:ncol(wide.df)]))
  new.df$counts <- colSums(wide.df[, 6:ncol(wide.df)])
  new.df$freq <- vegan::specnumber(wide.df[, 6:ncol(wide.df)], MARGIN = 2)
  new.df$pct <- (new.df$freq / nrow(wide.df)) * 100
  new <- new.df[new.df$pct >= 1, ]

  new.wide <- wide.df[, which(names(wide.df) %in% new[, 1])]
  bound <- cbind(wide.df[, 1:5], new.wide)

  Corr.df <- data.frame(cor(bound[, 6:ncol(bound)],
                            method = "spearman"))
  Corr.df[upper.tri(Corr.df, diag = "TRUE")] <- NA
  Corr.df$TAXA<- rownames(Corr.df)
  melted <- reshape2::melt(Corr.df, id.vars = c("TAXA"))
  sub.melted <- melted[c(melted$value >= value | melted$value <= -value),]
  sub.melted <- na.omit(sub.melted)
  return(sub.melted)
}
