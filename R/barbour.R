#==============================================================================
# Barbour et al. 1996 methodology
#==============================================================================

#'Barbour et al. 1996 method for testing metric sensitivity
#'
#'@param metrics.df = data frame of metric values for each station
#'@param ref.df = a data frame of only the upper.class values.
#'@param deg.df = a data frame of only the lower.class values.
#'@return A data frame of metrics scored using the Barbour et al. 1996 method
#'for testing metric sensitivity
#'@export
barbour <- function(metrics.df, ref.df, deg.df){
  Barbour <- data.frame(colnames(metrics.df[7:ncol(metrics.df)]))
  colnames(Barbour) <- "METRICS" #rename column #1 "Metrics"

  #Quartiles for ref.df and DEG Categorys
  if(ncol(ref.df) > 7){
    Barbour$Ref25 <- sapply(ref.df[, 7:ncol(ref.df)], quantile, 0.25, na.rm = TRUE)
    Barbour$Ref50 <- sapply(ref.df[, 7:ncol(ref.df)], quantile, 0.50, na.rm = TRUE)
    Barbour$Ref75 <- sapply(ref.df[, 7:ncol(ref.df)], quantile, 0.75, na.rm = TRUE)
    Barbour$DEG25 <- sapply(deg.df[, 7:ncol(deg.df)], quantile, 0.25, na.rm = TRUE)
    Barbour$DEG50 <- sapply(deg.df[, 7:ncol(deg.df)], quantile, 0.50, na.rm = TRUE)
    Barbour$DEG75 <- sapply(deg.df[, 7:ncol(deg.df)], quantile, 0.75, na.rm = TRUE)
  }
  
  if(ncol(ref.df) == 7){
    Barbour$Ref25 <- quantile(ref.df[, 7], 0.25, na.rm = TRUE)
    Barbour$Ref50 <- quantile(ref.df[, 7], 0.50, na.rm = TRUE)
    Barbour$Ref75 <- quantile(ref.df[, 7], 0.75, na.rm = TRUE)
    Barbour$DEG25 <- quantile(deg.df[, 7], 0.25, na.rm = TRUE)
    Barbour$DEG50 <- quantile(deg.df[, 7], 0.50, na.rm = TRUE)
    Barbour$DEG75 <- quantile(deg.df[, 7], 0.75, na.rm = TRUE)
  }

  #Barbour$MedianDiff <- abs(Barbour$Ref50 - Barbour$DEG50)
  #Barbour$DecreaseDiff <- Barbour$Ref25 - Barbour$DEG75
  #Barbour$IncreaseDiff <- Barbour$DEG25 - Barbour$Ref75

  #Round the metric quartiles to two decimal places
  Barbour[, -1] <- round(Barbour[, -1], 2) #the "-1" excludes column 1

  Barbour$DISTURBANCE <- ifelse(Barbour$Ref50 > Barbour$DEG50,
                                "DECREASE", "INCREASE")

  # Scoring from Barbour et al. 1996
  # 0 = The median of both ref.df and DEG plots overlap the interquartile range of the other category
  # 1 = One median (either ref.df or DEG) overlaps the interquartile range of the other category
  # 2 = The interquartiles overlap but neither median overlaps with the interquartile range of the other category
  # 3 = The interquartile ranges do not overlap
  Barbour$SENSITIVITY <- ifelse((Barbour$Ref50 <= Barbour$DEG75 &
                                   Barbour$Ref50>=Barbour$DEG25 &
                                   Barbour$DEG50<=Barbour$Ref75 &
                                   Barbour$DEG50>=Barbour$Ref25), 0,
                                ifelse((Barbour$Ref50<=Barbour$DEG75 &
                                          Barbour$Ref50>=Barbour$DEG25 |
                                          Barbour$DEG50<=Barbour$Ref75 &
                                          Barbour$DEG50>=Barbour$Ref25), 1,
                                       ifelse((Barbour$Ref25<=Barbour$DEG75 &
                                                 Barbour$Ref25>=Barbour$DEG25 |
                                                 Barbour$Ref75<=Barbour$DEG75 &
                                                 Barbour$Ref75>=Barbour$DEG25), 2, 3)))

  final.df <- Barbour[, c("METRICS", "DISTURBANCE", "SENSITIVITY")]
  final.df$METRICS <- as.character(final.df$METRICS)

  return(final.df)
}

#==============================================================================
#'Barbour Sensitivity Score
#'
#'@param Barbour.df = data frame of metrics scored using the Barbour et al. 1996 method
#'for testing metric sensitivity
#'@param score = a number 0 - 3 indicating which sensitivity score will be
#'represented by the new data frame
#'@return A data frame representing one sensitivity score from the
#' Barbour et al. 1996 method
#'@export
#'
barbour_sensitivity <- function(Barbour.df, score){
  b.sensitivity <- by(Barbour.df, Barbour.df$Barbour_Score, FUN = print)
  b.score <- if(score == 3){
    b.sensitivity$'3'
  }else{
    if(score == 2){
      b.sensitivity$'2'
    }else{
      if(score == 1){
        b.sensitivity$'1'
      }else{
        if(score == 0){
          b.sensitivity$'0'
        }
      }
    }
  }

  return(data.frame(b.score))
}
