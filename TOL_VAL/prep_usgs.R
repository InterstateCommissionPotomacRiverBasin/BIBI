# # Set working directory
# setwd("C:\\Users\\Zsmith\\Desktop\\Chessie_BIBI\\R _Script\\TOL_VAL")
# 
# # load BIBI package
# library(BIBI3)
# 
# #==============================================================================
# #==============================================================================
# # Prep USGS
# usgs <- read.csv("EPA_2010.csv")
# names(usgs) <- toupper(names(usgs))
# usgs$TAXON <- toupper(usgs$TAXON)
# usgs <- usgs[usgs$TAXON != "", ]
# 
# usgs$TAXON <- sapply(usgs$TAXON, function(x){
#   remove <- c(" SPP\\.", "GR\\.", "CF\\.", "UNDET\\.", "UNDETERMINED", "/", "\\?")
#   ifelse(grepl(paste(remove, collapse="|"), x),
#          gsub(paste(remove, collapse="|"), "", x), paste(x))
# })
# # Remove leading or trailing spaces from the TAXON column
# usgs$TAXON <- gsub("^\\s+|\\s+$", "", usgs$TAXON)
# usgs$TAXON <- gsub(" ","_", usgs$TAXON)
# 
# #===========================
# # Aggregate numeric columns
# usgs_numeric <- cbind(usgs[, "TAXON"], usgs[sapply(usgs, is.numeric)])
# names(usgs_numeric)[1] <- "TAXON"
# 
# agg_numeric <- aggregate(usgs_numeric[, 2:ncol(usgs_numeric)],
#                          by = list(usgs_numeric$TAXON), data = usgs_numeric,
#                          FUN = mean, na.rm = T)
# names(agg_numeric)[1] <- "TAXON"
# 
# #binary <- c("MEDIATE_DRAG", "EMERGE_SEASON_ALL_YEAR", "EMERGE_SYNCH", "EGGS_CEMENT",
# #            "EXIT_TEMPORARILY", "DIAPAUSE")
# #agg_numeric[, binary] <- round(agg_numeric[, binary], digits = 0)
# new <- data.frame(names(agg_numeric[, 2:ncol(agg_numeric)]))
# new$tmax <- apply(agg_numeric[,2:ncol(agg_numeric)], 2, max, na.rm =T)
# 
# new$tmin <- apply(agg_numeric[,2:ncol(agg_numeric)], 2, min, na.rm =T)
# 
# test <- new[new$tmax != new$tmin,]
# #===========================
# # Aggregate numeric columns
# usgs$TAXON <- as.factor(usgs$TAXON)
# usgs_factor <- data.frame(apply(usgs[sapply(usgs,is.factor)], 2, FUN = toupper))
# uf <- usgs_factor
# #cols.remove <- c("TRAITRECORD_ID", "FAMILY", "GENUS", "STUDY_LOCATION_STATE",
# #                 "STUDY_LOCATION_COUNTY", "STUDY_LOCATION_REGION", "STUDY_LATITUDE",
# #                 "STUDY_LONGITUDE", "STUDY_DATES", "DATA_ENTRY", "ADULT", "DATA_ENTRY_DATE")
# #uf <- uf[, -which(names(uf) %in% cols.remove)]
# uf2 <- apply(uf, 2, FUN = as.character)
# 
# agg <- aggregate(uf2 ~ TAXON, data = uf2, function(x) paste0(unique(x), collapse = ","))
# remove <- c(",,", ",,,", ",,,,", ",,,,,", ",,,,,,", ",,,,,,,", ",,,,,,,,", ",,,,,,,,,",
#             ",,,,,,,,,,")
# sub <- function(x)gsub(paste0(remove, collapse = "|"),",", x)
# agg2 <- apply(agg, 2, function(x) sub(x))
# 
# test <- apply(agg2, 2, function (x) ifelse(substr(x, 1, 1) %in% ",", substr(x, nchar(",") + 1, nchar(x)), x ))
# test2 <- apply(test, 2, function (x) ifelse(substr(x, nchar(x), nchar(x)) %in% ",", substr(x, 1, nchar(x)-1), x ))
# new2 <- data.frame(test2)
# 
# merged <- merge(agg_numeric, new2, by = "TAXON")
# 
# 
# 
