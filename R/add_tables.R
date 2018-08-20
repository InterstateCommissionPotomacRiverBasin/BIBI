
# setwd("//pike/data/Projects/Chessie_BIBI/BIBI_23Dec15")
# master <- read.csv("Master_Taxa_List_9_20_16.csv")
# master[master$TSN_DB %in% "", "TSN_DB"] <- NA
# master$TSN_R <- as.character(master$TSN_R)
# master$TSN_DB <- as.character(master$TSN_DB)
# master$X <- ifelse(is.na(master$TSN_DB), master$TSN_R, master$TSN_DB)
# master$TSN_FINAL <- gsub("BAY0", "BAY", master$TSN_FINAL)
# master$TSN_FINAL <- gsub("EMAP0", "EMAP", master$TSN_FINAL)
# master <- master[, !(names(master) %in% "X")]
# 
# trim <- function (x) gsub("^\\s+|\\s+$", "", x)
# taxa.col <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
#               "ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
#               "TRIBE", "GENUS", "SPECIES")
# master[, taxa.col] <- apply(master[, taxa.col], 2, trim)
# setwd("C:\\Users\\zsmith\\Desktop\\BIBI\\R\\bibi")
# devtools::use_data(master, overwrite = TRUE)



