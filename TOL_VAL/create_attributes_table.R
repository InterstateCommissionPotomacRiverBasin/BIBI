# Set working directory
setwd("C:\\Users\\zsmith\\Desktop\\ICPRB_PhotoVault\\Chessie_BIBI\\R _Script\\TOL_VAL")

# load BIBI package
library(BIBI)

#==============================================================================
# Prepare the BIBI master taxa list
#==============================================================================
# Import master taxa list

master2 <- read.csv("master_just_taxa_10_11_16.csv")
#master$SPECIES <- ifelse(master$GENUS != "" & master$SPECIES != "",
#                     apply(master[ , c("GENUS", "SPECIES")] , 1 ,
#                           paste , collapse = "_" ), "")

# use the fill_taxa function from the BIBI package to get the final_id
# of each row
new.master <- fill_taxa(master2)
# The final_id will be represented as Species. Create a final_id
# column on the original master list with new.master's Species column
master2$FINAL_ID <- new.master$SPECIES
test <- master2[duplicated(master2$TSN_FINAL),]
master3 <- unique(master2[, !(names(master2) %in% c("TSN_R", "TSN_FINAL", "TSN_DB",
                                                  "TSN_AGENCY_ID", "AGENCY_ID", "SYN",
                                                  "VALIDATION", "NOTES", "TAXON_LEVEL"))])
test2 <- master3[duplicated(master3$FINAL_ID), ]
#==============================================================================
old_att <- read.csv("old_attributes_to_merge.csv")
#names(old_att)[1] <- "TSN_R"
#old_att <- old_att[, !(names(old_att) %in% "FAMILY")]
#test <- old_att[duplicated(old_att$TSN_R),]
old_att <- old_att[, !(names(old_att) %in% "TSN")]
master2 <- merge(master2, old_att, by = "FINAL_ID", all.x = T)
test <- master2[duplicated(master2$TSN_FINAL), ]
#==============================================================================
fill_taxa2 <- function(taxon_list){
  taxon_list[taxon_list == ""] <- NA
  new <- data.frame(t(apply(taxon_list, 1, zoo::na.locf, na.rm = F)))
  return(new)
}
#==============================================================================
# Prepare NYSDEC master taxa / taxa attributes table
#==============================================================================
# Import NYSDEC table
nysdec <- read.csv("nysdec_attributes_2016.csv")
# All text to uppercasse
nysdec <- data.frame(apply(nysdec, 2, FUN = toupper))
remove <- c( "CF\\.", "/")
nysdec <- nysdec[!grepl(paste(remove, collapse = "|"), nysdec$GENSPECIES),]
# Split GENSPECIES column to just represent Genus
nysdec$GENUS <- gsub(" .*$", "", nysdec$GENSPECIES)
nysdec$GENUS <- sapply(nysdec$GENUS, function(x){
  remove <- c("SP.", "CF\\.", "UNDET.", "UNDETERMINED", "/", "\\?")
  ifelse(grepl(paste(remove, collapse="|"), x), "", paste(x))
})
# Remove all genera, groups, undeterminds, complexes, and uncertainties
nysdec$SPECIES <- sapply(nysdec$GENSPECIES, function(x){
  remove <- c("SP\\.", "CF\\.")
  ifelse(grepl(paste(remove, collapse="|"), x), "", paste(x))
})
remove <- c("UNDET\\.", "UNDETERMINED", "/")
nysdec <- nysdec[!(grepl(paste(remove, collapse="|"), nysdec$GENSPECIES)), ]
# Remove any text contained within parentheses
nysdec$SPECIES <- gsub("\ \\([^\\)]+\\)","", nysdec$SPECIES)
# Remove NR.
nysdec$SPECIES <- gsub("\ NR\\.","", nysdec$SPECIES)
# Remove GR.
nysdec$SPECIES <- gsub("\ GR\\.","", nysdec$SPECIES)
# Remove ?
nysdec$SPECIES <- gsub("\\?","", nysdec$SPECIES)
# Replace the space between genus and species with "_"
nysdec$SPECIES <- gsub(" ","_", nysdec$SPECIES)
# Reorder the dataframe
nysdec <- nysdec[, c("PHYLUM", "CLASS", "ORDER", "FAMILY",
                     "SUBFAMILY", "GENUS", "SPECIES", "GENSPECIES",
                     "TOLERANCE", "FEEDINGHAB", "NBI.P_TOLERANCE",
                     "NBI.N_TOLERANCE")]

nysdec.fill <- fill_taxa2(nysdec)


nysdec$FINAL_ID <- nysdec.fill$SPECIES

nysdec.final <- nysdec[, c("FINAL_ID", "TOLERANCE", "FEEDINGHAB",
                                "NBI.P_TOLERANCE","NBI.N_TOLERANCE")]

names(nysdec.final) <- c("FINAL_ID", "NYSDEC_TV", "NYSDEC_FFG",
                         "NYSDEC_NBI.P", "NYSDEC_NBI.N")

nysdec.final2 <- unique(nysdec.final)
ny <- nysdec.final2[, - which(names(nysdec.final2) %in% "NYSDEC_FFG")]
ny$NYSDEC_TV <- as.numeric(as.character(levels(ny$NYSDEC_TV)))[ny$NYSDEC_TV]
ny$NYSDEC_NBI.P <- as.numeric(as.character(levels(ny$NYSDEC_NBI.P)))[ny$NYSDEC_NBI.P]
ny$NYSDEC_NBI.N <- as.numeric(as.character(levels(ny$NYSDEC_NBI.N)))[ny$NYSDEC_NBI.N]
agg_ny <- aggregate(ny[,2:ncol(ny)] , by = list(ny$FINAL_ID), data = ny, FUN = mean, na.rm = T)
names(agg_ny)[1] <- "FINAL_ID"
agg_ny$NYSDEC_TV[is.nan(agg_ny$NYSDEC_TV)] <- ""
agg_ny$NYSDEC_NBI.P[is.nan(agg_ny$NYSDEC_NBI.P)] <- ""
agg_ny$NYSDEC_NBI.N[is.nan(agg_ny$NYSDEC_NBI.N)] <- ""
nysdec.final3 <- merge(agg_ny, unique(nysdec.final2[, c("FINAL_ID", "NYSDEC_FFG")]),
                  by = "FINAL_ID")
#colnames(nysdec.final3) <- paste("NYSDEC", colnames(nysdec.final3), sep = "_")

merged <- merge(master2, nysdec.final3, by = "FINAL_ID", all.x = T)


test <- nysdec.final3[duplicated(nysdec.final3$FINAL_ID),]
#to.print <- data.frame(merged[merged$TSN_R %in% NA, 1])
#names(to.print) <- "NYSDEC_TAXA"
#write.csv(to.print, "NYSDEC_TAXA_TO_ADD.csv")
#==============================================================================
# Prepare DC taxa attributes table
#==============================================================================
dc <- read.csv("DC_attributes.csv")
# All text to uppercasse
dc <- data.frame(apply(dc, 2, FUN = toupper))
names(dc) <- toupper(names(dc))

# Remove all genera, groups, undetermineds, complexes, and uncertainties
dc.cols <- c("SUBPHYLUM_CLASS", "ORDER", "FAMILY", "GENUS")
dc[,dc.cols] <- sapply(dc[, dc.cols], function(x){
  remove <- c("SP//.", "GR//.", "CF//.", "UNDET//.", "UNDETERMINED", "/", "\\?")
  ifelse(grepl(paste(remove, collapse="|"), x), "", paste(x))
})

dc.fill <- fill_taxa2(dc)

# Replace the space between genus and species with "_"
dc$FINAL_ID <- dc.fill$GENUS
dc.col <- c("FINAL_ID", "DC_TV", "DC_FFG", "DC_HABIT")
dc.final <- dc[, dc.col]

test <- dc.final[duplicated(dc.final$FINAL_ID),]

merged2 <- merge(merged, dc.final, by = "FINAL_ID", all.x = T)
#==============================================================================
# Prepare EPA RBP taxa attributes table
#==============================================================================
rbp <- read.csv("EPA_RBP.csv")
# All text to uppercasse
rbp <- data.frame(apply(rbp, 2, FUN = toupper))
names(rbp) <- toupper(names(rbp))
# Replace the space between genus and species with "_"
rbp$FINAL_ID <- gsub(" ","_", rbp$FINAL_ID)
rbp.col <- c("FINAL_ID", "SOUTHEAST_NC_TV", "UPPER_MIDWEST_WI_TV", "MIDWEST_OH_TV",
             "NORTHWEST_ID_TV", "MID_ATLANTIC_MACS_TV", "P_FFG", "S_FFG",
             "P_HABIT", "S_HABIT")
rbp$P_FFG <- as.character(ifelse(rbp$P_FFG %in% "GC", "CG", as.character(rbp$P_FFG)))
rbp$S_FFG <- as.character(ifelse(rbp$S_FFG %in% "GC", "CG", as.character(rbp$S_FFG)))
rbp$P_FFG <- as.character(ifelse(rbp$P_FFG %in% "FC", "CF", as.character(rbp$P_FFG)))
rbp$S_FFG <- as.character(ifelse(rbp$S_FFG %in% "FC", "CF", as.character(rbp$S_FFG)))
rbp$P_FFG <- as.character(ifelse(rbp$P_FFG %in% "FG", "CF", as.character(rbp$P_FFG)))
rbp$S_FFG <- as.character(ifelse(rbp$S_FFG %in% "FG", "CF", as.character(rbp$S_FFG)))
rbp.final <- rbp[, rbp.col]
colnames(rbp.final) <- paste("RBP", colnames(rbp.final), sep = "_")
colnames(rbp.final)[1] <- "FINAL_ID"

test <- rbp.final[duplicated(rbp.final$FINAL_ID),]


merged3 <- merge(merged2, rbp.final, by = "FINAL_ID", all.x = T)
test <- merged3[duplicated(merged3$TSN_FINAL),]
#==============================================================================
# WVDEP Traits
wvdep <- read.csv("WVDEP_Traits_2016.csv")
wvdep$FINAL_ID <- toupper(wvdep$FINAL_ID)
names(wvdep) <- toupper(names(wvdep))
# Remove leading or trailing spaces from the TAXON column
wvdep$FINAL_ID <- gsub("^\\s+|\\s+$", "", wvdep$FINAL_ID)
wvdep$FINAL_ID <- gsub(" ","_", wvdep$FINAL_ID)
# Remove taxa names that contain "/" (complex) or "()" (essentially a duplicate)
wvdep <- wvdep[!(grepl("/", wvdep$FINAL_ID)), ]
wvdep <- wvdep[!(grepl(")", wvdep$FINAL_ID)), ]

merged4 <- merge(merged3, wvdep, by = "FINAL_ID", all.x = T)
#==============================================================================
# PADEP Traits
padep <- read.csv("PADEP_Modified.csv")
# Remove leading or trailing spaces from the TAXON column
padep$FINAL_ID <- gsub("^\\s+|\\s+$", "", padep$FINAL_ID)

merged5 <- merge(merged4, padep, by = "FINAL_ID", all.x = T)

test <- merged5[duplicated(merged5$TSN_FINAL),]

#==============================================================================
# Prepare EDAS taxa attributes table
#==============================================================================
edas <- read.csv("EDAS_attributes.csv")
# All text to uppercasse
edas <- data.frame(apply(edas, 2, FUN = toupper))
# Remove all genera, groups, undetermineds, complexes, and uncertainties
edas$FINAL_ID <- sapply(edas$FINAL_ID, function(x){
  remove <- c("SP//.", "GR//.", "CF//.", "UNDET//.", "UNDETERMINED", "/", "\\?")
  ifelse(grepl(paste(remove, collapse="|"), x), "", paste(x))
})
edas <- edas[edas$FINAL_ID != "", ]
# Replace the space between genus and species with "_"
edas$FINAL_ID <- gsub(" ","_", edas$FINAL_ID)

edas.col <- c("FINAL_ID", "EDAS_TV_FINAL", "LIMESTONE_IBI_TV", "WAB_TV",
              "RBP_TV", "MACS_TV",
              "PA_TV", "MD_TV", "NC_TV", "DE_TV", "WV_TV",
              "ITIS_TV", "VA_TV_MAIS", "EDAS_TV", "EDAS_FAM_TV",
              "KY_TV", "FFG_EDAS_FINAL", "EDAS_TROPHIC_GROUP", "KY_FFG", "EDAS_FFG",
              "ITIS_FFG", "MAIS_DICTIONARY_FFG", "WV_FFG", "HABIT_EDAS_FINAL",
              "KY_HABIT", "EDAS_HABIT", "WV_HABIT")

edas.final <- edas[, edas.col]
colnames(edas.final) <- paste("EDAS", colnames(edas.final), sep = "_")
colnames(edas.final)[1] <- "FINAL_ID"


#merged4 <- merge(merged3, edas.final, by = "FINAL_ID", all.x = T)
#==============================================================================
# Prepare EPA taxa attributes table from Greg Pond
# NRSA values
#==============================================================================
epa <- read.csv("epa_taxa_attributes.csv")
# All text to uppercasse
epa <- data.frame(apply(epa, 2, FUN = toupper))
epa.col <- c("FINAL_ID", "EPA_FFG", "EPA_HABIT", "EPA_TV",
             "FFG_WSA", "HABIT_WSA", "WSA_TV")
epa.final <- epa[, epa.col]
merged6 <- merge(merged5, epa.final, by = "FINAL_ID", all.x = T)

test <- merged5[!(merged6$TSN_R %in% NA),]

check <- data.frame(apply(test, 2, FUN = max))

#==============================================================================
#==============================================================================
#==============================================================================
# Prep USGS
usgs <- read.csv("EPA_2010.csv")
names(usgs) <- toupper(names(usgs))
usgs$TAXON <- toupper(usgs$TAXON)
usgs <- usgs[usgs$TAXON != "", ]

usgs$TAXON <- sapply(usgs$TAXON, function(x){
  remove <- c(" SPP\\.", "GR\\.", "CF\\.", "UNDET\\.", "UNDETERMINED", "/", "\\?")
  ifelse(grepl(paste(remove, collapse="|"), x),
         gsub(paste(remove, collapse="|"), "", x), paste(x))
})
# Remove leading or trailing spaces from the TAXON column
usgs$TAXON <- gsub("^\\s+|\\s+$", "", usgs$TAXON)
usgs$TAXON <- gsub(" ","_", usgs$TAXON)

#===========================
# Aggregate numeric columns
usgs_numeric <- cbind(usgs[, "TAXON"], usgs[sapply(usgs, is.numeric)])
names(usgs_numeric)[1] <- "TAXON"

agg_numeric <- aggregate(usgs_numeric[, 2:ncol(usgs_numeric)],
                         by = list(usgs_numeric$TAXON), data = usgs_numeric,
                         FUN = mean, na.rm = T)
names(agg_numeric)[1] <- "TAXON"

#binary <- c("MEDIATE_DRAG", "EMERGE_SEASON_ALL_YEAR", "EMERGE_SYNCH", "EGGS_CEMENT",
#            "EXIT_TEMPORARILY", "DIAPAUSE")
#agg_numeric[, binary] <- round(agg_numeric[, binary], digits = 0)
new <- data.frame(names(agg_numeric[, 2:ncol(agg_numeric)]))
new$tmax <- apply(agg_numeric[,2:ncol(agg_numeric)], 2, max, na.rm =T)

new$tmin <- apply(agg_numeric[,2:ncol(agg_numeric)], 2, min, na.rm =T)

test <- new[new$tmax != new$tmin,]
#===========================
# Aggregate non-numeric columns
usgs$TAXON <- as.factor(usgs$TAXON)
usgs_factor <- data.frame(apply(usgs[sapply(usgs,is.factor)], 2, FUN = toupper))
uf <- usgs_factor
#cols.remove <- c("TRAITRECORD_ID", "FAMILY", "GENUS", "STUDY_LOCATION_STATE",
#                 "STUDY_LOCATION_COUNTY", "STUDY_LOCATION_REGION", "STUDY_LATITUDE",
#                 "STUDY_LONGITUDE", "STUDY_DATES", "DATA_ENTRY", "ADULT", "DATA_ENTRY_DATE")
#uf <- uf[, -which(names(uf) %in% cols.remove)]
uf2 <- apply(uf, 2, FUN = as.character)

agg <- aggregate(uf2 ~ TAXON, data = uf2, function(x) paste0(unique(x), collapse = ","))
remove <- c(",,", ",,,", ",,,,", ",,,,,", ",,,,,,", ",,,,,,,", ",,,,,,,,", ",,,,,,,,,",
            ",,,,,,,,,,")
sub <- function(x)gsub(paste0(remove, collapse = "|"),",", x)
agg2 <- apply(agg, 2, function(x) sub(x))

test <- apply(agg2, 2, function (x) ifelse(substr(x, 1, 1) %in% ",", substr(x, nchar(",") + 1, nchar(x)), x ))
test2 <- apply(test, 2, function (x) ifelse(substr(x, nchar(x), nchar(x)) %in% ",", substr(x, 1, nchar(x)-1), x ))
new2 <- data.frame(test2)
new2 <- new2[, !(names(new2) %in% c("TAXON.1", "ORDER", "FAMILY", "GENUS"))]

usgs_final <- merge(agg_numeric, new2, by = "TAXON")
names(usgs_final)[1] <- "FINAL_ID"

merged7 <- merge(merged6, usgs_final, by = "FINAL_ID", all.x = T)

test <- merged7[duplicated(merged7$FINAL_ID),]
test <- merged7[duplicated(merged7$TSN_FINAL),]
#==============================================================================
# Prep USGS 2007
usgs_fam <- read.csv("USGS_2007_FAMILY_TOLERANCE.csv")
usgs_fam$FAMILY <- toupper(usgs_fam$FAMILY)
names(usgs_fam)[1] <- "FINAL_ID"
usgs_gen <- read.csv("USGS_2007_GENUS_TOLERANCE.csv")
usgs_gen$GENUS <- toupper(usgs_gen$GENUS)
names(usgs_gen)[1] <- "FINAL_ID"

merged8 <- merge(merged7, usgs_fam, by = "FINAL_ID", all.x = T, all.y = T)

merged9 <- merge(merged8, usgs_gen, by = "FINAL_ID", all.x = T, all.y = T)


#==============================================================================
# Only taxa from the BIBI data base
final.df <- merged7[complete.cases(merged7$TSN_R), ]

ord <- final.df[final.df$ORDER == final.df$FINAL_ID,]
fam <- final.df[final.df$FAMILY == final.df$FINAL_ID,]
gen <- final.df[final.df$GENUS == final.df$FINAL_ID,]

final.df$test <- apply(final.df, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )

test1 <- final.df[, c("FINAL_ID", "test")]
#==============================================================================
nysdec.num <- c("NYSDEC_TOLERANCE", "NYSDEC_NBI.P", "NYSDEC_NBI.N")
as.numeric.factor <- function(x) {as.numeric(as.character(unlist(x)))}
nysdec.final[, nysdec.num] <- apply(nysdec.final[, nysdec.num], 2, FUN = as.numeric.factor)
check <- apply(nysdec.final[, nysdec.num], 2, FUN = max, na.rm = T)
check <- apply(nysdec.final[, nysdec.num], 2, FUN = min, na.rm = T)


#==============================================================================
# Mean tolerance values
TV2 <- final.df[, 1:17]
TV3 <- final.df[, grepl("_TV", names(final.df))]
TV3 <- TV3[, !(names(TV3) %in% c("RBP_UPPER_MIDWEST_WI_TV", "RBP_NORTHWEST_ID_TV",
                                 "EDAS_LIMESTONE_IBI_TV", "EDAS_EDAS_TV",
                                 "EDAS_MACS_TV", "EDAS_VA_TV_MAIS", "WSA_TV"))]
as.numeric.factor <- function(x) {as.numeric(as.character(unlist(x)))}
TV4 <- apply(TV3, 2, FUN = as.numeric.factor)


TV <- cbind(TV2, TV4)
TV$BIBI_TV <- rowMeans(TV[, 18:ncol(TV)], na.rm = TRUE)
TV$BIBI_TV <- ifelse(!is.na(TV$BIBI_TV), round(TV$BIBI_TV, 0), NA)

write.csv(table(TV$BIBI_TV), "TV_COUNT.csv")

write.csv(apply(TV[, 18:ncol(TV)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(TV)) * 100)), "TV_TABLE.csv")

TV.df <- data.frame(ALL = apply(TV[, 18:ncol(TV)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(TV)) * 100)))

tv_ord <- TV[TV$FINAL_ID %in% ord$FINAL_ID,]
#tv_fam <- tv_fam[, !(names(tv_fam) %in% "BIBI_TV")]
table(tv_ord$BIBI_TV)
TV.df$ORD <- apply(tv_ord[, 18:ncol(tv_ord)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(tv_ord)) * 100))

tv_fam <- TV[TV$FINAL_ID %in% fam$FINAL_ID,]
#tv_fam <- tv_fam[, !(names(tv_fam) %in% "BIBI_TV")]
table(tv_fam$BIBI_TV)
TV.df$FAM <- apply(tv_fam[, 18:ncol(tv_fam)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(tv_fam)) * 100))

tv_gen <- TV[TV$FINAL_ID %in% gen$FINAL_ID,]
#tv_fam <- tv_fam[, !(names(tv_fam) %in% "BIBI_TV")]
table(tv_gen$BIBI_TV)
TV.df$GEN <- apply(tv_gen[, 18:ncol(tv_gen)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(tv_gen)) * 100))
TV.df <- round(TV.df)
TV.df <- TV.df[order(TV.df$ALL),]
write.csv(TV.df, "TV_table_10_3_16.csv")
#==============================================================================
# Most frequent FFG

FFG2 <- final.df[, 1:17]
FFG3 <- final.df[, grepl("_FFG", names(final.df))]
FFG3 <- FFG3[, !(names(FFG3) %in% c("RBP_S_FFG", "EDAS_EDAS_FFG",
                                  "EDAS_MAIS_DICTIONARY_FFG"))]


FFG <- cbind(FFG2, FFG3)
FFG[, 18:23] <- apply(FFG[, 18:23], 2, FUN = as.character)
FFG[FFG == ""] <- NA
FFG$CG <- apply(FFG[, 18:23], 1, function(x) sum(grepl("CG", x)))
FFG$CF <- apply(FFG[, 18:23], 1, function(x) sum(grepl("CF", x)))
FFG$SH <- apply(FFG[, 18:23], 1, function(x) sum(grepl("SH", x)))
FFG$SC <- apply(FFG[, 18:23], 1, function(x) sum(grepl("SC", x)))
FFG$PR <- apply(FFG[, 18:23], 1, function(x) sum(grepl("PR", x)))
FFG$PC <- apply(FFG[, 18:23], 1, function(x) sum(grepl("PC", x)))
FFG$PI <- apply(FFG[, 18:23], 1, function(x) sum(grepl("PH", x)) + sum(grepl("HB", x)) + sum(grepl("PI", x)))
FFG$OM <- apply(FFG[, 18:23], 1, function(x) sum(grepl("OM", x)))
FFG$PA <- apply(FFG[, 18:23], 1, function(x) sum(grepl("PA", x)) + sum(grepl("PS", x)))

FFG$NUM_MAX <- apply(FFG[, 24:32], 1, function(x) sum(x==max(x)))
FFG$SUM <- apply(FFG[, 24:32], 1, function(x) sum(x))

FFG$BIBI_FFG <- ifelse(FFG$NUM_MAX > 1 & FFG$SUM == 0, NA,
                      ifelse(FFG$NUM_MAX > 1 & FFG$SUM > 0, 
                             apply(FFG[, 24:32], 1, function(x) paste0(names(which(x == max(x))), collapse = ",")),
                             ifelse(FFG$NUM_MAX == 1, 
                                    apply(FFG[, 24:32], 1, function(x) names(which.max(x))), "ERROR")))
ffg_correct <- read.csv("FFG_CORRECTIONS.csv", colClasses =  c("character", "character"))
FFG$BIBI_FFG <- ifelse(FFG$NUM_MAX > 1 & FFG$FINAL_ID %in% ffg_correct$FINAL_ID,
                       ffg_correct$BIBI_FFG, FFG$BIBI_FFG)
test <- FFG[FFG$NUM_MAX > 1 & FFG$SUM > 0, ]

write.csv(table(FFG$BIBI_FFG), "FFG_COUNT.csv")
(nrow(test) / nrow(FFG)) * 100

write.csv(apply(FFG[, 18:ncol(FFG)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(FFG)) * 100)), "FFG_PCT.csv")
ffg.df <- data.frame(ALL = apply(FFG[, 18:ncol(FFG)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(FFG)) * 100)))
ffg_ord <- FFG[FFG$FINAL_ID %in% ord$FINAL_ID,]
#ffg_fam <- ffg_fam[, !(names(ffg_fam) %in% "BIBI_FFG")]
table(ffg_ord$BIBI_FFG)
ffg.df$ORD <- apply(ffg_ord[, 18:ncol(ffg_ord)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(ffg_ord)) * 100))

ffg_fam <- FFG[FFG$FINAL_ID %in% fam$FINAL_ID,]
#ffg_fam <- ffg_fam[, !(names(ffg_fam) %in% "BIBI_FFG")]
table(ffg_fam$BIBI_FFG)
ffg.df$FAM <- apply(ffg_fam[, 18:ncol(ffg_fam)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(ffg_fam)) * 100))

ffg_gen <- FFG[FFG$FINAL_ID %in% gen$FINAL_ID,]
#ffg_fam <- ffg_fam[, !(names(ffg_fam) %in% "BIBI_FFG")]
table(ffg_gen$BIBI_FFG)
ffg.df$GEN <- apply(ffg_gen[, 18:ncol(ffg_gen)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(ffg_gen)) * 100))

write.csv(ffg.df, "FFG_table_10_3_16.csv")
#==============================================================================
# Most frequent Habit

HABIT2 <- final.df[, 1:17]
HABIT3 <- final.df[, grepl("HABIT", names(final.df))]
#as.numeric.factor <- function(x) {as.numeric(unlist(x))}
#HABIT4 <- apply(HABIT3, 2, FUN = as.character)
HABIT4 <- HABIT3[, !(names(HABIT3) %in% c("HABIT_PRIM", "HABIT_SEC",
                                          "HABIT_COMMENTS", "RBP_S_HABIT",
                                          "EDAS_EDAS_HABIT", "HABIT_WSA"))]

HABIT <- cbind(HABIT2, HABIT4)
HABIT[, 18:22] <- apply(HABIT[, 18:22], 2, FUN = as.character)
HABIT[HABIT == ""] <- NA
HABIT$CB <- apply(HABIT[, 18:22], 1, function(x) sum(grepl("CB", x)))
HABIT$SW <- apply(HABIT[, 18:22], 1, function(x) sum(grepl("SW", x)))
HABIT$CN <- apply(HABIT[, 18:22], 1, function(x) sum(grepl("CN", x)))
HABIT$BU <- apply(HABIT[, 18:22], 1, function(x) sum(grepl("BU", x)))
HABIT$SP <- apply(HABIT[, 18:22], 1, function(x) sum(grepl("SP", x)))
HABIT$DV <- apply(HABIT[, 18:22], 1, function(x) sum(grepl("DV", x)))
HABIT$SK <- apply(HABIT[, 18:22], 1, function(x) sum(grepl("SK", x)))
HABIT$SE <- apply(HABIT[, 18:22], 1, function(x) sum(grepl("SE", x)))
HABIT$PL <- apply(HABIT[, 18:22], 1, function(x) sum(grepl("PL", x)))


HABIT$NUM_MAX <- apply(HABIT[,23:31], 1, function(x) sum(x==max(x)))
HABIT$SUM <- apply(HABIT[,23:31], 1, function(x) sum(x))

HABIT$BIBI_HABIT <- ifelse(HABIT$NUM_MAX > 1 & HABIT$SUM == 0, NA,
                      ifelse(HABIT$NUM_MAX > 1 & HABIT$SUM > 0, 
                             apply(HABIT[,23:31], 1, function(x) paste0(names(which(x == max(x))), collapse = ",")),
                             ifelse(HABIT$NUM_MAX == 1, 
                                    apply(HABIT[,23:31], 1, function(x) names(which.max(x))), "ERROR")))
habit_correct <- read.csv("HABIT_CORRECTIONS.csv", colClasses =  c("character", "character"))
HABIT$BIBI_HABIT <- ifelse(HABIT$NUM_MAX > 1 & HABIT$FINAL_ID %in% habit_correct$FINAL_ID,
                       habit_correct$BIBI_HABIT, HABIT$BIBI_HABIT)
test <- HABIT[HABIT$NUM_MAX > 1 & HABIT$SUM > 0, ]
write.csv(table(HABIT$BIBI_HABIT), "HABIT_COUNT.csv")
(nrow(test) / nrow(HABIT)) * 100

write.csv(apply(HABIT[, 18:ncol(HABIT)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(HABIT)) * 100)), "HABIT_PCT.csv")
habit.df <- data.frame(ALL = apply(HABIT[, 18:ncol(HABIT)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(HABIT)) * 100)))
habit_ord <- HABIT[HABIT$FINAL_ID %in% ord$FINAL_ID,]
#habit_fam <- habit_fam[, !(names(habit_fam) %in% "BIBI_HABIT")]
table(habit_ord$BIBI_HABIT)
habit.df$ORD <- apply(habit_ord[, 18:ncol(habit_ord)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(habit_ord)) * 100))

habit_fam <- HABIT[HABIT$FINAL_ID %in% fam$FINAL_ID,]
#habit_fam <- habit_fam[, !(names(habit_fam) %in% "BIBI_HABIT")]
table(habit_fam$BIBI_HABIT)
habit.df$FAM <- apply(habit_fam[, 18:ncol(habit_fam)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(habit_fam)) * 100))

habit_gen <- HABIT[HABIT$FINAL_ID %in% gen$FINAL_ID,]
#habit_fam <- habit_fam[, !(names(habit_fam) %in% "BIBI_HABIT")]
table(habit_gen$BIBI_HABIT)
habit.df$GEN <- apply(habit_gen[, 18:ncol(habit_gen)], 2, function(x) 100 - ((sum(is.na(x)) / nrow(habit_gen)) * 100))
write.csv(habit.df, "Habit_table_10_3_16.csv")
#==============================================================================
#==============================================================================
#FINAL STEPS
final.df$BIBI_TV <- TV$BIBI_TV
final.df$BIBI_FFG <- FFG$BIBI_FFG
final.df$BIBI_HABIT <- HABIT$BIBI_HABIT
final.df <- within(final.df, BIBI_HABIT[FINAL_ID %in% "ISONYCHIIDAE"] <- "SW")
final.df2 <- final.df
final.df2[is.na(final.df2)] <- ""
final.df2[final.df2 == ""] <- NA
test <- final.df2[duplicated(final.df2$TSN_FINAL),]
test2 <- final.df2[final.df2$FINAL_ID %in% "ISONYCHIIDAE", "BIBI_HABIT"]
test3 <- unique(final.df2[, c("FINAL_ID", "TSN_R", "PHYLUM", "SUBPHYLUM",
                       "CLASS", "SUBCLASS", "ORDER", "SUBORDER", "FAMILY",
                       "SUBFAMILY", "TRIBE", "GENUS", "SPECIES")])
test3 <- test3[duplicated(test3), ]

write.csv(final.df2, "Master_Taxa_List_10_11_16.csv", row.names = FALSE)
#==============================================================================
#==============================================================================




# Extract columns with max and min = 1
# This format should make it easy to create new functions for these metrics

#===========================
# Aggregate numeric columns
usgs_numeric <- cbind(usgs[, "TAXON"], usgs[sapply(usgs, is.numeric)])
names(usgs_numeric)[1] <- "TAXON"

agg_numeric <- aggregate(usgs_numeric[, 2:ncol(usgs_numeric)],
                         by = list(usgs_numeric$TAXON), data = usgs_numeric,
                         FUN = mean, na.rm = T)
names(agg_numeric)[1] <- "TAXON"

#binary <- c("MEDIATE_DRAG", "EMERGE_SEASON_ALL_YEAR", "EMERGE_SYNCH", "EGGS_CEMENT",
#            "EXIT_TEMPORARILY", "DIAPAUSE")
#agg_numeric[, binary] <- round(agg_numeric[, binary], digits = 0)
new <- data.frame(names(agg_numeric[, 2:ncol(agg_numeric)]))
names(new) <- "METRICS"
new$METRICS <- as.character(new$METRICS)
new$tmax <- apply(agg_numeric[,2:ncol(agg_numeric)], 2, max, na.rm =T)

new$tmin <- apply(agg_numeric[,2:ncol(agg_numeric)], 2, min, na.rm =T)

names.list <- ifelse(new$tmax == 1 & new$tmin == 1, new$METRICS, NA)

nl <- na.omit(names.list)

P_A2 <- final.df[, 1:17]
P_A3 <- final.df[, names(final.df) %in% nl]
P_A <- cbind(P_A2, P_A3)
