#==============================================================================
# METRICS FOR ASSESSING TAXA BY PREDEFINED GROUPS
# Functional Feeding Groups, Habitat, Composition, etc.
#==============================================================================
#'Data Frame of a specific list of taxa
#'
#'@param NameList = uninque list of taxa.
#'@param Taxa.df = Wide data frame format of taxonomic counts.
#'@return Data frame of a specific list of taxa. Used in habit and functional
#'feeding group (FFG) fuctions to extract only the taxa specified in the object
#'NameList from the wide data frame of taxa (i.e. Family level or Genus level)
#'@export

group_taxa <- function(NameList, Taxa.df){
  ID <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER")
  taxa.list <- as.character(unlist(NameList))
  #taxa.list[which(c(1, diff(taxa.list)) != 0)]
  #idx <- match(taxa.list, names(Taxa.df))
  #idx <- idx[! is.na(idx)]
  taxa_list.df <- data.frame(Taxa.df[, names(Taxa.df) %in% c(ID, taxa.list)])
  taxa_list.df <- taxa_list.df[, !(names(taxa_list.df) %in% "UNIDENTIFIED")]
  taxa_list.df[is.na(taxa_list.df)] <- 0 #NA = zero
  final_taxa.df <- if(ncol(taxa_list.df) < 6) {
    0
  } else {
    if(ncol(taxa_list.df) == 6) {
      taxa_list.df[, 6]
    } else {
      if(ncol(taxa_list.df) > 6)
        rowSums(taxa_list.df[, 6:ncol(taxa_list.df)])
    }
  }
  return(final_taxa.df)
}

#==============================================================================
#'Vector of taxa richness for a specific list of taxa
#'
#'@param NameList = uninque list of taxa.
#'@param Taxa.df = Wide data frame format of taxonomic counts.
#'@return A vector of taxa richness for a specific list of taxa representing
#'each sampling event. NameList from the wide data frame of taxa
#' (i.e. Family level or Genus level)
#'@export

group_rich <- function(NameList, Taxa.df){
  ID <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER")
  taxa.list <- as.character(unlist(NameList))
  taxa_list.df <- data.frame(Taxa.df[, names(Taxa.df) %in% c(ID, taxa.list)])
  taxa_list.df <- taxa_list.df[, !(names(taxa_list.df) %in% "UNIDENTIFIED")]
  taxa_list.df[is.na(taxa_list.df)] <- 0 #NA = zero
  if(ncol(taxa_list.df) < 6) {
    final_taxa.df <- 0
  } else {
    final_taxa.df <- vegan::specnumber(taxa_list.df[, 6:ncol(taxa_list.df)])
  }
  return(final_taxa.df)
}
#==============================================================================
#Only taxa present in the data are used to create each data frame
#Therefore, absent taxa do not have a column of counts for each station
#The metrics require a value for each variable called on
#This function asks if a taxon is present in the data frame and thus present in at least one of your stations
#If the taxon is not present a temporary column of zeros is formed to represent the taxon for a specific metric
#'Temporary vector with zero values
#'

#'@param Taxon = The name of the taxon necessary for calculating a metric but
#'is not found in the data set.
#'@param Level = Taxanomic level (e.g. Class, Order, Family, Genus, etc.).
#'@return Creates a temporary vector containing zero values the length of the
#'data frame specified by the object Level.  Prevents errors created by missing taxon names
#'necessary for metric calculations. For example, if no Trichoptera taxa were observed
#'in the data set then there would be no column named "Trichoptera" in the wide order
#'level data frame.  When metrics, such as pct_ept or pct_trichoptera, search for a column
#'named "Trichoptera" a null value is returned and the function fails. Therefore,
#'this function temporarly fills the missing columns with zeros and the metric
#'can be calculated.
#'@export

blank_col <- function(Taxon, Level){
  if(Taxon %in% colnames(Level)) Level[, Taxon] else rep(0, nrow(Level))
}

#==============================================================================
#'The percent of a group
#'
#'@param Taxa.df = Wide data frame format of taxonomic counts.
#'@param Info = Taxonomic attributes table.
#'@param Group = The taxonomic group to be assessed
#'@param Group_Level = The specific level or levels of the group to be assessed.
#'@param Level = The taxonomic level used during the assessment.
#'@return The percentage of taxa representing a predefined group. Typically,
#'this function is used to assess functional feeding group and habits.
#'@export

pct_attribute <- function(Taxa.df, Info, Group, Group_Level, Level = "FAMILY"){
  #split.taxa <- split(Info[, Level], Info[, Group])
  #name.list <- split.taxa[Group_Level]
  new.group <- c(Group_Level)
  grep.taxa <- Info[grepl(paste(new.group,collapse="|"), Info[, Group]), ]
  name.list <- as.list(unique(grep.taxa$FINAL_ID))
  group.taxa <- group_taxa(name.list, Taxa.df)
  if(sum(group.taxa) == 0){
    final.vec <- 0
  }else{
    final.vec <- (group.taxa / rowSums(Taxa.df[, 6:ncol(Taxa.df)])) * 100
  }
  return(final.vec)
}

#==============================================================================
#'The richness of a group
#'
#'@param taxa.wide = Wide data frame format of taxonomic counts.
#'@param master = Taxonomic attributes table.
#'@param attribute.column = The name of the column that contains the
#'attribute of interest.
#'@param attribute.interest = The specific attribute of interest
#'@param rank = The taxonomic level used during the assessment.
#'@return The richness of taxa representing a predefined group. Typically,
#'this function is used to assess functional feeding group and habits.
#'@export

rich_attribute <- function(taxa.wide, master = BIBI::master, attribute.column,
                           attribute.interest, rank = "FAMILY"){
  new.group <- c(attribute.interest)
  grep.taxa <- master[grepl(paste(new.group, collapse="|"), master[, attribute.column]), ]
  name.list <- as.list(grep.taxa$FINAL_ID)
  ID <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER")
  taxa.list <- as.character(unlist(name.list))
  group.rich <- group_rich(name.list, taxa.wide)
  final_taxa.vec <- group.rich
  #taxa_list.df <- data.frame(taxa.wide[, names(taxa.wide) %in% c(ID, taxa.list)])
  #taxa_list.df[is.na(taxa_list.df)] <- 0 #NA = zero
  #if(ncol(taxa_list.df) < 6) {
  #  final_taxa.wide <- 0
  #} else {
  #  final_taxa.wide <- vegan::specnumber(taxa_list.df[, 6:ncol(taxa_list.df)])
  #}
  return(final_taxa.vec)
}

#==============================================================================
#'The percent of the most dominant group
#'
#'@param Long = Long data frame format of taxonomic counts.
#'@param Info = Taxonomic attributes table.
#'@param Group = The taxonomic group to be assessed
#'@param Level = The taxonomic level used during the assessment.
#'@return The percentage of taxa representing by the most dominant (abundant)
#'group. Typically, this function is used to assess the functional feeding
#'groups and habits.
#'@export

pct_dom1_group <- function(Long, Info, Group, Level){
  taxa.info <- Info[, c("FINAL_ID", Group)]
  merged <- merge(Long, taxa.info, by.x = Level, by.y = "FINAL_ID",
                   all.x = TRUE)
  wide.df <- wide(merged, Group)
  final.df <- pct_dom(wide.df, 1)
  return(final.df)
}
