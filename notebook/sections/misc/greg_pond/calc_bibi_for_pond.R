library(tidyverse)
library(mmir)
data(master)
#master <- master %>% 
#  rename_all(tolower) %>% 
#  mutate_if(is.factor, as.character) %>% 
#  mutate_if(is.character, tolower)
#------------------------------------------------------------------------------
rsc_200.df <- readxl::read_excel("D:/ZSmith/Projects/Chessie_BIBI/greg_pond/rsc_taxa_10_16_17.xlsx",
                                 sheet = "200_count") %>% 
  mutate(type = "200")

rsc_120.df <- readxl::read_excel("D:/ZSmith/Projects/Chessie_BIBI/greg_pond/rsc_taxa_10_16_17.xlsx",
                                 sheet = "120_count") %>% 
  mutate(type = "120")

rsc.df <- bind_rows(rsc_200.df, rsc_120.df) %>% 
  rename(sample_id = "Sample ID",
         orginal_id = "Final ID",
         count = "Count",
         tsn = "Taxonomic Serial Number") %>% 
  select(sample_id, type, orginal_id, count, tsn) %>% 
  mutate(tsn = as.character(tsn),
         bioregion = "sep", 
         unique_id = paste(sample_id, type, sep = "_")) %>% 
  select(unique_id, everything())
rm(rsc_200.df, rsc_120.df)
#------------------------------------------------------------------------------
check.dubs <- master %>% 
  select(FINAL_ID, TSN_FINAL, PHYLUM:SPECIES, BIBI_TV, BIBI_FFG, BIBI_HABIT, BECK_CLASS, ASPT) %>% 
  distinct() %>% 
  group_by(TSN_FINAL) %>% 
  summarize(n = n())
rsc.df <- master %>% 
  select(FINAL_ID, TSN_FINAL, PHYLUM:SPECIES, BIBI_TV, BIBI_FFG, BIBI_HABIT, BECK_CLASS, ASPT) %>% 
  distinct() %>% 
  right_join(rsc.df, by = c("TSN_FINAL" = "tsn")) %>% 
  BIBI::clean_taxa()
#------------------------------------------------------------------------------
check.gen.pct <- rsc.df %>% 
  group_by(unique_id) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  group_by(unique_id, GENUS, total) %>% 
  summarize(count = sum(count)) %>% 
  mutate(percent = count / total * 100) %>% 
  filter(GENUS == "UNIDENTIFIED",
         percent > 10)
#------------------------------------------------------------------------------
taxa.rank <- rlang::quo(GENUS)

att.df <- master %>% 
  select(FINAL_ID, BIBI_TV, BIBI_FFG, BIBI_HABIT, ASPT) %>% 
  distinct()

test3 <- rsc.df %>% 
  select(unique_id) %>% 
  distinct() %>% 
  mutate(CT_CAECIDOTEA = taxa_pct(rsc.df, unique_id, count, GENUS, "CAECIDOTEA"),
         PCT_COLLECT = taxa_pct(rsc.df, att.df, unique_id, count, BIBI_FFG, c("CG", "CF"),
                                att.df, "FINAL_ID")),
         PCT_COLLECT_NO_INTOL = taxa_pct(rsc.df, unique_id, count, BIBI_FFG, c("CG", "CF"),
                                         exclusion.col = BIBI_TV, exclusion.vec = 1:7))
  
  
  mutate(PCT_COLLECT = taxa_pct(rsc.df, unique_id, count, BIBI_FFG, c("CG", "CF")),
         PCT_COLLECT_NO_INTOL = taxa_pct(rsc.df, unique_id, count, BIBI_FFG, c("CG", "CF"),
                                          exclusion.col = BIBI_TV, exclusion.vec = 1:7),
         test = taxa_rich(rsc.df, unique_id, BIBI_FFG, !!taxa.rank, "CF", 
                          exclusion.col = BIBI_TV, exclusion.vec = 1:7),
         test2 = taxa_rich(rsc.df, unique_id, BIBI_FFG, !!taxa.rank, "CF"),
         PCT_EPT_RICH_NO_TOL = taxa_pct_rich(rsc.df, unique_id, ORDER, !!taxa.rank,
                                             taxon = c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA"),
                                             exclusion.col = BIBI_TV, exclusion.vec = 7:10),
         PCT_EPT_RICH_NONE = taxa_pct_rich(rsc.df, unique_id, ORDER, !!taxa.rank,
                                             taxon = c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA"),
                                             exclusion.col = BIBI_TV, exclusion.vec = c(0:10, NA)),
         PCT_EPT_RICH_NO_INTOL = taxa_pct_rich(rsc.df, unique_id, ORDER, !!taxa.rank,
                                           taxon = c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA"),
                                           exclusion.col = BIBI_TV, exclusion.vec = 1:3)
         )



test <- rsc.df %>% 
  select(unique_id) %>% 
  distinct() %>% 
  mutate(GOLD = taxa_pct(rsc.df, unique_id, count, CLASS, c("GASTROPODA", "OLIGOCHAETA")) +
           taxa_pct(rsc.df, unique_id, count, ORDER, "DIPTERA"),
         MARGALEFS = taxa_div(rsc.df, unique_id, count, NULL, !!taxa.rank, job = "margalef"),
         PCT_EPT_RICH = taxa_pct_rich(rsc.df, unique_id, ORDER, !!taxa.rank,
                                      taxon = c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA")),
         RICH_EPHEMEROPTERA = taxa_rich(rsc.df, unique_id, ORDER, !!taxa.rank, "EPHEMEROPTERA"),
         PCT_COLLECT = taxa_pct(rsc.df, unique_id, count, BIBI_FFG, c("CG", "CF")),
         PCT_PREDATOR = taxa_pct(rsc.df, unique_id, count, BIBI_FFG, "PR"),
         RICH_FILTER = taxa_rich(rsc.df, unique_id, BIBI_FFG, !!taxa.rank, "CF"),
         RICH_CLIMB = taxa_rich(rsc.df, unique_id, BIBI_HABIT, !!taxa.rank, "CB"),
         #ASPT_MOD = ,
         #HBI = ,
         PCT_CAECIDOTEA = taxa_pct(rsc.df, unique_id, count, GENUS, "CAECIDOTEA"),
         PCT_EPT_NO_HYDRO = taxa_pct(rsc.df, unique_id, count, ORDER, c("EPHEMEROPTER", "PLECOPTERA", "TRICHOPTERA")) - 
           taxa_pct(rsc.df, unique_id, count, FAMILY, "HYDROPSYCHIDAE"),
         PCT_DIPTERA = taxa_pct(rsc.df, unique_id, count, ORDER, "DIPTERA"),
         HURLBERTS_PIE = taxa_div(rsc.df, unique_id, count, NULL, !!taxa.rank, job = "gini_simpson"),
         PCT_GATHER = taxa_pct(rsc.df, unique_id, count, BIBI_FFG, "CG"),
         PCT_BURROW = taxa_pct(rsc.df, unique_id, count, BIBI_HABIT, "BU"),
         RICH_BURROW = taxa_rich(rsc.df, unique_id, BIBI_HABIT, !!taxa.rank, "BU"),
         PCT_DOM1 = mmir::pct_dom(rsc.df, unique_id, count, !!taxa.rank, 1),
         PCT_TOLERANT_5_10 = taxa_pct(rsc.df, unique_id, count, BIBI_TV, 5:10),
         PCT_TOLERANT_7_10 = taxa_pct(rsc.df, unique_id, count, BIBI_TV, 7:10),
         RICH_TOL = taxa_rich(rsc.df, unique_id, BIBI_TV, !!taxa.rank, 7:10),
         PCT_COTE = taxa_pct(rsc.df, unique_id, count, ORDER, c("COLEOPTERA", "ODONATA", 
                                                                "TRICHOPTERA", "EPHEMEROPTER")),
         PCT_EPT_RICH_NO_TOL = taxa_pct_rich(rsc.df, unique_id, ORDER, !!taxa.rank,
                                             taxon = c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA"),
                                             exclusion.col = BIBI_TV, exclusion.vec = 7:10),
         RICH_TRICHOPTERA = taxa_rich(rsc.df, unique_id, ORDER, !!taxa.rank, "TRICHOPTERA"),
         RICH_COLLECT = taxa_rich(rsc.df, unique_id, BIBI_FFG, !!taxa.rank, c("CG", "CF")),
         PCT_CLING = taxa_pct(rsc.df, unique_id, count, BIBI_HABIT, "CL"),
         RICH_CLING = taxa_rich(rsc.df, unique_id, BIBI_FFG, !!taxa.rank, "CL"),
         PCT_OLIGO_CHIRO = taxa_pct(rsc.df, unique_id, count, CLASS, "OLIGOCHAETA") +
           taxa_pct(rsc.df, unique_id, count, FAMILY, "CHIRONOMIDAE"),
         RICH_PREDATOR = taxa_rich(rsc.df, unique_id, BIBI_FFG, !!taxa.rank, "PR"),
         PCT_MOD_TOL_4_6 = taxa_pct(rsc.df, unique_id, count, BIBI_TV, 4:6)
         )

#------------------------------------------------------------------------------

rsc.sub <- rsc.df %>% 
  filter(unique_id == "3201702_200")
test2 <- rsc.sub %>% 
  select(unique_id) %>% 
  distinct() %>% 
  arrange(unique_id) %>% 
  mutate(GOLD = taxa_pct(rsc.sub, unique_id, count, CLASS, c("GASTROPODA", "OLIGOCHAETA")) +
           taxa_pct(rsc.sub, unique_id, count, ORDER, "DIPTERA"),
         MARGALEFS = taxa_div(rsc.sub, unique_id, count, NULL, !!taxa.rank, job = "margalef"),
         PCT_EPT_RICH = taxa_pct_rich(rsc.sub, unique_id, ORDER, !!taxa.rank,
                                      taxon = c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA")),
         RICH_EPHEMEROPTERA = taxa_rich(rsc.sub, unique_id, ORDER, !!taxa.rank, "EPHEMEROPTERA"),
         PCT_COLLECT = taxa_pct(rsc.sub, unique_id, count, BIBI_FFG, c("CG", "CF")),
         PCT_PREDATOR = taxa_pct(rsc.sub, unique_id, count, BIBI_FFG, "PR"),
         RICH_FILTER = taxa_rich(rsc.sub, unique_id, BIBI_FFG, !!taxa.rank, "CF"),
         RICH_CLIMB = taxa_rich(rsc.sub, unique_id, BIBI_HABIT, !!taxa.rank, "CB"),
         #ASPT_MOD = ,
         #HBI = ,
         PCT_CAECIDOTEA = taxa_pct(rsc.sub, unique_id, count, GENUS, "CAECIDOTEA"),
         PCT_EPT_NO_HYDRO = taxa_pct(rsc.sub, unique_id, count, ORDER, c("EPHEMEROPTER", "PLECOPTERA", "TRICHOPTERA")) - 
           taxa_pct(rsc.sub, unique_id, count, FAMILY, "HYDROPSYCHIDAE"),
         PCT_DIPTERA = taxa_pct(rsc.sub, unique_id, count, ORDER, "DIPTERA"),
         HURLBERTS_PIE = taxa_div(rsc.sub, unique_id, count, NULL, !!taxa.rank, job = "gini_simpson"),
         PCT_GATHER = taxa_pct(rsc.sub, unique_id, count, BIBI_FFG, "CG"),
         PCT_BURROW = taxa_pct(rsc.sub, unique_id, count, BIBI_HABIT, "BU"),
         RICH_BURROW = taxa_rich(rsc.sub, unique_id, BIBI_HABIT, !!taxa.rank, "BU"),
         PCT_DOM1 = mmir::pct_dom(rsc.sub, unique_id, count, !!taxa.rank, 1),
         PCT_TOLERANT_5_10 = taxa_pct(rsc.sub, unique_id, count, BIBI_TV, 5:10),
         PCT_TOLERANT_7_10 = taxa_pct(rsc.sub, unique_id, count, BIBI_TV, 7:10),
         RICH_TOL = taxa_rich(rsc.sub, unique_id, BIBI_TV, !!taxa.rank, 7:10),
         PCT_COTE = taxa_pct(rsc.sub, unique_id, count, ORDER, c("COLEOPTERA", "ODONATA", 
                                                                "TRICHOPTERA", "EPHEMEROPTER")),
         #PCT_EPT_RICH_NO_TOL =
         RICH_TRICHOPTERA = taxa_rich(rsc.sub, unique_id, ORDER, !!taxa.rank, "TRICHOPTERA"),
         RICH_COLLECT = taxa_rich(rsc.sub, unique_id, BIBI_FFG, !!taxa.rank, c("CG", "CF")),
         PCT_CLING = taxa_pct(rsc.sub, unique_id, count, BIBI_HABIT, "CL"),
         RICH_CLING = taxa_rich(rsc.sub, unique_id, BIBI_FFG, !!taxa.rank, "CL"),
         PCT_OLIGO_CHIRO = taxa_pct(rsc.sub, unique_id, count, CLASS, "OLIGOCHAETA") +
           taxa_pct(rsc.sub, unique_id, count, FAMILY, "CHIRONOMIDAE"),
         RICH_PREDATOR = taxa_rich(rsc.sub, unique_id, BIBI_FFG, !!taxa.rank, "PR"),
         PCT_MOD_TOL_4_6 = taxa_pct(rsc.sub, unique_id, count, BIBI_TV, 4:6)
  )


identical(test[1,], test2)





# Import BIBI thresholds
thresh.df <- readxl::read_excel("H:/Projects/Chessie_BIBI/report/FINAL_May25_2017/2017_Data/Metric_Thresholds/metric_thresholds.xlsx",
                   sheet = "thresholds") %>% 
  filter(!taxonomic_resolution %in% "Order",
         spatial_resolution %in% c("COAST", "SEP"))
unique(thresh.df$metric)
