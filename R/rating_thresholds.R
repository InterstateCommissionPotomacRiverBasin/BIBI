# 
# library(readxl)
# file.dir <- "H:/Projects/Chessie_BIBI/report/FINAL_May25_2017/2017_Data/Scores_Ratings/BIBI_Scores_Ratings_06292017.xlsx"
# rating.list <- lapply(excel_sheets(file.dir)[-1], function(sheet.i) {
#   read_excel(sheet.i, path = file.dir) %>% 
#     select(SPATIAL, HALF_REF_10, REF_10, REF_25, REF_50) %>% 
#     distinct() %>% 
#     mutate(res = sheet.i)
# })
# 
# rating.thresh <- bind_rows(rating.list) %>% 
#   clean_df() %>% 
#   mutate(
#     spatial_resolution = case_when(
#     stringr::str_detect(res, "ches") ~ "chesapeake_wide",
#     stringr::str_detect(res, "bioregion") ~ "bioregion",
#     stringr::str_detect(res, "region") ~ "region",
#     TRUE ~ "ERROR"
#   ),
#   taxonomic_resolution = case_when(
#     stringr::str_detect(res, "order") ~ "order",
#     stringr::str_detect(res, "family") ~ "family",
#     stringr::str_detect(res, "genus") ~ "genus",
#     TRUE ~ "ERROR"
#   )
#   ) %>% 
#   select(-spatial_resolution, -res) %>% 
#   rename(spatial_resolution = spatial) %>% 
#   select(spatial_resolution, taxonomic_resolution, everything())
# 
# data.table::fwrite(rating.thresh, "rating_threshold_06292017.csv")
