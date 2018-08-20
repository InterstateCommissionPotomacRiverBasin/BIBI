
check_range <- function(Wide.df, Param, Low, High){
  if(any(Wide.df[, Param] > High | Wide.df[, Param] < Low)) {
    cat("Check the following EVENT_ID(s) for errors associated with", Param,":",
        unique(Wide.df[(Wide.df[, Param] > High | Wide.df[, Param] < Low), c("EVENT_ID")]))
    param.name <- as.name(Param)
    param.name <- unique(Wide.df[(Wide.df[, Param] > High | Wide.df[, Param] < Low), c("EVENT_ID")])
  }else{
    cat(Param, "Checked")
  }
}

check_missing <- function(Wide.df, Param){
  test <- if(any(is.na(Wide.df[, Param]))){
    cat("No", Param, "value for the following EVENT_ID(s):",
        unique(Wide.df[is.na(Wide.df[, Param]) , c("EVENT_ID")]))
  }else{
    cat("No missing values associated with", Param)
  }
}
#Habitat=======================================================================
check_habitat <- function(Wide.df){
  #Check for missing values
  check_missing(Data, "BANKS")
  check_missing(Data, "CH_ALT")
  check_missing(Data, "HAB_HETERO")
  check_missing(Data, "INSTR_COND")
  check_missing(Data, "RIP_ZONE")
  check_missing(Data, "EMBED")
  check_missing(Data, "SUM")

  #Check for values our of range
  check_range(Data, "BANKS", 0, 20)
  check_range(Data, "CH_ALT", 0, 20)
  check_range(Data, "HAB_HETERO", 0, 20)
  check_range(Data, "INSTR_COND", 0, 20)
  check_range(Data, "RIP_ZONE", 0, 20)
  check_range(Data, "EMBED", 0, 20)
  check_range(Data, "SUM", 0, 120)
}

#Water Quality=======================================================================
check_wq <- function(){
  #Check for missing values
  check_missing(Data, "SPCOND")
  check_missing(Data, "PH")
  check_missing(Data, "DO")

  #Check for values out of "normal" range
  check_range(Data, "SPCOND", 30, 2000)
  check_range(Data, "PH", 2, 12)
  check_range(Data, "DO", 0, 20)
}

#Taxa Values=======================================================================
check_taxa <- function(Wide.df){
  check_missing(Data, "REPORTING_VALUE")
  check_range(Data, "REPORTING_VALUE", 0, 200)
}


