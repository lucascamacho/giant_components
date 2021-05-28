
#------------------------------------------------------------------------#
#---------------------------- clean_simdata() ----------------------------
#------------------------------------------------------------------------#

# Author: Kate
# Date: 04/2021
# Function that summarises the results from simulations of the coevolutionary
# model.
# input: simdata (results of simulation), structure (structure data)
# output: organised simdata 

#------------------------------------------------------------------------

clean_simdata <- function(simdata, structure){
  
  simdata <- data.table(simdata)
  comp_df <- simdata[, setdiff(names(simdata), "sp_id"), with = FALSE]
  comp_df <- comp_df[, lapply(.SD, mean, na.rm = TRUE), by = .(network_id, m)]
  names(comp_df)[1] <- "id"
  comp_df$id <- gsub(".txt.*", "", comp_df$id)
  all_df <- merge(x = comp_df, y = structure, by = "id", all.x = TRUE)
  all_df$m <- as.factor(all_df$m)
  
  return(all_df)
  
}

