#' Seconds to HH:MM:SS
#' 
#' Converts seconds into HH:MM:SS time format.
#' 
#' @param start_time Start time from Sys.time()
#' @param end_time End time from Sys.time()
#' 
secondsToPeriodTime <- function(start_time, end_time = Sys.time()){
  
  # Calculate time elapsed in seconds
  te <- difftime(
      time1 = end_time,
      time2 = start_time,
      units = c("secs")
      )
  
  # Convert seconds to period time
  te_dt <- lubridate::seconds_to_period(
    as.integer(round(te,0))
  )
  
  # Create formatted time
  sprintf(
    fmt = '%02d:%02d:%02d', 
    te_dt@hour, 
    lubridate::minute(te_dt), 
    lubridate::second(te_dt)
  )
  
}
