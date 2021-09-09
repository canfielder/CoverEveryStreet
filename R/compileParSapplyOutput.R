#' Compile parallel processing output matrix
#' 
#' The following function compiles the output from **parSapply** function from 
#' the **Parallel** package(1xk matrix of simple features dataframes) into a
#'single dataframe.
#' 
#' @param network.mtx Output from parSapply, 1xk matrix comprised of simple 
#'  features dataframes for each cluster subset of the street network 
#'  (matrix)
#' @param k Number of cores to use for parSapply process (integer)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#'  
#' @return Street network with a new column of listing the segments connected
#

compileParSapplyOutput <- function(
  network.mtx, 
  k, 
  start_time = NA, 
  verbose = 0
  ) {
  
  # Establish starting time
  if (is.na(start_time)){
    start_time <- Sys.time()
  }

  segment_label <- "Compile Cluster Outputs"
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label} - Start:\t{elapsed_time}"))
  }
  
  # Convert each matrix column to dataframe and store in list
  network_compile.list <- vector(mode = "list", length = k)
  
  for (i in seq(1,k,1)){
    network_compile.list[[i]] <- dplyr::as_tibble(network.mtx[,i])
  }
  
  # Bind all subset network dataframes in the list and remove cluster label
  network <- network_compile.list %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(-.cluster)
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("{segment_label} - Complete:\t{elapsed_time}"))
  }
  
  # Return 
  network
}
