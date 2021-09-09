#' Create LINESTRING Shapefile Object
#' 
#' Converts two points into a LINESTRING shapefile object.
#' 
#' @param from_lon Longitude of From point. (numeric)
#' @param from_lat Latitude of From point. (numeric)
#' @param to_lon Longitude of To point. (numeric)
#' @param to_lat Latitude of To point. (numeric)
#'
#' @return LINESTRING shapefile object
#' 
makeLinestring <- function (from_lon, from_lat, to_lon, to_lat) {
  
  # Define to/from point vector
  to_from_vector <- c(from_lon, to_lon, from_lat, to_lat)
  
  # Convert point vector into LINESTRING object
  sf::st_linestring(
    matrix(
      data = to_from_vector, 
      nrow = 2, 
      ncol =  2
    )
  )
}


#' Create LINESTRING Object - Dataframe Process
#' 
#' Creates a LINESTRING shapefile object within a dataframe. This 
#' requires the To and From latitude/longitude points to be columns
#' within the dataframe.
#' 
#' @param network Dataframe with from_lon, from_lat, to_lon, and to_lat columns
#'   (dataframe)
#' @param orthogonal Boolean flag indicating whether to create LINESTRING of 
#' the edge or a LINSETRING of the orthogonal projection of the edge
#'
#' @return Shapefile Dataframe with LINESTRING object as new column
#'
makeLinestringVector <- function(network, orthogonal = TRUE) { 
  # Import library for parallel processing. Use of pipe throws and error if 
  # the library is not defined or the package is not included with the pipe, 
  # which is very awkward
  library(dplyr)
  
  # Select correct geographic points, depending on if we are constructing
  # a normal or orthogonal projection
  
  if(!orthogonal){
    network_geo <- network %>% 
      select(from_lon, from_lat, to_lon, to_lat)
  } else{
    network_geo <- network %>% 
      select(from_lon_orth, from_lat_orth, to_lon_orth, to_lat_orth) %>% 
      rename_with(~gsub("_orth", "", .x))
  }
  # The following commented out code was the old method. The new method should
  # be equivalent, but in case it isn't, this code is maintained so I can 
  # switch back.
  
  # network_geo %>% 
  #   purrr::pmap(makeLinestring) %>% 
  #   sf::st_as_sfc(crs = 4326) %>% 
  #   {
  #     tibble(
  #       edge_id = network$edge_id, 
  #       geometry = .)
  #   }%>% 
  #   dplyr::left_join(
  #     y = network,
  #     by = "edge_id"
  #   ) %>% 
  #   sf::st_sf()
  
  # Create LINESTRING
  ls <- network_geo %>% 
    purrr::pmap(makeLinestring) %>% 
    sf::st_as_sfc(crs = 4326) 
  
  # Combine LINESTRING with edge_id to create dataframe
  ls_df <- tibble(
    edge_id = network$edge_id, 
    geometry = ls
  )
  
  # Append LINESTRING to network dataframe
  network %>%
    dplyr::left_join(
      y = ls_df,
      by = "edge_id"
    ) %>% 
    dplyr::select(geometry, dplyr::everything())

  }


#' Execute makeLinestringVector function in parallel
#' 
#' Executes the function makeLinestringVector using parallel processing
#' 
#' @param network Street network (dataframe)
#' @param k Number of cores to use for parallel process (integer)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#'  
#' @return Street network with each edge converted into a LINESTRING 
#'  (simple features dataframe)
#' 
makeLinestringVectorPar <- function(
  network, 
  orthogonal = FALSE, 
  k          = 6, 
  start_time = Sys.time(), 
  verbose    = 0
) {
  
  # CLUSTERS ##################################################################
  # Split network into sub-networks based on clusters
  network.list <- networkClusters(network, k)
  
  # ROUND 1 - PARALLEL PROCESSING #############################################
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("Parallel Processing - Start:\t\t\t{elapsed_time}"))
  }
  
  # Define dataset variables and required sub-functions for export to cluster
  cluster_data <- c("network", "verbose")
  cluser_functions <- c("makeLinestring", "makeLinestringVector")
  cluster_variables <- c(cluster_data, cluser_functions)
  
  # Create cluster
  cl <- parallel::makeCluster(k)
  
  # Export required variables to cluster
  parallel::clusterExport(
    cl      = cl, 
    varlist = cluster_variables, 
    envir   = environment()
  )
  
  # Execute parallel implementation of segment labeling
  network.mtx <- parallel::parSapply(
    cl         = cl, 
    X          = network.list, 
    FUN        = makeLinestringVector, 
    orthogonal = orthogonal
  )
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("Parallel Processing - Complete:\t\t\t{elapsed_time}"))
  }
  
  # PROCESSING ################################################################
  # Convert each matrix column to dataframe and store in list
  network_compile.list <- vector(mode = "list", length = k)
  for (i in seq(1,k,1)){
    network_compile.list[[i]] <- dplyr::as_tibble(network.mtx[,i])
  }
  
  # Bind all subset network dataframes in the list
  network <- dplyr::bind_rows(network_compile.list)
  
  # Remove cluster label
  network <- network %>% 
    dplyr::select(-.cluster) %>% 
    dplyr::select(geometry, dplyr::everything())
  
  network
}


#' Calculate the length of a LINESTRING in a dataframe
#' 
#' @param network Street network with LINESTRING column (simple feature 
#'  dataframe)
#'  
#' @return Street network with a column for the LINESTRING length (simple 
#' feature dataframe)
#
lineStringLength <- function(network){
  # Import library for parallel processing. Use of pipe throws and error if 
  # the library is not defined or the package is not included with the pipe, 
  # which is very awkward
  library(dplyr)
  
  network %>% 
    dplyr::mutate(
      edge_length = as.vector(sf::st_length(geometry))
    )
}


#' Calculate the length of a LINESTRING in a dataframe - parallel
#' 
#' @param network Street network with LINESTRING column (simple feature 
#'  dataframe)
#' @param k Number of cores to use for parallel process (integer)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer)
#'  
#' @return Street network with a column for the LINESTRING length (simple 
#' feature dataframe)
#
lineStringLengthPar <- function(
  network, 
  k = 6, 
  start_time = Sys.time(),
  verbose = 0){
  
  # CLUSTERS ##################################################################
  # Split network into sub-networks based on clusters
  network.list <- networkClusters(network, k)
  
  # ROUND 1 - PARALLEL PROCESSING #############################################
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("Parallel Processing - Start:\t\t\t{elapsed_time}"))
  }
  
  # Define dataset variables and required sub-functions for export to cluster
  cluster_data <- c("network", "verbose")
  cluser_functions <- c("lineStringLength")
  cluster_variables <- c(cluster_data, cluser_functions)
  
  # Create cluster
  cl <- parallel::makeCluster(k)
  
  # Export required variables to cluster
  parallel::clusterExport(
    cl      = cl, 
    varlist = cluster_variables, 
    envir   = environment()
  )
  
  # Execute parallel implementation of segment labeling
  network.mtx <- parallel::parSapply(
    cl         = cl, 
    X          = network.list, 
    FUN        = lineStringLength
  )
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("Parallel Processing - Complete:\t\t\t{elapsed_time}"))
  }
  
  # PROCESSING ################################################################
  # Convert each matrix column to dataframe and store in list
  network_compile.list <- vector(mode = "list", length = k)
  for (i in seq(1,k,1)){
    network_compile.list[[i]] <- dplyr::as_tibble(network.mtx[,i])
  }
  
  # Bind all subset network dataframes in the list
  network <- dplyr::bind_rows(network_compile.list)
  
  # Remove cluster label
  network <- network %>% 
    dplyr::select(-.cluster)
  
  network
}