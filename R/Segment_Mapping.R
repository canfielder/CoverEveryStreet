#' Segment Mapping
#' 
#' The following functions are related to generating the segment map dataset.
#' The segment pad dataset lists the connected segments for every individual
#' segment.
#' 
#' Determine the connected segments to each segment
#' 
#' The following function determines which segments are connected to the
#' provided segment.
#' 
#' @param network Street network (dataframe)
#' @param input_segment_id Segment identifier (character)
#' 
#' @return Segments connected to the input segment (list)
#'
connectedSegments <- function(network, input_segment_id) {
  
  # Library loaded for pipe use in parallel processes
  library(dplyr)
  
  # Determine node ids within the segment
  ## To Node
  segment_to_ids <- network %>% 
    dplyr::filter(segment_id == input_segment_id) %>% 
    dplyr::pull(to_id)
  
  ## From Node
  segment_from_ids <- network %>% 
    dplyr::filter(segment_id == input_segment_id) %>% 
    dplyr::pull(from_id)
  
  ## Compile
  segment_node_ids <- unique(c(segment_to_ids, segment_from_ids))
  
  # Reduce to segment ids which also use those node ids
  network %>% 
    dplyr::filter(
      (to_id %in% segment_node_ids) |
        (from_id %in% segment_node_ids) 
    ) %>% 
    dplyr::filter(
      segment_id != input_segment_id
    ) %>% 
    dplyr::pull(segment_id)
}


#' Determine the connected segments to each segment across a complete
#'  dataset
#' 
#' The following function executes the connectedSegments across a 
#' complete street network
#' 
#' @param network Street network (dataframe)
#' 
#' @return Street network with a new column of listing the segments connected
#'  each rows segment (dataframe)
#'
connectedSegmentsApply <- function(network){
  
  # Library loaded for pipe use in parallel processes
  library(dplyr)
  
  network <- network %>% 
    dplyr::mutate(
      connected_segment_ids = purrr::pmap(
        .l = list(
          input_segment_id = segment_id,
          network          = list(network)
        ),
        .f = connectedSegments
      )
    ) %>% 
    dplyr::select(segment_id, connected_segment_ids) %>% 
    dplyr::distinct() %>% 
    tibble::deframe()
}


#' Determine the connected segments to each segment across a complete
#'  dataset through parallel processing
#' 
#' The following function executes the **connectedSegmentsApply** utilizing
#' parallel processing to increase performance speed.
#' 
#' @param network Street network (dataframe)
#'  each rows segment (dataframe)
#' @param k Number of cores to use for parallel process (integer)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#'  
#' @return Street network with a new column of listing the segments connected


connectedSegmentsApplyPar <- function(
  network, k, start_time = Sys.time(), verbose = 0
){
  library(dplyr)
  # DEFINE CLUSTERS ###########################################################
  # Split network into sub-networks based on clusters
  network.list <- networkClusters(network, k)
  
  # PARALLEL PROCESSNG ########################################################
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("Parallel Processing - Start:\t{elapsed_time}"))
  }
  
  # Define dataset variables and required sub-functions for export to cluster
  cluster_data <- c("network.list")
  cluser_functions <- c("connectedSegments", "connectedSegmentsApply")
  cluster_variables <- c(cluster_data, cluser_functions)
  
  # Create Cluster
  cl <- parallel::makeCluster(k)
  
  # Export required variables to cluster
  parallel::clusterExport(
    cl      = cl, 
    varlist = cluster_variables, 
    envir   = environment()
  )
  
  # Execute parallel implementation of segment labeling
  network.mtx <- parallel::parSapply(
    cl  = cl, 
    X   = network.list, 
    FUN = connectedSegmentsApply
  )
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("Parallel Processing - Complete:\t{elapsed_time}"))
  }
  
  # PROCESS OUTPUT ############################################################
  # Compile parallel subset lists into a single list --------------------------
  
  # Initialize output list
  segment_map <- list()
  
  if (verbose >= 1){
    print(glue::glue("Compile Segment Maps: Start"))
  }
  
  for (i in seq(1,k,1)){
    # Define subset list
    tmp.list <- network.mtx[[i]]
    
    # Extract names from each list
    segment_ids <- names(segment_map)
    tmp.list.names <- names(tmp.list)
    
    # Reduce iteration list to only segment ids not already captured in the
    # segment map
    new_segments.bool <- !c(tmp.list.names %in% segment_ids)
    tmp.list <- tmp.list[new_segments.bool]
    
    # Add iteration list to final output list
    segment_map <- c(segment_map, tmp.list)
    
    if (verbose >= 2){
      print(glue::glue(strrep("_", 80)))
      print(glue::glue("Iteration:\t\t\t\t{i}"))
      print(glue::glue("Entries - iteration:\t\t{length(tmp.list.names)}"))
      print(glue::glue("New segments:\t\t\t{sum(new_segments.bool)}"))
      print(glue::glue("Entries - total:\t\t\t{length(segment_map)}"))
    }
  }
  
  if (verbose >= 1){
    print(glue::glue("Compile Segment Maps: Complete"))
  }
  
  
  segment_map
  
}
