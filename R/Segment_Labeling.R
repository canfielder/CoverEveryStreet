#' Label Segments
#' 
#' The following functions are related to defining segments in a street 
#' network. Segments are considered a collection of edges (straight
#' line segments) which meet the following criteria:
#' 
#'  - Are between two street intersections (nodes with three or more
#'  connected edges)
#'  - Are between a street intersection and a dead end (node with 
#'  only used once as either a To or From node)
#'  - Are in a closed loop (like a track or closed walking path) 
#'  
#'  
#' Identify Segments Between Nodes
#' 
#' Assigns a unique identifier to segments in a street network.
#' 
#' @param network Street network (dataframe)
#' @param file_date Date for intersection and parallel file names. Must be in
#'  the form YYYY-MM-DD (character)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#'
#' @return Street network with unique identifiers assigned to edges to define
#'  segments (dataframe)
#' 
labelSegments <- function(
  network, file_date, start_time = Sys.time(), verbose = 0
  ){
  # Import library for parallel processing. Use of pipe throws and error if 
  # the library is not defined or the package is not included with the pipe, 
  # which is very awkward
  library(dplyr)
  start_time <- Sys.time()
  
  # SETUP######################################################################
  # Load compiled vector of nodes that will stop the walking function
  halt_walk_nodes_init <- haltSegmentWalkNodes(
    network   = network, 
    file_date = file_date, 
    verbose = verbose
    )
  
  # Determine previously used segment ids and the number of cases
  # where segment ID is NA
  # Determine status of segment ids used, how many edge ids have been 
  # assigned a segment id, and how many edges still require a segment id.
  total_obs <- length(network$edge_id)
  utilized_segment_ids <- unique(purrr::discard(network$segment_id, is.na))
  n_segment_ids <- length(utilized_segment_ids)
  n_edges_labeled <- sum(!is.na(network$segment_id))
  n_edges_unlabeled <- total_obs - n_edges_labeled
  
  if (verbose >= 2){
    m0 <- glue::glue(strrep("_", 80))
    m1 <- glue::glue("Status Counts:")
    m2 <- glue::glue("Total Observations:\t\t{total_obs}")
    m3 <- glue::glue("Unique Segment IDs:\t\t{n_segment_ids}")
    m4 <- glue::glue("Edges w/ a Segment ID:\t\t{n_edges_labeled}")
    m5 <- glue::glue("Edges w/o a Segment ID:\t\t{n_edges_unlabeled}")
    print(glue::glue("{m0}\n{m1}\n{m2}\n{m3}\n{m4}\n{m5}\n{m0}"))
  }
  
  # Iteration Loop ############################################################
  # Initialize loop counter
  n <- 1
  while (n_edges_unlabeled > 0){
    if (verbose >= 2){
      hypen_print <- strrep("_", 35)
      print(glue::glue("{hypen_print} NEW LOOP {hypen_print}"))
    }
    
    # Set Up Walk -------------------------------------------------------------
    # Select starting point edge
    starting_edge <- network %>% 
      filter(is.na(segment_id)) %>% 
      dplyr::slice_sample(n=1) %>% 
      dplyr::pull(edge_id)
    
    if (verbose >= 2){
      print(glue::glue("Edge ID - Start:\t\t{starting_edge}"))
    }
    # Initialize vector to capture all edges within a single segment
    edges_in_segment <- vector(mode = "character", length = 0)
    
    # Start Walk - To Direction -----------------------------------------------
    edges_in_segment<- walkSegment(
      network = network,
      edge = starting_edge,
      edges_in_segment = edges_in_segment,
      halt_walk_nodes = halt_walk_nodes_init,
      direction = "to",
      verbose = verbose
    )
    
    # Start Walk - From Direction ---------------------------------------------
    edges_in_segment<- walkSegment(
      network = network,
      edge = starting_edge,
      edges_in_segment = edges_in_segment,
      halt_walk_nodes = halt_walk_nodes_init,
      direction = "from",
      verbose = verbose
    )
    
    # Update Network Table ---------------------------------------------------
    ## Generate Random Segment ID
    segment_id_rdm <- randomID(unique(utilized_segment_ids))
    
    if (verbose >= 2){
      print(glue::glue("Segment ID:\t\t\t\t{segment_id_rdm}"))
    }
    
    ## Reduce to unique values in case of double counting
    edges_in_segment <- unique(edges_in_segment)
    
    ## Label edges in network table with segment id
    network <- network %>% 
      mutate(
        segment_id = dplyr::if_else(
          condition = (edge_id %in% edges_in_segment),
          true =  segment_id_rdm, 
          false = segment_id
        )
      )
    
    # Recalculate while loop parameters ---------------------------------------
    # Determine previously used segment ids and the number of cases
    # where segment ID is NA
    # Determine status of segment ids used, how many edge ids have been 
    # assigned a segment id, and how many edges still require a segment id.
    utilized_segment_ids <- unique(purrr::discard(network$segment_id, is.na))
    n_segment_ids <- length(utilized_segment_ids)
    n_edges_labeled <- sum(!is.na(network$segment_id))
    n_edges_unlabeled <- total_obs - n_edges_labeled
    
    if (verbose >= 2){
      m0 <- glue::glue(strrep("_", 80))
      m1 <- glue::glue("Status Counts:")
      m2 <- glue::glue("Total Observations:\t\t\t\t{total_obs}")
      m3 <- glue::glue("Unique Segment IDs:\t\t\t\t{n_segment_ids}")
      m4 <- glue::glue("Edges w/ a Segment ID:\t\t\t{n_edges_labeled}")
      m5 <- glue::glue("Edges w/o a Segment ID:\t\t\t{n_edges_unlabeled}")
      print(glue::glue("{m0}\n{m1}\n{m2}\n{m3}\n{m4}\n{m5}\n{m0}"))
    }
    
    # Print Updates -----------------------------------------------------------
    if ((n %% 1e2 == 0) & verbose >= 1) {
      end_time <- Sys.time()
      elapsed_time <- secondsToPeriodTime(start_time, end_time)
      print(glue::glue("Elapsed Time:\t\t\t\t{elapsed_time}"))
      
      percent_rem <- round(100 * (total_obs - n_edges_unlabeled)/total_obs,1)
      print(glue::glue("Percent Edges Remaining:\\tt{percent_rem}%"))
    }
    
    # Safety Valve -----------------------------------------------------------
    if ((n > 1e5) & verbose >= 1){
      stop(glue::glue("Maximum iterations reached:\t\t{n}"))
    }
    
    # Increase iteration counter
    n <- n + 1
    
  }
  
  if (verbose >= 1){
    m0 <- glue::glue(strrep("_", 80))
    m1 <- glue::glue("Status Counts:")
    m2 <- glue::glue("Total Observations:\t\t{total_obs}")
    m3 <- glue::glue("Unique Segment IDs:\t\t{n_segment_ids}")
    m4 <- glue::glue("Edges w/ a Segment ID:\t\t{n_edges_labeled}")
    m5 <- glue::glue("Edges w/o a Segment ID:\t\t{n_edges_unlabeled}")
    print(glue::glue("{m0}\n{m1}\n{m2}\n{m3}\n{m4}\n{m5}\n{m0}"))
  }
  
  network
}


#' Walk along a segment
#' 
#' Walks along a segment, from edge to edge, determining if each subsequent
#' edge remains part of the segment
#' 
#' @param network Street network (dataframe)
#' @param edge Unique identifier of an edge in the street network (character)
#' @param edges_in_segment Edges which have been identified as part of the 
#'  segment which is being walked (character vector)
#' @param halt_walk_nodes Nodes which are at the end of a segment 
#'  (character vector)
#' @param direction Either "to" or "from". Indicates which direction along the
#'  segment to walk (character )
#' @param verbose Flag indicating the level of print outputs the function 
#'  produces. The larger the number the more granular the outputs (numeric)
#'  
#' @return Vector of edge identifiers for all the edges associated with 
#'  the walked segment
#' 
walkSegment <- function(
  network,
  edge,
  edges_in_segment,
  halt_walk_nodes,
  direction,
  verbose = 0
) {
  # Add edge to vector of edges within segment
  edges_in_segment <- c(edges_in_segment, edge)
  
  # Select Node to walk to 
  if (direction == "to"){
    node <- network %>% 
      filter(edge_id == edge) %>% 
      pull(to_id)
  } else {
    node <- network %>% 
      dplyr::filter(edge_id == edge) %>% 
      dplyr::pull(from_id)
  }
  
  if (verbose >= 3) {
    print(glue::glue("Direction:\t\t\t\t{direction}"))
    print(glue::glue("Selected Node:\t\t\t\t{node}"))
  }
  
  # Check if walk should be halted
  HALT_WALK <- node %in% halt_walk_nodes
  
  if (verbose >= 3) {
    print(glue::glue("Halt Walk:\t\t\t\t{HALT_WALK}"))
  }
  
  # Recursive choice - act on whether walk should stop
  if (HALT_WALK){
    
    edges_in_segment
    
  } else {
    
    # Add node to halt walk nodes. This is to prevent getting trapped in a 
    # loop, like a track or closed-loop walking path.
    halt_walk_nodes <- c(halt_walk_nodes, node)
    
    
    # Select next edge
    if (direction == "to"){
      
      next_edge <- network %>% 
        dplyr::filter(from_id == node) %>% 
        dplyr::pull(edge_id)
      
    } else {
      
      next_edge <- network %>% 
        dplyr::filter(to_id == node) %>% 
        dplyr::pull(edge_id)
      
    }
    
    if (verbose >= 3) {
      print(glue::glue("Edge ID - Next:\t\t\t{next_edge}"))
    }
    
    # Verify that next_edge is a single length vector. If it is not, there is 
    # an error in the properly defining the halt nodes
    if (length(next_edge) > 1 ) {
      msg1 <- glue::glue("Next edge selected does not have length of one.")
      msg2 <- glue::glue("Next Edge Length:\t\t{length(next_edge)}")
      msg3 <- glue::glue("Current Edge:\t\t{edge}")
      msg4 <- glue::glue("Problem Node:\t\t{node}")
      stop_msg <- glue::glue("{msg1}\n{msg2}\n{msg3}\n{msg4}")
      stop(stop_msg)
      
    }
    
    # Walk to next edge
    walkSegment(
      network          = network,
      edge             = next_edge,
      edges_in_segment = edges_in_segment,
      halt_walk_nodes  = halt_walk_nodes,
      direction        = direction,
      verbose          = verbose
    ) 
    
  }
}

#' Execute labelSegment function in parallel
#' 
#' Executes the function labelSegment using parallel processing
#' 
#' @param network Street network (dataframe)
#' @param file_date Date for intersection and parallel file names. Must be in
#'  the form YYYY-MM-DD (character)
#' @param k Number of cores to use for parallel process (integer)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#'  
#' @return Street network with unique identifiers assigned to edges to define
#'  segments (dataframe)
#' 
labelSegmentsPar <- function(
  network, 
  file_date, 
  k = 6, 
  start_time = Sys.time(),
  verbose = 0
  ) {
  # SETUP ####################################################################
  # import
  halt_walk_nodes <- haltSegmentWalkNodes(
    network   = network, 
    file_date = file_date, 
    verbose = verbose
  )
  
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
  cluster_data <- c("network.list", "start_time", "file_date", "verbose")
  cluser_functions <- c(
    "randomID", "walkSegment", "secondsToPeriodTime", "haltSegmentWalkNodes"
    )
  cluster_variables <- c(cluster_data, cluser_functions)
  
  # Create cluster
  cl <- makeCluster(k)
  
  # Export required variables to cluster
  clusterExport(
    cl      = cl, 
    varlist = cluster_variables, 
    envir   = environment()
  )
  
  # Execute parallel implementation of segment labeling
  network.mtx <- parallel::parSapply(
    cl  = cl, 
    X   = network.list, 
    FUN = labelSegments, 
    file_date = file_date,
    verbose = verbose, 
    start_time = start_time
  )
  
  # Stop cluster
  stopCluster(cl)
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("Parallel Processing - Complete:\t\t\t{elapsed_time}"))
  }
  
  # PROCESSING ################################################################
  # Convert each matrix column to dataframe and store in list
  network_compile.list <- vector(mode = "list", length = k)
  for (i in seq(1,k,1)){
    network_compile.list[[i]] <- as.data.frame(network.mtx[,i])
  }
  
  # Bind all subset network dataframes in the list
  network <- bind_rows(network_compile.list)
  
  # Remove cluster label
  network <- network %>% select(-.cluster)
  
  # Define all node/segment pairs
  ## To Nodes
  to_node_w_seg <- network %>% 
    select(to_id, segment_id) %>% 
    rename(id = to_id)
  
  ## From nodes
  from_node_w_seg <- network %>% 
    select(from_id, segment_id)  %>% 
    rename(id = from_id)
  
  ## All Nodes
  node_seg_pairs <- bind_rows(to_node_w_seg, from_node_w_seg)
  
  # Identify segments that have been assigned two segment ids due to being 
  # split during parallel processing.
  multi_labeled_segments <- node_seg_pairs %>% 
    group_by(id) %>% 
    mutate(
      unique_segments = n_distinct(segment_id)
    ) %>% 
    ungroup() %>% 
    mutate(
      is_halt_node = dplyr::if_else(
        condition = id %in% halt_walk_nodes, 
        true = TRUE, false = FALSE
      ),
      is_split_segment =
        dplyr::if_else(
          condition = (
            (unique_segments == 2) &
              !is_halt_node
          ),
          true = TRUE, false = FALSE
        )
    ) %>% 
    dplyr::filter(is_split_segment) %>% 
    dplyr::pull(segment_id) %>% 
    unique()
  
  # Remove all segment ids associated with segments labeled multiple times
  network <- network %>% 
    mutate(
      segment_id = if_else(
        condition = (segment_id %in% multi_labeled_segments), 
        true = NA_character_, 
        false = segment_id
      )
    )
  
  if (verbose >= 1){
    # Record / print reduction stats
    n_total <- length(network$edge_id)
    n_segments <- sum(!is.na(network$segment_id))
    per_reduction <- round(100 * ((n_total - n_segments)/n_total), 1)
    print(glue::glue("Segment Labeling Error Rate\t\t\t{per_reduction}%"))
  }
  
  # ROUND 2 - STANDARD PROCESSING #############################################
  # Round 2 fixes errors introduced by the parallel method.
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("Standard Processing - Start:\t\t\t{elapsed_time}"))
  }
  
  network <- labelSegments(
    network = network, 
    file_date = file_date, 
    verbose = verbose
  )
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("Standard Processing - Complete:\t{elapsed_time}"))
  }
  
  # RETURN ####################################################################
  
  network
  
}