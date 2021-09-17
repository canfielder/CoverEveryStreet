#' This script contains functions related to the identification of
#' streets represented as two parallel lines where the segment
#' division on both lines is not compatible.
#' 
#' Parallel Segment Identification Processing
#' 
#' The following are the dataframe processing steps required before identifying
#' parallel streets with incompatible segment division.
#' 
#' @param network Dataframe which represents a street network (dataframe)
#' @param orth_dist distance of orthogonal projection, in meters (numeric)
#' 
#' @return Network with additional properties calculated for parallel street
#'  processing (dataframe)

preProcessingParallelIdentification <- function(network, orth_dist = 15){
  # Calculate edge center point
  network <- network %>% 
    dplyr::mutate(
      center_lat = (to_lat + from_lat)/2,
      center_lon = (to_lon + from_lon)/2
    )
  
  # Calculate edge angle to normal x,y axis. Angle is reported in units of pi 
  # radians. Invert all negative angles.
  network <- network %>% 
    dplyr::mutate(
      theta = atan2(
        y = to_lat - from_lat,
        x = to_lon - from_lon
      ) / pi ,
      theta = if_else(
        theta < 0, 
        theta + 1,
        theta
      )
    )
  
  # Calculate angle of orthogonal projection to the edge. This is done by shifting
  # the edge angle clockwise by by 1/2.
  network <- network %>% 
    dplyr::mutate(
      theta_orth = theta + 1/2,
      theta_orth = if_else(
        theta_orth >= 2, 
        theta_orth - 2, 
        theta_orth
      )
    )
  
  # Determine points of orthogonal projection. Because our data is representing
  # a real world area, a 2D rendering is only a approximation. This affects the 
  # implementation of length, which we need to determine the start and end 
  # points of the orthogonal projection. Depending on the angle of the orthogonal 
  # projection, length in the 2D space could be off by as much as a factor of 
  # 2 PI. Therefore, in order for the length calculation to be conservative, we 
  # will use a 2 PI factor to adjust the 2D length.
  
  # Calculate distance measures
  dist_adj_factor <- 2 * pi
  dist_adj <- orth_dist * dist_adj_factor
  m_to_deg = 111000
  orth_deg <- dist_adj/m_to_deg
  
  network <- network %>% 
    dplyr::mutate(
      orth_x        = orth_deg * cos(theta_orth*pi),
      orth_y        = orth_deg * sin(theta_orth*pi),
      to_lat_orth   = (center_lat + orth_y),
      to_lon_orth   = (center_lon + orth_x),
      from_lat_orth = (center_lat - orth_y),
      from_lon_orth = (center_lon - orth_x),
    )
  
  network
}



#' Export Parallel Edge Results
#' 
#' The following function saves the parallel edge check results depending on 
#' if a parallel edge was detected or not.
#' 
#' @param network Dataframe which represents a street network (dataframe)
#' @param input_edge_id Network street edge to which is being compared to other
#'  edges to detect if they are parallel (char)
#' @param parallel_detected Flag indicating a parallel match was found (bool)
#' @param parallel_edges_detected Vector of parallel edges detected (char vector)
#' @param parallel_edges_detected Vector of parallel segments detected 
#'  (char vector)
#'
#' @return Network with parallel check results incorporated (dataframe)
#' 
exportParallelEdgeResults <- function(
  network, 
  input_edge_id, 
  parallel_detected = FALSE, 
  parallel_edges_detected = NA, 
  parallel_segments_detected = NA
) {
  if(parallel_detected) {
    
    # Parallel Segments Detected
    
    network %>% 
      mutate(
        parallel_segments = dplyr::if_else(
          condition = (edge_id == input_edge_id),
          true      = list(parallel_segments_detected),
          false     = parallel_segments
        ),
        n_parallel_segments = dplyr::if_else(
          condition = (edge_id == input_edge_id),
          true      = as.double(length(parallel_segments_detected)),
          false     = n_parallel_segments
        ),
        parallel_edges = dplyr::if_else(
          condition = (edge_id == input_edge_id),
          true      = list(parallel_edges_detected),
          false     = parallel_edges
        )
      )
    
  } else {
    
    # Parallel Segments Not Detected
    
    network %>% 
      mutate(
        n_parallel_segments = dplyr::if_else(
          condition = (edge_id == input_edge_id),
          true      = 0, 
          false     = n_parallel_segments
        )
      )
  }
}


#' Perform pre-processing steps for network dataframes
#' 
#' The following function performs all of the necessary pre-processing steps
#' on the network dataframes, generating normal and orthogonal 
#' simple features dataframes.
#' 
#' @param network Street network (dataframe)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#' @param network_complete Street network of the complete area. This is for 
#'  parallel processing. If not employing parallel processing, the **network**
#'  input covers the complete area (simple features dataframe) 
#'
#' @return Network with parallel segments identified (simple feature network)
#
networkLinestringProcessing <- function(
  network, 
  parallel = FALSE,
  start_time = NA,
  k = 6,
  verbose = 0
) {
  
  # Establish starting time
  if (is.na(start_time)){
    start_time <- Sys.time()
  }
  
  if (verbose >= 2) {
    if (parallel){
      print(glue::glue("Parallel Processing: ON"))
    } else{
      print(glue::glue("Parallel Processing: OFF"))
    }
  }
  
  # Orthogonal ----------------------------------------------------------------
  if (verbose >= 2) {
    print(glue::glue("Create LINESTRING - Orthogonal:\t\t\tStart"))
  }
  
  # Pre-processing
  network_orth_input <- preProcessingParallelIdentification(
    network   = network,
    orth_dist = 15
  )
  
  if (parallel){
    # Create LINESTRING - Orthogonal
    network_orth <- makeLinestringVectorPar(
      network    = network_orth_input,
      orthogonal = TRUE,
      k          = k,
      start_time = start_time,
      verbose    = verbose
    ) 
    } else {
      network_orth <- makeLinestringVector(
        network    = network_orth_input,
        orthogonal = TRUE
      ) 
    }
  
  if (verbose >= 2) {
    print(glue::glue("Create LINESTRING - Orthogonal:\t\t\tComplete"))
  }
  
  # Normal --------------------------------------------------------------------
  # Pre-processing
  network <- preProcessingParallelIdentification(
    network   = network,
    orth_dist = 15
  )
  
  if (verbose >= 2) {
    print(glue::glue("Create LINESTRING - Normal:\t\t\tStart"))
  }
  
  # Create LINESTRING - Normal
  if (parallel){
    network <- makeLinestringVectorPar(
      network,
      orthogonal = FALSE,
      k          = k,
      start_time = start_time,
      verbose    = verbose
    )
  } else {
    network <- makeLinestringVector(
      network    = network,
      orthogonal = FALSE
      )
    }
  
  if (verbose >= 2) {
    print(glue::glue("Create LINESTRING - Normal:\t\t\tComplete"))
  }
  
  if (verbose >= 2) {
    print(glue::glue("Calculate LINESTRING Length - Normal:\t\tStart"))
  }
  
  # Calculate LINESTRING length
  if (parallel){
    network <- lineStringLengthPar(
      network    = network,
      k          = k,
      start_time = start_time,
      verbose    = verbose
    )
  } else {
    network <- lineStringLength(
      network    = network
    )
    }
  
  if (verbose >= 2) {
    print(glue::glue("Calculate LINESTRING Length - Normal:\t\tComplete"))
  }

  
  # Export --------------------------------------------------------------------
  output.list <- vector(mode = "list", length = 2)
  names(output.list) <- c("normal", "orthogonal")
  
  # Assign dataframes to export list
  output.list[["normal"]] <- network
  output.list[["orthogonal"]] <- network_orth
  
  output.list
}


#' Perform pre-processing steps for network dataframes - Complete Network
#'
#' @param network Street network (dataframe)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#'
#' @return Network with parallel segments identified (simple feature network)
#
networkLinestringProcessingComplete <- function(
  network, 
  parallel = FALSE,
  start_time = NA,
  k = 6,
  verbose = 0
) {
  if (verbose >= 2) {
    print(glue::glue("Generate LINESTRINGS Features - Complete Network: Start"))
  }
  
  
  # Establish starting time
  if (is.na(start_time)){
    start_time <- Sys.time()
  }
    
    # Pre-processing
  network <- preProcessingParallelIdentification(
      network   = network,
      orth_dist = 15
    )
  
  # Determine if process should be executed in parallel or not.
  if (parallel) {
    if (verbose >= 2) {
      print(glue::glue("Parallel Processing: ON"))
    }
   
    if (verbose >= 2) {
      print(glue::glue("Create LINESTRING - Complete:\t\t\tStart"))
    }
      
      # Create LINESTRING - Normal
    network <- makeLinestringVectorPar(
      network    = network,
      orthogonal = FALSE,
      k          = k,
      start_time = start_time,
      verbose    = verbose
      )
      
    if (verbose >= 2) {
      print(glue::glue("Create LINESTRING - Complete:\t\t\tComplete"))
    }
    
    if (verbose >= 2) {
      print(glue::glue("Calculate LINESTRING Length - Complete:\t\tStart"))
    }
      
      # Calculate LINESTRING length
    network <- lineStringLengthPar(
      network    = network,
      k          = k,
      start_time = start_time,
      verbose    = verbose
    )
      
    if (verbose >= 2) {
      print(glue::glue("Calculate LINESTRING Length - Complete:\t\tComplete"))
    }
  
  } else {
    if (verbose >= 2) {
      print(glue::glue("Parallel Processing: OFF"))
    }
    
    if (verbose >= 2) {
      print(glue::glue("Create LINESTRING - Complete:\t\t\tStart"))
    }
    
    # Create LINESTRING - Normal
    network <- makeLinestringVector(
      network    = network,
      orthogonal = FALSE
    )
    
    if (verbose >= 2) {
      print(glue::glue("Create LINESTRING - Complete:\t\t\tComplete"))
    }
    
    if (verbose >= 2) {
      print(glue::glue("Calculate LINESTRING Length - Complete:\t\tStart"))
    }
    
    # Calculate LINESTRING length
    network <- lineStringLength(
      network    = network
    )
    
    if (verbose >= 2) {
      print(glue::glue("Calculate LINESTRING Length - Compare:\t\tComplete"))
    }
  }
  
  # Export --------------------------------------------------------------------
  if (verbose >= 2) {
    print(glue::glue("Generate LINESTRINGS Features - Complete Network: Complete"))
  }
  
  network

}


#' Compare edges to identify parallel cases
#' 
#' The following function compares each edge in a supplied network against all
#' other qualifying edges, in order to identify which edges are parallel to 
#' one another.
#' 
#' @param network_cluster Subset of the complete street network, reduced via a
#'  kmeans clustering method(simple features dataframe)
#' @param network_complete Street network of the complete area. (simple 
#'  features dataframe) 
#' @param dist Maximum distance an edge in the **network_complete** network
#'  can be from an edge in the **network_cluster** network to be kept, in 
#'  meters (numeric)
#'
#' @return Network reduced to only edges within a defined distance of the 
#'  cluster subset network (simple features dataframe)
#' 
cropComparisonNetwork <- function(
  network_cluster, 
  network_complete, 
  dist = 100,
  verbose = 0
  ) {
  if (verbose >= 1) {
    time_elapsed <- secondsToPeriodTime(start_time)
    print(glue::glue("Extract simple feature geometries:\t\t\t\t\t\t{time_elapsed}"))
  }
  # Extract simple feature geometries
  cluster_sf_geo <- network_cluster %>% dplyr::pull(geometry)
  network_sf_geo <- network_complete %>% dplyr::pull(geometry)
  
  # Determine which edges in the complete network are within x distance of the cluster
  # Allows for reduction in size of comparison network
  
  if (verbose >= 1) {
    time_elapsed <- secondsToPeriodTime(start_time)
    print(glue::glue("Identify Edges within defined distance of cluster network:\t\t{time_elapsed}"))
  }
  
  edge_to_edge_test <- sf::st_is_within_distance(
    x      = network_sf_geo, 
    y      = cluster_sf_geo, 
    sparse = FALSE,
    dist   = dist
  )
  
  # Determine number of network edges
  n_edges_network <- dim(edge_to_edge_test)[1]
  
  # If any edge in the network is within x distance of any edge in the 
  # cluster subset, select True.
  
  if (verbose >= 1) {
    time_elapsed <- secondsToPeriodTime(start_time)
    print(glue::glue("Reduce matrix to yes/no vector:\t\t\t{time_elapsed}"))
  }
  
  keep_edge <- apply(
    X      = edge_to_edge_test, 
    MARGIN = 1, 
    FUN    = any
  )
  
  # Assign output vector to comparison network dataframe
  network_complete$near_cluster <- keep_edge
  
  # reduce to only edges within x distance from any cluster edge
  network_complete %>% 
    dplyr::filter(near_cluster)
}


#' Compare edges to identify parallel cases
#' 
#' The following function compares each edge in a supplied network against all
#' other qualifying edges, in order to identify which edges are parallel to 
#' one another.
#' 
#' @param network Street network (simple features dataframe)
#' @param network_orth Network comprised of orthogonal projections of all edges
#'  within the provided **network** ( simple features dataframe)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#' @param network_complete Street network of the complete area. This is for 
#'  parallel processing. If not employing parallel processing, the **network**
#'  input covers the complete area (simple features dataframe) 
#'
#' @return Network with parallel segments identified (simple feature network)
#' 
compareParallelEdges <- function(
  network,
  network_orth,
  network_compare = NA, 
  k = 6,
  start_time = NA, 
  verbose = 0 
) {
  library(dplyr)
  
  # Establish starting time
  if (is.na(start_time)){
    start_time <- Sys.time()
  }
  
  if (verbose >= 1) {
    time_elapsed <- secondsToPeriodTime(start_time)
    print(glue::glue("Compare Parallel Segments - Start:\t\t\t\t\t{time_elapsed}"))
    n_sub <- length(network$edge_id)
    n_total <-length(network_compare$edge_id)
    print(glue::glue("Cluster Sub-network Size:\t\t\t\t\t\t{n_sub}"))
    print(glue::glue("Comparison Network Size:\t\t\t\t\t\t{n_total}"))
    print(glue::glue(""))
    
  }
  
  # PROCESSING ################################################################
  # Define complete network for comparison
  if (!(is.data.frame(network_compare))) {
    network_compare <- network
  } 
  
  # Initialize Parallel Separation columns
  network <- network %>%
    dplyr::mutate(
      parallel_edges      = list(c("<EMPTY>")),
      parallel_segments   = list(c("<EMPTY>")),
      n_parallel_segments = -1,
    )
  
  # LOOP ######################################################################
  # Initialize Inputs
  ## Subset radius
  radius_m <- 10
  safety_buffer <- 0.10
  radius_adj <- radius_m * 2 * pi * (1 + safety_buffer)
  
  ## Angle
  angle_buffer <- 0.015
  
  ## Initial while loop value check
  n_remaining <- network %>%
    dplyr::filter(n_parallel_segments == -1) %>%
    dplyr::pull(edge_id) %>%
    length()
  
  ## Iteration counter
  iter_counter <- 0
  
  while (n_remaining > 0) {
    iter_counter <- iter_counter + 1
    
    # Runaway While Loop escape valve
    if(iter_counter %% 1e4 == 0) {
      print(glue::glue('Maximum iterations reached.'))
      break
    }
    
    # Time Updates
    if (verbose >= 1) {
      if(iter_counter %% 5e1 == 0) {
        time_elapsed <- secondsToPeriodTime(start_time)
        print(glue::glue(" "))
        print(glue::glue(strrep("*", 80)))
        print(glue::glue('Iteration:\t\t\t\t\t\t\t\t{iter_counter}'))
        print(glue::glue('Time Elapsed:\t\t\t\t\t\t\t\t{time_elapsed}'))
        print(glue::glue(strrep("*", 80)))
        print(glue::glue(" "))
      }
    }
    
    # Define comparison network for this iteration
    network_compare_iter <-  network_compare
    
    # Setup ---------------------------------------------------------------------
    if (verbose >= 2){
      print(glue::glue(""))
      print(glue::glue(strrep("#", 80)))
      print(glue::glue('Iteration:\t\t\t\t\t\t\t\t{iter_counter}'))
    }
    
    # Select single random row that has not been processed
    iter.sf <- network %>%
      dplyr::filter(n_parallel_segments == -1) %>%
      dplyr::slice_sample(n = 1)
    
    # Extract Inputs
    edge_id_iter    = iter.sf %>% dplyr::pull(edge_id)
    center_lat      = iter.sf %>% dplyr::pull(center_lat)
    center_lon      = iter.sf %>% dplyr::pull(center_lon)
    segment_id_iter = iter.sf %>% dplyr::pull(segment_id)
    theta_iter      = iter.sf %>% dplyr::pull(theta)
    
    if (verbose >= 3){
      print(glue::glue('Edge ID:\t\t\t\t\t\t{edge_id_iter}'))
      print(glue::glue('Number of Edges - Total:\t\t\t\t{length(network$edge_id)}'))
    }
    
    # REDUCE DATASET BEFORE SHAPEFILE PROCESSES #################################
    # Subset - Radius ---------------------------------------------------------
    network_compare_iter <- subset_lat_lon(
      network    = network_compare_iter,
      lat_center = center_lat,
      lon_center = center_lon,
      radius     = radius_adj
    )
    
    n_subset <- length(network_compare_iter$edge_id)
    
    if (verbose >= 3){
      print(glue::glue('Number of Edges - Subset - Radius:\t\t\t{n_subset}'))
    }
    
    if (n_subset == 0) {
      network <- exportParallelEdgeResults(
        network                    = network,
        input_edge_id              = edge_id_iter,
        parallel_detected          = FALSE,
        parallel_edges_detected    = NA,
        parallel_segments_detected = NA
      )
      
      # Iteration  actions
      n_remaining <- network %>%
        dplyr::filter(n_parallel_segments == -1) %>%
        dplyr::pull(edge_id) %>%
        length()
      next
    }
    
    # Subset - Angle ----------------------------------------------------------
    # Define max/min angle values
    theta_max <- theta_iter + angle_buffer
    theta_min <- theta_iter - angle_buffer
    
    # Execute subset
    network_compare_iter <- network_compare_iter %>%
      filter(
        between(
          x     = theta,
          left  = theta_min,
          right = theta_max
        )
      )
    
    n_subset <- length(network_compare$edge_id)
    
    if (verbose >= 3){
      print(glue::glue('Number of Edges - Subset - Angle:\t\t\t{n_subset}'))
      if (n_subset != 0){
        print(glue::glue(strrep("<>", 40)))
      }
    }
    
    if (n_subset == 0) {
      network <- exportParallelEdgeResults(
        network                    = network,
        input_edge_id              = edge_id_iter,
        parallel_detected          = FALSE,
        parallel_edges_detected    = NA,
        parallel_segments_detected = NA
      )
      
      # Iteration  actions
      n_remaining <- network %>%
        dplyr::filter(n_parallel_segments == -1) %>%
        dplyr::pull(edge_id) %>%
        length()
      next
    }
    
    # Subset - Segment ---------------------------------------------------------
    # Remove all instances of the selected segment from the network
    ## A segment won't be parallel with itself (assumption)
    network_compare_iter <- network_compare_iter %>%
      dplyr::filter(segment_id != segment_id_iter)
    
    n_subset <- length(network_compare_iter$edge_id)
    
    if (verbose >= 3){
      print(glue::glue('Number of Edges - Subset - Segment ID:\t\t\t{n_subset}'))
    }
    
    if (n_subset == 0) {
      network <- exportParallelEdgeResults(
        network                    = network,
        input_edge_id              = edge_id_iter,
        parallel_detected          = FALSE,
        parallel_edges_detected    = NA,
        parallel_segments_detected = NA
      )
      
      # Iteration  actions
      n_remaining <- network %>%
        dplyr::filter(n_parallel_segments == -1) %>%
        dplyr::pull(edge_id) %>%
        length()
      next
    }
    
    # PARALLEL INTERSECTION CHECK #############################################
    # Initialize results vector
    parallel_edges_iter <- c()
    parallel_segments_iter <- c()
    
    # Counter for printing
    j <- 1
    
    for (comparison_edge_idx in seq_along(network_compare_iter$edge_id)) {
      
      # Select orthogonal LINESTRING of main edge to check with comparison edge
      # intersection
      ls_orth <- network_orth %>%
        dplyr::filter(edge_id == edge_id_iter) %>%
        dplyr::pull(geometry)
      
      # Select Row to compare
      row_comp <- network_compare_iter %>% 
        dplyr::slice(comparison_edge_idx)
      
      # Extract parameters for comparison
      ls_comp         <- row_comp %>%  dplyr::pull(geometry)
      edge_id_comp    <- row_comp %>%  dplyr::pull(edge_id)
      segment_id_comp <- row_comp %>%  dplyr::pull(segment_id)
      
      # Check if the orthogonal projection and edge intersect
      ls_intersect <- sf::st_intersects(
        ls_orth, 
        ls_comp,
        sparse = FALSE
      )[[1]]
      
      how_close <- NA
      
      # Document edge and segment if they did intersect
      if (ls_intersect) {
        parallel_edges_iter    <- c(parallel_edges_iter, edge_id_comp)
        parallel_segments_iter <- c(parallel_segments_iter, segment_id_comp)
        
      } else {
        how_close <- round(sf::st_distance(ls_orth, ls_comp),2)
      }
      
      if (verbose >= 4) {
        print(glue::glue("Edge ID:\t\t\t\t\t\t{edge_id_iter}"))
        print(glue::glue("Edge ID - Comparison:\t\t\t\t\t{edge_id_comp}"))
        print(glue::glue("Segment ID:\t\t\t\t\t{segment_id_iter}"))
        print(glue::glue("Segment ID - Comparison:\t\t\t\t{segment_id_comp}"))
        print(glue::glue("Intersect:\t\t\t\t\t\t{ls_intersect}"))
        
        if (!is.na(how_close)){
          print(glue::glue("Distance:\t\t\t\t\t\t{how_close}"))
        }
        
        if (j < length(network_compare_iter$edge_id)) {
          print(glue::glue(strrep("- ", 40)))
        }
      }
      
      j <- j + 1
    }
    
    # Check if any edges were found to be parallel
    if (length(parallel_segments_iter) > 0){
      
      network <- exportParallelEdgeResults(
        network                    = network,
        input_edge_id              = edge_id_iter,
        parallel_detected          = TRUE,
        parallel_edges_detected    = parallel_edges_iter,
        parallel_segments_detected = parallel_segments_iter
      )
      
    } else {
      
      network <- exportParallelEdgeResults(
        network                    = network,
        input_edge_id              = edge_id_iter,
        parallel_detected          = FALSE,
        parallel_edges_detected    = NA,
        parallel_segments_detected = NA
      )
      
    }
    
    # Iteration Actions
    n_remaining <- network %>%
      dplyr::filter(n_parallel_segments == -1) %>%
      dplyr::pull(edge_id) %>%
      length()
    
  }
  
  # REMOVE SIMPLE FEATURES ####################################################
  # Drop simple features column. It causes problems when merging all the subset
  # networks together. It is easier and more reliable to just regenerate the 
  # simple feature geometries later when needed.s
  network <- network %>% 
    dplyr::select(-geometry)
  
  # Return ####################################################################
  if (verbose >= 1) {
    time_elapsed <- secondsToPeriodTime(start_time)
    print(glue::glue(""))
    print(glue::glue(strrep("_", 80)))
    print(glue::glue("Compare Parallel Segments - Complete:\t\t\t{time_elapsed}"))
  }
  
  network
}


#' Identify parallel segments - single cluster
#' 
#' The following function identifies all cases of parallel edges within a 
#' street network, subsetted to a single cluster of the entire network.
#' 
#' @param cluster_label Cluster label
#' @param network_w_cluster Street network which has been assigned a cluster label
#'  in the column **.cluster** (dataframe)
#' @param crop_compare_network Flag indicating whether to reduce the size of the 
#' comparison network with the function **cropComparisonNetwork** (logical)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#'  
#' @return Street network with a new column of listing the segments connected
#' 
identifyParallelSegmentsCluster <- function(
  cluster_label,
  network_w_cluster,
  crop_compare_network = FALSE, 
  start_time = NA, 
  verbose = 0 
  ) {
  library(dplyr)
  
  # Establish starting time
  if (is.na(start_time)){
    start_time <- Sys.time()
  }
  
  # PROCESSING ################################################################
  segment_label_a <- "Cluster Sub-network processing"
  segment_label_b <- "Cluster simple features processing"
  segment_label_c <- "Reduce comparison network - latitude/longitude"
  segment_label_d <- "Reduce comparison network - simple features"
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label_a} - {cluster_label} - Start:\t\t\t\t{elapsed_time}"))
  }
  
  if (verbose >= 2){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label_b} - {cluster_label} - Start:\t\t\t\t{elapsed_time}"))
  }
  
  print(glue::glue("Norm/orth processing"))
  # Simple features processing
  network_norm_orth <- networkLinestringProcessing(
    network_w_cluster, 
    parallel = FALSE
  )
  
  if (verbose >= 2){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label_b} - {cluster_label} - Complete:\t\t\t{elapsed_time}"))
  }

  # Define complete network
  network_complete <- network_norm_orth[["normal"]] %>% 
    dplyr::select(-.cluster)
  
  # Reduce to cluster label sub-network
  network <- network_norm_orth[["normal"]] %>% 
    dplyr::filter(.cluster == cluster_label)
  
  network_orth <- network_norm_orth[["orthogonal"]] %>% 
    dplyr::filter(.cluster == cluster_label)
  
  # Reduce Comparison Network -------------------------------------------------
  ## Latitude / Longitude
  if (verbose >= 2){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label_c} - {cluster_label} - Start:\t\t{elapsed_time}"))
  }
  
  # Define latitude/longitude bounding box of cluster area
  geom_cluster <- network %>% dplyr::pull(geometry)
  bbox_cluster <- sf::st_bbox(geom_cluster)
  
  # Determine center latitude/longitude point
  lon_center <- mean(c(bbox_cluster[1] , bbox_cluster[3]))
  lat_center <- mean(c(bbox_cluster[2], bbox_cluster[4]))
  
  # Determine distance from bbox corner to center
  bbox_corner_to_center <- sqrt(
    (lon_center - bbox_cluster[1])^2 + (lat_center - bbox_cluster[2])^2
    )
  
  # Determine latitude/longitude subset radius
  deg_to_m_conv <- 1.11e5
  dist_adj_factor <- 1.25
  subset_radius = (bbox_corner_to_center * dist_adj_factor * deg_to_m_conv)
  
  # Subset complete network to within radius of bbox center point
  network_compare <- subset_lat_lon(
    network    = network_complete, 
    lat_center = lat_center, 
    lon_center = lon_center, 
    radius     = subset_radius
  )
  
  if (verbose >= 2){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label_c} - {cluster_label} - Complete:\t\t{elapsed_time}"))
  }
  
  
  # Reduce Complete Network to Radius around Cluster
  if (crop_compare_network) {
    if (verbose >= 2){
      end_time <- Sys.time()
      elapsed_time <- secondsToPeriodTime(start_time)
      print(glue::glue("{segment_label_d} - {cluster_label} - Start:\t\t{elapsed_time}"))
      print(glue::glue("Number of Cluster Edges - {cluster_label}:\t\t\t\t\t\t{length(network$edge_id)}"))
      print(glue::glue("Number of Comparison Edges - {cluster_label}:\t\t\t\t\t\t{length(network_compare$edge_id)}"))
    }
    
    network_compare <- cropComparisonNetwork(
      network_cluster  = network,
      network_complete = network_compare,
      dist             = 100, 
      verbose          = verbose
      )
    
    if (verbose >= 2){
      end_time <- Sys.time()
      elapsed_time <- secondsToPeriodTime(start_time)
      print(glue::glue("{segment_label_d} - {cluster_label} - Complete:\t{elapsed_time}"))
    }
  } else{
    if (verbose >= 2){
      end_time <- Sys.time()
      elapsed_time <- secondsToPeriodTime(start_time)
      print(glue::glue("{segment_label_d} - {cluster_label} - Step Skipped:\t\t{elapsed_time}"))
    }
  }

  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label_a} - {cluster_label} - Complete:\t\t\t\t{elapsed_time}"))
  }
  
  # COMPARE || EDGES ##########################################################
  compareParallelEdges(
    network         = network,
    network_orth    = network_orth,
    network_compare = network_compare, 
    start_time      = start_time, 
    verbose         = verbose
  )
  
}


#' Identify parallel segments through parallel processing
#' 
#' The following function executes the **identifyParallelSegmentsCluster** utilizing
#' parallel processing to increase performance speed..
#' 
#' @param network Street network (simple features dataframe)
#' @param k Number of cores to use for parallel process (integer)
#' @param crop_compare_network Flag indicating whether to reduce the size of the 
#' comparison network with the function **cropComparisonNetwork** (logical)
#' @param start_time Time stamp to base all print outputs off of. Must be of
#'  the form "YYYY-MM-DD HH:MM:SS" (double)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#' @log_progress Flag indicating if a log file should be kept of the parallel 
#'  processing progress (logical)
#'  
#' @return Street network with a new column of listing the segments connected
#
identifyParallelSegmentsPar <- function(
  network, 
  k = 6, 
  crop_compare_network = FALSE,
  start_time = NA, 
  verbose = 0,
  log_progress = FALSE
  ) {

  # Establish starting time
  if (is.na(start_time)){
    start_time <- Sys.time()
  }
  
  # CLUSTER LABELING ##########################################################
  segment_label <- "Network Processing"
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label} - Start:\t\t{elapsed_time}"))
  }
  
  # Generate cluster subset datasets
  network.list <- networkClusters(network, k)
  
  # Extract cluster label names
  cluster_labels <- names(network.list)
  
  # Recombine cluster subsets
  network_w_cluster <- dplyr::bind_rows(network.list)
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label} - Complete:\t\t{elapsed_time}"))
  }
  
  # PARALLEL CLUSTER PROCESSNG ################################################
  segment_label <- "Parallel Cluster Processing"
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time)
    print(glue::glue("{segment_label} - Start:\t{elapsed_time}"))
  }
  
  # Define dataset variables and required sub-functions for export to cluster
  cluster_data <- c("network.list", "k", "start_time")
  cluster_functions <- c(
    "exportParallelEdgeResults", "preProcessingParallelIdentification", 
    "makeLinestring", "makeLinestringVector", "networkClusters",
    "secondsToPeriodTime", "compareParallelEdges", 
    "networkLinestringProcessing", "lineStringLength", "subset_lat_lon", 
    "line_length", "cropComparisonNetwork"
  )
  cluster_variables <- c(cluster_data, cluster_functions)
  
  # Create Cluster
  file_dt <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  file_name <- glue::glue("identifyParallelSegmentsPar_{file_dt}.log")
  file_path <- file.path("./..", "parallel_processing_logs", file_name)
  # print(file_path)
  if (log_progress) {
    cl <- parallel::makeCluster(k,  outfile = file_path)
  } else {
    cl <- parallel::makeCluster(k)
  }

  
  # Export required variables to cluster
  parallel::clusterExport(
    cl      = cl, 
    varlist = cluster_variables, 
    envir   = environment()
  )
  
  # Execute parallel implementation of segment labeling
  network.mtx <- parallel::parSapply(
    cl                   = cl, 
    X                    = cluster_labels, 
    FUN                  = identifyParallelSegmentsCluster,
    crop_compare_network = crop_compare_network,
    network_w_cluster    = network_w_cluster,
    start_time           = start_time,
    verbose              = verbose
  )
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  if (verbose >= 1){
    end_time <- Sys.time()
    elapsed_time <- secondsToPeriodTime(start_time, end_time)
    print(glue::glue("{segment_label} - Complete:\t{elapsed_time}"))
  }
  
  # PROCESS CLUSTER OUTPUTS ###################################################
  # Compile network matrix dataframes
  network <- compileParSapplyOutput(
    network.mtx = network.mtx, 
    k           = k, 
    start_time  = start_time, 
    verbose     = verbose
  )
  # Return ####################################################################
  network
  
}