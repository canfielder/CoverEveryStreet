#' Define network nodes which halt segment walks
#' 
#' The following function generates a vector of all network nodes that stop 
#' a segment walk loop. These nodes include:
#' 
#' - Intersection: Node is connected to more than two edges
#' - Parallel Split: Node is defined to split a long segment which is parallel
#'   to multiple shorter segments in order to match the short segment connection 
#'   point.
#' - Dead Ends: Nodes are only connected to a single edge
#' 
#' @param network Street network (dataframe)
#' @param file_date Date for intersection and parallel file names. Must be in
#'  the form YYYY-MM-DD (character)
#' @param verbose Flag indicating the granularity of print outputs the function
#'  generates. The larger the number, the more granular (integer) 
#'
#' @return Network nodes which will halt a segment walk (character vector)
#
haltSegmentWalkNodes <- function(network, file_date, verbose = 0) {
  
  # Define file paths
  file_dir <- "./../data"
  
  # Define intersection node file path
  intersection_file_name <- glue::glue("intersection_nodes_{file_date}.RDS")
  intersection_path <- file.path(file_dir, intersection_file_name)
  
  # Define parallel separation node file path
  parallel_file_name <- glue::glue("parallel_separation_nodes_{file_date}.RDS")
  parallel_path <- file.path(file_dir, parallel_file_name)
  
  # Import intersection nodes
  intersection_nodes <- readr::read_rds(intersection_path)
  
  # Import parallel separation nodes
  if (file.exists(parallel_path)) {
    parallel_separation_nodes <- readr::read_rds(parallel_path)
  } else{
    # Initialize Empty Vector
    parallel_separation_nodes <- vector(mode = "character", length = 0)
  }

  # Define Dead End Nodes
  ## Defined, not loaded, because depending on if the street network is
  ## the complete network or just a subset will determine change which
  ## nodes are dead ends
  to_nodes <-  network %>% 
    dplyr::select(to_id) %>% 
    dplyr::rename(id = to_id)
  
  from_nodes <-  network %>% 
    dplyr::select(from_id) %>% 
    dplyr::rename(id = from_id)
  
  dead_end_nodes <- dplyr::bind_rows(
    from_nodes,
    to_nodes
  ) %>% 
    dplyr::count(id) %>% 
    dplyr::filter(n == 1) %>% 
    dplyr::pull(id)
  
  if (verbose >= 1){
    print(glue::glue("Halt Walk Nodes Loaded."))
  }
  
  if (verbose >= 2){
    n_inter <- length(intersection_nodes)
    n_par <- length(parallel_separation_nodes)
    n_dead_ends <- length(dead_end_nodes)
    
    print(glue::glue("Number of Each Node Type:"))
    print(glue::glue("Intersections:\t\t\t\t\t{n_inter}"))
    print(glue::glue("Parallel Separation:\t\t\t\t{n_par}"))
    print(glue::glue("Dead End:\t\t\t\t\t\t{n_dead_ends}"))
  }
  
  # Compiled into single vector
  c(intersection_nodes, parallel_separation_nodes, dead_end_nodes)
}


