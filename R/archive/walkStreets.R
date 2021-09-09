#' Identify Segments Between Nodes
#' 
#' Assigns a unique identifier to all lines between two nodes. 
#' 
#' @param network.df Dataframe which represents a street network (dataframe)
#' @param input_edge_id Network street edge id to start process at (char)
#' @param segment_id_rdm Randomly generated identifier to apply to segment 
#'  (char)
#' @param intersection_nodes Vector of to and from nodes identified as 
#'  intersections (character vector)
#' @param n_recurse Counter for keeping track of how many recursive iterations 
#'  have occurred (int)
#'
#' @return Network streets dataframe with segment_id applied
#' 
walkStreets <- function(
  network.df, input_edge_id, segment_id_rdm, intersection_nodes, n_recurse
  ) {

  recurse_count <- n_recurse + 1

  # Recursion Pit Check
  if (recurse_count > 100) {
    print(glue('Recurse limit reached.'))
    print(glue('Edge ID at Limit: {input_edge_id}'))
    exit()
  }
  
  # print(glue(' Edge: {input_edge_id}'))
  # print(glue('Number of Edge Id(s): {length(input_edge_id)}'))
  
  # ACTION --------------------------------------------------------------------
  # Filter network data to selected edge
  row.df <- network.df %>% 
    filter(
      edge_id == input_edge_id
    ) 
  
  # print(glue('Size of row dataframe: {length(row.df$edge_id)}'))
  
  # Segment ID - Before Processing
  segment_id_exist <- row.df %>% pull(segment_id)
  
  # Assign segment id to selected row
  row.df <- row.df %>% 
    mutate(segment_id = segment_id_rdm)
  
  # Select TO node
  to_node <- row.df %>% pull(to_id)
  # print(glue('To node: {to_node}'))
  
  # Replace sampled row in main dataset
  network.df <- replaceDataframeRow(network.df, row.df)
  
  # RECURSION CHECKS ----------------------------------------------------------
  
  # Does Street Continue
  # Determine if the to_node of the current edge continues as the 
  # from node of a different edge
  bool_dead_end <- !(to_node %in% network.df$from_id)
  # print(glue('Dead End:\t\t{bool_dead_end}'))
  
  # Intersection Check
  # Determine if the TO node is an intersection
  bool_intersection <- to_node %in% intersection_nodes
  # print(glue('Intersection:\t\t{bool_intersection}'))
  
  # Loop Check
  # This check verifies that the segment we are walking is not a connected 
  # loop, like a track. If we detect it is a loop, we need to break this 
  # process once the loop is complete. 
  bool_loop <- segment_id_exist == segment_id_rdm
  # print(glue('Loop:\t\t{bool_loop}'))
  # print(glue('Segment - Previous:\t\t{segment_id_exist}'))
  # print(glue('Segment - Current:\t\t{segment_id_rdm}'))
  
  # RECURSION PATH ------------------------------------------------------------
  
  if (any(bool_dead_end, bool_intersection, bool_loop)) {
    # print(glue('Return dataframe.'))
    # Return update dataframe
    network.df
    
  } else {
    # print('Continue Walk.')
    # Extract edge_id of subsequent street
    input_edge_id <- network.df %>% 
      filter(from_id == to_node) %>% 
      pull(edge_id)
    
    # Walk
    walkStreets(
      network.df, input_edge_id, segment_id_rdm, intersection_nodes, n_recurse
      )
  }
}



