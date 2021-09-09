#' Edge Properties
#' The following functions are related to determining if a segment which 
#' runs parallel to at least two other segments should be split into 
#' multiple segments.

#' Determine the slope of a linear line
#' 
#' @param x1 X coordinate for point 1 (numeric)
#' @param y1 Y coordinate for point 1 (numeric)
#' @param x2 X coordinate for point 2 (numeric)
#' @param y2 Y coordinate for point 2 (numeric)
#' 
#' @return Slope of a line between points one and two (numeric)
#' 
lineSlope <- function(x_1, y_1, x_2, y_2) {
  (y_2 - y_1) / (x_2 - x_1)
}


#' Determine the y-intercept of a linear line
#' 
#' @param x X coordinate for a point along the line (numeric)
#' @param y Y coordinate for a point along the line  (numeric)
#' @param slope Slope of the line (numeric)
#' 
#' @return Y-intercept of a line through the point x,y with 
#'  the provided slope
#' 
yIntercept <- function(x, y, slope){
  # y = mx + b >>> b = y - mx
  y - slope * x
}


#' Determine the distance between two points
#' 
#' Calculates the length of the line between point 1 and point 2. This is based
#' one a 2D coordinate system, utilizing the Pythagorean Theorem.
#' 
#' @param x1 X coordinate for point 1 (numeric)
#' @param y1 Y coordinate for point 1 (numeric)
#' @param x2 X coordinate for point 2 (numeric)
#' @param y2 Y coordinate for point 2 (numeric)
#' 
#' @return Distance between points one and two (numeric)
#' 
distPointToPoint <- function(x1, y1, x2, y2){
  sqrt((y2 - y1)^2 + (x2 - x1)^2)
}


#' Determine the distance between a line and a point
#' 
#' Calculates the distance between a set point and point along a line, 
#' given the x coordinate value of a point along the line.
#' 
#' @param x_line X coordinate for a point along the line (numeric)
#' @param x_point X coordinate for the point
#' @param y_point Y coordinate for the point
#' @param slope Slope of the line (numeric)
#' @param y_intercept Y intercept of the line
#' 
#' @return Distance between a point along the line and a set point (numeric)
#' 
distLineToPoint <- function(x_line, x_point, y_point, slope, y_intercept){
  
  # Determine y value along line given x
  y_line <- slope * x_line + y_intercept
  
  # Determine x_line between point and line
  distPointToPoint(x_l, y_line, x_point, y_point)
}


#' Determine the minimum distance between a line and a point
#' 
#' Calculates the minimum distance between a point and a line.
#' 
#' @param x1 X coordinate at one end of the line (numeric)
#' @param x2 X coordinate at one end of the line (numeric)
#' @param slope Slope of the line (numeric)
#' @param y_intercept Y intercept of the line (numeric)
#' @param length Y Length of the line (numeric)
#' @param x_point X coordinate at the point (numeric)
#' @param y_point Y coordinate at the point (numeric)
#' @param seq_length Number of points along the line to test (numeric)
#' 
#' @return Minimum distance between a point and a line (numeric)
# 
minimumDistance <- function(x1, x2, slope, y_intercept, length, x_point, y_point, seq_length = 100){
  # Create sequence of x values
  x_sequence <- seq(
    from = x1, to = x2 , by = length/seq_length)
  
  # Optimize distance calculation from points along edge to intersection point
  optimize(
    f = distLineToPoint, 
    x_sequence, 
    x_point, 
    y_point, 
    slope,
    y_intercept
  )
}


#' Generate a table of all possible pairs from a list of values
#' 
#' @param input_list List of values (list)
#' 
#' @return Table of all possible pairs of the values provided in input_list 
#'  (dataframe)
# 
listPairs <- function(input_list){
  as.data.frame(
    x = t(combn(unlist(input_list), 2)) 
  ) %>% 
    rename(
      item_1 = V1,
      item_2 = V2
    )
}


#' Determine the shared point between two segments
#' 
#' @param segment_id_1 Segment ID for segment 1
#' @param segment_id_2 Segment ID for segment 2
#' @param node_table Table with properties of all To and From nodes. This 
#'  parameters must be wrapped in a list (dataframe)
#' 
#' @return Table of all possible pairs of the values provided in input_list 
#'  (dataframe)
# 
sharedShortSegmentPoint <- function(segment_id_1, segment_id_2, node_table) {
  
  # Extract nodes related to segment 1
  seg_1_nodes <- node_table %>% 
    filter(
      (segment_id == segment_id_1)
    ) %>% 
    pull(node_id)
  
  # Extract nodes related to segment 1
  seg_2_nodes <- node_table %>% 
    filter(
      (segment_id == segment_id_2)
    ) %>% 
    pull(node_id)
  
  # Determine intersection
  intersect_node <- intersect(seg_1_nodes, seg_2_nodes)
  
  # Determine longitude/latitude of intersection point
  short_intersect_lon <- node_table %>% 
    filter(node_id == intersect_node) %>% 
    pull(node_lon) %>% 
    unique()
  
  short_intersect_lat <- node_table %>% 
    filter(node_id == intersect_node) %>% 
    pull(node_lat) %>% 
    unique()
  
  # Export lon/lat vector
  short_intersect_point <- c(short_intersect_lon, short_intersect_lat)
  names(short_intersect_point) <- c('short_intersect_lon', 'short_intersect_lat')
  short_intersect_point
}


#' Verify that two segments are adjacent
#' 
#' @param segment_id_1 Segment ID for segment 1
#' @param segment_id_2 Segment ID for segment 2
#' @param segment_map Mapping list of the adjacent segments to each segment. 
#'  This input must be wrapped in a list. (vector)
#' 
#' @return Results if segment 2 is adjacent to segment 1 (bool)
# 
verifyAdjacentSegments <- function(segment_id_1, segment_id_2){
  
  # Generate vector of adjacent segments
  segment_1_adj <- segment_map[segment_id_1][[1]]
  
  # Check if second segment is in vector of segment 1 adjacent segments
  segment_id_2 %in% segment_1_adj
  
}