#' Replace Dataframe Rows
#' 
#' Replaces rows in a dataframe. The dataframe structure remains the same, only
#' the values in the replaced rows may change.
#' 
#' @param network Street network (dataframe)
#' @param k Number of clusters to break network into (integer)
#'
#' @return Network split into separate dataframes based on cluster label (list)
#
networkClusters <- function(network, k){
  # Extract kmeans inputs
  points <- network %>% 
    select(from_lon, from_lat)
  
  # Generate cluster labels and assign to network
  network <- points %>%  
    stats::kmeans(centers = k, nstart = 1) %>% 
    broom::augment(points) %>% 
    dplyr::select(.cluster) %>% 
    dplyr::bind_cols(network)
  
  # Split network into sub-networks based on clusters
  split(
    x = network , 
    f = network$.cluster
  )
}