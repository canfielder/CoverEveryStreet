#' Subset Street Networks
#' 
#' The following functions are related to reducing a street networks sized
#' based on a specific method.

#' Line Length
#'
#' Generate the length of a straight line. This can be done using a 
#' two diimensional assumption, or 3D projection
#' 
#' @param lat_1 Latitude of Point 1 (numeric)
#' @param lon_1 Longitude of Point 1 (numeric)
#' @param lat_2 Latitude of Point 1 (numeric)
#' @param lon_2 Longitude of Point 2 (numeric)
#' @param proj Boolean indicating whether to calculate the line length as a 
#'  3D projection (True) or 2D line (False)
#' @param crs Projection code (numeric)
#' @param m_to_deg_conv Assumption for convering latitude/longitude degrees 
#'  to meters (numeric).
#'
#' @return Length of the line, in meters, defined by Point 1 and 2.
line_length <- function(
  lat_1, lon_1, lat_2, lon_2, proj = FALSE, 
  crs = 4326, m_to_deg_conv = 111000){
  if (proj){
      #` Convert two points to a line and calculate length, in meters
      sf::st_sfc(
        sf::st_linestring(
          rbind(
            c(lat_1, lon_1), 
            c(lat_2, lon_2)
          )
        ), 
        crs = crs
      ) %>% 
        sf::st_length() %>% 
        as.vector() 
    
  } else {
    # Calculate the length of the vector between points 1 and 2, assuming a 2D projection
    deg_lat <- abs(lat_2 - lat_1)
    deg_lon <- abs(lon_2 - lon_1)
    deg_vect <- sqrt(deg_lat^2 + deg_lon^2)
    
    # Convert degrees to meters
    deg_vect * m_to_deg_conv
  }
}


#' Line Length - Projection
#'
#' Generate the length of a shapefile LINESTRING.
#' 
#' @param lat_1 Latitude of Point 1 (numeric)
#' @param lon_1 Longitude of Point 1 (numeric)
#' @param lat_2 Latitude of Point 1 (numeric)
#' @param lon_2 Longitude of Point 2 (numeric)
#'
#' @return Length of the line, in meters, defined by Point 1 and 2.
#' 
line_length_proj <- function(lat_1, lon_1, lat_2, lon_2, crs = 4326) {
  #` Convert two points to a line and calculate length, in meters
  sf::st_sfc(
    sf::st_linestring(
      rbind(
        c(lat_1, lon_1), 
        c(lat_2, lon_2)
      )
    ), 
    crs = crs
  ) %>% 
    sf::st_length() %>% 
    as.vector()
}


#' Subset Network
#'
#' Reduce a street network dataframe to all strings within a radius of
#' a selected point (latitude/longitude). This function requires a dataframe
#' with the following columns: to_lat, to_lon, from_lat, from_lon.
#'
#' @param network Street network (dataframe)
#' @param lat_center Latitude of the center of the subsetted network (numeric)
#' @param lon_center Longitude of the center of the subsetted network (numeric)
#' @param radius Radius for subsetting network, in meters (numeric)
#' 
#' @return All streets within the input radius
subset_lat_lon <- function(
  network, 
  lat_center,
  lon_center, 
  radius = 50
  ){
  
  network %>%
    mutate(
      to_dist = purrr::map2_dbl(
        to_lat, to_lon,
        ~line_length(
          lat_1 = lat_center,
          lat_2 =.x,
          lon_1 = lon_center,
          lon_2 = .y
        )
      ),
      from_dist = purrr::map2_dbl(
        from_lat, from_lon,
        ~line_length(
          lat_1 = lat_center,
          lat_2 =.x,
          lon_1 = lon_center,
          lon_2 = .y
        )
      )
    ) %>%     
    filter(
      (to_dist   <= radius) | 
      (from_dist <= radius)
    ) %>% 
    select(-c(to_dist, from_dist))
  
}