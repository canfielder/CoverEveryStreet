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
#' @param df Dataframe with from_lon, from_lat, to_lon, and to_lat columns.
#'   (dataframe)
#'
#' @return Shapefile Dataframe with LINESTRING object as new column
#'
makeLinestringVector <- function(df) { 
  df %>% 
  select(from_lon, from_lat, to_lon, to_lat) %>% 
    purrr::pmap(makeLinestring) %>% 
    sf::st_as_sfc(crs = 4326) %>% 
    {
      tibble(
        edge_id = df$edge_id, 
        sf_linestring = .)
    }%>% 
    left_join(
      y = df,
      by = "edge_id"
    ) %>% 
    sf::st_sf()
}