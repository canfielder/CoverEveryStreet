#' Visualization Functions
#' 
#' 
#' ggplot Segements
#' 
#' Generates a ggplot object using the geom_segement method.
#' 
#' @param df Dataframe which contains the following columns, at a minimum:
#'  to_lat
#'  to_lon
#'  from_lat
#'  from_lon
#'  highway
#'  tiger.name_base
#'  edge_id
#'  to_id
#'  from_id
#' @param color_feature Column feature which will provide the color aesthetic (char)
#'
#' @return ggplot object of street network
#'
plotStreets <- function(df, color_feature) {
  
  # Convert color feature to column variable
  color_feature <- rlang::sym(color_feature)
  
  # plot
  ggplot(
    data = df
    ) +
    geom_segment(
      mapping = aes(
        x      = from_lon,
        xend   = to_lon,
        y      = from_lat,
        yend   = to_lat,
        color  = !!color_feature, 
        text_1 = highway, 
        text_2 = tiger.name_base,
        text_3 = edge_id,
        text_4 = segment_id,
        text_5 = from_id,
        text_6 = to_id
      )
    ) + 
      coord_map() +
      theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "white")
        )
  
  }