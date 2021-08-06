#' Replace Dataframe Rows
#' 
#' Replaces rows in a dataframe. The dataframe structure remains the same, only
#' the values in the replaced rows may change.
#' 
#' @param df Original dataframe. (dataframe)
#' @param df_row Subsetted rows of original dataframe. (dataframe)
#'
#' @return Updated dataframe.
#

replaceDataframeRow <- function(df, df_row) {
  
  # Replace row with same edge_id
  df[match(df_row$edge_id, df$edge_id), ] <- df_row
  
  df
}