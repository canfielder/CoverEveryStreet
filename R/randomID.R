#' Random ID Generator
#'
#' Generate a random alphanumeric id.
#' 
#' @param utilized_ids Previously utilized ids. (character vector)
#' @param n_char Character length of the id. (integer)
#'
#' @return Unique alphanumeric random id of n_char length
#'

randomID <- function(utilized_ids, n_char=10) {
  
  id <- stringi::stri_rand_strings(
    length = n_char, 
    n = 1, 
    pattern ='[a-z0-9]' 
  )
  
  if (id %in% utilized_ids) {
    randomID(utilized_ids, n_char=n_char) 
  
    } else {
    
      id
  }
}