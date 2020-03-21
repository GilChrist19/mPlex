###############################################################################
#     ________ __ __  _______ 
#    / ____/ //_//  |/  / __ \
#   / /   / ,<  / /|_/ / /_/ /
#  / /___/ /| |/ /  / / _, _/ 
#  \____/_/ |_/_/  /_/_/ |_|  
#    
###############################################################################
########################################################################
# Delete files in a directory
########################################################################
eraseDirectory <- function(directory){
  # check directory exists
  if(!dir.exists(directory)){
    cat("no such directory exists\n")
    return(NULL)
  }
  dirFiles = list.files(path = directory)
  # begin deleting contents
  if(length(dirFiles)>0){
    for(i in dirFiles){
      cat("removing file: ",file.path(directory, i),"\n", sep = "")
      #file.remove(file.path(directory, i))
      unlink(x = file.path(directory, i), recursive = TRUE)
    }
  }
  # end deleting contents
}# end function


###############################################################################
# PLOTTING UTILITIES
###############################################################################
#' Utility to Imitate ggplot2 Colors
#'
#' Sample at equally spaced intervals along the color wheel
#'
#' @param n number of colors
#' @param alpha transparency
#'
#' @export
ggCol_utility <- function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}

