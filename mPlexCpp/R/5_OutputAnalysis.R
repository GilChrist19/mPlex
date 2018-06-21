###############################################################################
#                            ____  __          ______          
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
#                                               /_/   /_/      
###############################################################################
###############################################################################
# PATCH SPLITTER
###############################################################################

#' Split Output by Patch
#'
#' Split each run output into multiple files by patch.
#'
#' @usage splitOutput(directory, numCores)
#'
#' @param directory Directory where output was written to; must not end in path seperator
#' @param numCores How many cores to use for reading/writing files
#'
#' @return *.csv files for each patch
#' @export
splitOutput <- function(directory, numCores=1){
  # get all files in directory
  dirFiles = list.files(path = directory, pattern = ".*\\.csv$")
  
  for(file in dirFiles){
    cat("processing ",file,"\n",sep="")
    fileIn = data.table::fread(input = file.path(directory, file), sep = ",",
                               header = TRUE, verbose = FALSE, showProgress = FALSE,
                               data.table = TRUE,nThread = numCores, logical01 = TRUE)
    
    # for each file, get all the patches and split into multiple files
    for(patch in unique(fileIn$Patch)){
      data.table::fwrite(x = fileIn[Patch==patch],
                         file = file.path(directory, sub(pattern = ".csv",
                                                         replacement = file.path("_Patch",patch,".csv",fsep = ""),
                                                         x = file, fixed = TRUE)), 
                         logical01 = TRUE, showProgress = FALSE, verbose = FALSE,
                         nThread = numCores)
    }
    cat("removing ",file,"\n",sep="")
    file.remove(file.path(directory, file))
  } # end file loop
} # end function

###############################################################################
# PATCH SPLITTER
###############################################################################
# NEED ANALYSIS FUNCTIONS FOR 3 DRIVES - 2 FUNCTIONS HERE








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



# NEED PLOT FUNCTION






