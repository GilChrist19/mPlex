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

#' Plot mPlex
#'
#' Plots the analyzed output of mPlex.
#'
#' @usage Plot_mPlex(directory, whichPatches=NULL,totalPop=TRUE)
#'
#' @param directory Path to directory of analyzed output files
#' @param whichPatches Vector of patches to plot, must be less than 15
#' @param totalPop Boolean, to plot the total population or not.
#'
#' @details This function plots output from the analyze function. totPop is not currently 
#' used. If there are several repetitions available, it only plots the first one.
#'
#' @export
Plot_mPlex <- function(directory, whichPatches = NULL, totalPop = TRUE){
  
  #keep old plot parameters to resetB later
  oldPar <- par(no.readonly = TRUE)
  
  #Get the data to plot, remove aggKey
  dirFiles <- list.files(path = directory, pattern = ".*\\.csv$")
  dirFiles <- dirFiles[-grep(pattern = "AggKey", x = dirFiles, fixed = TRUE, useBytes = TRUE)]
  
  #test file to get sizes
  columnNames <- scan(file = file.path(directory, dirFiles[1]),
                      what = character(), sep = ",", quiet = TRUE, nlines = 1)
  testFile <- matrix(data = scan(file = file.path(directory, dirFiles[1]), what = integer(),
                                            sep = ",", skip = 1, quiet = TRUE),
                                ncol = length(columnNames), byrow = TRUE)
  
  # get genotypes, num patches, etc
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch_[0-9]+", text = dirFiles)))

  # check if user chose specific patches
  if(length(patches)>15 && is.null(whichPatches)){
    stop("There are more than 15 patches in the simulation.
         Please select less than 15 to plot.")
  }
  if(!is.null(whichPatches)){
    # make sure not too many. I dont' know what that number is, but at some point
    #  the plotting function breaks down.
    if(length(whichPatches)>15){
      stop("Please select less than 15 patches.")
    }
    # select the patches, if they select less than the total number of patches.
    if(length(whichPatches)<=length(patches)){
      patches = patches[whichPatches]
    }
  }
  
  genotypes <- columnNames[-1]
  numPatches <- length(patches)
  numGen <- length(genotypes)
  
  # create data list holder, and then fill it
  maleData=femaleData=array(data = 0L, dim = c(nrow(testFile),numGen+1,numPatches),
                            dimnames = list(NULL, columnNames, patches))
                   
                   
  for(patch in patches){
    names = grep(pattern = patch, x = dirFiles, fixed = TRUE, value = TRUE)
    femaleData[,,patch] = matrix(data = scan(file = file.path(directory, dirFiles[1]), what = integer(),
                                            sep = ",", skip = 1, quiet = TRUE),
                                ncol = numGen+1, byrow = TRUE)
    maleData[,,patch] = matrix(data = scan(file = file.path(directory, dirFiles[2]), what = integer(),
                                            sep = ",", skip = 1, quiet = TRUE),
                                ncol = numGen+1, byrow = TRUE)
  }
  
  
  # not sure how to implement this now, so just not using it :-)
  #if(!totalPop){numGen <- numGen-1}
  
  col <- ggCol_utility(n = numGen)
  
  #setup plot layout
  lmatrix <- matrix(data = 1:(numPatches*3), nrow = numPatches, ncol = 3, byrow = TRUE)
  if(numPatches>1){
    #fill in rest of plot labels
    lmatrix[2:numPatches, c(1,2)] <- matrix(data = 4:(3+2*(numPatches-1)),
                                            ncol = 2, byrow = TRUE)
    #legend gets whole right side
    lmatrix[,3] <- 3
  }
  
  layout(lmatrix, widths = c(3,3,1))
  
  #plot first patch and the legend
  #male
  par(mar = c(2,3,3,1), las = 1, font.lab = 2, font.axis = 2, font.main = 2, cex.main = 1.75)
  matplot(femaleData[,1+(1:numGen), patches[1]], type = "l", lty = 1,
          main = "Female Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(femaleData[,1+(1:numGen), patches[1]])),
          xlim = c(0, dim(maleData)[1]), yaxs = "i", xaxs = "i",
          col = col, panel.first=grid())
  title(ylab = "Population", line = 2)
  box(lwd = 2)

  #female
  par(mar = c(2,2,3,1), las = 1)
  matplot(maleData[,1+(1:numGen),patches[1]], type = "l", lty = 1,
          main = "Male Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(maleData[,1+(1:numGen), patches[1]])),
          xlim = c(0, dim(maleData)[1]), yaxs = "i", xaxs = "i",
          col = col, panel.first=grid())
  mtext(patches[1], side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
  box(lwd = 2)
  
  #legend
  par(mar = c(0,0,0,0), font=2)
  plot.new()
  legend(x = "left", legend = genotypes[1:numGen] , col = col,
         bty = "n", bg = "transparent",lty = 1, lwd=3,cex = 1)
  
  
  ##rest of the patches
  if(numPatches>1){
    for(patch in patches[-1]){
      par(mar = c(2,3,1,1), las = 1)
      matplot(femaleData[,1+(1:numGen), patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(femaleData[,1+(1:numGen), patch])),
              xlim = c(0, dim(maleData)[1]), yaxs = "i", xaxs = "i",
              col = col, panel.first=grid())
      title(ylab = "Population", line = 2)
      box(lwd = 2)
      
      par(mar = c(2,2,1,1))
      matplot(maleData[,1+(1:numGen),patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(maleData[,1+(1:numGen), patch])),
              xlim = c(0, dim(maleData)[1]), yaxs = "i", xaxs = "i",
              col = col, panel.first=grid())
      mtext(patch, side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
      box(lwd = 2)
    }#end patch loop
  }#end if
  
  #reset par()
  par(oldPar)
}





