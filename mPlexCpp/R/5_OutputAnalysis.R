###############################################################################
#                            ____  __          ______          
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
#                                               /_/   /_/      
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
# GENOTYPE GROUPINGS
###############################################################################

#' Generate genotypes of interest for mPlex-mLoci or DaisyDrive
#'
#' This function generates a *.csv of genotypes of interest and the grouping scheme, ie, 
#' if one locus isn't interesting and should be considered the same. This file will be used 
#' for data analysis.
#'
#' @param outputFile Name of file to output. Must end in .csv
#' @param genotypes A list of each locus containing the genotypes of interest at that locus. Default is all genotypes
#' @param collapse A vector of each locus containing TRUE/FALSE. If TRUE, the genotypes of interest at that locus are collapsed and the output is the sum of all of them.
#'
#' @export
genOI_mLoci_Daisy <- function(outputFile, genotypes, collapse){
  
  #check that the collapse length is equal to genotype length
  if(length(genotypes) != length(collapse)){
    stop("collapse must be specified for each loci.
         length(collapse) == length(genotypes)")
  }
  
  # check for null genotypes
  for(i in 1:length(genotypes)){
    #if null, look at all possible genotypes at that locus
    if( is.null(genotypes[[i]]) ){
      genotypes[[i]] <- c("HH","HR","HS","HW","RR","RS","RW","SS","SW","WW")
    }
  }
  
  
  # generate all combinations of genotypes of interest
  holdGens <- expand.grid(genotypes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  holdGens <- do.call(what = paste0, args = holdGens)
  
  
  #do collapse if there is some
  for(i in 1:length(genotypes)){
    #if collapse is true, collapse the genotypes so all are searched for as one
    if(collapse[i]){
      genotypes[[i]] <- file.path("(", paste0(genotypes[[i]],collapse = "|"), ")", fsep = "")
    }
  }
  #expand all combinations of alleles at each site
  gOI <- expand.grid(genotypes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  #bind all combinations into complete genotypes
  gOI <- do.call(what = paste0, args = gOI)
  
  
  # create output vector
  outDF <- data.frame("Key"=holdGens, "Group"=0, stringsAsFactors = FALSE)
  
  # group output into groups
  for(i in 1:length(gOI)){
    outDF[grep(pattern = gOI[i], x = holdGens),"Group"] <- i
  }
  
  # write output
  write.csv(x = outDF, file = outputFile, row.names = FALSE)
  
}

#' Generate genotypes of interest for mPlex-oLocus
#'
#' This function generates a *.csv of genotypes of interest and the grouping scheme, ie, 
#' if one locus isn't interesting and should be considered the same. This file will be used 
#' for data analysis.
#'
#' @param outputFile Name of file to output. Must end in .csv
#' @param alleles A list of lists that contain the genotypes of interest at each locus. Default is all genotypes
#' @param collapse A list of lists containing TRUE/FALSE for each locus. If TRUE, the genotypes of interest at that locus are collapsed and the output is the sum of all of them.
#'
#' @export
genOI_oLocus <- function(outputFile, alleles, collapse){
  
  # safety checks
  #check that the number of loci is equal to the genotype length
  if(length(alleles)!=2){
    stop("There are 2 alleles in this simulation
         list(list(locus_1, locus_2), list(locus_1, locus_2))")
  }
  #check that the collapse length is equal to genotype length
  if(length(alleles) != length(collapse) || lengths(alleles) != lengths(collapse)){
    stop("collapse must be specified for each locus in each allele.
         length(collapse) == length(alleles)
         lengths(collapse) == lengths(alleles)")
  }
  
  
  # check for null genotypes
  for(outer in 1:2){
    for(inner in 1:length(collapse[[1]])){
      #if null, look at all possible genotypes at that locus
      if( is.null(alleles[[outer]][[inner]]) ){
        alleles[[outer]][[inner]] <- c("H", "R", "S", "W")
      }
    }#end loop over each loci
  }#end loop over each allele
  
  
  
  # generate all combinations of genotypes of interest
  holdAlleles <- vector(mode = "list", length = 2)
  
  #do collapse if there is some
  for(outer in 1:2){
    #collapse possible alleles
    holdAlleles[[outer]] <- do.call(what = paste0,
                                    args = expand.grid(alleles[[outer]],
                                                       KEEP.OUT.ATTRS = FALSE,
                                                       stringsAsFactors = FALSE))
    
    # loop over alleles at each locus
    for(inner in 1:length(collapse[[1]])){
      #if collapse is true, collapse the genotypes so all are searched for as one
      if(collapse[[outer]][inner]){
        alleles[[outer]][[inner]] <- file.path("(", paste0(alleles[[outer]][[inner]],collapse = "|"), ")", fsep = "")
      }
    }#end loop over each loci
    
    #expand and paste all possible loci combinations in each allele
    alleles[[outer]] <- do.call(what = paste0,
                                args = expand.grid(alleles[[outer]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
    
  }#end loop over each allele
  
  #expand/bind all combinations of alleles at each site
  gOI <- do.call(what = paste0,
                 args = expand.grid(alleles, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  #expand/bind all possible alleles
  holdAlleles <- do.call(what = paste0,
                         args = expand.grid(holdAlleles,
                                            KEEP.OUT.ATTRS = FALSE,
                                            stringsAsFactors = FALSE))
  
  
  # create output dataframe
  outDF <- data.frame("Key"=holdAlleles, "Group"=0, stringsAsFactors = FALSE)
  
  # group output into groups
  for(i in 1:length(gOI)){
    outDF[grep(pattern = gOI[i], x = holdAlleles),"Group"] <- i
  }
  
  
  # write output
  write.csv(x = outDF, file = outputFile, row.names = FALSE)
  
}

###############################################################################
# AGGREGATE
###############################################################################

#' Aggregate mPlex output
#'
#' This function takes all the files in a directory and aggregates them based on 
#' an "aggKey" file supplied in the write folder. This key is generated by "genOI*" 
#' functions for the appropriate experiment. Output is organized by time and group of 
#' interested and saved by experiment number, sex, and patch.
#' File structure as shown below: \cr
#'  * someTopDirectory
#'    * readDirectory
#'      * F_Run_001_Patch_000.csv
#'      * F_Run_002_Patch_000.csv
#'      * M_Run_001_Patch_000.csv
#'      * M_Run_002_Patch_000.csv
#'    * writeDirectory1
#'      * AggKey1.csv
#'    * writeDirectory2
#'      * AggKey2.csv
#'    * ...
#'
#' @param readDirectory Directory with simulation output, all male/female/patch/experiment number files.
#' @param writeDirectory Directory to save analyzed data. Must have aggregation key in it!
#' @param simTime Simulation time
#'
#' @export
SimAggregation <- function(readDirectory, writeDirectory, simTime){
  
  # list all files
  readFiles <- list(list.files(path = readDirectory, pattern = 'M_', full.names = TRUE),
                    list.files(path = readDirectory, pattern = 'F_', full.names = TRUE))
  
  # get file with largest size for buffer
  #  Definitely female, so only check them.
  largeFile <- readFiles[[2]][which.max(x = file.size(readFiles[[2]]))]
  
  # generate write file names
  writeFiles <- vector(mode = "list", length = 2)
  writeDirectory <- path.expand(path = writeDirectory)
  
  for(i in 1:2){
    hold <- strsplit(x = readFiles[[i]], split = "/", fixed = TRUE, useBytes = TRUE)
    writeFiles[[i]] <- file.path(writeDirectory,
                                 unlist(x = hold)[seq.int(from = 0, to = length(hold)*length(hold[[1]]), by = length(hold[[1]]))])
  }
  
  # read in genotype collapse key
  genKey <- read.csv(file = list.files(path = writeDirectory, pattern = "AggKey", full.names = TRUE),
                     header = TRUE, stringsAsFactors = FALSE)
  
  # c++ for analysis
  mPlexCpp:::simAgg(readFiles_ = readFiles, writeFiles_ = writeFiles,
                    largeFile_ = largeFile, simTime_ = simTime, genKey_ = genKey)
  
} # end function


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





