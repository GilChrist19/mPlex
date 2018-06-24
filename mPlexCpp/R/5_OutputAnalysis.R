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
                               data.table = TRUE,nThread = numCores,
                               logical01 = TRUE, key = "Patch")
    
    # for each file, get all the patches and split into multiple files
    for(patch in unique(fileIn$Patch)){
      
      fileName = sub(pattern = ".csv",
                     replacement = file.path("_Patch",
                                             formatC(x = patch,width = 6,format = "d",flag = "0"),
                                             ".csv",fsep = ""),
                     x = file, fixed = TRUE)
      
      data.table::fwrite(x = fileIn[J(patch)][ ,Patch:=NULL],
                         file = file.path(directory,fileName), 
                         logical01 = TRUE, showProgress = FALSE, verbose = FALSE,
                         nThread = numCores)
    }
    cat("removing ",file,"\n",sep="")
    file.remove(file.path(directory, file))
  } # end file loop
} # end function

###############################################################################
# PATCH AGGREGATION
###############################################################################

#' Analyze output for mPlex-mLoci or DaisyDrive
#'
#' This function takes all the files in a directory and analyzes the population by
#' genotype of interest. It saves output by patch. The files are organzed
#' by time, then each genotype and total population.Data are analyzed by matching the genotypes of interest.
#'
#' @param readDirectory Directory where output was written to; should not end in path seperator
#' @param saveDirectory Directory to save analyzed data. Don't save to readDirectory
#' @param genotypes A list of each locus containing the genotypes of interest at that locus. Default is all genotypes
#' @param collapse A vector of each locus containing TRUE/FALSE. If TRUE, the genotypes of interest at that locus are collapsed and the output is the sum of all of them.
#' @param numCores How many cores to read file in with
#'
#' @export
AnalyzeOutput_mLoci_Daisy <- function(readDirectory, saveDirectory,
                                      genotypes, collapse, numCores){
  
  # Check save directory
  if(saveDirectory == readDirectory){
    stop("Please create a save directory in a new directory.")
  }
  
  # get list of all files, then unique patches
  dirFiles = list.files(path = readDirectory, pattern = ".*\\.csv$")
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch[0-9]+", text = dirFiles)))
  
  # import one file:get simTime, check genotypes for safety checks
  testFile <- data.table::fread(input = file.path(readDirectory, dirFiles[1]), sep = ",",
                                header = TRUE, verbose = FALSE, showProgress = FALSE,
                                data.table = TRUE, nThread = numCores,
                                logical01 = TRUE, drop = c("Age","Mate"))
  
  simTime <- data.table::uniqueN(testFile$Time)
  
  
  #safety checks
  #check that the number of loci is equal to the genotype length
  if(nchar(testFile$Genotype[1])/2 != length(genotypes)){
    stop("genotypes must be the length of loci in the critters to analyze.
         list(c(locus_1),c(locus_2),etc)
         NULL -> all genotypes")
  }
  #check that the collapse length is equal to genotype length
  if(length(genotypes) != length(collapse)){
    stop("collapse must be specified for each loci.
         length(collapse) == length(genotypes)")
  }
  
  
  #do collapse if there is some
  for(i in 1:length(genotypes)){
    #if null, look at all possible genotypes at that locus
    if( is.null(genotypes[[i]]) ){
      genotypes[[i]] <- c("HH","HR","HS","HW","RR","RS","RW","SS","SW","WW")
    }
    #if collapse is true, collapse the genotypes so all are searched for as one
    if(collapse[i]){
      genotypes[[i]] <- file.path("(", paste0(genotypes[[i]],collapse = "|"), ")", fsep = "")
    }
  }
  #expand all combinations of alleles at each site
  gOI <- expand.grid(genotypes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  #bind all combinations into complete genotypes
  gOI <- do.call(what = paste0, args = gOI)
  
  
  
  
  #create matrices to store information, set time
  mMatrix = fMatrix = matrix(data = 0, nrow = simTime, ncol = length(gOI)+2,
                             dimnames = list(NULL, c("Time", gOI, "Total Pop.")))
  
  mMatrix[ ,1] = fMatrix[ ,1] = 0:(simTime-1)
  
  
  # initialize progress bar
  pb = txtProgressBar(min = 0,max = length(patches),style = 3)
  pbVal = 0
  
  #loop over each patch
  for(patch in patches){
    #read in male/female files for this run and patch
    # fixed and perl about same speed with PCRE_use_JIT = TRUE, fixed over 
    #  2x faster if PCRE_use_JIT = FALSE
    names = grep(pattern = patch, x = dirFiles, fixed = TRUE, value = TRUE)
    
    # female file, drop patch number, key the time
    fFile = data.table::fread(input = file.path(readDirectory, names[1]), sep = ",",
                              header = TRUE, verbose = FALSE, showProgress = FALSE,
                              data.table = TRUE, nThread = numCores,
                              logical01 = TRUE, drop = c("Age","Mate"), key = "Time")
    
    # male file, drop patch number, key the time
    mFile = data.table::fread(input = file.path(readDirectory, names[2]), sep = ",",
                              header = TRUE, verbose = FALSE, showProgress = FALSE,
                              data.table = TRUE, nThread = numCores,
                              logical01 = TRUE, drop = "Age", key = "Time")
    
    # at <150,000 rows, subsetting on data frames is faster
    #  anything above that though, the key in data.table is constant time, and 
    #  makes a huge difference
    
    
    #loop over simulation time
    for(loopTime in 1:simTime){
      # because my output is zero indexed, but R is 1 indexed. 
      zeroCountTime = loopTime - 1
      #subset time objects for ease of reading
      mTimeObj <- mFile[J(zeroCountTime), Genotype]
      fTimeObj <- fFile[J(zeroCountTime), Genotype]
      #loop over genotypes of interest
      for(gen in gOI){
        #match genotype pattens, store how many were found
        mMatrix[loopTime, gen] <- sum(grepl(pattern = gen, x = mTimeObj, useBytes = TRUE))
        fMatrix[loopTime, gen] <- sum(grepl(pattern = gen, x = fTimeObj, useBytes = TRUE))
      }#end gOI loop
      
      
      # set total population
      #  fread reads empties as NA values, so if empty, leave zero, else fill with length
      #  This is becase we3 may not have grabbed the genotypes that exist
      if(!anyNA(mTimeObj)){
        mMatrix[loopTime, "Total Pop."] <- length(mTimeObj)
      }
      if(!anyNA(fTimeObj)){
        fMatrix[loopTime, "Total Pop."] <- length(fTimeObj)
      }
    } # end time loop
    
    
    # write output for each patch
    #  female
    fileName <- sub(pattern = ".csv", replacement = "_Aggregate.csv", x = names[1], fixed = TRUE)
    write.csv(x = fMatrix, file = file.path(saveDirectory,fileName), row.names = FALSE)
    
    
    # male
    fileName <- sub(pattern = ".csv", replacement = "_Aggregate.csv", x = names[2], fixed = TRUE)
    write.csv(x = mMatrix, file = file.path(saveDirectory,fileName), row.names = FALSE)
    
    
    # some indication that it's working
    pbVal = pbVal +1
    setTxtProgressBar(pb = pb, value = pbVal)
    
    # reset total populations to default value, in case they change
    mMatrix[loopTime, "Total Pop."] = fMatrix[loopTime, "Total Pop."] = 0
    
  } # end patch loop
} # end function

#' Analyze output for mPlex-oLocus
#'
#' @details This function takes all the files in a directory and analyzes the population by
#' genotype of interest. It saves output by patch, corresponding to the male and 
#' female of each patch. Data are analyzed by matching the genotypes of interest.
#'
#' @param readDirectory directory where output was written to; should not end in path seperator
#' @param saveDirectory directory to save analyzed data. Default is readDirectory
#' @param alleles A list of lists that contain the genotypes of interest at each locus. Default is all genotypes
#' @param collapse A list of lists containing TRUE/FALSE for each locus. If TRUE, the genotypes of interest at that locus are collapsed and the output is the sum of all of them.
#' @param numCores How many cores to read files in with, default is 1
#'
#' @export
AnalyzeOutput_oLocus <- function(readDirectory, saveDirectory,
                                 alleles, collapse, numCores=1){
  
  # Check save directory
  if(saveDirectory == readDirectory){
    stop("Please create a save directory in a new directory.")
  }
  
  # get list of all files, then unique patches
  dirFiles = list.files(path = readDirectory, pattern = ".*\\.csv$")
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch[0-9]+", text = dirFiles)))
  
  # import one file:get simTime, check genotypes for safety checks
  testFile <- data.table::fread(input = file.path(readDirectory, dirFiles[1]), sep = ",",
                                header = TRUE, verbose = FALSE, showProgress = FALSE,
                                data.table = TRUE, nThread = numCores,
                                logical01 = TRUE, drop = c("Age","Mate"))
  
  simTime <- data.table::uniqueN(testFile$Time)
  
  #safety checks
  #check that the number of loci is equal to the genotype length
  if(length(alleles)!=2){
    stop("There are 2 alleles in this simulation
         list(list(locus_1, locus_2), list(locus_1, locus_2))")
  }
  if(any(nchar(testFile$Genotype[1])/2 != lengths(alleles))){
    stop("Each allele list must be the length of loci in the allele.
         list(list(locus_1, locus_2), list(locus_1, locus_2))
         NULL -> all possible alleles")
  }
  #check that the collapse length is equal to genotype length
  if(length(alleles) != length(collapse) || lengths(alleles) != lengths(collapse)){
    stop("collapse must be specified for each locus in each allele.
         length(collapse) == length(alleles)
         lengths(collapse) == lengths(alleles)")
  }
  
  
  #do collapse if there is some
  for(outer in 1:2){
    for(inner in 1:length(collapse[[1]])){
      #if null, look at all possible genotypes at that locus
      if( is.null(alleles[[outer]][[inner]]) ){
        alleles[[outer]][[inner]] <- c("H", "R", "S", "W")
      }
      #if collapse is true, collapse the genotypes so all are searched for as one
      if(collapse[[outer]][inner]){
        alleles[[outer]][[inner]] <- file.path("(", paste0(alleles[[outer]][[inner]],collapse = "|"), ")", fsep = "")
      }
    }#end loop over each loci
    
    #expand and paste all possible loci combinations in each allele
    alleles[[outer]] <- do.call(what = paste0,
                                args = expand.grid(alleles[[outer]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
    
  }#end loop over each allele
  #expand all combinations of alleles at each site
  gOI <- expand.grid(alleles, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  #bind all combinations into complete genotypes
  gOI <- do.call(what = paste0, args = gOI)
  
  
  #create matrices to store information, set time
  mMatrix = fMatrix = matrix(data = 0, nrow = simTime, ncol = length(gOI)+2,
                             dimnames = list(NULL, c("Time", gOI, "Total Pop.")))
  
  mMatrix[ ,1] = fMatrix[ ,1] = 0:(simTime-1)
  
  
  # initialize progress bar
  pb = txtProgressBar(min = 0,max = length(patches),style = 3)
  pbVal = 0
  
  
  #loop over each patch
  for(patch in patches){
    #read in male/female files for this run and patch
    # fixed and perl about same speed with PCRE_use_JIT = TRUE, fixed over 
    #  2x faster if PCRE_use_JIT = FALSE
    names = grep(pattern = patch, x = dirFiles, fixed = TRUE, value = TRUE)
    
    # female file, drop patch number, key the time
    fFile = data.table::fread(input = file.path(readDirectory, names[1]), sep = ",",
                              header = TRUE, verbose = FALSE, showProgress = FALSE,
                              data.table = TRUE, nThread = numCores,
                              logical01 = TRUE, drop = c("Age","Mate"), key = "Time")
    
    # male file, drop patch number, key the time
    mFile = data.table::fread(input = file.path(readDirectory, names[2]), sep = ",",
                              header = TRUE, verbose = FALSE, showProgress = FALSE,
                              data.table = TRUE, nThread = numCores,
                              logical01 = TRUE, drop = "Age", key = "Time")
    
    # at <150,000 rows, subsetting on data frames is faster
    #  anything above that though, the key in data.table is constant time, and 
    #  makes a huge difference
    
    
    #loop over simulation time
    for(loopTime in 1:simTime){
      # because my output is zero indexed, but R is 1 indexed. 
      zeroCountTime = loopTime - 1
      #subset time objects for ease of reading
      mTimeObj <- mFile[J(zeroCountTime), Genotype]
      fTimeObj <- fFile[J(zeroCountTime), Genotype]
      #loop over genotypes of interest
      for(gen in gOI){
        #match genotype pattens, store how many were found
        mMatrix[loopTime, gen] <- sum(grepl(pattern = gen, x = mTimeObj, useBytes = TRUE))
        fMatrix[loopTime, gen] <- sum(grepl(pattern = gen, x = fTimeObj, useBytes = TRUE))
      }#end gOI loop
      
      
      # set total population
      #  fread reads empties as NA values, so if empty, leave zero, else fill with length
      #  This is becase we3 may not have grabbed the genotypes that exist
      if(!anyNA(mTimeObj)){
        mMatrix[loopTime, "Total Pop."] <- length(mTimeObj)
      }
      if(!anyNA(fTimeObj)){
        fMatrix[loopTime, "Total Pop."] <- length(fTimeObj)
      }
    } # end time loop
    
    
    # write output for each patch
    #  female
    fileName <- sub(pattern = ".csv", replacement = "_Aggregate.csv", x = names[1], fixed = TRUE)
    write.csv(x = fMatrix, file = file.path(saveDirectory,fileName), row.names = FALSE)
    
    
    # male
    fileName <- sub(pattern = ".csv", replacement = "_Aggregate.csv", x = names[2], fixed = TRUE)
    write.csv(x = mMatrix, file = file.path(saveDirectory,fileName), row.names = FALSE)
    
    
    # some indication that it's working
    pbVal = pbVal +1
    setTxtProgressBar(pb = pb, value = pbVal)
    
    # reset total populations to default value, in case they change
    mMatrix[loopTime, "Total Pop."] = fMatrix[loopTime, "Total Pop."] = 0
    
  } # end patch loop
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
#' @details This function plots output from the analyze function. Setting totalPop
#' to FALSE keeps it from plotting the total population. The function breaks down
#'  if you plot too many patches, but it looks bad long before that.
#'
#' @export
Plot_mPlex <- function(directory, whichPatches = NULL, totalPop = TRUE){
  
  #keep old plot parameters to reset later
  oldPar <- par(no.readonly = TRUE)
  
  #Get the data to plot
  dirFiles = list.files(path = directory, pattern = ".*\\.csv$")
  
  #test file to get sizes
  columnNames <- scan(file = file.path(directory, dirFiles[1]), what = character(),sep = ",", quiet = TRUE, nlines = 1)
  testFile <- matrix(data = scan(file = file.path(directory, dirFiles[1]), what = integer(),
                                            sep = ",", skip = 1, quiet = TRUE),
                                ncol = length(columnNames), byrow = TRUE)
  
  # get genotypes, num patches, etc
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch[0-9]+", text = dirFiles)))

  # check if user chose specific patches
  if(!is.null(whichPatches)){
    # make sure not too many. I dont' know what that number is, but at some point
    #  the plotting function breaks down.
    if(length(whichPatches)>14){
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
  
  
  
  if(!totalPop){numGen <- numGen-1}
  
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
          col = col)
  title(ylab = "Population", line = 2)
  box(lwd = 2)
  grid()
  
  #female
  par(mar = c(2,2,3,1), las = 1)
  matplot(maleData[,1+(1:numGen),patches[1]], type = "l", lty = 1,
          main = "Male Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(maleData[,1+(1:numGen), patches[1]])),
          xlim = c(0, dim(maleData)[1]), yaxs = "i", xaxs = "i",
          col = col)
  mtext(patches[1], side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
  box(lwd = 2)
  grid()
  
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
              col = col)
      title(ylab = "Population", line = 2)
      box(lwd = 2)
      grid()
      
      par(mar = c(2,2,1,1))
      matplot(maleData[,1+(1:numGen),patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(maleData[,1+(1:numGen), patch])),
              xlim = c(0, dim(maleData)[1]), yaxs = "i", xaxs = "i",
              col = col)
      mtext(patch, side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
      box(lwd = 2)
      grid()
    }#end patch loop
  }#end if
  
  #reset par()
  par(oldPar)
}





