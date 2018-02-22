###############################################################################
#                                    ____  _
#                          _ __ ___ |  _ \| | _____  __
#                         | '_ ` _ \| |_) | |/ _ \ \/ /
#                         | | | | | |  __/| |  __/>  <
#                         |_| |_| |_|_|   |_|\___/_/\_\
#
###############################################################################
###############################################################################
# Release Functions
###############################################################################

#' Make List of Mosquito Releases
#'
#' Sets up a release schedule for a single patch, returns a list to be used in
#' \code{\link{oneDay_Releases_Patch}}
#'
#' @param releaseStart Day releases start
#' @param releaseEnd Day releases end
#' @param releaseInterval Interval between releases
#' @param genMos List of genotypes for new Mosquitoes
#' @param numMos Integer number of Mosquitoes to create
#' @param minAge Integer specifying the minimum age of Mosquitoes
#' @param maxAge Integer specifying the maximum age of Mosquiotes
#' @param ageDist Distribution for ages of Mosquitoes. Must be length(maxAge-minAge+1)
#'
#' @details This function creates a release schedule and the genotypes and ages
#' of mosquitoes to be released at that time. The release schedule is determined
#' by releaseStart, releaseEnd, and releaseInterval. The mosquitoes to be released
#' are determined by genMos (a list of genotypes), numMos (a vector of how many of
#' each genotype to release), and minAge/maxAge/ageDist. The last three determine
#' the age distribution of the mosquitoes to be released.
#'
#' @return List of release dates and the population to be released on that day.
#'
#' @examples
#' # to setup for 3 patches but only release in the first with a defined release schedule:
#'
#' patchReleases = replicate(n = 3,
#'                           expr = list(maleReleases = NULL,femaleReleases = NULL,larvaeReleases=NULL),
#'                           simplify = FALSE)
#'
#' set.seed(42)
#' patchReleases[[1]]$femaleReleases = Release_basicRepeatedReleases(releaseStart = 5,
#' releaseEnd = 30,
#' releaseInterval = 5,
#' genMos = c("A", "B", "C"),
#' numMos = c(10, 20, 30),
#' minAge = 1,
#' maxAge = 10,
#' ageDist = rep(x = 1, times = 10-1+1)/(10-1+1)) #uniform distribution of ages
#'
#' set.seed(42)
#' patchReleases[[1]]$maleReleases = Release_basicRepeatedReleases(releaseStart = 50,
#' releaseEnd = 60,
#' releaseInterval = 1,
#' genMos = c("A", "B", "C"),
#' numMos = c(10, 20, 30),
#' minAge = 5,
#' maxAge = 10,
#' ageDist = c(0,0,0,0,0,1)) #all mosqutioes will be 10
#'
#' set.seed(42)
#' patchReleases[[1]]$larvaeReleases = Release_basicRepeatedReleases(releaseStart = 1,
#' releaseEnd = 5,
#' releaseInterval = 1,
#' genMos = c("A", "B", "C"),
#' numMos = c(10, 20, 30),
#' minAge = 5,
#' maxAge = 10,
#' ageDist = c(0,1,0,1,0,0)) #half of mosquitoes are 6, half are 8
#'
#' @export
Release_basicRepeatedReleases <- function(releaseStart, releaseEnd, releaseInterval, genMos, numMos, minAge, maxAge, ageDist){
  #genMos is a list of genotypes to relaese
  #numMos is a vector of the number of mosquitoes you want to make, corresponding
  #  to the genotypes of genMos

  #minAge, maxAge are the min/max age range. To get a single age, must be length
  # 2 with a c(1,0) vector in ageDist
  #ageDist - probabilities to sample from for age range. must be length
  #  minAge:maxAge

  # check timing of releases
  if(releaseInterval > (releaseEnd - releaseStart)){
    stop("interval between releases cannot be greater than time between start and end of releases")
  }

  # Initialize release times. Initialize return list
  releaseTimes = seq(from=releaseStart,to = releaseEnd,by = floor(releaseInterval))
  releaseList = vector(mode="list",length=length(releaseTimes))

  # Initialize genotypes and ages
  genotypeVector = rep.int(x = genMos, times = numMos)
  ageVector <- sample(x = minAge:maxAge, size = sum(numMos), replace = TRUE, prob = ageDist)

  # check for male/female/larvae. Fill appropriate list.
  for(tx in 1:length(releaseTimes)){
    releaseList[[tx]]$genVec = genotypeVector
    releaseList[[tx]]$ageVec = ageVector
    releaseList[[tx]]$tRelease = releaseTimes[[tx]]
  }

  return(releaseList)
}

###############################################################################
# Network Initialization Functions
###############################################################################

#' Create new Mosquito objects
#'
#' Creates new Mosquito objects specifically for simulation setup
#'
#' @usage CreateMosquitoes_Distribution_Genotype(numMos, minAge, maxAge, ageDist, aTypes)
#'
#' @param numMos Integer number of Mosquitoes to create
#' @param minAge Integer specifying the minimum age of Mosquitoes
#' @param maxAge Integer specifying the maximum age of Mosquiotes
#' @param ageDist Distribution for ages of Mosquitoes. Must be length(maxAge-minAge+1)
#' @param aTypes Nested list containing the alleles and their probabilities at each locus
#'
#' @details This function creates new mosquitoes. It is similar to
#'  \code{\link{CreateMosquitoes_Defined_Genotype}}, but intended for simulation setup.
#' It takes a list of alleles and their frequency in the population, then creates
#' genotypes based on that distribution. Then it samples the age distribtuion, and
#' creates new Mosquitoes with the proper distributions of genotypes and age.
#'
#' @return List of Mosquito objects length(sum(numMos))
#'
#' @examples
#' set.seed(42)
#' numMos <- 10
#' minAge <- 1
#' maxAge <- 10
#' ageDist <- rep(x = 1, times = 10-1+1)/(10-1+1) #uniform distribution
#'
#' alleloTypes <- vector(mode = "list", length = 2) #2 loci
#' alleloTypes[[1]]$alleles <- c("W","H")
#' alleloTypes[[1]]$probs <- c(1,0)
#' alleloTypes[[2]]$alleles <- c("W","H")
#' alleloTypes[[2]]$probs <- c(1,0)
#'
#' CreateMosquitoes_Distribution_Genotype(numMos, minAge, maxAge, ageDist, AlleloTypes)
#'
#' @export
CreateMosquitoes_Distribution_Genotype <- function(numMos, minAge, maxAge, ageDist, aTypes){

  population <- vector(mode = "list", length = numMos)
  genotypes <- vector(mode = "list", length = length(aTypes))

  for(i in 1:numMos){

    #generate genotypes from distribution
    for(locus in 1:length(genotypes)){
      hold <- sample(x = aTypes[[locus]]$alleles,
                                   size = 2, replace = T,
                                   prob = aTypes[[locus]]$probs)
      genotypes[[locus]] <- sort(x = hold)
    }

    #generate age
    holdAge <- sample(x = minAge:maxAge, size = 1, replace = FALSE, prob = ageDist)

    #create new mosquito
    population[[i]] <- Mosquito(genotype = paste0(unlist(genotypes), collapse = ""),
                                    age = holdAge)
  }

  return(population)
}

###############################################################################
# Post-processing of Output
###############################################################################

#' Split Output by Patch
#'
#' Split each run output into multiple files by patch.
#'
#' @usage splitOutput(directory)
#'
#' @param directory directory where output was written to; must not end in path seperator
#' @param multiCore Set whether or not to use multiple cores for writing output
#'
#' @return *.csv files for each patch and each run
#' @export
splitOutput <- function(directory, multiCore=FALSE){
  dirFiles = list.files(path = directory, pattern = ".*\\.csv$")

  if(multiCore){
    nThread <- getDTthreads()
  } else{
    nThread <- 1
  }

  for(file in dirFiles){
    cat("processing ",file,"\n",sep="")
    fileIn = data.table::fread(input = file.path(directory, file))
    # for each file, get all the patches and split into multiple files
    for(patch in unique(fileIn$Patch)){
      patchIn = fileIn[Patch==patch]
      patchName = gsub(pattern = ".csv",replacement = paste0("_Patch",patch,".csv"),x = file)
      data.table::fwrite(x = patchIn, file = file.path(directory, patchName), nThread = nThread)
    }
    cat("removing ",file,"\n",sep="")
    file.remove(file.path(directory, file))
  }
}

#' Analyze output for mPlex-mLoci or DaisyDrive
#'
#' This function takes all the files in a directory and analyzes the population by
#' genotype of interest. It saves output by run, and inside each run it contains
#' arrays corresponding to the male and female of each patch. The arrays are organzed
#' by time, then each genotype and total population, and the array depth is each
#' patch. Data are analyzed by matching the genotypes of interest.
#'
#' @param readDirectory directory where output was written to; should not end in path seperator
#' @param saveDirectory directory to save analyzed data. Default is readDirectory
#' @param genotypes A list of each locus containing the genotypes of interest at that locus. Default is all genotypes
#' @param collapse A vector of each locus containing TRUE/FALSE. If TRUE, the genotypes of interest at that locus are collapsed and the output is the sum of all of them.
#'
#' @return A *.rds object containing list(metaData=character, maleData=array, femaleData=array)
#' @export
AnalyzeOutput_mLoci_Daisy <- function(readDirectory, saveDirectory=NULL, filename, genotypes, collapse){

  #get list of all files, unique runs, and unique patches
  dirFiles = list.files(path = readDirectory, pattern = ".*\\.csv$")
  runID = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Run[0-9]+", text = dirFiles)))
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch[0-9]+", text = dirFiles)))
  patches = patches[order(as.integer(substring(text = patches, first = 6)))]

  #import one file:get simTime, check genotypes for safety checks
  testFile <- data.table::fread(input = file.path(readDirectory, dirFiles[1]))
  simTime <- data.table::uniqueN(x = testFile$Time)

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
  if(pcre_config()["JIT"]){options(PCRE_use_JIT = TRUE)}

  #do collapse if there is some
  for(i in 1:length(genotypes)){
    #if null, look at all possible genotypes at that locus
    if( is.null(genotypes[[i]]) ){
      genotypes[[i]] <- c("HH","HR","HS","HW","RR","RS","RW","SS","SW","WW")
    }
    #if collapse is true, collapse the genotypes so all are searched for as one
    if(collapse[i]){
      genotypes[[i]] <- paste0("(",paste0(genotypes[[i]], collapse = "|") , ")", collapse = "")
    }
  }
  #expand all combinations of alleles at each site
  gOI <- expand.grid(genotypes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  #bind all combinations into complete genotypes
  gOI <- do.call(what = paste0, args = gOI)


  #create arrays to store information
  mArray = fArray = array(data = 0, dim = c(simTime, length(gOI)+2, length(patches)),
                          dimnames = list(NULL, c("Time", gOI, "Total Pop."), patches) )
  mArray[,1,] = fArray[,1,] = 1:simTime
  note <- "THIS IS A NOTE ABOUTE THE DATA. Make it reproducible."

  #loop over each run
  for(run in runID){
    #loop over each patch
    for(patch in patches){
      #read in male/female files for this run and patch
      mName = grep(pattern = paste("ADM", run, patch, sep = ".*"),
                   x = dirFiles,ignore.case = FALSE, perl = TRUE,
                   value = TRUE, useBytes = TRUE)[1]
      mFile = data.table::fread(input = file.path(readDirectory, mName))
      fName = grep(pattern = paste("ADF", run, patch, sep = ".*"),
                   x = dirFiles, ignore.case = FALSE, perl = TRUE,
                   value = TRUE, useBytes = TRUE)[1]
      fFile = data.table::fread(input = file.path(readDirectory, fName))

      #loop over simulation time
      for(loopTime in 1:simTime){
        #subset time objects for ease of reading
        mTimeObj <- mFile[Time==loopTime, Genotype]
        fTimeObj <- fFile[Time==loopTime, Genotype]
        #loop over genotypes of interest
        for(gen in gOI){
          #match genotype pattens, store how many were found
          mArray[loopTime, gen, patch] <- length(grep(pattern = gen, x = mTimeObj,
                                                      ignore.case = FALSE, perl = TRUE,
                                                      value = TRUE, useBytes = TRUE))
          fArray[loopTime, gen, patch] <- length(grep(pattern = gen, x = fTimeObj,
                                                      ignore.case = FALSE, perl = TRUE,
                                                      value = TRUE, useBytes = TRUE))
        }#end gOI loop

        #set total population
        if(!all(grepl(pattern = "NULL", x = mTimeObj, fixed = TRUE))){
          mArray[loopTime, "Total Pop.", patch] <- length(mTimeObj)
        }
        if(!all(grepl(pattern = "NULL", x = fTimeObj, fixed = TRUE))){
          fArray[loopTime, "Total Pop.", patch] <- length(fTimeObj)
        }

      }#end time loop
    }#end patch loop

    #save output for each run.
    if(is.null(saveDirectory)){saveDirectory <- readDirectory}
    fileName <- paste0(format(x = Sys.Date(), "%Y%m%d"), "_", run, "_",
                       filename, ".rds")

    saveRDS(object = list(metaData=note, maleData=mArray, femaleData=fArray),
            file = file.path(saveDirectory,fileName),
            compress = "gzip")

  }#end run loop
}#end function

#' Analyze output for mPlex-oLocus
#'
#' @details This function takes all the files in a directory and analyzes the population by
#' genotype of interest. It saves output by run, and inside each run it contains
#' arrays corresponding to the male and female of each patch. The arrays are organzed
#' by time, then each genotype and total population, and the array depth is each
#' patch. Data are analyzed by matching the genotypes of interest.
#'
#' @param readDirectory directory where output was written to; should not end in path seperator
#' @param saveDirectory directory to save analyzed data. Default is readDirectory
#' @param alleles A list of lists that contain the genotypes of interest at each locus. Default is all genotypes
#' @param collapse A list of lists containing TRUE/FALSE for each locus. If TRUE, the genotypes of interest at that locus are collapsed and the output is the sum of all of them.
#'
#' @return A *.rds object containing list(metaData=character, maleData=array, femaleData=array)
#' @export
AnalyzeOutput_oLocus <- function(readDirectory, saveDirectory=NULL, filename, alleles, collapse){

  #must give in order: "H", "R", "S", "W"
  #alleles <- list(list(c("H","W"),"W", "W"),list(c("H","W"),NULL, "W"))
  #collapse <- list(c(F,T,F), c(T,F,F))

  #get list of all files, unique runs, and unique patches
  dirFiles = list.files(path = readDirectory, pattern = ".*\\.csv$")
  runID = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Run[0-9]+", text = dirFiles)))
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch[0-9]+", text = dirFiles)))
  patches = patches[order(as.integer(substring(text = patches, first = 6)))]

  #import one file:get simTime, check genotypes for safety checks
  testFile <- data.table::fread(input = file.path(readDirectory, dirFiles[1]))
  simTime <- data.table::uniqueN(x = testFile$Time)

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
    for(inner in 1:length(alleles[[1]])){
      #if null, look at all possible genotypes at that locus
      if( is.null(alleles[[outer]][[inner]]) ){
        alleles[[outer]][[inner]] <- c("H", "R", "S", "W")
      }
      #if collapse is true, collapse the genotypes so all are searched for as one
      if(collapse[[outer]][inner]){
        alleles[[outer]][[inner]] <- paste0("(",paste0(alleles[[outer]][[inner]], collapse = "|") , ")", collapse = "")
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


  #create arrays to store information
  mArray = fArray = array(data = 0, dim = c(simTime, length(gOI)+2, length(patches)),
                          dimnames = list(NULL, c("Time", gOI, "Total Pop."), patches) )
  mArray[,1,] = fArray[,1,] = 1:simTime
  note <- "THIS IS A NOTE ABOUTE THE DATA. Make it reproducible."


  #loop over each run
  for(run in runID){
    #loop over each patch
    for(patch in patches){
      #read in male/female files for this run and patch
      mName = grep(pattern = paste("ADM", run, patch, sep = ".*"),
                   x = dirFiles,ignore.case = FALSE, perl = TRUE,
                   value = TRUE, useBytes = TRUE)[1]
      mFile = data.table::fread(input = file.path(readDirectory, mName))
      fName = grep(pattern = paste("ADF", run, patch, sep = ".*"),
                   x = dirFiles, ignore.case = FALSE, perl = TRUE,
                   value = TRUE, useBytes = TRUE)[1]
      fFile = data.table::fread(input = file.path(readDirectory, fName))

      #loop over simulation time
      for(loopTime in 1:simTime){
        #subset time objects for ease of reading
        mTimeObj <- mFile[Time==loopTime, Genotype]
        fTimeObj <- fFile[Time==loopTime, Genotype]
        #loop over genotypes of interest
        for(gen in gOI){
          #match genotype pattens, store how many were found
          mArray[loopTime, gen, patch] <- length(grep(pattern = gen, x = mTimeObj,
                                                      ignore.case = FALSE, perl = TRUE,
                                                      value = TRUE, useBytes = TRUE))
          fArray[loopTime, gen, patch] <- length(grep(pattern = gen, x = fTimeObj,
                                                      ignore.case = FALSE, perl = TRUE,
                                                      value = TRUE, useBytes = TRUE))
        }#end gOI loop

        #set total population
        if(!all(grepl(pattern = "NULL", x = mTimeObj, fixed = TRUE))){
          mArray[loopTime, "Total Pop.", patch] <- length(mTimeObj)
        }
        if(!all(grepl(pattern = "NULL", x = fTimeObj, fixed = TRUE))){
          fArray[loopTime, "Total Pop.", patch] <- length(fTimeObj)
        }

      }#end time loop
    }#end patch loop

    #save output for each run.
    if(is.null(saveDirectory)){saveDirectory <- readDirectory}
    fileName <- paste0(format(x = Sys.Date(), "%Y%m%d"), "_", run, "_",
                       filename, ".rds")

    saveRDS(object = list(metaData=note, maleData=mArray, femaleData=fArray),
            file = file.path(saveDirectory,fileName),
            compress = "gzip")

  }#end run loop
}#end function

###############################################################################
# Random Others
###############################################################################

#' R Dirichlet
#'
#' Calculates a Dirichlet distribution.
#'
#' @usage JaredDirichlet(n, alpha)
#'
#' @param n Number of draws to perform
#' @param alpha vector length(n) as a shape parameter
#'
#' @details This function makes draws from a Dirichlet dristribution. It is slow,
#' but entirely written in R.
#'
#' @return matrix(nrow=n, ncol=length(alpha))
#'
#' @examples
#' set.seed(42)
#' JaredDirichlet(n=4, alpha = c(0.1,0.2,0.3,0.4))
#'
#' @export
JaredDirichlet <- function(n=1,alpha){
  Gam <- matrix(0,n,length(alpha))
  for(i in 1:length(alpha)) {Gam[,i] <- rgamma(n,shape=alpha[i])}
  Gam/rowSums(Gam)
}

#' calc_parameter
#'
#' Gets the genotype-dependent parameter of choice
#'
#' @usage get_parameter(paramList, genotype)
#'
#' @param paramList sublist from referenceList, created by Make_reference functions.
#' @param genotype Genotype of interest
#'
#' @details This function gets whatever genotype-dependent factor desired. The inputs
#' are whatever parameter list of interest, then the genotype it is applied to.
#'
#' @return double
#'
#' @examples
#' set.seed(42)
#' JaredDirichlet(n=4, alpha = c(0.1,0.2,0.3,0.4))
#'
#' @export
calc_parameter <- function(paramList, genotype){

  #default return value
  holdValue = 1

  #counter for locus in loop
  counter=1
  for(i in 1:length(paramList)){

    #get the genotype at that locus, then pull value from parameters
    locus <- substr(x = genotype, start = counter, stop = counter+1) #use character vector and subset? Would be better.
    test <- paramList[[i]][[locus]]

    #if the value exists, use it.
    if(!is.null(test)){
      holdValue <- holdValue*test
    }

    #increment counter
    counter <- counter + 2
  }

  return(holdValue)
}

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
#' @usage Plot_mPlex(file, totalPop=TRUE)
#'
#' @param file Path to .gzip file from analyze function
#' @param totalPop Boolean, to plot the total population or not.
#'
#' @details This function plots output from the analyze function. Setting totalPop
#' to FALSE keeps it from plotting the total population.
#'
#' @export
Plot_mPlex <- function(file = NULL, totalPop = TRUE){

  #keep old plot parameters to reset later
  oldPar <- par(no.readonly = TRUE)

  #Get the data to plot
  Data <- readRDS(file = file)

  #Get genotypes and number of patches
  genotypes <- dimnames(Data$maleData)[[2]][-1]
  patches <- dimnames(Data$maleData)[[3]]
  numPatches <- length(patches)
  numGen <- length(genotypes)

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
  matplot(Data$femaleData[,1+(1:numGen), patches[1]], type = "l", lty = 1,
          main = "Female Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(Data$femaleData[,1+(1:numGen), patches[1]])),
          xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
          col = col)
  title(ylab = "Population", line = 2)
  box(lwd = 2)
  grid()

  #female
  par(mar = c(2,2,3,1), las = 1)
  matplot(Data$maleData[,1+(1:numGen),patches[1]], type = "l", lty = 1,
          main = "Male Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(Data$maleData[,1+(1:numGen), patches[1]])),
          xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
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
      matplot(Data$femaleData[,1+(1:numGen), patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(Data$femaleData[,1+(1:numGen), patch])),
              xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
              col = col)
      title(ylab = "Population", line = 2)
      box(lwd = 2)
      grid()

      par(mar = c(2,2,1,1))
      matplot(Data$maleData[,1+(1:numGen),patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(Data$maleData[,1+(1:numGen), patch])),
              xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
              col = col)
      mtext(patch, side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
      box(lwd = 2)
      grid()
    }#end patch loop
  }#end if

  #reset par()
  par(oldPar)
}














