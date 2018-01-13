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

CreateMosquitoes_Defined_Genotype <- function(genMos, numMos, minAge, maxAge, ageDist){
  #genMos is a list of genotypes to relaese
  #numMos is a vector of the number of mosquitoes you want to make, corresponding
  #  to the genotypes of genMos

  #minAge, maxAge are the min/max age range. To get a single age, must be length
  # 2 with a c(1,0) vector in ageDist
  #ageDist - probabilities to sample from for age range. must be length
  #  minAge:maxAge

  #return list
  population <- vector(mode = "list", length = sum(numMos))

  #external counter
  count = 1L

  #loop over each genotype
  for(gen in 1:length(genMos)){
    #loop over number of mosquitoes of that genotype
    for(num in 1:numMos[gen]){
      #generate age
      holdAge <- sample(x = minAge:maxAge, size = 1, replace = FALSE, prob = ageDist)

      #create new mosquito
      population[[count]] <- Mosquito$new(genotype = genMos[gen], age = holdAge)

      count = count + 1L
    }
  }

  return(population)
}


Release_basicRepeatedReleases <- function(releaseStart, releaseEnd, releaseInterval, releaseVector, sex="M"){

  # check timing of releases
  if(releaseInterval > (releaseEnd - releaseStart)){
    stop("interval between releases cannot be greater than time between start and end of releases")
  }

  # name and check releaseVector. Initialize release times. Initialize return list
  releaseTimes = seq(from=releaseStart,to = releaseEnd,by = floor(releaseInterval))
  releaseList = vector(mode="list",length=length(releaseTimes))

  # check for male/female/larvae. Fill appropriate list.
  if(sex=="M"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$nuM = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
    }
  } else if(sex=="F"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$nuF = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
    }
  } else if(sex=="L"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$larvae = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
    }
  } else {
    stop(paste0("expected character in 'M','F','L' in argument 'sex', got: ",sex))
  }
  return(releaseList)
}

###############################################################################
# Network Initialization Functions
###############################################################################
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
    population[[i]] <- Mosquito$new(genotype = paste0(unlist(genotypes), collapse = ""),
                                    age = holdAge)
  }

  return(population)
}

CreateMosquitoes_Eggs <- function(genMos, numMos){
  #genMos is a list of genotypes to relaese
  #numMos is a vector of the number of mosquitoes you want to make, corresponding
  #  to the genotypes of genMos


  #return list
  population <- vector(mode = "list", length = sum(numMos))

  #external counter
  count = 1L

  #loop over each genotype
  for(gen in 1:length(genMos)){
    #skip if there are no mosquitoes of this genotype
    if(numMos[gen]==0){next}
    #loop over number of mosquitoes of that genotype
    for(num in 1:numMos[gen]){
      #create new mosquito
      population[[count]] <- Mosquito$new(genotype = genMos[gen], age = 0)

      count = count + 1L
    }
  }

  #list of new mosquitoes
  return(population)

}

###############################################################################
# Post-processing of Output
###############################################################################
#' Split Output by Patch
#'
#' Split each run output into multiple files by patch.
#'
#' @param directory directory where output was written to; must not end in path seperator
#'
#' @export
splitOutput <- function(directory){
  dirFiles = list.files(path = directory)

  # for each file read it in
  for(file in dirFiles){
    cat("processing ",file,"\n",sep="")
    fileIn = read.csv(file.path(directory, file))
    # for each file, get all the patches and split into multiple files
    for(patch in unique(fileIn$Patch)){
      patchIn = fileIn[fileIn$Patch==patch,]
      patchName = gsub(pattern = ".csv",replacement = paste0("_Patch",patch,".csv"),x = file)
      write.csv(x = patchIn,file = file.path(directory, patchName),row.names = FALSE)
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
#' @param collapse A list of each locus containing TRUE/FALSE. If TRUE, the genotypes of interest at that locus are collapsed and the output is the sum of all of them.
#'
#' @return A *.rds object containing list(metaData=character, maleData=array, femaleData=array)
#' @export
AnalyzeOutput_mLoci_Daisy <- function(readDirectory, saveDirectory=NULL, genotypes, collapse){

  #get list of all files, unique runs, and unique patches
  dirFiles = list.files(path = readDirectory, pattern = ".*\\.csv$")
  runID = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Run[0-9]+", text = dirFiles)))
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch[0-9]+", text = dirFiles)))

  #import one file:get simTime, check genotypes for safety checks
  testFile <- read.csv(file = file.path(readDirectory, dirFiles[1]),
                       header = TRUE, stringsAsFactors = FALSE)
  simTime <- unique(testFile$Time)

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
      genotypes[[i]] <- "(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)"
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
  mArray = fArray = array(data = 0, dim = c(length(simTime), length(gOI)+2, length(patches)),
                          dimnames = list(NULL, c("Time", gOI, "Total Pop."), patches) )
  mArray[,1,] = fArray[,1,] = simTime
  note <- "THIS IS A NOTE ABOUTE THE DATA. Make it reproducible."


  #loop over each run
  for(run in runID){
    #loop over each patch
    for(patch in patches){
      #read in male/female files for this run and patch
      mName = grep(pattern = paste("ADM", run, patch, sep = ".*"),
                   x = dirFiles, ignore.case = FALSE, value = TRUE)
      mFile = read.csv(file = file.path(readDirectory, mName),
                       header = TRUE, stringsAsFactors = FALSE)
      fName = grep(pattern = paste("ADF", run, patch, sep = ".*"),
                   x = dirFiles, ignore.case = FALSE, value = TRUE)
      fFile = read.csv(file = file.path(readDirectory, fName),
                       header = TRUE, stringsAsFactors = FALSE)

      #loop over simulation time
      for(loopTime in simTime){
        #subset time objects for ease of reading
        mTimeObj <- mFile$Genotype[mFile$Time == loopTime]
        fTimeObj <- fFile$Genotype[fFile$Time == loopTime]
        #loop over genotypes of interest
        for(gen in gOI){
          #match genotype pattens, store how many were found
          mArray[loopTime, gen, patch] <- length(grep(pattern = gen, x = mTimeObj, ignore.case = FALSE))
          fArray[loopTime, gen, patch] <- length(grep(pattern = gen, x = fTimeObj, ignore.case = FALSE))
        }#end gOI loop

        #set total population
        mArray[loopTime, "Total Pop.", patch] <- length(mFile$Genotype[mFile$Time == loopTime])
        fArray[loopTime, "Total Pop.", patch] <- length(fFile$Genotype[fFile$Time == loopTime])
      }#end time loop
    }#end patch loop

    #save output for each run.
    if(is.null(saveDirectory)){saveDirectory <- readDirectory}
    fileName <- paste0(format(x = Sys.Date(), "%Y%m%d"), "_", run, "_",
                       paste0(gOI, collapse = "_"), ".rds")

    saveRDS(object = list(metaData=note, maleData=mArray, femaleData=fArray),
            file = file.path(saveDirectory,fileName),
            compress = "gzip")

  }#end run loop
}#end function


#' Analyze output for mPlex-oLocus
#'
#' This function takes all the files in a directory and analyzes the population by
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
AnalyzeOutput_oLocus <- function(readDirectory, saveDirectory=NULL, alleles, collapse){

  #must give in order: "H", "R", "S", "W"
  #alleles <- list(list(c("H","W"),"W", "W"),list(c("H","W"),NULL, "W"))
  #collapse <- list(c(F,T,F), c(T,F,F))

  #get list of all files, unique runs, and unique patches
  dirFiles = list.files(path = readDirectory, pattern = ".*\\.csv$")
  runID = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Run[0-9]+", text = dirFiles)))
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch[0-9]+", text = dirFiles)))

  #import one file:get simTime, check genotypes for safety checks
  testFile <- read.csv(file = file.path(readDirectory, dirFiles[1]),
                       header = TRUE, stringsAsFactors = FALSE)
  simTime <- unique(testFile$Time)

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
        alleles[[outer]][[inner]] <- "(H|R|S|W)"
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
  mArray = fArray = array(data = 0, dim = c(length(simTime), length(gOI)+2, length(patches)),
                          dimnames = list(NULL, c("Time", gOI, "Total Pop."), patches) )
  mArray[,1,] = fArray[,1,] = simTime
  note <- "THIS IS A NOTE ABOUTE THE DATA. Make it reproducible."


  #loop over each run
  for(run in runID){
    #loop over each patch
    for(patch in patches){
      #read in male/female files for this run and patch
      mName = grep(pattern = paste("ADM", run, patch, sep = ".*"),
                   x = dirFiles, ignore.case = FALSE, value = TRUE)
      mFile = read.csv(file = file.path(readDirectory, mName),
                       header = TRUE, stringsAsFactors = FALSE)
      fName = grep(pattern = paste("ADF", run, patch, sep = ".*"),
                   x = dirFiles, ignore.case = FALSE, value = TRUE)
      fFile = read.csv(file = file.path(readDirectory, fName),
                       header = TRUE, stringsAsFactors = FALSE)

      #loop over simulation time
      for(loopTime in simTime){
        #subset time objects for ease of reading
        mTimeObj <- mFile$Genotype[mFile$Time == loopTime]
        fTimeObj <- fFile$Genotype[fFile$Time == loopTime]
        #loop over genotypes of interest
        for(gen in gOI){
          #match genotype pattens, store how many were found
          mArray[loopTime, gen, patch] <- length(grep(pattern = gen, x = mTimeObj, ignore.case = FALSE))
          fArray[loopTime, gen, patch] <- length(grep(pattern = gen, x = fTimeObj, ignore.case = FALSE))
        }#end gOI loop

        #set total population
        mArray[loopTime, "Total Pop.", patch] <- length(mFile$Genotype[mFile$Time == loopTime])
        fArray[loopTime, "Total Pop.", patch] <- length(fFile$Genotype[fFile$Time == loopTime])
      }#end time loop
    }#end patch loop

    #save output for each run.
    if(is.null(saveDirectory)){saveDirectory <- readDirectory}
    fileName <- paste0(format(x = Sys.Date(), "%Y%m%d"), "_", run, "_",
                       paste0(gOI, collapse = "_"), ".rds")

    saveRDS(object = list(metaData=note, maleData=mArray, femaleData=fArray),
            file = file.path(saveDirectory,fileName),
            compress = "gzip")

  }#end run loop
  }#end function

###############################################################################
# Random Others
###############################################################################

JaredDirichlet <- function(n=1,alpha){
  Gam <- matrix(0,n,length(alpha))
  for(i in 1:length(alpha)) {Gam[,i] <- rgamma(n,shape=alpha[i])}
  Gam/rowSums(Gam)
}

















