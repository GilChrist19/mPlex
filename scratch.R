###############################################################################
#                                    ____  _
#                          _ __ ___ |  _ \| | _____  __
#                         | '_ ` _ \| |_) | |/ _ \ \/ /
#                         | | | | | |  __/| |  __/>  <
#                         |_| |_| |_|_|   |_|\___/_/\_\
#
###############################################################################

#Things for parameters
numMos <- 10
age <- 10
alleloTypes <- vector(mode = "list", length = 3)

alleloTypes[[1]]$alleles <- c("W","H")
alleloTypes[[1]]$probs <- c(1,0)
alleloTypes[[2]]$alleles <- c("W","H")
alleloTypes[[2]]$probs <- c(1,0)
alleloTypes[[3]]$alleles <- c("W","H")
alleloTypes[[3]]$probs <- c(1,0)




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


for(critter in test){

  critter <- NULL

}





test <- CreateMosquitoes_Defined_Genotype(numMos = c(2,2,2), genMos = c("AABBCC", "DDDDDD", "DSDF"),minAge = 10, maxAge = 11,ageDist = c(1,0))
listedTest <- replicate(n = 10, expr = test, simplify = FALSE)
listedTest2 <- replicate(n = 10, expr = listedTest, simplify = FALSE)


numNodes = 10
hold <- vector(mode = "list", length = numNodes)
for(i in 1:numNodes){
  for(j in 1:numNodes){
    hold[[i]] <- c(hold[[i]], listedTest2[[i]][[j]])
  }
}



eggsHist <- rlnorm(n = 10000, meanlog = log(x = 1), sdlog = log(x = 2))
larvaHist <- rlnorm(n = 1000000, meanlog = log(x = 1+14+1), sdlog = log(x = 1.3))
larvaHist <- rpois(n = 10000, lambda = 15)
larvaHist <- rgamma(n = 10000, shape = 15, scale = 1)


hold <- hist(x = larvaHist, breaks = 50, plot = FALSE)

hold$counts <- hold$counts/sum(hold$counts)

plot(hold, xlim = c(0, max(larvaHist)))

reference <- MakeReferenceMultiplex_oLocus(H=c(0.9, 0.9,0.9),R=c(0,0,0), S=c(0,0,0)/3, d=c(0,0,0))






###############################################################################
# Original one-locus multiplex Function
###############################################################################

#given the alleles, this will create a vector of all unique genotypes, and
#  fill a cube as traditionally used in MGDrivE

#assume the inputs are
#this tells me there are 3 sites, with 2 potential resistant alleles, and background mutation
##

Multiplex_FULL <- function(H=c(0.9, .6), R=0.0001, S=R/3, d=0.0001){
  #     H=homing alleles, given as a vector of homing rate. The length of that
  #         vector is the multiplex number
  #     R=resistant allele, damaging, from NHEJ
  #     S=resistant allele, not damaging, from NHEJ. defined as 1/3 the rate of R
  #     d=background mutations

  gtype <- c("W", "H", "R", "S")

  library("gtools")
  library("stringi")
  library("matrixStats")

  #experimental
  library("slam")
  library(tractor.base)


  ###############################################################################
  #generate all genotypes, set up vectors and matrices

  #gets all allele possibilities.
  # expand.grid is faster, but requires a list.
  openGtype <- gtools::permutations(n = length(gtype),r =length(H) ,v = gtype, set = FALSE, repeats.allowed = TRUE )

  #convert to list, openGtype is not always the same dimension
  # ncol() is always equal to length(H)
  openGtype <- lapply(X = seq_len(length.out = ncol(openGtype)), FUN = function(i){openGtype[,i]})

  #stitch the alleles together
  alleleType <- do.call(what = stringi::stri_join, args = openGtype)

  #get all possible combinations of allelles, always some rows by 2 columns
  openAllele <- gtools::combinations(n = length(alleleType), r = 2, v = alleleType, set = FALSE, repeats.allowed = TRUE)

  #final allele combinations
  genotypes <- do.call(what = stringi::stri_join, args = list(openAllele[,1], openAllele[,2]))

  #create a thing to see if this works
  #these work for H=2
  setNames(object = numeric(length(genotypes)), nm = genotypes)
  matrix(data = 0, nrow = length(genotypes), ncol = length(genotypes), dimnames = list(genotypes, genotypes))
  tMatrix <- array(data = 0, dim = c(length(genotypes),length(genotypes),length(genotypes)), dimnames = list(genotypes,genotypes,genotypes))


  ###############################################################################
  #setup all probability matrices

  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(H), dimnames = list(gtype, NULL))

  homingProbs[1, ] <- 1-d-H #chance to stay W is 1-homing-background mutation
  homingProbs[2, ] <- H*(1-R-S) #chance to become homing is H*1-H*R-H*S
  #homing probs, - NHEJ probs
  homingProbs[3, ] <- H*R   #NHEJ caused resistance, detrimentalt allele
  homingProbs[4, ] <- d+H*S #good resistant allele, from NHEJ and background mutation rate

  #set up lists to hold probabilities
  homProbsList <- vector(mode = "list", length = length(H))
  noHomProbsList <- vector(mode = "list", length = length(H))

  #fill the lists
  for(i in 1:length(H)){
    homProbsList[[i]]$W <- homingProbs[,i]
    homProbsList[[i]]$H <- setNames(object = c(1-d,d), nm = c("H", "S"))
    homProbsList[[i]]$R <- 1
    homProbsList[[i]]$S <- 1

    noHomProbsList[[i]]$W <- setNames(object = c(1-d,d), nm = c("W", "S"))
    noHomProbsList[[i]]$H <- setNames(object = c(1-d,d), nm = c("H", "S"))
    noHomProbsList[[i]]$R <- 1
    noHomProbsList[[i]]$S <- 1
  }

  ###############################################################################
  #fill transition matrix

  #check if homing allele is present
  hCheck <- grepl(pattern = "H", x = genotypes, ignore.case = FALSE, perl = FALSE, fixed = TRUE)
  nHom <- length(H)

  #loop over all matings, male and female
  for (fi in 1:length(genotypes)){
    for (mi in 1:length(genotypes)){

      #female and male genotype for this mating
      fSplit <- strsplit(genotypes[fi], "")[[1]]
      mSplit <- strsplit(genotypes[mi], "")[[1]]

      #split into parental alleles
      fGList <- list(allele1 = fSplit[1:nHom], allele2 = fSplit[(nHom+1):(2*nHom)] )
      mGList <- list(allele1 = mSplit[1:nHom], allele2 = mSplit[(nHom+1):(2*nHom)] )

      #set up allele and probability lists
      fAllele1 <- vector(mode = "list", length = nHom)
      fProbs1 <- vector(mode = "list", length = nHom)

      fAllele2 <- vector(mode = "list", length = nHom)
      fProbs2 <- vector(mode = "list", length = nHom)

      mAllele1 <- vector(mode = "list", length = nHom)
      mProbs1 <- vector(mode = "list", length = nHom)

      mAllele2 <- vector(mode = "list", length = nHom)
      mProbs2 <- vector(mode = "list", length = nHom)

      #fill female lists




      if(hCheck[fi]){
        #homing in females
        #lists of allele letter and probabilities
        for (i in 1:nHom){

          #female allele1
          if(fGList[[1]][i]=="W"){
            fAllele1[[i]] <- c("W", "H", "R", "S")
            fProbs1[[i]] <- homProbsList[[i]]$W
          } else if(fGList[[1]][i]=="H"){
            fAllele1[[i]] <- c("H", "S")
            fProbs1[[i]] <- homProbsList[[i]]$H
          } else if(fGList[[1]][i]=="R"){
            fAllele1[[i]] <- "R"
            fProbs1[[i]] <- homProbsList[[i]]$R
          } else if(fGList[[1]][i]=="S"){
            fAllele1[[i]] <- "S"
            fProbs1[[i]] <- homProbsList[[i]]$S
          }

          #female allele2
          if(fGList[[2]][i]=="W"){
            fAllele2[[i]] <- c("W", "H", "R", "S")
            fProbs2[[i]] <- homProbsList[[i]]$W
          } else if(fGList[[2]][i]=="H"){
            fAllele2[[i]] <- c("H", "S")
            fProbs2[[i]] <- homProbsList[[i]]$H
          } else if(fGList[[2]][i]=="R"){
            fAllele2[[i]] <- "R"
            fProbs2[[i]] <- homProbsList[[i]]$R
          } else if(fGList[[2]][i]=="S"){
            fAllele2[[i]] <- "S"
            fProbs2[[i]] <- homProbsList[[i]]$S
          }

        }

      } else {
        #no homing in females
        #lists of allele letter and probabilities
        for (i in 1:nHom){

          #female allele1
          if(fGList[[1]][i]=="W"){
            fAllele1[[i]] <- c("W", "S")
            fProbs1[[i]] <- noHomProbsList[[i]]$W
          } else if(fGList[[1]][i]=="H"){
            fAllele1[[i]] <- c("H", "S")
            fProbs1[[i]] <- noHomProbsList[[i]]$H
          } else if(fGList[[1]][i]=="R"){
            fAllele1[[i]] <- "R"
            fProbs1[[i]] <- noHomProbsList[[i]]$R
          } else if(fGList[[1]][i]=="S"){
            fAllele1[[i]] <- "S"
            fProbs1[[i]] <- noHomProbsList[[i]]$S
          }

          #female allele2
          if(fGList[[2]][i]=="W"){
            fAllele2[[i]] <- c("W", "S")
            fProbs2[[i]] <- noHomProbsList[[i]]$W
          } else if(fGList[[2]][i]=="H"){
            fAllele2[[i]] <- c("H", "S")
            fProbs2[[i]] <- noHomProbsList[[i]]$H
          } else if(fGList[[2]][i]=="R"){
            fAllele2[[i]] <- "R"
            fProbs2[[i]] <- noHomProbsList[[i]]$R
          } else if(fGList[[2]][i]=="S"){
            fAllele2[[i]] <- "S"
            fProbs2[[i]] <- noHomProbsList[[i]]$S
          }

        }

      }

      #fill male lists
      if(hCheck[mi]){
        #homing in males
        #lists om allele letter and probabilities
        for (i in 1:nHom){

          #male allele1
          if(mGList[[1]][i]=="W"){
            mAllele1[[i]] <- c("W", "H", "R", "S")
            mProbs1[[i]] <- homProbsList[[i]]$W
          } else if(mGList[[1]][i]=="H"){
            mAllele1[[i]] <- c("H", "S")
            mProbs1[[i]] <- homProbsList[[i]]$H
          } else if(mGList[[1]][i]=="R"){
            mAllele1[[i]] <- "R"
            mProbs1[[i]] <- homProbsList[[i]]$R
          } else if(mGList[[1]][i]=="S"){
            mAllele1[[i]] <- "S"
            mProbs1[[i]] <- homProbsList[[i]]$S
          }

          #male allele2
          if(mGList[[2]][i]=="W"){
            mAllele2[[i]] <- c("W", "H", "R", "S")
            mProbs2[[i]] <- homProbsList[[i]]$W
          } else if(mGList[[2]][i]=="H"){
            mAllele2[[i]] <- c("H", "S")
            mProbs2[[i]] <- homProbsList[[i]]$H
          } else if(mGList[[2]][i]=="R"){
            mAllele2[[i]] <- "R"
            mProbs2[[i]] <- homProbsList[[i]]$R
          } else if(mGList[[2]][i]=="S"){
            mAllele2[[i]] <- "S"
            mProbs2[[i]] <- homProbsList[[i]]$S
          }

        }

      } else {
        #no homing in males
        #lists om allele letter and probabilities
        for (i in 1:nHom){

          #male allele1
          if(mGList[[1]][i]=="W"){
            mAllele1[[i]] <- c("W", "S")
            mProbs1[[i]] <- noHomProbsList[[i]]$W
          } else if(mGList[[1]][i]=="H"){
            mAllele1[[i]] <- c("H", "S")
            mProbs1[[i]] <- noHomProbsList[[i]]$H
          } else if(mGList[[1]][i]=="R"){
            mAllele1[[i]] <- "R"
            mProbs1[[i]] <- noHomProbsList[[i]]$R
          } else if(mGList[[1]][i]=="S"){
            mAllele1[[i]] <- "S"
            mProbs1[[i]] <- noHomProbsList[[i]]$S
          }

          #male allele2
          if(mGList[[2]][i]=="W"){
            mAllele2[[i]] <- c("W", "S")
            mProbs2[[i]] <- noHomProbsList[[i]]$W
          } else if(mGList[[2]][i]=="H"){
            mAllele2[[i]] <- c("H", "S")
            mProbs2[[i]] <- noHomProbsList[[i]]$H
          } else if(mGList[[2]][i]=="R"){
            mAllele2[[i]] <- "R"
            mProbs2[[i]] <- noHomProbsList[[i]]$R
          } else if(mGList[[2]][i]=="S"){
            mAllele2[[i]] <- "S"
            mProbs2[[i]] <- noHomProbsList[[i]]$S
          }

        }

      }



      #get all combinations of those lists, separately letters and numbers, but in same order
      fKid1 <- expand.grid(fAllele1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      fPKid1 <- expand.grid(fProbs1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      fKid2 <- expand.grid(fAllele2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      fPKid2 <- expand.grid(fProbs2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      mKid1 <- expand.grid(mAllele1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      mPKid1 <- expand.grid(mProbs1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      mKid2 <- expand.grid(mAllele2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      mPKid2 <- expand.grid(mProbs2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      #add letters together, multiply probabilities
      fKid1Collapse <- do.call(what = stringi::stri_join, args = fKid1)
      fPKid1Collapse <- matrixStats::rowProds(x = as.matrix(fPKid1), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")
      fKid2Collapse <- do.call(what = stringi::stri_join, args = fKid2)
      fPKid2Collapse <- matrixStats::rowProds(x = as.matrix(fPKid2), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")

      mKid1Collapse <- do.call(what = stringi::stri_join, args = mKid1)
      mPKid1Collapse <- matrixStats::rowProds(x = as.matrix(mPKid1), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")
      mKid2Collapse <- do.call(what = stringi::stri_join, args = mKid2)
      mPKid2Collapse <- matrixStats::rowProds(x = as.matrix(mPKid2), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")

      #since these are 2 alleles in same parent, ofspring can't get both, so combine to 1 long list of all possible alleles
      # reorder allele labels, expand.grid calls them backwords from my definition
      fPoss <- stringi::stri_reverse(str = c(fKid1Collapse, fKid2Collapse))
      fProbsPoss <- c(fPKid1Collapse, fPKid2Collapse)

      mPoss <- stringi::stri_reverse(str = c(mKid1Collapse, mKid2Collapse))
      mProbsPoss <- c(mPKid1Collapse, mPKid2Collapse)

      #all combinations of male and female alleles
      offspringgenotypes <- expand.grid(fPoss, mPoss, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      offspringprobabilities <- expand.grid(fProbsPoss, mProbsPoss, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      #recombine into genotypes and their probabilities, 2 lists in same order.
      # collapse is 2 to 1 because expand.grid orders vectors backwords from combinations and permutations above
      offspringgenocollapse <- do.call(what = stringi::stri_join, args = list(offspringgenotypes[,2], offspringgenotypes[,1]))

      offspringprobcollapse <- matrixStats::rowProds(x = as.matrix(offspringprobabilities), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")/4



      #normalize probabilities
      #offspringprobcollapse <- offspringprobcollapse/4

      #reverse allele permutations that don't match allele combinations
      index <- which((offspringgenocollapse %in% genotypes)==FALSE)

      offspringgenocollapse[index] <- do.call(what = stringi::stri_join, args = list(offspringgenotypes[index,1], offspringgenotypes[index,2]))

      #set names on vector for better access
      nIndex <- match(x = offspringgenocollapse, table = genotypes)


     # offspringprobcollapse <- setNames(object = offspringprobcollapse, nm = offspringgenocollapse)



      #tMatrix[fi,mi, nIndex ] <- tMatrix[fi,mi, nIndex ] + offspringprobcollapse[nIndex]



      #store in matrix
      for (i in 1:length(offspringgenocollapse)){

        tMatrix[fi,mi, nIndex[i] ] <- tMatrix[fi,mi, nIndex[i] ] + offspringprobcollapse[i]

      }

    }
  }

  return(tMatrix)

}

###############################################################################
# Figure out how to analyze
###############################################################################
genotypes <- list(c("HH", "HR", "HS", "HW"), c("WW"), c("WW"))
genotypes <- list(c("HH"), c("WW"), c("WW"))
genotypes <- list(c("HH"), NULL, c("WW"))
collapse <- c(T,F,F)

AnalyzeOutput_mLoci_Daisy <- function(readDirectory, saveDirectory=NULL, genotypes, collapse){

  #get list of all files, unique runs, and unique patches
  dirFiles = list.files(path = readDirectory)
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
  for(i in 1:numLoci){
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
      fName = grep(pattern = paste("AF1", run, patch, sep = ".*"),
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
    fileName <- paste0(run, "_", paste0(gOI, collapse = "_"),
                       "_", format(x = Sys.Date(), "%Y%m%d"), ".rds")

    saveRDS(object = list(metaData=note, maleData=mArray, femaleData=fArray),
            file = file.path(saveDirectory,fileName),
            compress = "gzip")

  }#end run loop
}








#' Aggregate Female Output by Genotype
#'
#' Aggregate over male mate genotype to convert female matrix output into vector output.
#'
#' @param directory directory where output was written to; must not end in path seperator
#' @param genotypes character vector of possible genotypes; found in \code{driveCube$genotypesID}
#' @param remove boolean flag to remove original (unaggregated) file
#'
#' @export
aggregateFemales <- function(directory, genotypes, remove = FALSE){
  dirFiles = list.files(path = directory)

  femaleFiles = dirFiles[grep(pattern = "AF1",x = dirFiles)]

  for(file in femaleFiles){
    cat("processing ",file,"\n",sep="")
    thisPatch = regmatches(x = file, m = regexpr(pattern = "Patch[0-9]+", text = file))
    thisOutput = read.csv(file = file.path(directory, file))

    timePatch = thisOutput[,c(1,2)]
    thisOutput = thisOutput[,-c(1,2)]
    thisOutputNames = names(thisOutput)

    sizeGenotype = nchar(genotypes)[1]
    subNames = substr(thisOutputNames,1,sizeGenotype)

    aggregateOut = vapply(X = genotypes,FUN = function(x){
      cols = which(subNames==x)
      rowSums(thisOutput[,cols])
    },FUN.VALUE = numeric(nrow(thisOutput)))

    fileName = sub(pattern = "_",replacement = "_Aggregate_",x = file)
    cat("writing ",fileName,"\n",sep="")
    write.table(x = cbind(timePatch,aggregateOut),file=file.path(directory, fileName),
                row.names = FALSE,col.names = TRUE,sep = ",",quote = FALSE)
    if(remove){
      cat("removing ",file,"\n",sep="")
      file.remove(file.path(directory, file))
    }
  }
}

#' Retrieve Output
#'
#' Read in output from directory. The resulting object will be a nested list; outermost nesting dimension indexes runID, within runID elements are split by sex
#' and innermost nesting is over patches.
#'
#' @param directory directory where output was written to; must not end in path seperator
#' @param genotypes character vector of possible genotypes; found in \code{driveCube$genotypesID}
#'
#' @export
retrieveOutput <- function(directory, genotypes){
  dirFiles = list.files(path = directory)

  runID = strsplit(x = dirFiles,split = "_")
  runID = unique(x = grep( pattern = "Run", x = unlist(x = runID), value = TRUE))

  output = setNames(object = vector(mode = "list",length = length(runID)), runID)
  for(run in runID){
    cat("processing ",run,"\n",sep="")

    runFiles = dirFiles[grep(pattern = run,x = dirFiles)]
    patches = regmatches(x = runFiles, m = regexpr(pattern = "Patch[0-9]+", text = runFiles))
    patches = as.integer(gsub(pattern = "Patch",replacement = "",x = patches))
    patches = sort(unique(patches))

    output[[run]]$F = setNames(object = vector(mode = "list",length = length(patches)),
                               nm = paste0("Patch",patches))
    output[[run]]$M = setNames(object = vector(mode = "list",length = length(patches)),
                               nm = paste0("Patch",patches))

    # retrieve males for this run
    males = runFiles[grep(pattern = "ADM",x = runFiles)]
    for(male in males){
      thisPatch = regmatches(x = male, m = regexpr(pattern = "Patch[0-9]+", text = male))
      thisOutput = read.csv(file = file.path(directory, male))
      thisOutput = as.matrix(thisOutput)
      rownames(thisOutput) = NULL
      output[[run]]$M[[thisPatch]] = thisOutput[,-c(1,2)]
    }

    # retrieve females for this run
    females = runFiles[grep(pattern = "AF1_Aggregate",x = runFiles)]
    for(female in females){
      thisPatch = regmatches(x = female, m = regexpr(pattern = "Patch[0-9]+", text = female))
      thisOutput = read.csv(file = file.path(directory, female))
      thisOutput = thisOutput[,-c(1,2)]
      output[[run]]$F[[thisPatch]] = as.matrix(thisOutput)
    }

  }
  return(output)
}
