###############################################################################
#                            ____  __          ______          
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
#                                               /_/   /_/      
###############################################################################
###############################################################################
#                        _____         _   _____ _ _
#                       |_   _|__  ___| |_|  ___(_) | ___
#                         | |/ _ \/ __| __| |_  | | |/ _ \
#                         | |  __/\__ \ |_|  _| | | |  __/
#                         |_|\___||___/\__|_|   |_|_|\___|
#
###############################################################################
###############################################################################
# Clean environment and source files
###############################################################################
rm(list=ls());gc()
library(mPlexCpp)






###############################################################################
# Setup Parameters for Network
###############################################################################

numPatch <- 5
set.seed(10)
migration <- matrix(data = runif(numPatch*numPatch), nrow = numPatch, ncol = numPatch)
migration <- migration/rowSums(migration)

patchPops = rep(1000L,numPatch)

directory1 = "~/Desktop/OUTPUT/mPlex/"
directory2 <- "~/Desktop/OUTPUT/MGDrivEHOLD/"

#setup alleles to initiate patches
alleloTypes <- vector(mode = "list", length = 1L) #1 locus
alleloTypes[[1]]$alleles <- c("W")
alleloTypes[[1]]$probs <- c(1L)

AllAlleles <- replicate(n = numPatch, expr = alleloTypes, simplify = FALSE)


# reproductionReference <- MakeReference_DaisyDrive(H = c(0.98, 0.5),
#                                                   R = c(0.0001, 0.0001),
#                                                   S = c(0.0003, 0.004),
#                                                   d = c(0, 0), eta = c("TIME"=4))

# reproductionReference <- MakeReference_Multiplex_oLocus(H = c(0.98, 0.5),
#                                                          R = c(0.0001, 0.0001),
#                                                          S = c(0.0003, 0.004),
#                                                          d = c(0, 0))

# This sets up a basic CRISPR drive, with perfect homing and no resistance or backgorund mutation
reproductionReference <- MakeReference_Multiplex_mLoci(H = c(0.97),
                                                       R = c(0.03),
                                                       S = c(0),
                                                       d = c(0.00))


###############################################################################
# Release Setup
###############################################################################

# create Release List
patchReleases = replicate(n = numPatch,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL),
                          simplify = FALSE)


# Create release object to pass to patches
holdRel <- Release_basicRepeatedReleases(releaseStart = 100L,
                                         releaseEnd = 110L,
                                         releaseInterval = 1,
                                         genMos = c("HH"),
                                         numMos = c(25L),
                                         minAge = 16L,
                                         maxAge = 24L,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)


holdRel2 <- Release_basicRepeatedReleases(releaseStart = 600L,
                                          releaseEnd = 610L,
                                          releaseInterval = 2L,
                                          genMos = c("RR"),
                                          numMos = c(10L),
                                          minAge = 16L,
                                          maxAge = 24L,
                                          ageDist = rep(x = 1, times = 24-16+1)/9)


patchReleases[[1]]$maleReleases <- c(holdRel, holdRel2)


###############################################################################
# Calculate parameters and initialize network
###############################################################################
netPar = NetworkParameters(nPatch = numPatch,
                           simTime = 1000,
                           alleloTypes = AllAlleles,
                           AdPopEQ = patchPops,
                           runID = 1L,
                           dayGrowthRate = 1.1,
                           beta = 32L, tEgg = 1, tLarva = 10, tPupa = 1)

migrationBatch <- basicBatchMigration(numPatches = numPatch)

sTime <- Sys.time()
for(i in 1:5){
  
  
  mPlex_oneRun(seed = 10,
               networkParameters = netPar,
               reproductionReference = reproductionReference,
               patchReleases = patchReleases,
               migrationMale = migration,
               migrationFemale = migration,
               migrationBatch = migrationBatch,
               outputDirectory = directory1,
               reproductionType = "mPlex_mLoci",
               verbose = TRUE)
  

}
eTime <- Sys.time()
print(x = difftime(time1 = eTime, time2 = sTime))

# aggregate by genotype.
AnalyzeOutput_mLoci_Daisy(readDirectory = directory1,
                          saveDirectory = directory2,
                          genotypes = list(NULL),
                          collapse = c(FALSE),
                          numCores = 1)

# plot for example
Plot_mPlex(directory = directory2, whichPatches = NULL, totalPop = TRUE)

# see profiling if done

cat("NumSamp  PercentSamp  CumPercentSamp  NumSampTree  PercentSampTree  Function")
system(sprintf("google-pprof --text --cum --lines  /bin/ls %sprofile.log", "~/Desktop/OUTPUT/"), intern = TRUE)




detach("package:mPlexCpp", unload=TRUE)
  


###############################################################################
# DON'T RUN THIS STUFF. SHOULD WORK, ISN'T ORGANIZED, MOSTLY FOR TESTING
###############################################################################

# repetitions wrapper - no reinitializing memory between reps. 
sTime <- Sys.time()
mPlex_runRepetitions(seed = 10,
                     numReps = 5, 
                     networkParameters = netPar,
                     reproductionReference = reproductionReference,
                     patchReleases = patchReleases,
                     migrationMale = migration,
                     migrationFemale = migration,
                     migrationBatch = migrationBatch,
                     outputDirectory = directory1,
                     reproductionType = "mPlex_mLoci",
                     verbose = FALSE)
eTime <- Sys.time()
print(x = difftime(time1 = eTime, time2 = sTime))

#splitOutput(directory = directory, numCores = 1)
AnalyzeOutput_oLocus(readDirectory = "~/Desktop/HOLD/MGDrivE (copy 1)/",
                     saveDirectory = "~/Desktop/HOLD/MGDrivE (copy 1) (copy 1)/",
                     alleles = list(list(c(NULL)),list(c(NULL))),
                     collapse = list(c(F),c(F)),
                     numCores = 1)

AnalyzeOutput_mLoci_Daisy(readDirectory = "~/Desktop/HOLD/MGDrivE/",
                          saveDirectory = "~/Desktop/HOLD/mPlex/",
                          genotypes = list(NULL),
                          collapse = c(FALSE),
                          numCores = 1)


Plot_mPlex(directory = "~/Desktop/HOLD/MGDrivE (copy 1) (copy 1)/", whichPatches = NULL, totalPop = TRUE)

detach("package:mPlexCpp", unload=TRUE)




###############################################################################
# write split function
###############################################################################
# unsure what was here. It didn't make sense, so I deleted it. 
#  probably just scratch space anyway

###############################################################################
# When profiling, use this to see it
###############################################################################

cat("NumSamp  PercentSamp  CumPercentSamp  NumSampTree  PercentSampTree  Function")
system(sprintf("google-pprof --text --cum --lines  /bin/ls %sprofile.log", directory_det_c), intern = TRUE)



cat("NumSamp  PercentSamp  CumPercentSamp  NumSampTree  PercentSampTree  Function")
system(sprintf("google-pprof --text --cum --lines  /bin/ls %sprofile.log", "~/Desktop/OUTPUT/"), intern = TRUE)


mPlexCpp::genOI_mLoci_Daisy(outputFile = "~/Desktop/OUTPUT/mPlex/Aggregate/0_AggKey.csv", genotypes = list(NULL), collapse = c(FALSE))



readDirectory <- "~/Desktop/OUTPUT/mPlex/experimentTest"
writeDirectory <- "~/Desktop/OUTPUT/mPlex/Aggregate"
simTime=1000

testFunc <- function(readDirectory, writeDirectory, simTime){
  
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

  
  
  
  microbenchmark::microbenchmark(  mPlexCpp:::testFunc(readFiles_ = readFiles, writeFiles_ = writeFiles, largeFile_ = largeFile, simTime_ = simTime, genKey_ = genKey),
                                   times = 10)
  
  
  
  
  mPlexCpp:::testFunc(readFiles_ = readFiles, writeFiles_ = writeFiles, largeFile_ = largeFile, simTime_ = simTime, genKey_ = genKey)
  mPlexCpp:::testFunc2(readFiles_ = readFiles, writeFiles_ = writeFiles, largeFile_ = largeFile, simTime_ = simTime, genKey_ = genKey)
  
  
  
  
  
  Unit: seconds
  expr
  mPlexCpp:::testFunc(readFiles_ = readFiles, writeFiles_ = writeFiles,      largeFile_ = largeFile, simTime_ = simTime, genKey_ = genKey)
  mPlexCpp:::testFunc2(readFiles_ = readFiles, writeFiles_ = writeFiles,      largeFile_ = largeFile, simTime_ = simTime, genKey_ = genKey)
  min       lq     mean   median       uq      max neval
  8.803576 8.891493 8.951752 8.979497 9.021589 9.077743    10
  8.866992 8.900587 8.948064 8.953239 8.987247 9.054045    10
  

  
  
  
}



genKey <- read.csv(file = "~/Desktop/testFile1.csv", header = TRUE, stringsAsFactors = FALSE)
mPlexCpp:::testFunc(genKey_ = genKey)


genotypes <- list(c("HH","HW"),NULL, c("HH","HW"))
collapse <- c(FALSE, TRUE, FALSE)
outputFile <- "~/Desktop/testFile1.csv"





alleles = list(list(c(NULL),c("H")),list(c(NULL),c("W","R","B")))
collapse = list(c(F,T),c(F,T))





genOI_oLocus <- function(outputFile, genotypes, collapse){
  
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
  
  # set cores for data.table
  #  There is lots of internal multithreading now, tis a problem.
  oCores <- data.table::setDTthreads(threads = numCores, restore_after_fork = FALSE)
  
  # Check save directory
  if(saveDirectory == readDirectory){
    stop("Please create a save directory in a new directory.")
  }
  
  # get list of all files, then unique patches
  dirFiles = list.files(path = readDirectory, pattern = ".*\\.csv$")
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch_[0-9]+", text = dirFiles)))
  
  # import one file:get simTime, check genotypes for safety checks
  testFile <- data.table::fread(input = file.path(readDirectory, dirFiles[1]), sep = ",",
                                header = TRUE, verbose = FALSE, showProgress = FALSE,
                                data.table = TRUE, stringsAsFactors = FALSE,
                                logical01 = FALSE, drop = c("Age","Mate"))
  
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
                              data.table = TRUE, stringsAsFactors = FALSE,
                              logical01 = FALSE, drop = c("Age","Mate"), key = "Time")
    
    # male file, drop patch number, key the time
    mFile = data.table::fread(input = file.path(readDirectory, names[2]), sep = ",",
                              header = TRUE, verbose = FALSE, showProgress = FALSE,
                              data.table = TRUE, stringsAsFactors = FALSE,
                              logical01 = FALSE, drop = "Age", key = "Time")
    
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
  
  # reset number of cores (if set for some other reason)
  setDTthreads(threads = oCores)
  
  } # end function














