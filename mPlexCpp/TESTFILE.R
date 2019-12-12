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
# Setup Directories
###############################################################################
topDirectory <- "~/Desktop/OUTPUT/mPlex"
simDir = "~/Desktop/OUTPUT/mPlex/simDir"
aggDir <- "~/Desktop/OUTPUT/mPlex/aggDir"

if(!dir.exists(paths = topDirectory)){
  dir.create(path = topDirectory)
} else {
    unlink(x = list.files(topDirectory, full.names = TRUE), recursive = TRUE, force = TRUE)
}
for(i in c(simDir, aggDir)){dir.create(path = i)}


###############################################################################
# Setup Parameters for Network
###############################################################################

numPatch <- 5
set.seed(10)
migration <- matrix(data = runif(numPatch*numPatch), nrow = numPatch, ncol = numPatch)
migration <- migration/rowSums(migration)

patchPops = rep(1000L,numPatch)

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
simTime <- 1000
netPar = NetworkParameters(nPatch = numPatch,
                           simTime = simTime,
                           alleloTypes = AllAlleles,
                           AdPopEQ = patchPops,
                           runID = 1L,
                           dayGrowthRate = 1.1,
                           beta = 32L, tEgg = 1, tLarva = 10, tPupa = 2)

migrationBatch <- basicBatchMigration(numPatches = numPatch)

mPlex_oneRun(seed = 10,
             networkParameters = netPar,
             reproductionReference = reproductionReference,
             patchReleases = patchReleases,
             migrationMale = migration,
             migrationFemale = migration,
             migrationBatch = migrationBatch,
             outputDirectory = simDir,
             reproductionType = "mPlex_mLoci",
             verbose = TRUE)


# setup aggregation key
#  this example sets all genotypes as different
genOI_mLoci_Daisy(outputFile = file.path(aggDir, "0_AggKey.csv"), genotypes = list(NULL), collapse = c(FALSE))


# aggregate experiment by aggregation key
SimAggregation(readDirectory = simDir, writeDirectory = aggDir, simTime = simTime)


# plot for example
Plot_mPlex(directory = aggDir, whichPatches = NULL, totalPop = TRUE)


# see profiling if done
cat("NumSamp  PercentSamp  CumPercentSamp  NumSampTree  PercentSampTree  Function")
system(sprintf("google-pprof --text --cum --lines  /bin/ls %sprofile.log", "~/Desktop/OUTPUT/"), intern = TRUE)


detach("package:mPlexCpp", unload=TRUE)
  


###############################################################################
# DON'T RUN THIS STUFF. SHOULD WORK, ISN'T ORGANIZED, MOSTLY FOR TESTING
###############################################################################

# repetitions wrapper - no reinitializing memory between reps. 
mPlex_runRepetitions(seed = 10,
                     numReps = 5, 
                     networkParameters = netPar,
                     reproductionReference = reproductionReference,
                     patchReleases = patchReleases,
                     migrationMale = migration,
                     migrationFemale = migration,
                     migrationBatch = migrationBatch,
                     outputDirectory = simDir,
                     reproductionType = "mPlex_mLoci",
                     verbose = FALSE)




# setup aggregation key
#  this example sets all genotypes as different
genOI_mLoci_Daisy(outputFile = file.path(aggDir, "0_AggKey.csv"), genotypes = list(NULL), collapse = c(FALSE))
# example for oLocus
#genOI_oLocus(outputFile = file.path(aggDir, "0_AggKey.csv"), alleles = list(list(c(NULL)),list(c(NULL))), collapse = list(c(F),c(F)))


# aggregate experiment by aggregation key
SimAggregation(readDirectory = simDir, writeDirectory = aggDir, simTime = simTime)


# plot for example
Plot_mPlex(directory = aggDir, whichPatches = NULL, totalPop = TRUE)



detach("package:mPlexCpp", unload=TRUE)



###############################################################################
# SCRATCH SPACE
###############################################################################
# Use this for profiling.
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

  # c++ for analysis
  mPlexCpp:::simAgg(readFiles_ = readFiles, writeFiles_ = writeFiles,
                    largeFile_ = largeFile, simTime_ = simTime, genKey_ = genKey)
  
}





