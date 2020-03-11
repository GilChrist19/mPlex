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

numPatch <- 2
set.seed(10)
migration <- matrix(data = runif(numPatch*numPatch), nrow = numPatch, ncol = numPatch)
migration <- migration/rowSums(migration)
migration <- diag(numPatch)

patchPops = rep(1000L,numPatch)

#setup alleles to initiate patches
alleloTypes <- vector(mode = "list", length = 1L) #1 locus
alleloTypes[[1]]$alleles <- c("W")
alleloTypes[[1]]$probs <- c(1L)

AllAlleles <- replicate(n = numPatch, expr = alleloTypes, simplify = FALSE)




# alleloTypes <- vector(mode = "list", length = 15L) #3 loci
# alleloTypes[[1]]$alleles <- c("W","R")
# alleloTypes[[1]]$probs <- c(1,2)
# alleloTypes[[2]]$alleles <- c("W","R")
# alleloTypes[[2]]$probs <- c(0,1)
# alleloTypes[[3]]$alleles <- c("W","H")
# alleloTypes[[3]]$probs <- c(1,0)
# alleloTypes[[4]]$alleles <- c("W","H")
# alleloTypes[[4]]$probs <- c(1,1)
# alleloTypes[[5]]$alleles <- c("W","H")
# alleloTypes[[5]]$probs <- c(1,0)
# alleloTypes[[6]]$alleles <- c("W","R")
# alleloTypes[[6]]$probs <- c(1,2)
# alleloTypes[[7]]$alleles <- c("W","R")
# alleloTypes[[7]]$probs <- c(0,1)
# alleloTypes[[8]]$alleles <- c("W","H")
# alleloTypes[[8]]$probs <- c(1,1)
# alleloTypes[[9]]$alleles <- c("W","H")
# alleloTypes[[9]]$probs <- c(1,2)
# alleloTypes[[10]]$alleles <- c("W","H")
# alleloTypes[[10]]$probs <- c(1,1)
# alleloTypes[[11]]$alleles <- c("W","R")
# alleloTypes[[11]]$probs <- c(1,2)
# alleloTypes[[12]]$alleles <- c("W","R")
# alleloTypes[[12]]$probs <- c(0,1)
# alleloTypes[[13]]$alleles <- c("W","H")
# alleloTypes[[13]]$probs <- c(1,0)
# alleloTypes[[14]]$alleles <- c("W","H")
# alleloTypes[[14]]$probs <- c(1,1)
# alleloTypes[[15]]$alleles <- c("W","H")
# alleloTypes[[15]]$probs <- c(1,0)
# alleloTypes[[16]]$alleles <- c("W","R")
# alleloTypes[[16]]$probs <- c(1,2)
# alleloTypes[[17]]$alleles <- c("W","R")
# alleloTypes[[17]]$probs <- c(0,1)
# alleloTypes[[18]]$alleles <- c("W","H")
# alleloTypes[[18]]$probs <- c(1,1)
# alleloTypes[[19]]$alleles <- c("W","H")
# alleloTypes[[19]]$probs <- c(1,2)
# alleloTypes[[20]]$alleles <- c("W","H")
# alleloTypes[[20]]$probs <- c(1,1)
# 
# AllAlleles <- replicate(n = numPatch, expr = alleloTypes, simplify = FALSE)




# reproductionReference <- MakeReference_DaisyDrive(cRateM = c(0),
#                                                   hRateM = c(0),
#                                                   rRateM = c(0),
#                                                   dM = c(0), eta = c("TIME"=4))
# 
# reproductionReference <- MakeReference_Multiplex_oLocus(cRateM = c(0.98),
#                                                         hRateM = c(0.1),
#                                                         rRateM = c(0.3),
#                                                         dM = c(0))

# This sets up a basic CRISPR drive, with perfect homing and no resistance or backgorund mutation
reproductionReference <- MakeReference_Multiplex_mLoci(cRateM = c(0),
                                                       hRateM = c(0.03),
                                                       rRateM = c(0.03),
                                                       cRateF = c(1),
                                                       dM = c(0))


###############################################################################
# Release Setup
###############################################################################

# create Release List
patchReleases = replicate(n = numPatch,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL,
                                      matedFemaleReleases = NULL),
                          simplify = FALSE)

# Create release object to pass to patches
# holdRel <- basicRepeatedReleases(releaseStart = 100L,
#                                          releaseEnd = 110L,
#                                          releaseInterval = 5,
#                                          genMos = c("HH"),
#                                          numMos = c(25L),
#                                          minAge = 16L,
#                                          maxAge = 24L,
#                                          ageDist = rep(x = 1, times = 24-16+1)/9)
# 
# 
# holdRel2 <- basicRepeatedReleases(releaseStart = 600L,
#                                           releaseEnd = 610L,
#                                           releaseInterval = 2L,
#                                           genMos = c("RR"),
#                                           numMos = c(10L),
#                                           minAge = 16L,
#                                           maxAge = 24L,
#                                           ageDist = rep(x = 1, times = 24-16+1)/9)
# 
# 
# patchReleases[[2]]$maleReleases <- c(holdRel, holdRel2)




patchReleases[[1]]$matedFemaleReleases <- basicRepeatedReleases(releaseStart = 10,
                                                                releaseEnd = 20,
                                                                releaseInterval = 1,
                                                                genMos = c("RR"),
                                                                genMosMate = c("HH"),
                                                                numMos = c(100L),
                                                                minAge = 16L,
                                                                maxAge = 24L,
                                                                ageDist = rep(x = 1, times = 24-16+1)/9)



patchReleases[[2]]$maleReleases <- basicRepeatedReleases(releaseStart = 10,
                                          releaseEnd = 20,
                                          releaseInterval = 1,
                                          genMos = c("RR"),
                                          numMos = c(100L),
                                          minAge = 16L,
                                          maxAge = 24L,
                                          ageDist = rep(x = 1, times = 24-16+1)/9)

###############################################################################
# Calculate parameters and initialize network
###############################################################################
simTime <- 100
netPar = NetworkParameters(nPatch = numPatch,
                           simTime = simTime,
                           sampTime = 2,
                           AdPopEQ = patchPops,
                           runID = 1L,
                           dayGrowthRate = 1.1,
                           beta = 32L, tEgg = 1, tLarva = 10, tPupa = 2)

migrationBatch <- basicBatchMigration(numPatches = numPatch)

startTime <- Sys.time()
runMPlex(seed = 10,
         numThreads = 1,
         numReps = 1, 
         networkParameters = netPar,
         reproductionReference = reproductionReference,
         initAlleles = AllAlleles,
         patchReleases = patchReleases,
         migrationMale = migration,
         migrationFemale = migration,
         migrationBatch = migrationBatch,
         outputDirectory = simDir,
         reproductionType = "mPlex_mLoci",
         verbose = TRUE)
difftime(time1 = Sys.time(), time2 = startTime)

# setup aggregation key
#  this example sets all genotypes as different
genOI_mLoci_Daisy(outputFile = file.path(aggDir, "0_AggKey.csv"), genotypes = list(NULL), collapse = c(FALSE))


# aggregate experiment by aggregation key
simAggregation(readDirectory = simDir, writeDirectory = aggDir, simTime = netPar$simTime, sampTime = netPar$sampTime)


# plot for example
plotmPlexSingle(directory = aggDir, whichPatches = NULL, totalPop = TRUE)
plotmPlexMult(directory = aggDir,whichPatches = NULL, totalPop = TRUE, nonZeroGen = FALSE)


# see profiling if done
cat("NumSamp  PercentSamp  CumPercentSamp  NumSampTree  PercentSampTree  Function")
system(sprintf("google-pprof --text --cum --lines  /bin/ls %sprofile.log", "~/Desktop/OUTPUT/"), intern = TRUE)


detach("package:mPlexCpp", unload=TRUE)
  


###############################################################################
# DON'T RUN THIS STUFF. SHOULD WORK, ISN'T ORGANIZED, MOSTLY FOR TESTING
###############################################################################

# repetitions wrapper - no reinitializing memory between reps. 
runMPlex(seed = 10,
         numReps = 2, 
         numThreads = 2,
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





