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

numPatch <- 4
set.seed(10)
# migration <- matrix(data = runif(numPatch*numPatch), nrow = numPatch, ncol = numPatch)
# migration <- migration/rowSums(migration)
migration <- diag(nrow = numPatch)


patchPops = rep(1000L,numPatch)

#setup alleles to initiate patches
reference <- list('eta'=numeric(0), 'phi'=numeric(0), 'omega'=numeric(0),
                  'xiF'=numeric(0), 'xiM'=numeric(0), 's'=numeric(0))


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
                           AdPopEQ = patchPops,
                           runID = 1L,
                           dayGrowthRate = 1.1,
                           beta = 32L, tEgg = 1, tLarva = 10, tPupa = 1)

migrationBatch <- basicBatchMigration(numPatches = numPatch)

samplingScheme <- list("samplingDays"=c(5,5,5,5,5),
                       "samplingCoverage"=c(0.2,0.1,0.1,0.1,0.1))


startTime <- Sys.time()
runMPlex(seed = 10,
         numThreads = 4,
         networkParameters = netPar,
         reproductionReference = reference,
         patchReleases = patchReleases,
         migrationMale = migration,
         migrationFemale = migration,
         migrationBatch = migrationBatch,
         samplingParameters = samplingScheme,
         outputDirectory = simDir,
         verbose = TRUE)

endTime <- Sys.time()
print(difftime(time1 = endTime, time2 = startTime))



# test read things in simDir
files <- list.files(path = simDir, full.names = TRUE, pattern = "E_")
allFiles <- list()
for(i in 1:length(files)){
  allFiles[[i]] <- as.matrix(read.csv(file = files[i], header = TRUE))
}

myID <- lapply(X = allFiles, FUN = function(x){x[ ,3]})

#myID[[1]] %in% myID[[2]]

for( i in 1:3){
  print(sum(myID[[i]] %in% myID[[i+1]]))
  print(sum(myID[[i+1]] %in% myID[[i]]))
  cat("\n")
}







# # setup aggregation key
# #  this example sets all genotypes as different
# genOI_mLoci_Daisy(outputFile = file.path(aggDir, "0_AggKey.csv"), genotypes = list(NULL), collapse = c(FALSE))
# 
# 
# # aggregate experiment by aggregation key
# SimAggregation(readDirectory = simDir, writeDirectory = aggDir, simTime = simTime)
# 
# 
# # plot for example
# Plot_mPlex(directory = aggDir, whichPatches = NULL, totalPop = TRUE)


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
         numThreads = 1,
         networkParameters = netPar,
         reproductionReference = reference,
         patchReleases = patchReleases,
         migrationMale = migration,
         migrationFemale = migration,
         migrationBatch = migrationBatch,
         samplingParameters = samplingScheme,
         outputDirectory = simDir,
         verbose = TRUE)



# 
# # setup aggregation key
# #  this example sets all genotypes as different
# genOI_mLoci_Daisy(outputFile = file.path(aggDir, "0_AggKey.csv"), genotypes = list(NULL), collapse = c(FALSE))
# # example for oLocus
# #genOI_oLocus(outputFile = file.path(aggDir, "0_AggKey.csv"), alleles = list(list(c(NULL)),list(c(NULL))), collapse = list(c(F),c(F)))
# 
# 
# # aggregate experiment by aggregation key
# SimAggregation(readDirectory = simDir, writeDirectory = aggDir, simTime = simTime)
# 
# 
# # plot for example
# Plot_mPlex(directory = aggDir, whichPatches = NULL, totalPop = TRUE)



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





