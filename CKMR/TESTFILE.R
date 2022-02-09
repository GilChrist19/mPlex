###############################################################################
#     ________ __ __  _______ 
#    / ____/ //_//  |/  / __ \
#   / /   / ,<  / /|_/ / /_/ /
#  / /___/ /| |/ /  / / _, _/ 
#  \____/_/ |_/_/  /_/_/ |_|  
#    
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
library(CKMR)



###############################################################################
# Setup Directories
###############################################################################
topDirectory <- "~/Desktop/OUTPUT/mPlex"
simDir <- file.path(topDirectory,"simDir")
aggDir <- file.path(topDirectory,"aggDir")

if(!dir.exists(paths = topDirectory)){
  dir.create(path = topDirectory)
} else {
    unlink(x = list.files(topDirectory, full.names = TRUE), recursive = TRUE, force = TRUE)
}
for(i in c(simDir, aggDir)){dir.create(path = i)}


###############################################################################
# Setup Parameters for Network
###############################################################################
simTime <- 100
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
holdRel <- basicRepeatedReleases(releaseStart = 100L,
                                         releaseEnd = 110L,
                                         releaseInterval = 1,
                                         genMos = c("HH"),
                                         numMos = c(25L),
                                         minAge = 16L,
                                         maxAge = 24L,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)


holdRel2 <- basicRepeatedReleases(releaseStart = 600L,
                                          releaseEnd = 610L,
                                          releaseInterval = 2L,
                                          genMos = c("RR"),
                                          numMos = c(10L),
                                          minAge = 16L,
                                          maxAge = 24L,
                                          ageDist = rep(x = 1, times = 24-16+1)/9)


patchReleases[[1]]$maleReleases <- c(holdRel, holdRel2)





###############################################################################
# Sampling Setup
###############################################################################
# sampling is different for every patch, every life stage, and every time point
# final object is a list with 2 arrays in it, one for when to sample, one for 
#  how many to sample

# sampling time
#  Boolean array: simTime x lifeStages x numPatches
#  lifeStages is always 5 - egg, larva, pupa, male, female
# Basic, every patch has every life stage (5 of them) sampled every day
sampDay <- array(data = TRUE, dim = c(simTime, 5, numPatch))

# test, don't sample any life stage in patch 1
#  patch 2 gets forgotten half way through
sampDay[ , ,1] <- FALSE
sampDay[50:simTime, ,2] <- FALSE

# sampling coverage
#  double array: simTime x lifeStages x numPatches
#  lifeStages is always 5 - egg, larva, pupa, male, female
# Example, every patch and every stage at 10% every day
sampCov <- array(data = 0.1, dim = c(simTime, 5, numPatch))

# put in list for simulation
#  the names in this list are fixed, do not change them
samplingScheme <- list("samplingDays"=sampDay, "samplingCoverage"=sampCov)


###############################################################################
# Calculate parameters and initialize network
###############################################################################
netPar = NetworkParameters(nPatch = numPatch,
                           simTime = simTime,
                           AdPopEQ = patchPops,
                           runID = 1L,
                           dayGrowthRate = 1.1,
                           beta = 32L, tEgg = 1, tLarva = 10, tPupa = 1)

migrationBatch <- basicBatchMigration(numPatches = numPatch)



startTime <- Sys.time()
runCKMR(seed = 10,
         numThreads = 2,
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






# see profiling if done
cat("NumSamp  PercentSamp  CumPercentSamp  NumSampTree  PercentSampTree  Function")
system(sprintf("google-pprof --text --cum --lines  /bin/ls %sprofile.log", "~/Desktop/OUTPUT/"), intern = TRUE)


detach("package:CKMR", unload=TRUE)
  


###############################################################################
# DON'T RUN THIS STUFF. SHOULD WORK, ISN'T ORGANIZED, MOSTLY FOR TESTING
###############################################################################

# repetitions wrapper - no reinitializing memory between reps. 
runCKMR(seed = 10,
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





detach("package:CKMR", unload=TRUE)



###############################################################################
# SCRATCH SPACE
###############################################################################
# Use this for profiling.
cat("NumSamp  PercentSamp  CumPercentSamp  NumSampTree  PercentSampTree  Function")
system(sprintf("google-pprof --text --cum --lines  /bin/ls %sprofile.log", "~/Desktop/OUTPUT/"), intern = TRUE)









