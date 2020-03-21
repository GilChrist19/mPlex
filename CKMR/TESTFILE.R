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
runCKMR(seed = 10,
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









