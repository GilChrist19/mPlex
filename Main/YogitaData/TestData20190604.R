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
# 20190618
#  Updated for new data for Yogita
#  Sampling everyday now
#  Test migration at 10% and 20%
#  Keeping patches at 2, she asked for more.
#
# 20190916
#  New test data for Yogita
#  10 x 10 patches landscape
#  migration rate of 10% daily, split equally over all other patches
#  16 adults per patch (8 male, 8 female)
#  10,000 day simulation time
#  Sampling
#    mother and larvae only
#    daily sampling
#    10% of population at each node
#    SPARSE OVER NODES - unsure if I can do this, gonna see
#  1 repetition\
#
# 20190917
#  added a bunch of shit to sample randomly over populations
#
# 20190930
#  30x30 patch landscape - 900 nodes total
#  Rest of the params are the same
#
# 20200219
#  new landscape, setup using 72% lifetime stay probs, exp. rate from MGDrivE, 
#  and 16.6062 meters between patches. 
#  Update network parameters to match threshold stuff: 2/5/1 stage times
#  Run on 3 pop sizes: 6,12,18
#
# 20200302
#  Same landscape and parameters as before.
#  Run on 2 new pop sizes: 24, 30
#
# 20200312
#  Same landscape
#  Run on 12,18,24,30 pop sizes
#  Sample adult males and larval offspring for a set
#  Sample only larvae for a set
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

# numPatch <- 16 #900
# #set.seed(10)
# migration <- matrix(data = 0.10/(numPatch-1), nrow = numPatch, ncol = numPatch)
# diag(x = migration) <- 1-0.10

migration <- as.matrix(x = read.csv(file = "~/Downloads/20200312_FO/gridMoveMat.csv",header = FALSE))
numPatch <- nrow(migration)


patchPops = rep(12,numPatch)

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
# holdRel <- Release_basicRepeatedReleases(releaseStart = 100L,
#                                          releaseEnd = 110L,
#                                          releaseInterval = 1,
#                                          genMos = c("HH"),
#                                          numMos = c(25L),
#                                          minAge = 16L,
#                                          maxAge = 24L,
#                                          ageDist = rep(x = 1, times = 24-16+1)/9)
#
#
# holdRel2 <- Release_basicRepeatedReleases(releaseStart = 600L,
#                                           releaseEnd = 610L,
#                                           releaseInterval = 2L,
#                                           genMos = c("RR"),
#                                           numMos = c(10L),
#                                           minAge = 16L,
#                                           maxAge = 24L,
#                                           ageDist = rep(x = 1, times = 24-16+1)/9)
#
#
# patchReleases[[1]]$maleReleases <- c(holdRel, holdRel2)


###############################################################################
# Calculate parameters and initialize network
###############################################################################
simTime <- 2000
netPar = NetworkParameters(nPatch = numPatch,
                           simTime = simTime,
                           AdPopEQ = patchPops,
                           runID = 1L,
                           dayGrowthRate = 1.175,
                           beta = 20, tEgg = 2, tLarva = 5, tPupa = 1,
                           muAd = 0.09)

migrationBatch <- basicBatchMigration(batchProbs = 0, numPatches = numPatch)

# sample larvae and females, 10% of the populations
samplingScheme <- list("samplingDays"=c(1,1,1,1,1),
                       "samplingCoverage"=c(0,0.1,0,0,0))


startTimes <- Sys.time()
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
         verbose = FALSE)
difftime(Sys.time(), startTimes)

"2 cores"
12
  8.07min
18
  8.54min
24
  8.87min
30
  9.24min

"4 cores"
12
  6.3min
18
  6.68min
24
  6.98min
30
  7.33min

detach("package:mPlexCpp", unload=TRUE)



####################
# Sparse sample the patches
####################
# Yogita wants sparse sampling over patches
# I don't innately do that, so, here I figure out how to do that
# Let's sample 10% of the patches at each time step. 

set.seed(10)

# 1 remove empty files
unlink(x = list.files(path = simDir, pattern = "(E|P|M)_", full.names = TRUE),
       recursive = TRUE, force = TRUE)

# 2 read in all files
fFiles <- list.files(path = simDir, pattern = "F_", full.names = TRUE)
lFiles <- list.files(path = simDir, pattern = "L_", full.names = TRUE)

expList <- list('fList'=list(),'lList'=list())

# loop over patches and read it all in
for(nPatch in 1:numPatch){
  # read in and store female
  expList$fList[[nPatch]] <- matrix(data = scan(file = fFiles[nPatch], what = integer(),sep = ",", skip = 1, quiet = TRUE),
                          ncol = 6, byrow = TRUE, dimnames = list(NULL,c('Time','Age','myID','momID','dadID','Mate'))
                          )
  
  # read in and store larvae
  expList$lList[[nPatch]] <- matrix(data = scan(file = lFiles[nPatch], what = integer(),sep = ",", skip = 1, quiet = TRUE),
                            ncol = 5, byrow = TRUE, dimnames = list(NULL,c('Time','Age','myID','momID','dadID'))
                            )
}

# sparse sampling and setup output
#  lets sample 10% of the populations every day - this is exact sampling for right now
#  output in 2 files, one for females and one for larvae
nSamp <- 3

retMatList <- list('femMat' = matrix(data = 0, nrow = 0, ncol = 7,
                                     dimnames = list(NULL,c('Patch','Time','Age','myID','momID','dadID','Mate'))),
                   'larvaeMat' = matrix(data = 0, nrow = 0, ncol = 6,
                                        dimnames = list(NULL,c('Patch','Time','Age','myID','momID','dadID')))
                   )

# loop over every day
for(dayTime in 1:simTime){
  # draw a sample of the patches to use for today
  #  patches are numbered 0-(nPatch - 1), but R is 1 indexed
  whichPatches <- sort.int(x = sample(x = 1:numPatch, size = nSamp, replace = FALSE))
  
  ##########
  # subset for sampling
  ##########
  # run over larvae and females
  for(stage in 1:2){
    # subset total list, then for each sample
    #  get index of time matching current time, return those rows
    # DON'T DROP DIMENSIONS - NEED LATER
    sampledObs <- lapply(X = expList[[stage]][whichPatches], function(x){
      x[which(x[ ,"Time"] == dayTime), ,drop=FALSE]
    })
    
    # get the number of rows in each list element, this corresponds to the sample for that day
    #  unlist it, use it to replicate the patch that it came from
    #  now the path labels are the same lenght as the number of observations
    # This is because observations were stochastic, so there could be many or none 
    #  on a given day. 
    patchSampled <- rep(x = whichPatches, unlist(lapply(X = sampledObs, FUN = NROW)))
    
    # bind all the observations into 1 object
    obsMat <- Reduce(f = rbind, x = sampledObs)
    
    # bind the patch label to the sample observations
    retMatList[[stage]] <- rbind(retMatList[[stage]], cbind("Patch"=patchSampled, obsMat))
  } # end loop over stage
} # end loop over day

# save output
write.table(x = retMatList$femMat, file = file.path(topDirectory,"femaleObservations.csv"),
            sep = ",", row.names = FALSE, col.names = TRUE)
write.table(x = retMatList$larvaeMat, file = file.path(topDirectory,"larvaeObservations.csv"),
            sep = ",", row.names = FALSE, col.names = TRUE)




###############################################################################
# DON'T RUN THIS STUFF. SHOULD WORK, ISN'T ORGANIZED, MOSTLY FOR TESTING
###############################################################################

# repetitions wrapper - no reinitializing memory between reps.
mPlex_runRepetitions(seed = 10,
                     numReps = 5,
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

detach("package:mPlexCpp", unload=TRUE)



###############################################################################
# SCRATCH SPACE
###############################################################################
# Use this for profiling.
cat("NumSamp  PercentSamp  CumPercentSamp  NumSampTree  PercentSampTree  Function")
system(sprintf("google-pprof --text --cum --lines  /bin/ls %sprofile.log", "~/Desktop/OUTPUT/"), intern = TRUE)



