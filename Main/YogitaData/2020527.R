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
# 20200421
#  25x25 landscape
#  16.6062m distance between nodes on the grid
#  72% lifetime stay probability
#  exp. rate from MGDrivE
#  Expected Displacement: 30.978m
#  Expected Cumulative Distance: 34.394m
#  Mosquitoes per patch: 40
#  Sampling Time: daily
#  Sampling Places: 130 133 136 139 142 145 205 208 211 214 217 220 280 283 286 
#                   289 292 295 355 358 361 364 367 370 430 433 436 439 442 445
#                   505 508 511 514 517 520
#  Sampling Fraction: 35/7% of female pop, 23/7% of male pop, no aquatic
# then, sampling Fraction: 35% of female pop, 23% of male pop, no aquatic, once a week





###############################################################################
# Clean environment and source files
###############################################################################
rm(list=ls());gc()
library(CKMR)



###############################################################################
# Setup Directories
###############################################################################
topDirectory <- "~/Desktop/OUTPUT/mPlex/CKMR"
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

migration <- as.matrix(x = read.csv(file = "~/Desktop/OUTPUT/mPlex/gridMoveMat.csv",header = FALSE))
numPatch <- nrow(migration)

simTime <- 2100
patchPops = rep(40,numPatch)

#setup alleles to initiate patches
reference <- list('eta'=numeric(0), 'phi'=numeric(0), 'omega'=numeric(0),
                  'xiF'=numeric(0), 'xiM'=numeric(0), 's'=numeric(0))


###############################################################################
# Release Setup
###############################################################################

# create Release List
#  there are no releases, this is just null
patchReleases = replicate(n = numPatch,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL),
                          simplify = FALSE)


###############################################################################
# Sampling Setup
###############################################################################
# sampling is different for every patch, and every life stage
# final object is a list with 2 matrices in it, one for when to sample, one for 
#  how many to sample

# Basic, every patch has every life stage (5 of them), never sampled
#  to never sample, must put time to simTime + 1, otherwise modulo throws errors
sampDay <- matrix(data = simTime + 1, nrow = numPatch, ncol = 5)

# coverage
#  all aquatic stages get 0% coverage
#  femalesget 35/7% coverage
#  males get 23/7% coverage
sampCov <- matrix(data = c(0,0,0,0.23,0.35), nrow = numPatch, ncol = 5, byrow = TRUE)

# set specific patches to get sampled every week
#  only sample male and female
sampPlaces <- c(130,133,136,139,142,145,205,208,211,214,217,220,280,283,286,289,
                292,295,355,358,361,364,367,370,430,433,436,439,442,445,505,508,
                511,514,517,520)

sampDay[sampPlaces,4:5] <- 7


# list to pass to mPlex
samplingScheme <- list("samplingDays"=sampDay, "samplingCoverage"=sampCov)


###############################################################################
# Calculate parameters and initialize network
###############################################################################

netPar = NetworkParameters(nPatch = numPatch,
                           simTime = simTime,
                           AdPopEQ = patchPops,
                           runID = 1L,
                           dayGrowthRate = 1.175,
                           beta = 20, tEgg = 2, tLarva = 5, tPupa = 1,
                           muAd = 0.09)

migrationBatch <- basicBatchMigration(batchProbs = 0, numPatches = numPatch)

#startTimes <- Sys.time()
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
         verbose = FALSE)
#difftime(Sys.time(), startTimes)

# 2100 days on laptop took 5 minutes


##########
# remove empty files
##########
# combineFiles.R fails when there are empty files, so remove them
allFiles <- list.files(path = simDir, full.names = TRUE)
emptyFiles <- which(file.size(allFiles) < 100)
file.remove(allFiles[emptyFiles])


# bash command to delete the empty files
# turns out, the combine script fails on empty files. sigh.
# 
# find . -name "*.tif" -type 'f' -size -160k -delete
# # link: https://superuser.com/questions/644272/how-do-i-delete-all-files-smaller-than-a-certain-size-in-all-subfolders


##########
# Combine Files
##########
# see combineFiles.R



##########
# Remove initial population
##########
Remove any days with "0" parents
this is to get rid of "0" parents, the ones who started the population so there was 
 no structure before that.
I will check male and female files, then remove both.
Or, screw it, remove first 100 days, set maxVal = 100.
 

-F: set column delimiter
$5/$6: column numbers
print if they equal 0
NR: number record, it is the line number, so keep header
then print first column from those columns.
then, bye eye, get max value

awk -F "," '($5 == "0") && ($6 == "0")' 000_F.csv | awk -F "," '{print $1}'

awk -F "," '(NR==1) || ($1 > 100)' 000_F.csv > newFemFile.csv










detach("package:mPlexCpp", unload=TRUE)


