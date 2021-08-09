###############################################################################
#                            ____  __          ______
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/
#                                               /_/   /_/
###############################################################################
###############################################################################
# 20210805
#  Update to 20210722
#  Keep the sampling parameters, only need the 1000 pop, add males back in though.
#
###############################################################################
# Clean environment and source files
###############################################################################
rm(list=ls());gc()
library(CKMR)
source("~/Desktop/mPlex/Main/YogitaData/combineFiles.R")

set.seed(seed = 10)
simTime <- 190
numThreads <- 1
###############################################################################
# Setup Directories
###############################################################################
topDirectory <- "~/Desktop/OUTPUT/CKMR"

if(!dir.exists(paths = topDirectory)){
  dir.create(path = topDirectory)
} else {
  unlink(x = list.files(topDirectory, full.names = TRUE), recursive = TRUE, force = TRUE)
}



###############################################################################
# Setup Parameters for Network
###############################################################################
# single patch
migration <- as.matrix(x = 1)
numPatch <- nrow(migration)

# batch migration
#  set to 0
migrationBatch <- basicBatchMigration(batchProbs = 0, numPatches = numPatch)

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

# male/female get sampled daily
sampDay[ ,c(2,4,5)] <- 1

# coverage
#  sample females 10% per week, converted to daily
#  sample larvae 10% per week, converted to daily
sampCov <- matrix(data = c(0,0.10,0,0.10,0.10)/7, nrow = numPatch, ncol = 5, byrow = TRUE)

# list to pass to mPlex
samplingScheme <- list("samplingDays"=sampDay, "samplingCoverage"=sampCov)


###############################################################################
# Setup Sweep
###############################################################################
# Setup desired parameter ranges here
#  This then creates a dataframe with every combination

# sweep over cube changing parameters
# nRep: indices of reps to do
# nPop: population sizes to test
paramCombo <- as.matrix(expand.grid('nRep' = 1:15,
                                    'nPop' = c(1000) ))
numPC <- NROW(paramCombo)


########################################
# Loop over parameters
########################################
for(i in 1:numPC){
  
  ####################
  # Setup Folder
  ####################
  # width = 7 handles up to a million population size
  simDir <- file.path(topDirectory,
                      formatC(x = paramCombo[i,'nPop'], width = 7, format = "d", flag = "0"),
                      formatC(x = paramCombo[i,'nRep'], width = 3, format = "d", flag = "0"))
  dir.create(path = simDir, recursive = TRUE)

  ####################
  # Network Parameters
  ####################
  netPar = NetworkParameters(nPatch = numPatch,
                             simTime = simTime,
                             AdPopEQ = paramCombo[i,'nPop'],
                             runID = 1L,
                             dayGrowthRate = 1.175,
                             beta = 20, tEgg = 2, tLarva = 5, tPupa = 1,
                             muAd = 0.09)
  
  ####################
  # Run Sim!
  ####################
  # force evaluation of the sampling
  seed <- sample(x = (-.Machine$integer.max):.Machine$integer.max, size = 1, replace = FALSE)
  runCKMR(seed = seed,
          numThreads = numThreads,
          networkParameters = netPar,
          reproductionReference = reference,
          patchReleases = patchReleases,
          migrationMale = migration,
          migrationFemale = migration,
          migrationBatch = migrationBatch,
          samplingParameters = samplingScheme,
          outputDirectory = simDir,
          verbose = FALSE)

  ####################
  # remove empty files
  ####################
  # combineFiles.R fails when there are empty files, so remove them
  allFiles <- list.files(path = simDir, full.names = TRUE)
  emptyFiles <- which(file.size(allFiles) < 100)
  file.remove(allFiles[emptyFiles])
  
  ####################
  # Combine Files
  ####################
  # no output
  combineFiles(mainDir=simDir, workIndicator=25, fPattern=c("L","F","M"))

  ####################
  # Remove initial population
  ####################
  # Remove any days with "0" parents
  # this is to get rid of "0" parents, the ones who started the population so there was 
  #  no structure before that.
  # I will check male and female files, then remove both.
  # Or, screw it, remove first 100 days, set maxVal = 100.
  #  
  # 
  # -F: set column delimiter
  # $5/$6: column numbers
  # print if they equal 0
  # NR: number record, it is the line number, so keep header
  # then print first column from those columns.
  # then, bye eye, get max value
  # 
  # This one takes off anyone with 0 parents
  # awk -F "," '($5 == "0") && ($6 == "0")' 000_F.csv | awk -F "," '{print $1}'
  # 
  # This one takes off first 100 days
  # awk -F "," '(NR==1) || ($1 > 100)' 000_F.csv > newFemFile.csv
  
  readFiles <- c('/000_F.csv', '/000_L.csv', '/000_M.csv')
  writeFiles <- c('/cut_F.csv', '/cut_L.csv', '/cut_M.csv')
  for(cFile in 1:length(writeFiles)){
    
    # build command, 
    cmd <- file.path("-F ',' '(NR==1) || ($1 > 100)' ", simDir, readFiles[cFile], 
                     ' > ', simDir, writeFiles[cFile], fsep = '' )
    
    # Pass down to cmdline
    system2(command = 'awk', args = cmd)
    
  } # end loop
  
  ####################
  # Cleanup
  ####################
  # Only need the final cut pops
  # May as well remove the rest
  allFiles <- list.files(path = simDir, full.names = TRUE)
  trashFiles <- grep(pattern = 'cut', x = allFiles, fixed = TRUE, value = TRUE, invert = TRUE)
  file.remove(trashFiles)
  
} # end parameter sweep


detach("package:CKMR", unload=TRUE)

