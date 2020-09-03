###############################################################################
#   _    __          _ __     
#  | |  / /__  _____(_) /_  __
#  | | / / _ \/ ___/ / / / / /
#  | |/ /  __/ /  / / / /_/ / 
#  |___/\___/_/  /_/_/\__, /  
#                    /____/   
#
# Marshall Lab - CKMR Project
# Jared Bennett
# jared_bennett@berkeley.edu
###############################################################################
###############################################################################
# 20200903
#  Created this script to run simulations for the initial Verily sequencing analysis.
#  This was copied from earlier simulation files, in the "YogitaData" directory
#  It has been re-arranged, lots of path stuff has been moved to the top, scratch 
#  space is fully managed, output is properly stored as zipped files. It loops 
#  over different migration patterns, but could be expanded for other things.
#
#  Need to run over 2 adult mortalities:
#   0.09
#   0.125
#  The kernel is defined on a daily basis, so this mortality change doesn't affect it.
#
###############################################################################
# Clean Environment
###############################################################################
rm(list=ls());gc()
#library(CKMR)


simTime <- 190 # 100 day burn-in and 3 months of data
patchPops <- 30
nCore <- 2



###############################################################################
# Setup Directories
###############################################################################
workDir <- "~/Desktop/mPlex/Main/Verily"
outDir <- "~/Desktop/OUTPUT/mPlex"
simDir <- file.path(outDir,"scratch")
kernels <- c("c2_kernel1","c2_kernel_exp80")
saveDir <- setNames(object = file.path(outDir, paste0(gsub(pattern = "-", replacement = "",
                                                           x = Sys.Date(), fixed = TRUE),
                                                      "_", kernels)),
                    nm = kernels)

if(!dir.exists(paths = outDir)){
  dir.create(path = outDir)
} else {
  unlink(x = list.files(outDir, full.names = TRUE), recursive = TRUE, force = TRUE)
}
for(i in c(simDir, saveDir)) dir.create(i)


# source aux functions
source(file.path(workDir,"combineFiles.R"))


###############################################################################
# Load Network
###############################################################################
# load kernels
for(k in kernels){
  load(file.path(workDir,paste0(k,".rds")))
}

# load node information for later
load(file.path(workDir,"c2_centroids_info.rds"))

# number of patches in the simulation
numPatch <- NROW(c2_centroids_info)


###############################################################################
# Sampling Setup
###############################################################################
# sampling can be different for every patch, and every life stage
# final object is a list with 2 matrices in it, one for when to sample, one for 
#  how many to sample

# coverage
#  all aquatic stages get 0% coverage
#  Females: need to catch ~440 per week, over all 28 traps
#           This is 15.7 mosquitoes per trapping, or the entire expected female population
#  Males: are about 2/3 as effective as females
sampCov <- matrix(data = c(0,0,0,0.66,1), nrow = numPatch, ncol = 5, byrow = TRUE)

# sampling time
# Basic, every patch has every life stage (5 of them), never sampled
#  to never sample, must put time to simTime + 1, otherwise modulo throws errors
#  Get males/females weekly, from designated trap sites
sampDay <- matrix(data = simTime + 1, nrow = numPatch, ncol = 5)
sampDay[as.logical(c2_centroids_info$trap),4:5] <- 7


# list to pass to mPlex
samplingScheme <- list("samplingDays"=sampDay, "samplingCoverage"=sampCov)


###############################################################################
# Random Setup
###############################################################################
# These are basically unused in this simulation

# create Release List
#  there are no releases, this is just null
patchReleases = replicate(n = numPatch,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL),
                          simplify = FALSE)

# setup alleles to initiate patches
reference <- list('eta'=numeric(0), 'phi'=numeric(0), 'omega'=numeric(0),
                  'xiF'=numeric(0), 'xiM'=numeric(0), 's'=numeric(0))

# we don't use batch migration here
migrationBatch <- CKMR::basicBatchMigration(batchProbs = 0, numPatches = numPatch)


###############################################################################
# Calculate parameters and Run Sims
###############################################################################
# finalize parameters
netPar = CKMR::NetworkParameters(nPatch = numPatch,
                                 simTime = simTime,
                                 AdPopEQ = rep.int(x = patchPops, times = numPatch),
                                 runID = 1L,
                                 dayGrowthRate = 1.175,
                                 beta = 20, tEgg = 2, tLarva = 5, tPupa = 1,
                                 muAd = 0.125)

# loop over kernels
for(k in kernels){
  ##########
  # run simulation
  ##########
  CKMR::runCKMR(seed = 10,
                 numThreads = nCore,
                 networkParameters = netPar,
                 reproductionReference = reference,
                 patchReleases = patchReleases,
                 migrationMale = get(k),
                 migrationFemale = get(k),
                 migrationBatch = migrationBatch,
                 samplingParameters = samplingScheme,
                 outputDirectory = simDir,
                 verbose = FALSE)
  
  
  ########################################
  # Analysis and Cleanup
  ########################################
  ##########
  # Remove Empty Files
  ##########
  #  combineFiles.R fails when there are empty files
  allFiles <- list.files(path = simDir, full.names = TRUE)
  emptyFiles <- which(file.size(allFiles) < 100)
  file.remove(allFiles[emptyFiles])
  
  
  ##########
  # Combine Files
  ##########
  #  no output - writes files in place
  combineFiles(mainDir=simDir, workIndicator=999999, fPattern=c("M","F"))
  
  
  ##########
  # Remove Initial Population
  ##########
  #  Remove any days with "0" parents, the ones who started the population, because 
  #  there is no structure before that.
  #  Explanation:
  #   -F: set column delimiter
  #   $1: column numbers
  #   print if they are greater than line 100
  
  # read files are hardCoded in the combineFiles.R script
  readFiles <- file.path(simDir, c('000_F.csv', '000_M.csv'))
  writeFiles <- file.path(saveDir[[k]], c('cut_F.csv', 'cut_M.csv'))
    
  # loop over files, remove first 100 days, store in the save directory
  for(cFile in 1:2){
    # build command, 
    cmd <- file.path("-F ',' '(NR==1) || ($1 > 100)' ", readFiles[cFile], 
                     ' > ', writeFiles[cFile], fsep = '' )
  
    # Pass down to cmdline
    system2(command = 'awk', args = cmd)
  } # end loop

  
  ##########
  # Store
  ##########
  # build terminal command
  # -C change to the following directory
  # -c output to this place (sent to stdout here)
  # -f read from this file
  # | lbzip2 pipe stdin to lbzip2
  # -9 lbzip2 level 9
  # -n number of cores to use
  # > output to here
  cmd <- file.path('-C ', outDir, ' --remove-files -cf - ', basename(saveDir[[k]]),
                   ' | lbzip2 -9 -n ', nCore, ' > ', saveDir[[k]], '.tar.bz2', fsep = '')

  # store
  system2(command = 'tar', args = cmd )

  
  ##########
  # Clean
  ##########
  # clear scratch for next run
  unlink(x = list.files(path = simDir, full.names = TRUE))
  
} # end loop over kernels

# remove scratch space
unlink(x = simDir, recursive = TRUE)


