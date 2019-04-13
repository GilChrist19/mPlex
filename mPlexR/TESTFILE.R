###############################################################################
#                                    ____  _
#                          _ __ ___ |  _ \| | _____  __
#                         | '_ ` _ \| |_) | |/ _ \ \/ /
#                         | | | | | |  __/| |  __/>  <
#                         |_| |_| |_|_|   |_|\___/_/\_\
#
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
library(mPlexR)

###############################################################################
# Setup Parameters for Network
###############################################################################

migration = diag(3) #matrix(data = c(0.99, 0, 0.05, 0.05, 0, 0.99, 0.05, 0.5, 0.03,0.03,0.99,0.04, 0.02,0.04,0.04,0.99), nrow = 4, ncol = 4, byrow = TRUE) #migration matrix
N = nrow(migration) #number of patches
patchPops = rep(50L,N) #population of eachpatch
directory <- "~/Desktop/OUTPUT/mPlex/"

    #setup alleles to initiate patches
alleloTypes <- vector(mode = "list", length = 2L) #3 loci
alleloTypes[[1]]$alleles <- c("W")
alleloTypes[[1]]$probs <- c(1L)
alleloTypes[[2]]$alleles <- c("W","H")
alleloTypes[[2]]$probs <- c(1,0)
# alleloTypes[[3]]$alleles <- c("W","H")
# alleloTypes[[3]]$probs <- c(1,0)

AllAlleles <- replicate(n = N, expr = alleloTypes, simplify = FALSE)

###############################################################################
# reproduction Setup
###############################################################################
  #setup reference for offspring production.
  # This must match the reproductionType in network initialization
  #   "DaisyDrive" = MakeReference_DaisyDrive()
  #   "mPlex_oLocus" = MakeReference_Multiplex_oLocus()
  #   "mPlex_mLoci" = MakeReference_Multiplex_mLoci()

#these numbers are made up. Just need them all the same length, and that length
# must match the length of AlleloTypes
s_frac <- vector(mode = "list", length = length(alleloTypes))
s_frac[[1]] <- list("HHHH"=1)
reproductionReference <- MakeReference_DaisyDrive(H = c(0.98, 0.5),
                                                  R = c(0.0001, 0.0001),
                                                  S = c(0.0003, 0.004),
                                                  d = c(0, 0),
                                                  s_frac = s_frac)

###############################################################################
# Release Setup
###############################################################################

  # create Release List
patchReleases = replicate(n = N,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      larvaeReleases = NULL),
                          simplify = FALSE)


  # Create release object to pass to patches
holdRel <- Release_basicRepeatedReleases(releaseStart = 500L,
                                                                 releaseEnd = 510L,
                                                                 releaseInterval = 2L,
                                                                 genMos = c("HHHH"),
                                                                 numMos = c(50L),
                                                                 minAge = 16L,
                                                                 maxAge = 24L,
                                                                 ageDist = rep(x = 1, times = 24-16+1)/9)


holdRel2 <- Release_basicRepeatedReleases(releaseStart = 600L,
                                         releaseEnd = 610L,
                                         releaseInterval = 2L,
                                         genMos = c("SSSS"),
                                         numMos = c(10L),
                                         minAge = 16L,
                                         maxAge = 24L,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)

for(i in seq(1,N,1)){
  patchReleases[[i]]$maleReleases <- holdRel
}


###############################################################################
# Calculate parameters and initialize network
###############################################################################

    # calculate network parameters, auxiliary function
netPar = Network.Parameters(nPatch = N,simTime = 1000L,
                            alleloTypes = AllAlleles,
                            AdPopEQ = patchPops,
                            runID = 1L,
                            dayGrowthRate = 1.1,
                            beta = 32L)

    # initialize network!
network = Network$new(networkParameters = netPar,
                      patchReleases = patchReleases,
                      reproductionType = "DaisyDrive",
                      offspringReference = reproductionReference,
                      migrationMale = migration,
                      migrationFemale = migration,
                      directory = directory)


#Rprof(interval = 0.01, line.profiling = TRUE)
set.seed(seed = 42)
network$oneRun()
#summaryRprof(lines = "both")
network$reset()


###############################################################################
# Post-Run Stuff
###############################################################################
splitOutput(directory = directory)
AnalyzeOutput_mLoci_Daisy(readDirectory = directory,
                          saveDirectory = "~/Desktop/OUTPUT/mPlex/HOLD/",
                          fileName = "ProfileTest",
                          genotypes = list(c(NULL),c(NULL)),
                          collapse = c(F,T))

Plot_mPlex(file = "/home/jared/Desktop/OUTPUT/mPlex/HOLD/20190412_Run1_ProfileTest.rds", totalPop = T)












detach("package:mPlexR", unload=TRUE)
