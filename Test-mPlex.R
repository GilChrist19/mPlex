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

migration = diag(1) #migration matrix
N = nrow(migration) #number of patches
patchPops = rep(50,N) #population of eachpatch
directory <- "~/Desktop/HOLD"

    #setup alleles to initiate patches
alleloTypes <- vector(mode = "list", length = 3) #3 loci
alleloTypes[[1]]$alleles <- c("W","H")
alleloTypes[[1]]$probs <- c(1,0)
alleloTypes[[2]]$alleles <- c("W","H")
alleloTypes[[2]]$probs <- c(1,0)
alleloTypes[[3]]$alleles <- c("W","H")
alleloTypes[[3]]$probs <- c(1,0)

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
reproductionReference <- MakeReference_Multiplex_mLoci(H = c(0.9, 0.9, 0.9),
                                                       R = c(0.0001, 0.0001, 0.0001),
                                                       S = c(0.00003,0.00003,0.00003),
                                                       d = c(0.00001,0.00001,0.00001))

###############################################################################
# Release Setup
###############################################################################

  # create Release List
patchReleases = replicate(n = N,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      larvaeReleases = NULL),
                          simplify = FALSE)

  # Create list of mosquitoes for a release
releaseList <- CreateMosquitoes_Defined_Genotype(genMos = c("HHHHHH", "WHWWWW", "HHWWWW"),
                                                 numMos = c(100,100,100),
                                                 minAge = 16,
                                                 maxAge = 36,
                                                 ageDist = rep(x = 1, times = 36-16+1)/21)

  # Create release object to pass to patches
patchReleases[[1]]$larvaeReleases <- Release_basicRepeatedReleases(releaseStart = 5,
                                                                 releaseEnd = 10,
                                                                 releaseInterval = 5,
                                                                 releaseVector = releaseList,
                                                                 sex = "L")

###############################################################################
# Calculate parameters and initialize network
###############################################################################

    # calculate network parameters, auxiliary function
netPar = Network.Parameters(nPatch = N,simTime = 250,
                            alleloTypes = AllAlleles,
                            AdPopEQ = patchPops,
                            parallel = FALSE)

    # initialize network!
network = Network$new(networkParameters = netPar,
                      patchReleases = patchReleases,
                      reproductionType = "mPlex_mLoci",
                      offspringReference = reproductionReference,
                      migrationMale = migration,
                      migrationFemale = migration,
                      directory = directory)



    #reset network
network$oneRun()
network$reset()


###############################################################################
# Post-Run Stuff
###############################################################################
splitOutput(directory = directory)





















