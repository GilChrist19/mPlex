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
source("~/Desktop/mPlex/1_MosquitoClass.R")
source("~/Desktop/mPlex/2_PatchClass.R")
source("~/Desktop/mPlex/2_PatchSimulation.R")
source("~/Desktop/mPlex/3_NetworkClass.R")
source("~/Desktop/mPlex/Network_Parameters_Equilibrium.R")
source("~/Desktop/mPlex/Auxiliary_Functions.R")

###############################################################################
# Setup Parameters for Network
###############################################################################

migration = diag(4) #migration matrix
N = nrow(migration) #number of patches
patchPops = rep(10,N) #population of eachpatch
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
netPar = Network.Parameters(nPatch = N,simTime = 25,
                            alleloTypes = AllAlleles, tAdult = 21,
                            AdPopEQ = patchPops)

    # initialize network!
network = Network$new(networkParameters = netPar,
                      patchReleases = patchReleases,
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





















