#test this thing

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

migration = diag(1) #migration matrix
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
                          expr = list(maleReleases = NULL,femaleReleases = NULL),
                          simplify = FALSE)

  # Create list of mosquitoes for a release
releaseList <- CreateMosquitoes_Defined_Genotype(genMos = c("WWWWWW", "WHWWWW", "HHWWWW"),
                                                 numMos = c(10,10,10),
                                                 minAge = 16,
                                                 maxAge = 36,
                                                 ageDist = rep(x = 1, times = 36-16+1)/21)

  # Create release object to pass to patches
patchReleases[[1]]$maleReleases <- Release_basicRepeatedReleases(releaseStart = 10,
                                                                 releaseEnd = 100,
                                                                 releaseInterval = 10,
                                                                 releaseVector = releaseList,
                                                                 sex = "M")

###############################################################################
# Calculate parameters and initialize network
###############################################################################

    # calculate network parameters, auxiliary function
netPar = Network.Parameters(nPatch = N,simTime = 10,
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






















