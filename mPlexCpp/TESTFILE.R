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
###############################################################################
# Clean environment and source files
###############################################################################
rm(list=ls());gc()
library(mPlexCpp)






###############################################################################
# Setup Parameters for Network
###############################################################################

numPatch <- 2
migration <- matrix(data = runif(numPatch*numPatch), nrow = numPatch, ncol = numPatch)
migration <- migration/rowSums(migration)
patchPops = rep(20L,numPatch)
directory = "~/Desktop/HOLD/MGDrivE/"

#setup alleles to initiate patches
alleloTypes <- vector(mode = "list", length = 1L) #3 loci
alleloTypes[[1]]$alleles <- c("W")
alleloTypes[[1]]$probs <- c(1L)
# alleloTypes[[2]]$alleles <- c("W","H")
# alleloTypes[[2]]$probs <- c(1,0)
# alleloTypes[[3]]$alleles <- c("W","H")
# alleloTypes[[3]]$probs <- c(1,0)

AllAlleles <- replicate(n = numPatch, expr = alleloTypes, simplify = FALSE)





# reproductionReference <- MakeReference_DaisyDrive(H = c(0.98, 0.5),
#                                                   R = c(0.0001, 0.0001),
#                                                   S = c(0.0003, 0.004),
#                                                   d = c(0, 0), eta = c("TIME"=4))

reproductionReference <- MakeReference_Multiplex_oLocus()
# reproductionReference$mendelian[[1]]$W <- c(1.0,1.5,2.3)
# reproductionReference$mendelian[[1]]$H <- c(1.0,1.5,2.3)
# reproductionReference$mendelian[[1]]$R <- c(1.0,1.5,2.3)
# reproductionReference$mendelian[[1]]$S <- c(1.0,1.5,2.3)
# 
# reproductionReference$mendelian[[2]]$W <- c(1.0,1.5,2.3)
# reproductionReference$mendelian[[2]]$H <- c(1.0,1.5,2.3)
# reproductionReference$mendelian[[2]]$R <- c(1.0,1.5,2.3)
# reproductionReference$mendelian[[2]]$S <- c(1.0,1.5,2.3)

###############################################################################
# Release Setup
###############################################################################

# create Release List
patchReleases = replicate(n = numPatch,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      larvaeReleases = NULL),
                          simplify = FALSE)


# Create release object to pass to patches
holdRel <- Release_basicRepeatedReleases(releaseStart = 10L,
                                         releaseEnd = 20L,
                                         releaseInterval = 1,
                                         genMos = c("HH"),
                                         numMos = c(50L),
                                         minAge = 16L,
                                         maxAge = 24L,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)


holdRel2 <- Release_basicRepeatedReleases(releaseStart = 600L,
                                          releaseEnd = 610L,
                                          releaseInterval = 2L,
                                          genMos = c("SS"),
                                          numMos = c(10L),
                                          minAge = 16L,
                                          maxAge = 24L,
                                          ageDist = rep(x = 1, times = 24-16+1)/9)

# for(i in seq(1,numPatch,1)){
#   patchReleases[[i]]$maleReleases <- holdRel
# }

patchReleases[[1]]$maleReleases <- holdRel


###############################################################################
# Calculate parameters and initialize network
###############################################################################
netPar = NetworkParameters(nPatch = numPatch,
                           simTime = 100L,
                           alleloTypes = AllAlleles,
                           AdPopEQ = patchPops,
                           runID = 1L,
                           dayGrowthRate = 1.1,
                           beta = 32L)

migrationBatch <- basicBatchMigration(numPatches = numPatch)


mPlex_oneRun(seed = 10,
             networkParameters = netPar,
             reproductionReference = reproductionReference,
             patchReleases = patchReleases,
             migrationMale = migration,
             migrationFemale = migration,
             migrationBatch = migrationBatch,
             output_directory = directory,
             reproductionType = "mPlex_oLocus",
             verbose = TRUE)







detach("package:mPlexCpp", unload=TRUE)
