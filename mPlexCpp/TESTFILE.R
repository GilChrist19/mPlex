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

numPatch <- 10
set.seed(10)
migration <- matrix(data = runif(numPatch*numPatch), nrow = numPatch, ncol = numPatch)
migration <- migration/rowSums(migration)

patchPops = rep(100L,numPatch)

directory1 = "~/Desktop/HOLD/MGDrivE (copy 1)/"
directory2 <- "~/Desktop/HOLD/MGDrivE (copy 1) (copy 1)/"

#setup alleles to initiate patches
alleloTypes <- vector(mode = "list", length = 1L) #1 locus
alleloTypes[[1]]$alleles <- c("W")
alleloTypes[[1]]$probs <- c(1L)

AllAlleles <- replicate(n = numPatch, expr = alleloTypes, simplify = FALSE)


# reproductionReference <- MakeReference_DaisyDrive(H = c(0.98, 0.5),
#                                                   R = c(0.0001, 0.0001),
#                                                   S = c(0.0003, 0.004),
#                                                   d = c(0, 0), eta = c("TIME"=4))

# reproductionReference <- MakeReference_Multiplex_oLocus(H = c(0.98, 0.5),
#                                                          R = c(0.0001, 0.0001),
#                                                          S = c(0.0003, 0.004),
#                                                          d = c(0, 0))

# This sets up a basic CRISPR drive, with perfect homing and no resistance or backgorund mutation
reproductionReference <- MakeReference_Multiplex_mLoci(H = 1.0, R = 0, S = 0, d = 0,
                                                       eta = c("HH"=0.5))


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
netPar = NetworkParameters(nPatch = numPatch,
                           simTime = 2000L,
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
             output_directory = directory1,
             reproductionType = "mPlex_mLoci",
             verbose = TRUE)



# split the output by patch
splitOutput(directory = directory1, numCores = 1)

# aggregate by genotype.
AnalyzeOutput_mLoci_Daisy(readDirectory = directory1,
                          saveDirectory = directory2,
                          genotypes = list(NULL),
                          collapse = c(FALSE),
                          numCores = 1)

# plot for example
Plot_mPlex(directory = directory2, whichPatches = NULL, totalPop = FALSE)

detach("package:mPlexCpp", unload=TRUE)
  
  






###############################################################################
# DON'T RUN THIS STUFF. SHOULD WORK, ISN'T ORGANIZED, MOSTLY FOR TESTING
###############################################################################


# repetitions wrapper - no reinitializing memory between reps. 
dirVec <- paste0(directory, 1:4)

mPlex_runRepetitions(seed = 10,
                     networkParameters = netPar,
                     reproductionReference = reproductionReference,
                     patchReleases = patchReleases,
                     migrationMale = migration,
                     migrationFemale = migration,
                     migrationBatch = migrationBatch,
                     output_directory = dirVec,
                     reproductionType = "mPlex_mLoci",
                     verbose = TRUE)


splitOutput(directory = directory, numCores = 1)
AnalyzeOutput_oLocus(readDirectory = "~/Desktop/HOLD/MGDrivE (copy 1)/",
                     saveDirectory = "~/Desktop/HOLD/MGDrivE (copy 1) (copy 1)/",
                     alleles = list(list(c(NULL)),list(c(NULL))),
                     collapse = list(c(F),c(F)),
                     numCores = 1)

AnalyzeOutput_mLoci_Daisy(readDirectory = "~/Desktop/HOLD/MGDrivE/",
                          saveDirectory = "~/Desktop/HOLD/mPlex/",
                          genotypes = list(NULL),
                          collapse = c(FALSE),
                          numCores = 1)


Plot_mPlex(directory = "~/Desktop/HOLD/MGDrivE (copy 1) (copy 1)/", whichPatches = NULL, totalPop = TRUE)

detach("package:mPlexCpp", unload=TRUE)




















###############################################################################
# write split function
###############################################################################
library(data.table)














test <- matrix(data = 10L,nrow = 10000,ncol = 20,dimnames = list(NULL,LETTERS[1:20]))
write.csv(x = test,file = "~/Desktop/HOLD/myTest.csv",row.names = FALSE)




system.time(AnalyzeOutput_oLocus(readDirectory = "~/Desktop/HOLD/MGDrivE/",
                                 saveDirectory = "~/Desktop/HOLD/mPlex/",
                                 alleles = list(list(c(NULL)),list(c(NULL))),
                                 collapse = list(c(F),c(F)),
                                 numCores = 1))





Rprof(filename = "~/Desktop/HOLD/profiling.out", interval = 0.01, line.profiling = TRUE)
AnalyzeOutput_oLocus(readDirectory = "~/Desktop/HOLD/MGDrivE (copy 1)/",
                     saveDirectory = "~/Desktop/HOLD/MGDrivE (copy 1) (copy 1)/",
                     alleles = list(list(c(NULL)),list(c(NULL))),
                     collapse = list(c(F),c(F)),
                     numCores = 1)



summaryRprof(filename = "~/Desktop/HOLD/profiling.out", lines = "both")





readRDS(file = "~/Desktop/HOLD/mPlex/20180622_Run000001_MYTEST2.rds")

patches




mName = grep(pattern = paste("ADM", patches[1], sep = ".*"),
             x = dirFiles,ignore.case = FALSE, perl = TRUE,
             value = TRUE, useBytes = TRUE)[1]



names = grep(pattern = patches[1], x = dirFiles, fixed = TRUE, value = TRUE)






microbenchmark::microbenchmark(grep(pattern = paste("ADM", patches[1], sep = ".*"),
                                    x = dirFiles,ignore.case = FALSE, perl = TRUE,
                                    value = TRUE, useBytes = TRUE)[1],
                               grep(pattern = patches[1], x = dirFiles, fixed = TRUE, value = TRUE),
                               times = 1000)





###############################################################################
# When profiling, use this to see it
###############################################################################

cat("NumSamp  PercentSamp  CumPercentSamp  NumSampTree  PercentSampTree  Function")
system(sprintf("google-pprof --text --cum --lines  /bin/ls %sprofile.log", directory_det_c), intern = TRUE)





