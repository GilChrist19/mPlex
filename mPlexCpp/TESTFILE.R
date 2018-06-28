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

numPatch <- 1
migration <- matrix(data = runif(numPatch*numPatch), nrow = numPatch, ncol = numPatch)
migration <- migration/rowSums(migration)
patchPops = rep(1000L,numPatch)
directory = "~/Desktop/HOLD/MGDrivE (copy 1)/"

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

# reproductionReference <- MakeReference_Multiplex_oLocus(H = c(0.98, 0.5),
#                                                          R = c(0.0001, 0.0001),
#                                                          S = c(0.0003, 0.004),
#                                                          d = c(0, 0))
reproductionReference <- MakeReference_Multiplex_mLoci()

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
                                      eggReleases = NULL),
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
                           simTime = 500L,
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






system.time(mPlex_oneRun(seed = 10,
             networkParameters = netPar,
             reproductionReference = reproductionReference,
             patchReleases = patchReleases,
             migrationMale = migration,
             migrationFemale = migration,
             migrationBatch = migrationBatch,
             output_directory = directory,
             reproductionType = "mPlex_oLocus",
             verbose = TRUE))






  
  


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





