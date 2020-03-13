## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  hold = TRUE,
  fig.width = 7,
  fig.height = 10.5,
  eval = TRUE
)

## ---- eval=FALSE--------------------------------------------------------------
#  ####################
#  # Load libraries
#  ####################
#  library(mPlexCpp)
#  
#  ####################
#  # Output Folder
#  ####################
#  outFolder <- "mPlex"
#  
#  simDir <- file.path(outFolder, "simDir")
#  aggDir <- file.path(outFolder, "aggDir")
#  for(i in c(outFolder, simDir, aggDir)){ dir.create(i) }
#  
#  
#  ####################
#  # Landscape
#  ####################
#  # a 3-node network with 2% per day migration rate
#  #  Remember, rows need to sum to 1.
#  sitesNumber <- 3
#  moveMat <- matrix(data = c(0.98, 0.02, 0,
#                             0.02, 0.98, 0,
#                             0, 0, 1), nrow = sitesNumber, ncol = sitesNumber, byrow = TRUE)
#  
#  
#  ####################
#  # Inheritance pattern
#  ####################
#  # 1-locus CRISPR-like drive system, with no extra genotype-specific costs
#  #  97% cutting rate, 100% homing rate, no backgorund mutation
#  reproductionReference <- MakeReference_Multiplex_mLoci(cRateM = c(0.97),
#                                                         hRateM = c(1.0),
#                                                         rRateM = c(0),
#                                                         dM = c(0))
#  
#  
#  ####################
#  # Setup Initial genotype ratios
#  ####################
#  # 1 locus, start completely wild-type
#  aTypes <- vector(mode = "list", length = 1L) #1 locus
#  aTypes[[1]]$alleles <- c("W")
#  aTypes[[1]]$probs <- c(1)
#  
#  # replicate so each patch starts the same
#  #  This is optional. If a length-1 list is supplied, it is internally replicated
#  #  and all patches begin the same. Otherwise, the list must have length equal
#  #  to the number of patches.
#  AllAlleles <- replicate(n = sitesNumber, expr = aTypes, simplify = FALSE)
#  
#  
#  ####################
#  # Setup releases and batch migration
#  ####################
#  # create Release List
#  patchReleases = replicate(n = sitesNumber,
#                            expr = list(maleReleases = NULL,
#                                        femaleReleases = NULL,
#                                        eggReleases = NULL,
#                                        matedFemaleReleases = NULL),
#                            simplify = FALSE)
#  
#  # Create release object to pass to patches
#  patchReleases[[1]]$maleReleases <- basicRepeatedReleases(releaseStart = 50,
#                                           releaseEnd = 60,
#                                           releaseInterval = 5,
#                                           genMos = c("HH"),
#                                           numMos = c(25),
#                                           minAge = 16,
#                                           maxAge = 24,
#                                           ageDist = rep(x = 1, times = 24-16+1)/9)
#  
#  # default migration rate is 0, so no actual batch migration
#  migrationBatch <- basicBatchMigration(numPatches = sitesNumber)
#  
#  
#  ####################
#  # Parameter Setup
#  ####################
#  netPar = NetworkParameters(nPatch = sitesNumber,
#                             simTime = 365,
#                             sampTime = 2,
#                             AdPopEQ = 100,
#                             runID = 1L,
#                             dayGrowthRate = 1.175,
#                             beta = 20, tEgg = 5, tLarva = 6, tPupa = 4, muAd = 0.09)
#  
#  ####################
#  # Run Simulation
#  ####################
#  runMPlex(seed = 10,
#           numThreads = 1,
#           numReps = 5,
#           networkParameters = netPar,
#           reproductionReference = reproductionReference,
#           initAlleles = AllAlleles,
#           patchReleases = patchReleases,
#           migrationMale = moveMat,
#           migrationFemale = moveMat,
#           migrationBatch = migrationBatch,
#           outputDirectory = simDir,
#           reproductionType = "mPlex_mLoci",
#           verbose = FALSE)
#  
#  
#  ####################
#  # Post Analysis
#  ####################
#  # setup aggregation key
#  #  this example sets all genotypes as different
#  genOI_mLoci_Daisy(outputFile = file.path(aggDir, "0_AggKey.csv"),
#                    genotypes = list(NULL), collapse = c(FALSE))
#  
#  
#  # aggregate experiment by aggregation key
#  simAggregation(readDirectory = simDir, writeDirectory = aggDir,
#                 simTime = netPar$simTime, sampTime = netPar$sampTime)
#  
#  
#  # plot for example
#  plotmPlexSingle(directory = aggDir, whichPatches = NULL, nonZeroGen = FALSE)
#  plotmPlexMult(directory = aggDir,whichPatches = NULL, nonZeroGen = FALSE, lwd=0.35,alpha=0.75)

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(mPlexCpp)


## -----------------------------------------------------------------------------
####################
# Output Folder
####################
outFolder <- "mPlex"

simDir <- file.path(outFolder, "simDir")
aggDir <- file.path(outFolder, "aggDir")
for(i in c(outFolder, simDir, aggDir)){ dir.create(i) }


## -----------------------------------------------------------------------------
####################
# Landscape
####################
# a 3-node network with 2% per day migration rate
#  Remember, rows need to sum to 1.
sitesNumber <- 3
moveMat <- matrix(data = c(0.98, 0.02, 0,
                           0.02, 0.98, 0,
                           0, 0, 1), nrow = sitesNumber, ncol = sitesNumber, byrow = TRUE)

moveMat

## -----------------------------------------------------------------------------
####################
# Inheritance pattern
####################
# 1-locus CRISPR-like drive system, with no extra genotype-specific costs
#  97% cutting rate, 100% homing rate, no backgorund mutation
reproductionReference <- MakeReference_Multiplex_mLoci(cRateM = c(0.97),
                                                       hRateM = c(1.00),
                                                       rRateM = c(0),
                                                       dM = c(0))

## -----------------------------------------------------------------------------
####################
# Setup Initial genotype ratios
####################
# 1 locus, start completely wild-type
aTypes <- vector(mode = "list", length = 1L) #1 locus
aTypes[[1]]$alleles <- c("W")
aTypes[[1]]$probs <- c(1)

# replicate so each patch starts the same
#  This is optional. If a length-1 list is supplied, it is internally replicated 
#  and all patches begin the same. Otherwise, the list must have length equal 
#  to the number of patches.
AllAlleles <- replicate(n = sitesNumber, expr = aTypes, simplify = FALSE)

## -----------------------------------------------------------------------------
####################
# Setup releases and batch migration
####################
# create release List
#  mPlex pulls things out by name
patchReleases = replicate(n = sitesNumber,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL,
                                      matedFemaleReleases = NULL),
                          simplify = FALSE)

# Create release object to pass to patches
patchReleases[[1]]$maleReleases <- basicRepeatedReleases(releaseStart = 50,
                                         releaseEnd = 60,
                                         releaseInterval = 5,
                                         genMos = c("HH"),
                                         numMos = c(25),
                                         minAge = 16,
                                         maxAge = 24,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)

## -----------------------------------------------------------------------------
# look at male release structure
patchReleases[[1]]$maleReleases

## -----------------------------------------------------------------------------
# batch migration is disabled by setting the probability to 0
migrationBatch <- basicBatchMigration(batchProbs=0, numPatches=sitesNumber)

## -----------------------------------------------------------------------------
####################
# Parameter Setup
####################
netPar = NetworkParameters(nPatch = sitesNumber,
                           simTime = 365,
                           sampTime = 1,
                           AdPopEQ = 100,
                           runID = 1L,
                           dayGrowthRate = 1.175,
                           beta = 20, tEgg = 5, tLarva = 6, tPupa = 4, muAd = 0.09)

## -----------------------------------------------------------------------------
####################
# Run Simulation
####################
runMPlex(seed = 10,
         numThreads = 1,
         numReps = 5, 
         networkParameters = netPar,
         reproductionReference = reproductionReference,
         initAlleles = AllAlleles,
         patchReleases = patchReleases,
         migrationMale = moveMat,
         migrationFemale = moveMat,
         migrationBatch = migrationBatch,
         outputDirectory = simDir,
         reproductionType = "mPlex_mLoci",
         verbose = FALSE)

# list folders to show that they have been created
list.files(path = outFolder)
list.files(path = simDir)

## -----------------------------------------------------------------------------
# read in male and female files
fDataFrame <- read.csv(file = list.files(path = simDir, full.names = TRUE)[1], 
                       header = TRUE)
mDataFrame <- read.csv(file = tail(x = list.files(path = simDir, full.names = TRUE), 
                                   n = 1), 
                       header = TRUE)

# look at male file header
colnames(mDataFrame)

## -----------------------------------------------------------------------------
# look at female file header
colnames(fDataFrame)

## -----------------------------------------------------------------------------
head(x = fDataFrame, n = 5)

## -----------------------------------------------------------------------------
fDataFrame[50:55,]

## -----------------------------------------------------------------------------
####################
# Post Analysis
####################
# setup aggregation key
#  this example sets all genotypes as different
genOI_mLoci_Daisy(outputFile = file.path(aggDir, "0_AggKey.csv"), 
                  genotypes = list(NULL), collapse = c(FALSE))

list.files(path = aggDir)

## -----------------------------------------------------------------------------
aggKey <- read.csv(file = file.path(aggDir, "0_AggKey.csv"), header = TRUE)
aggKey

## -----------------------------------------------------------------------------
# aggregate experiment by aggregation key
simAggregation(readDirectory = simDir, writeDirectory = aggDir, 
               simTime = netPar$simTime, sampTime = netPar$sampTime)

list.files(path = aggDir)

## -----------------------------------------------------------------------------
# read in male and female files
fMat <- as.matrix(read.csv(file = list.files(path = aggDir, full.names = TRUE)[2], 
                       header = TRUE, check.names = FALSE))
mMat <- as.matrix(read.csv(file = tail(x = list.files(path = aggDir, full.names = TRUE), 
                                   n = 1), 
                       header = TRUE, check.names = FALSE))

# look at male and female file headers
colnames(mMat)
colnames(fMat)

## -----------------------------------------------------------------------------
head(x = mMat, n = 3)
head(x = fMat, n = 3)

## -----------------------------------------------------------------------------
# plot the first repetition
plotmPlexSingle(directory = aggDir, whichPatches = NULL, nonZeroGen = FALSE)

## -----------------------------------------------------------------------------
# plot all 5 repetitions together
plotmPlexMult(directory = aggDir,whichPatches = NULL, nonZeroGen = FALSE, lwd=0.35,alpha=0.75)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## ---- eval=TRUE---------------------------------------------------------------
####################
# Load libraries
####################
library(mPlexCpp)

####################
# Output Folder
####################
outFolder <- "mPlex"

simDir <- file.path(outFolder, "simDir")
aggDir <- file.path(outFolder, "aggDir")
for(i in c(outFolder, simDir, aggDir)){ dir.create(i) }

                        
####################
# Landscape
####################
# a 3-node network with 2% per day migration rate
#  Remember, rows need to sum to 1.
sitesNumber <- 3
moveMat <- matrix(data = c(0.98, 0.02, 0,
                           0.02, 0.98, 0,
                           0, 0, 1), nrow = sitesNumber, ncol = sitesNumber, byrow = TRUE)


####################
# Inheritance pattern
####################
# 1-locus CRISPR-like drive system
#  97% cutting rate, 90% homing rate, 33% resistant 1 allele rate
#  (ie, 1:2 R1:R2 generation), no backgorund mutation
reproductionReference <- MakeReference_Multiplex_mLoci(cRateM = c(0.97),
                                                       hRateM = c(0.90),
                                                       rRateM = c(0.33),
                                                       dM = c(0), omega = c("HH"=0.5, "RR" = 0.9))


####################
# Setup Initial genotype ratios
####################
# 1 locus, start completely wild-type in patches 1 and 2
#  start with a mixture in patch 3
aTypes1 <- vector(mode = "list", length = 1L) #1 locus
aTypes1[[1]]$alleles <- c("W")
aTypes1[[1]]$probs <- c(1L)
aTypes3 <- vector(mode = "list", length = 1L) #1 locus
aTypes3[[1]]$alleles <- c("W", "R", "S")
aTypes3[[1]]$probs <- c(0.8, 0.15, 0.05)

# replicate so each patch starts the same
AllAlleles <- list(aTypes1, aTypes1, aTypes3)


####################
# Setup releases and batch migration
####################
# create Release List
patchReleases = replicate(n = sitesNumber,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL,
                                      matedFemaleReleases = NULL),
                          simplify = FALSE)

# Create release object to pass to patches
rel1 <- basicRepeatedReleases(releaseStart = 50,
                                         releaseEnd = 60,
                                         releaseInterval = 5,
                                         genMos = c("HH"),
                                         numMos = c(25),
                                         minAge = 16,
                                         maxAge = 24,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)

rel2 <- basicRepeatedReleases(releaseStart = 100,
                                         releaseEnd = 110,
                                         releaseInterval = 5,
                                         genMos = c("RR"),
                                         numMos = c(25),
                                         minAge = 16,
                                         maxAge = 24,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)




patchReleases[[1]]$maleReleases <- c(rel1, rel2)

# default migration rate is 0, so no actual batch migration
migrationBatch <- basicBatchMigration(numPatches = sitesNumber)


####################
# Parameter Setup
####################
netPar = NetworkParameters(nPatch = sitesNumber,
                           simTime = 365*3,
                           sampTime = 2,
                           AdPopEQ = c(100, 200, 300),
                           runID = 1L,
                           dayGrowthRate = 1.175,
                           beta = 20, tEgg = 5, tLarva = 6, tPupa = 4, muAd = 0.09)

####################
# Run Simulation
####################
runMPlex(seed = 10,
         numThreads = 1,
         numReps = 5, 
         networkParameters = netPar,
         reproductionReference = reproductionReference,
         initAlleles = AllAlleles,
         patchReleases = patchReleases,
         migrationMale = moveMat,
         migrationFemale = moveMat,
         migrationBatch = migrationBatch,
         outputDirectory = simDir,
         reproductionType = "mPlex_mLoci",
         verbose = FALSE)


####################
# Post Analysis
####################
# setup aggregation key
#  this example sets all genotypes as different
genOI_mLoci_Daisy(outputFile = file.path(aggDir, "0_AggKey.csv"), 
                  genotypes = list(NULL), collapse = c(FALSE))


# aggregate experiment by aggregation key
simAggregation(readDirectory = simDir, writeDirectory = aggDir, 
               simTime = netPar$simTime, sampTime = netPar$sampTime)


# plot for example
plotmPlexSingle(directory = aggDir, whichPatches = NULL, nonZeroGen = FALSE)
plotmPlexMult(directory = aggDir,whichPatches = NULL, nonZeroGen = FALSE, lwd=0.35,alpha=0.75)

## -----------------------------------------------------------------------------
aggKey <- read.csv(file = file.path(aggDir, "0_AggKey.csv"), header = TRUE)
aggKey

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  hold = TRUE,
  fig.width = 7,
  fig.height = 7,
  eval = TRUE
)

## ---- eval=TRUE---------------------------------------------------------------
####################
# Load libraries
####################
library(mPlexCpp)

####################
# Output Folder
####################
outFolder <- "mPlex"

simDir <- file.path(outFolder, "simDir")
aggDir <- file.path(outFolder, paste0("aggDir", 1:6))
for(i in c(outFolder, simDir, aggDir)){ dir.create(i) }

                        
####################
# Landscape
####################
# a 1-node network
sitesNumber <- 1
moveMat <- matrix(data = 1, nrow = sitesNumber, ncol = sitesNumber, byrow = TRUE)


####################
# Inheritance pattern
####################
# generate adult lifetime fitness cost
# Remember, each possible genotype must be specified.
# lociGenos are the alleles at each locus,
lociGenos <- c('HH','HR','HW','RR','RW','WW')

# complete genoes are all combinations at all 3 loci, 216 long for 3 genotypes 
#  and 6 different allele combinations
completeGenos <- expand.grid(lociGenos,lociGenos,lociGenos)

# find genotypes with HH at the third locus
#  We are ignoring all the rest of things, all combinations of 2 at each locus, etc.
#  Because I don't want to write all the code I need to do that.
locus1 <- grep(pattern = "HH", x = completeGenos$Var1, fixed = TRUE, useBytes = TRUE)

# combine loci names into a single list
completeGenos <- do.call(what = paste0, args = completeGenos)

# create named vector of fitness cost
#  since the default fitness is 1, we don't need to specify any genotypes that 
#  don't have a cost. Thus, only adding the ones with a fitness cost on the 
#  first locus.
omega <- setNames(object = rep.int(x = 0.6, times = length(locus1)),
                  nm = completeGenos[locus1])

# 3-locus CRISPR-like drive system
#  first locus is high homing, second is medium, third is poor
reproductionReference <- MakeReference_Multiplex_mLoci(cRateM = c(0.95, 0.8, 0.4),
                                                       hRateM = c(0.90, 0.8, 0.6),
                                                       rRateM = c(1, 1, 1),
                                                       dM = c(0, 0, 0),
                                                       omega = omega)


####################
# Setup Initial genotype ratios
####################
# 3 locus, start completely wild-type
aTypes <- vector(mode = "list", length = 3L) #3 loci in sim
aTypes[[1]]$alleles <- c("W")
aTypes[[1]]$probs <- c(1L)
aTypes[[2]]$alleles <- c("W")
aTypes[[2]]$probs <- c(1L)
aTypes[[3]]$alleles <- c("W")
aTypes[[3]]$probs <- c(1L)


####################
# Setup releases and batch migration
####################
# create Release List
patchReleases = replicate(n = sitesNumber,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL,
                                      matedFemaleReleases = NULL),
                          simplify = FALSE)

# Create release object to pass to patches
rel1 <- basicRepeatedReleases(releaseStart = 50,
                                         releaseEnd = 60,
                                         releaseInterval = 5,
                                         genMos = c("HHWWWW"),
                                         numMos = 25,
                                         minAge = 16,
                                         maxAge = 24,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)

rel2 <- basicRepeatedReleases(releaseStart = 400,
                                         releaseEnd = 410,
                                         releaseInterval = 5,
                                         genMos = c("WWHHWW"),
                                         numMos = 25,
                                         minAge = 16,
                                         maxAge = 24,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)

rel3 <- basicRepeatedReleases(releaseStart = 750,
                                         releaseEnd = 760,
                                         releaseInterval = 5,
                                         genMos = c("WWWWHH"),
                                         numMos = 25,
                                         minAge = 16,
                                         maxAge = 24,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)

patchReleases[[1]]$maleReleases <- c(rel1, rel2, rel3)


# default migration rate is 0, so no actual batch migration
migrationBatch <- basicBatchMigration(numPatches = sitesNumber)


####################
# Parameter Setup
####################
netPar = NetworkParameters(nPatch = sitesNumber,
                           simTime = 365*6,
                           sampTime = 2,
                           AdPopEQ = 250,
                           runID = 1L,
                           dayGrowthRate = 1.175,
                           beta = 20, tEgg = 5, tLarva = 6, tPupa = 4, muAd = 0.09)

####################
# Run Simulation
####################
runMPlex(seed = 10,
         numThreads = 1,
         numReps = 1, 
         networkParameters = netPar,
         reproductionReference = reproductionReference,
         initAlleles = list(aTypes),
         patchReleases = patchReleases,
         migrationMale = moveMat,
         migrationFemale = moveMat,
         migrationBatch = migrationBatch,
         outputDirectory = simDir,
         reproductionType = "mPlex_mLoci",
         verbose = FALSE)


####################
# Aggregation Keys
####################
# Set all genotypes as different
genOI_mLoci_Daisy(outputFile = file.path(aggDir[1], "0_AggKey.csv"), 
                  genotypes = list(NULL, NULL, NULL),
                  collapse = c(FALSE, FALSE, FALSE))

# look at the first locus, ignore everything at the other two
genOI_mLoci_Daisy(outputFile = file.path(aggDir[2], "0_AggKey.csv"), 
                  genotypes = list(NULL, NULL, NULL),
                  collapse = c(FALSE, TRUE, TRUE))

# look at the second locus, ignore other two
genOI_mLoci_Daisy(outputFile = file.path(aggDir[3], "0_AggKey.csv"), 
                  genotypes = list(NULL, NULL, NULL),
                  collapse = c(TRUE, FALSE, TRUE))

# look at the third locus, ignore other two
genOI_mLoci_Daisy(outputFile = file.path(aggDir[4], "0_AggKey.csv"), 
                  genotypes = list(NULL, NULL, NULL),
                  collapse = c(TRUE, TRUE, FALSE))

# H allele at every locus, either 1 H or 2
genOI_mLoci_Daisy(outputFile = file.path(aggDir[5], "0_AggKey.csv"), 
                  genotypes = list(c("HH","HR","HW"), c("HH","HR","HW"), c("HH","HR","HW")),
                  collapse = c(TRUE, TRUE, TRUE))

# H allele at any locus
#  this will aggregate all of the other loci, so the "other" category will 
#  contain anything with an H allele at any locus
genOI_mLoci_Daisy(outputFile = file.path(aggDir[6], "0_AggKey.csv"), 
                  genotypes = list(c("RR","RW","WW"), c("RR","RW","WW"), c("RR","RW","WW")),
                  collapse = c(TRUE, TRUE, TRUE))


####################
# Aggregate
####################
# aggregate experiment by aggregation key
for(x in aggDir){
  simAggregation(readDirectory = simDir, writeDirectory = x, 
               simTime = netPar$simTime, sampTime = netPar$sampTime)
}

## -----------------------------------------------------------------------------
# all genotypes different
aggKey <- read.csv(file = file.path(aggDir[1], "0_AggKey.csv"), header = TRUE)
head(x = aggKey, n = 5)
tail(x = aggKey, n = 5)

plotmPlexSingle(directory = aggDir[1], whichPatches = NULL, nonZeroGen = FALSE)

## -----------------------------------------------------------------------------
# First locus
aggKey <- read.csv(file = file.path(aggDir[2], "0_AggKey.csv"), header = TRUE)
head(x = aggKey, n = 5)
tail(x = aggKey, n = 5)

plotmPlexSingle(directory = aggDir[2], whichPatches = NULL, nonZeroGen = FALSE)

## -----------------------------------------------------------------------------
# second locus
aggKey <- read.csv(file = file.path(aggDir[3], "0_AggKey.csv"), header = TRUE)
head(x = aggKey, n = 5)
tail(x = aggKey, n = 5)

plotmPlexSingle(directory = aggDir[3], whichPatches = NULL, nonZeroGen = FALSE)

## -----------------------------------------------------------------------------
# third locus
aggKey <- read.csv(file = file.path(aggDir[4], "0_AggKey.csv"), header = TRUE)
head(x = aggKey, n = 5)
tail(x = aggKey, n = 5)

plotmPlexSingle(directory = aggDir[4], whichPatches = NULL, nonZeroGen = FALSE)

## -----------------------------------------------------------------------------
# H at every locus
aggKey <- read.csv(file = file.path(aggDir[5], "0_AggKey.csv"), header = TRUE)
aggKey

plotmPlexSingle(directory = aggDir[5], whichPatches = NULL, nonZeroGen = FALSE)

## -----------------------------------------------------------------------------
# H at any locus
aggKey <- read.csv(file = file.path(aggDir[6], "0_AggKey.csv"), header = TRUE)
aggKey

plotmPlexSingle(directory = aggDir[6], whichPatches = NULL, nonZeroGen = FALSE)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

