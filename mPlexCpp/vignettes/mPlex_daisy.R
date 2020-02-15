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
# a 1-node network, ie, just one, well-mixed population
sitesNumber <- 1
moveMat <- as.matrix(1)


####################
# Inheritance pattern
####################
# a 1-locus Daisy behaves like a Mendelian drive, and we add no extra costs to 
#  it here, along with no background mutations
reproductionReference <- MakeReference_DaisyDrive(H = c(0),
                                                  R = c(0),
                                                  S = c(0),
                                                  d = c(0))


####################
# Setup Initial genotype ratios
####################
# 1 locus, start completely wild-type
alleloTypes <- vector(mode = "list", length = 1L) #1 locus
alleloTypes[[1]]$alleles <- c("W")
alleloTypes[[1]]$probs <- c(1)

# replicate so each patch starts the same
AllAlleles <- replicate(n = sitesNumber, expr = alleloTypes, simplify = FALSE)


####################
# Setup releases and batch migration
####################
# create Release List
patchReleases = replicate(n = sitesNumber,
                          expr = list(maleReleases = NULL,
                                      femaleReleases = NULL,
                                      eggReleases = NULL),
                          simplify = FALSE)

# Create release object to pass to patches
patchReleases[[1]]$maleReleases <- basicRepeatedReleases(releaseStart = 25,
                                         releaseEnd = 26,
                                         releaseInterval = 1,
                                         genMos = c("HH"),
                                         numMos = c(10),
                                         minAge = 16,
                                         maxAge = 24,
                                         ageDist = rep(x = 1, times = 24-16+1)/9)

# default migration rate is 0, so no actual batch migration
migrationBatch <- basicBatchMigration(numPatches = sitesNumber)


####################
# Parameter Setup
####################
netPar = NetworkParameters(nPatch = sitesNumber,
                           simTime = 365,
                           sampTime = 2,
                           AdPopEQ = 500,
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
         reproductionType = "DaisyDrive",
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

## ---- eval=TRUE, include = FALSE----------------------------------------------
aggKey <- read.csv(file = file.path(aggDir, "0_AggKey.csv"), header = TRUE, 
                   stringsAsFactors = FALSE)
aggKey

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

