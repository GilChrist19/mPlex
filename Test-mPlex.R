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

migration = diag(20) #migration matrix
N = nrow(migration) #number of patches
patchPops = rep(15L,N) #population of eachpatch
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


  # Create release object to pass to patches
patchReleases[[1]]$larvaeReleases <- Release_basicRepeatedReleases(releaseStart = 5,
                                                                 releaseEnd = 10,
                                                                 releaseInterval = 5,
                                                                 genMos = c("HHHHHH", "WHWWWW", "HHWWWW"),
                                                                 numMos = c(100,100,100),
                                                                 minAge = 16,
                                                                 maxAge = 36,
                                                                 ageDist = rep(x = 1, times = 36-16+1)/21)

###############################################################################
# Calculate parameters and initialize network
###############################################################################

    # calculate network parameters, auxiliary function
netPar = Network.Parameters(nPatch = N,simTime = 150L,
                            alleloTypes = AllAlleles,
                            AdPopEQ = patchPops,
                            parallel = FALSE,
                            runID = 1L)

    # initialize network!
network = Network$new(networkParameters = netPar,
                      patchReleases = patchReleases,
                      reproductionType = "mPlex_mLoci",
                      offspringReference = reproductionReference,
                      migrationMale = migration,
                      migrationFemale = migration,
                      directory = directory)



    #reset network
Rprof(interval = 0.01, line.profiling = TRUE)
network$oneRun()
summaryRprof(lines = "both")
network$reset()


###############################################################################
# Post-Run Stuff
###############################################################################
splitOutput(directory = directory)
AnalyzeOutput_mLoci_Daisy(readDirectory = directory,
                          saveDirectory = "~/Desktop/HOLD",
                          genotypes = list(NULL,NULL,NULL),
                          collapse = c(TRUE,TRUE,TRUE))

Run1 <- readRDS(file = "~/Desktop/HOLD/20180118_Run1_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")
Run2 <- readRDS(file = "~/Desktop/HOLD/20180117_Run2_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")
Run3 <- readRDS(file = "~/Desktop/HOLD/20180117_Run3_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")
Run4 <- readRDS(file = "~/Desktop/HOLD/20180117_Run4_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")

par(mfrow=c(3,1))
R1M <- apply(X = Run1$maleData, MARGIN = c(1,2), FUN = mean)[,3]
R1Sd <- apply(X = Run1$maleData, MARGIN = c(1,2), FUN = sd)[,3]

plot(Run1$maleData[,1,1], R1M,pch=18,main = "20 patches with AdEq=5",
     xlab="Time",ylab="Population",xlim=c(0, max(Run1$maleData[,1,1])),
     ylim=c(min(R1M-R1Sd),max(R1M+R1Sd)))
lines(rbind(Run1$maleData[,1,1],Run1$maleData[,1,1],NA),rbind(R1M-R1Sd,R1M+R1Sd,NA))

R1M <- apply(X = Run2$maleData, MARGIN = c(1,2), FUN = mean)[,3]
R1Sd <- apply(X = Run2$maleData, MARGIN = c(1,2), FUN = sd)[,3]

plot(Run2$maleData[,1,1], R1M,pch=18,main = "20 patches with AdEq=10",
     xlab="Time",ylab="Population",xlim=c(0, max(Run2$maleData[,1,1])),
     ylim=c(min(R1M-R1Sd),max(R1M+R1Sd)))
lines(rbind(Run2$maleData[,1,1],Run2$maleData[,1,1],NA),rbind(R1M-R1Sd,R1M+R1Sd,NA))

R1M <- apply(X = Run3$maleData, MARGIN = c(1,2), FUN = mean)[,3]
R1Sd <- apply(X = Run3$maleData, MARGIN = c(1,2), FUN = sd)[,3]

plot(Run2$maleData[,1,1], R1M,pch=18,main = "20 patches with AdEq=20",
     xlab="Time",ylab="Population",xlim=c(0, max(Run2$maleData[,1,1])),
     ylim=c(min(R1M-R1Sd),max(R1M+R1Sd)))
lines(rbind(Run2$maleData[,1,1],Run2$maleData[,1,1],NA),rbind(R1M-R1Sd,R1M+R1Sd,NA))

R1M <- apply(X = Run4$maleData, MARGIN = c(1,2), FUN = mean)[,3]
R1Sd <- apply(X = Run4$maleData, MARGIN = c(1,2), FUN = sd)[,3]

plot(Run2$maleData[,1,1], R1M,pch=18,main = "20 patches with AdEq=50",
     xlab="Time",ylab="Population",xlim=c(0, max(Run2$maleData[,1,1])),
     ylim=c(min(R1M-R1Sd),max(R1M+R1Sd)))
lines(rbind(Run2$maleData[,1,1],Run2$maleData[,1,1],NA),rbind(R1M-R1Sd,R1M+R1Sd,NA))


