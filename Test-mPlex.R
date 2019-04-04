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
library(mPlexRCpp)

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



    #reset network
# library(foreach)
# library(parallel)
# library(iterators)
# library(doParallel)
#
#
# cl = parallel::makeForkCluster(nnodes = parallel::detectCores()-4L)
# parallel::clusterSetRNGStream(cl=cl,iseed=NULL)
# doParallel::registerDoParallel(cl)

#Rprof(interval = 0.01, line.profiling = TRUE)
set.seed(seed = 42)
network$oneRun()
#summaryRprof(lines = "both")
network$reset()

# parallel::stopCluster(cl)
# rm(cl);gc()


###############################################################################
# Post-Run Stuff
###############################################################################
splitOutput(directory = directory)
AnalyzeOutput_mLoci_Daisy(readDirectory = directory,
                          saveDirectory = "~/Desktop/OUTPUT/mPlex/",
                          fileName = "ProfileTest",
                          genotypes = list(c(NULL),c(NULL)),
                          collapse = c(F,T))

Plot_mPlex(file = "/home/jared/Desktop/OUTPUT/mPlex/20190404_Run1_ProfileTest.rds", totalPop = T)











Run1 <- readRDS(file = "~/Desktop/HOLD/20180119_Run1_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")
Run2 <- readRDS(file = "~/Desktop/HOLD/20180119_Run2_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")
Run3 <- readRDS(file = "~/Desktop/HOLD/20180119_Run3_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")
Run4 <- readRDS(file = "~/Desktop/HOLD/20180119_Run4_HH_HR_HS_HW_RR_RS_RW_SS_SW_WW.rds")
test <- readRDS(file = "~/Desktop/HOLD/20180202_Run1_HH_HR_HS_HW_RR_RS_RW_SS_SW_WW.rds")



file <- "~/Desktop/HOLD/20180120_Run4_HH_HR_HS_HW_RR_RS_RW_SS_SW_WW.rds"
library(viridisLite)






RUN  <- Run4

#male
R1M <- apply(X = RUN$maleData, MARGIN = c(1,2), FUN = mean)[,3]
R1Sd <- apply(X = RUN$maleData, MARGIN = c(1,2), FUN = sd)[,3]

plot(RUN$maleData[,1,1], R1M,pch=18,main = "20 patches with AdEq=50",
     xlab="Time",ylab="Male Population",xlim=c(0, max(RUN$maleData[,1,1])),
     ylim=c(min(R1M-R1Sd),max(R1M+R1Sd)), col = "blue")
lines(rbind(RUN$maleData[,1,1],RUN$maleData[,1,1],NA),rbind(R1M-R1Sd,R1M+R1Sd,NA), col = "skyblue1")

#female
R1M <- apply(X = RUN$femaleData, MARGIN = c(1,2), FUN = mean)[,3]
R1Sd <- apply(X = RUN$femaleData, MARGIN = c(1,2), FUN = sd)[,3]

plot(RUN$femaleData[,1,1], R1M,pch=18,main = "20 patches with AdEq=50",
     xlab="Time",ylab="Female Population",xlim=c(0, max(RUN$femaleData[,1,1])),
     ylim=c(min(R1M-R1Sd),max(R1M+R1Sd)), col = "magenta1", new = TRUE)
lines(rbind(RUN$femaleData[,1,1],RUN$femaleData[,1,1],NA),rbind(R1M-R1Sd,R1M+R1Sd,NA), col = "lightpink")






detach("package:mPlexR", unload=TRUE)
detach("package:mPlexRCpp", unload=TRUE)
