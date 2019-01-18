###############################################################################
#                            ____  __          ______          
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
#                                               /_/   /_/      
###############################################################################
###############################################################################
#     ______                                 _                     _____           _       __ 
#    / ____/___  ____ ___  ____  ____ ______(_)________  ____     / ___/__________(_)___  / /_
#   / /   / __ \/ __ `__ \/ __ \/ __ `/ ___/ / ___/ __ \/ __ \    \__ \/ ___/ ___/ / __ \/ __/
#  / /___/ /_/ / / / / / / /_/ / /_/ / /  / (__  ) /_/ / / / /   ___/ / /__/ /  / / /_/ / /_  
#  \____/\____/_/ /_/ /_/ .___/\__,_/_/  /_/____/\____/_/ /_/   /____/\___/_/  /_/ .___/\__/  
#                      /_/                                                      /_/           
#
###############################################################################
###############################################################################
# Clean environment and source files
###############################################################################
rm(list=ls());gc()
library(MGDrivE)
library(MGDrivEv2)
#library(mPlexCpp)
library(parallel)

###############################################################################
# Modified quantiles function
###############################################################################
AnalyzeQuantilesMOD <- function(readDirectory, writeDirectory, numFiles, numPatches, FPop){
  
  #get files
  repFiles = list.dirs(path = readDirectory, full.names = TRUE, recursive = FALSE)[1:numFiles]
  patchFiles = lapply(X = repFiles, FUN = list.files, pattern = ".*\\.csv$")
  
  #subset females/males
  malePatches <- lapply(X = patchFiles, FUN = grep, pattern = "ADM", fixed = TRUE, value=TRUE)
  femalePatches <- lapply(X = patchFiles, FUN = grep, pattern = "AF1_Aggregate", fixed = TRUE, value=TRUE)
  
  #generate a list of all patches to run over
  patchList = unique(regmatches(x = patchFiles[[1]],
                                m = regexpr(pattern = "Patch[0-9]+",
                                            text = patchFiles[[1]],
                                            perl = TRUE)))
  
  #read in a file initially to get variables and setup return array
  testFile <- data.table::fread(input = file.path(repFiles[1], patchFiles[[1]][1]),
                                header = TRUE, verbose = FALSE, showProgress = FALSE,
                                logical01 = FALSE, sep = ",", drop = "Time")
  
  #bunch of constants that get used several times
  numReps <- length(repFiles)
  columnNames <- c("Time", names(testFile))
  numRow <- dim(testFile)[1]
  numCol <- dim(testFile)[2]+1
  
  #setup input data holder
  popDataMale <- array(data = 0, dim = c(numRow,  numReps, numCol-1))
  popDataFemale <- array(data = 0, dim = c(numRow, numReps, numCol-1))
  
  #setup output data holder
  outputDataMale <- array(data = 0, dim = c(numRow, numCol, 1),
                          dimnames = list(NULL, columnNames, NULL))
  outputDataFemale <- array(data = 0, dim = c(numRow, numCol, 1),
                            dimnames = list(NULL, columnNames, NULL))
  
  outputDataMale[,1,] <- 1:numRow
  outputDataFemale[,1,] <- 1:numRow
  
  #loop over all patches and do stats.
  for(patch in patchList){
    #get male and female files, all repetitions of this patch
    maleFiles <- vapply(X = malePatches,
                        FUN = grep, pattern = patch, fixed=TRUE, value=TRUE,
                        FUN.VALUE = character(length = 1L))
    
    femaleFiles <- vapply(X = femalePatches,
                          FUN = grep, pattern = patch, fixed=TRUE, value=TRUE,
                          FUN.VALUE = character(length = 1L))
    
    
    #Read in all repetitions for this patch
    for(repetition in 1:numReps){
      popDataMale[ ,repetition, ] <- as.matrix(data.table::fread(input = file.path(repFiles[repetition], maleFiles[repetition]),
                                                                 header = TRUE, verbose = FALSE, showProgress = FALSE,
                                                                 logical01 = FALSE, sep = ",", drop = "Time"))
      popDataFemale[ ,repetition, ] <- as.matrix(data.table::fread(input = file.path(repFiles[repetition], femaleFiles[repetition]),
                                                                   header = TRUE, verbose = FALSE, showProgress = FALSE,
                                                                   logical01 = FALSE, sep = ",", drop = "Time"))
    }
    
    
    #do mean
    for(whichCol in 1:(numCol-1)){
      outputDataMale[ ,whichCol+1,1] <- .rowMeans(x = popDataMale[ , ,whichCol],
                                                  m = numRow, n = numReps)
      outputDataFemale[ ,whichCol+1,1] <- .rowMeans(x = popDataFemale[ , ,whichCol],
                                                    m = numRow, n = numReps)
    }
    
    
    #write output
    maleFileName <- file.path(writeDirectory,
                              file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                        ":", formatC(x = as.integer(gsub(pattern = "Patch", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                        "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                        "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                        "_Mean_Male.csv", fsep = "")
    )
    femaleFileName <- file.path(writeDirectory,
                                file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                          ":", formatC(x = as.integer(gsub(pattern = "Patch", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                          "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                          "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                          "_Mean_Female.csv", fsep = "")
    )
    
    data.table::fwrite(x = as.data.frame(outputDataMale[ , ,1]),
                       file = maleFileName, col.names = TRUE, verbose = FALSE,
                       showProgress = FALSE, logical01 = FALSE, nThread = 1)
    data.table::fwrite(x = as.data.frame(outputDataFemale[ , ,1]),
                       file = femaleFileName, col.names = TRUE, verbose = FALSE,
                       showProgress = FALSE, logical01 = FALSE, nThread = 1)
    
    #do variance
    for(whichCol in 1:(numCol-1)){
      outputDataMale[ ,whichCol+1, ] <- apply(X = popDataMale[ , ,whichCol], MARGIN = 1, FUN = var)
      
      outputDataFemale[ ,whichCol+1, ] <- apply(X = popDataFemale[ , ,whichCol], MARGIN = 1, FUN = var)
      
    }#end loop to calculate quantiles
    
    
    
    #write output
    #file names
    maleFileName <- file.path(writeDirectory,
                              file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                        ":", formatC(x = as.integer(gsub(pattern = "Patch", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                        "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                        "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                        "_Variance_Male.csv", fsep = "")
    )
    femaleFileName <- file.path(writeDirectory,
                                file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                          ":", formatC(x = as.integer(gsub(pattern = "Patch", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                          "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                          "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                          "_Variance_Female.csv", fsep = "")
    )
    
    #write output
    data.table::fwrite(x = as.data.frame(outputDataMale[ , ,1]),
                       file = maleFileName, col.names = TRUE, verbose = FALSE,
                       showProgress = FALSE, logical01 = FALSE, nThread = 1)
    data.table::fwrite(x = as.data.frame(outputDataFemale[ , ,1]),
                       file = femaleFileName, col.names = TRUE, verbose = FALSE,
                       showProgress = FALSE, logical01 = FALSE, nThread = 1)
    
  }#end loop over patches
}#end function

AnalyzeQuantilesMOD2 <- function(readDirectory, writeDirectory, numFiles, numPatches, FPop){
  
  #get files
  repFiles = list.dirs(path = readDirectory, full.names = TRUE, recursive = FALSE)[1:numFiles]
  patchFiles = lapply(X = repFiles, FUN = list.files, pattern = ".*\\.csv$")
  
  #subset females/males
  malePatches <- lapply(X = patchFiles, FUN = grep, pattern = "ADM", fixed = TRUE, value=TRUE)
  femalePatches <- lapply(X = patchFiles, FUN = grep, pattern = "ADF", fixed = TRUE, value=TRUE)
  
  #generate a list of all patches to run over
  patchList = unique(regmatches(x = patchFiles[[1]],
                                m = regexpr(pattern = "Patch[0-9]+",
                                            text = patchFiles[[1]],
                                            perl = TRUE)))
  
  #read in a file initially to get variables and setup return array
  testFile <- data.table::fread(input = file.path(repFiles[1], patchFiles[[1]][1]),
                                header = TRUE, verbose = FALSE, showProgress = FALSE,
                                logical01 = FALSE, sep = ",", drop = c("Time", "Total Pop."))
  
  #bunch of constants that get used several times
  numReps <- length(repFiles)
  columnNames <- c("Time", names(testFile))
  numRow <- dim(testFile)[1]
  numCol <- dim(testFile)[2]+1
  
  #setup input data holder
  popDataMale <- array(data = 0, dim = c(numRow,  numReps, numCol-1))
  popDataFemale <- array(data = 0, dim = c(numRow, numReps, numCol-1))
  
  #setup output data holder
  outputDataMale <- array(data = 0, dim = c(numRow, numCol, 1),
                          dimnames = list(NULL, columnNames, NULL))
  outputDataFemale <- array(data = 0, dim = c(numRow, numCol, 1),
                            dimnames = list(NULL, columnNames, NULL))
  
  outputDataMale[,1,] <- 1:numRow
  outputDataFemale[,1,] <- 1:numRow
  
  #loop over all patches and do stats.
  for(patch in patchList){
    #get male and female files, all repetitions of this patch
    maleFiles <- vapply(X = malePatches,
                        FUN = grep, pattern = patch, fixed=TRUE, value=TRUE,
                        FUN.VALUE = character(length = 1L))
    
    femaleFiles <- vapply(X = femalePatches,
                          FUN = grep, pattern = patch, fixed=TRUE, value=TRUE,
                          FUN.VALUE = character(length = 1L))
    
    
    #Read in all repetitions for this patch
    for(repetition in 1:numReps){
      popDataMale[ ,repetition, ] <- as.matrix(data.table::fread(input = file.path(repFiles[repetition], maleFiles[repetition]),
                                                                 header = TRUE, verbose = FALSE, showProgress = FALSE,
                                                                 logical01 = FALSE, sep = ",", drop = c("Time", "Total Pop.")))
      popDataFemale[ ,repetition, ] <- as.matrix(data.table::fread(input = file.path(repFiles[repetition], femaleFiles[repetition]),
                                                                   header = TRUE, verbose = FALSE, showProgress = FALSE,
                                                                   logical01 = FALSE, sep = ",", drop = c("Time", "Total Pop.")))
    }
    
    
    #do mean
    for(whichCol in 1:(numCol-1)){
      outputDataMale[ ,whichCol+1,1] <- .rowMeans(x = popDataMale[ , ,whichCol],
                                                  m = numRow, n = numReps)
      outputDataFemale[ ,whichCol+1,1] <- .rowMeans(x = popDataFemale[ , ,whichCol],
                                                    m = numRow, n = numReps)
    }
    
    
    #write output
    maleFileName <- file.path(writeDirectory,
                              file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                        ":", formatC(x = as.integer(gsub(pattern = "Patch", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                        "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                        "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                        "_Mean_Male.csv", fsep = "")
    )
    femaleFileName <- file.path(writeDirectory,
                                file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                          ":", formatC(x = as.integer(gsub(pattern = "Patch", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                          "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                          "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                          "_Mean_Female.csv", fsep = "")
    )
    
    data.table::fwrite(x = as.data.frame(outputDataMale[ , ,1]),
                       file = maleFileName, col.names = TRUE, verbose = FALSE,
                       showProgress = FALSE, logical01 = FALSE, nThread = 1)
    data.table::fwrite(x = as.data.frame(outputDataFemale[ , ,1]),
                       file = femaleFileName, col.names = TRUE, verbose = FALSE,
                       showProgress = FALSE, logical01 = FALSE, nThread = 1)
    
    #do variance
    for(whichCol in 1:(numCol-1)){
      outputDataMale[ ,whichCol+1, ] <- apply(X = popDataMale[ , ,whichCol], MARGIN = 1, FUN = var)
      
      outputDataFemale[ ,whichCol+1, ] <- apply(X = popDataFemale[ , ,whichCol], MARGIN = 1, FUN = var)
      
    }#end loop to calculate quantiles
    
    
    
    #write output
    #file names
    maleFileName <- file.path(writeDirectory,
                              file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                        ":", formatC(x = as.integer(gsub(pattern = "Patch", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                        "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                        "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                        "_Variance_Male.csv", fsep = "")
    )
    femaleFileName <- file.path(writeDirectory,
                                file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                          ":", formatC(x = as.integer(gsub(pattern = "Patch", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                          "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                          "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                          "_Variance_Female.csv", fsep = "")
    )
    
    #write output
    data.table::fwrite(x = as.data.frame(outputDataMale[ , ,1]),
                       file = maleFileName, col.names = TRUE, verbose = FALSE,
                       showProgress = FALSE, logical01 = FALSE, nThread = 1)
    data.table::fwrite(x = as.data.frame(outputDataFemale[ , ,1]),
                       file = femaleFileName, col.names = TRUE, verbose = FALSE,
                       showProgress = FALSE, logical01 = FALSE, nThread = 1)
    
  }#end loop over patches
}#end function

eraseDirectoryMOD <- function(directory){
  # check directory exists
  if(!dir.exists(directory)){
    cat("no such directory exists\n")
    return(NULL)
  }
  dirFiles = list.files(path = directory)
  # begin deleting contents
  if(length(dirFiles)>0){
    for(i in dirFiles){
      cat("removing file: ",file.path(directory, i),"\n", sep = "")
      #file.remove(file.path(directory, i))
      unlink(x = file.path(directory, i), recursive = TRUE)
    }
  }
  # end deleting contents
}# end function

###############################################################################
# Setup folder structure and loop values
###############################################################################
DataDir = "~/Desktop/HOLD/ComparisonData"
DataDir2 = "~/Desktop/HOLD/ComparisonData2"
AnalysisDir = "~/Desktop/HOLD/ComparisonAnalysis"
MGDrivEAnalysisDir = "~/Desktop/HOLD/ComparisonAnalysis/MGDrivE"
mPlexAnalysisDir = "~/Desktop/HOLD/ComparisonAnalysis/mPlex"


if(!dir.exists(DataDir)){dir.create(DataDir)}else{eraseDirectoryMOD(DataDir)}
if(!dir.exists(DataDir2)){dir.create(DataDir2)}else{eraseDirectoryMOD(DataDir2)}
if(!dir.exists(AnalysisDir)){dir.create(AnalysisDir)}else{eraseDirectoryMOD(AnalysisDir)}
if(!dir.exists(MGDrivEAnalysisDir)){dir.create(MGDrivEAnalysisDir)}else{eraseDirectoryMOD(MGDrivEAnalysisDir)}
if(!dir.exists(mPlexAnalysisDir)){dir.create(mPlexAnalysisDir)}else{eraseDirectoryMOD(mPlexAnalysisDir)}

FPop <- c(10,50,100,500)
Movement <- list(matrix(data = 1,nrow = 1,ncol = 1),matrix(data = 1/9,nrow = 3,ncol = 3),matrix(data = 1/100,nrow = 10,ncol = 10))
numFileAnalysis <- c(10,25,50,100,252)




FPop <- c(10,20)
Movement <- list(matrix(data = 1,nrow = 1,ncol = 1))
numFileAnalysis <- c(10)

simulationTime=1000 # Number of "days" run in the simulation
repetitions=numFileAnalysis[length(numFileAnalysis)] #252
numCores <- 1#detectCores()/2 #should give number of real cpu cores

# setup cluster
cl=parallel::makePSOCKcluster(names=numCores, outfile = "~/Desktop/HOLD/error.out")
parallel::clusterEvalQ(cl=cl,expr={
  library(MGDrivE)
  library(MGDrivEv2)
})

for(Kernel in Movement){
  for(Population in FPop){
    ###############################################################################
    # Setup Parameters for Network
    ###############################################################################
    # set parameters
    numNodes=NROW(Kernel)
    patchPops=rep(2*Population,numNodes)
    bioParameters=list(betaK=8,tEgg=6,tLarva=11,tPupa=4,popGrowth=1.096,muAd=0.09)
    
    # set cube with random interesting things
    #eM = 0.984, eF = 0.984, rM = 0.01, bM = 0.006, rF = 0.01, bF = 0.006
    driveCube = MGDrivE::Cube_HomingDrive(cM = 0.9, cF = 0.9, dF = 0,
                                          chM = 0.9, crM = 0.9, chF = 0.9, crF = 0.9,
                                          dhF = 0, drF = 0,
                                          omega = c("HH"=0.6, "HB"=0.6,
                                                    "HR"=0.8, "BB"=0.6, "WH"=0.8,
                                                    "WB"=0.8, "RB"=0.8)
                                          )
    
    # movement matrix and batch migration
    migration <- Kernel
    batchMigration <- MGDrivE::basicBatchMigration(batchProbs = c(0.0), sexProbs = c(.5, .5), numPatches = numNodes)
    
    
    # releases, male and female, ignore the egg ones cause I don't care
    patchReleases=replicate(n=numNodes,expr={list(maleReleases=NULL,femaleReleases=NULL,eggReleases = NULL)},simplify=FALSE)
    releasesParameters=list(
      releasesStart=100,releasesNumber=10,releasesInterval=1,
      releaseProportion=Population
    )
    maleReleasesVector1=generateReleaseVector(driveCube=driveCube,releasesParameters=releasesParameters,sex="M")
    
    patchReleases[[1]]$maleReleases=maleReleasesVector1
    
    
    #numCores and repsPerCore have to divide to whole numbers. Basically, pay attention
    # and balance runs so they make sense?
    repsPerCore <- repetitions/numCores
    repStart <- seq(1, repetitions, repetitions/numCores)
    OutputList <- vector(mode = "list", length = numCores)
    
    # loop over number of cores
    for(i in 1:numCores){
      # vector of folder names for each core to output to
      folderName=paste0(DataDir,"/Rep_", formatC(x = repStart[i]:(repStart[i]+repsPerCore-1), width = 4, format = "d", flag = "0"), sep = "")
      # create the folders if they don't exist, clear them if they do
      for(j in 1:repsPerCore){
        if(!dir.exists(folderName[j])){dir.create(folderName[j])}else{eraseDirectoryMOD(folderName[j])}
      }
      # store in list, the vector of folder names and the number of which rep to start on
      OutputList[[i]]$folder=folderName
      OutputList[[i]]$i=repStart[i]
    } # end loop
    
    
    
    # setup parallel cluster and run!
    parallel::clusterExport(
      cl=cl,
      varlist=c("simulationTime","numNodes","bioParameters","patchPops","patchReleases","migration",
                "driveCube","batchMigration", "AnalyzeQuantilesMOD", "Population", "DataDir", "MGDrivEAnalysisDir")
    )

    
    
    parallel::parLapply(cl=cl,X=OutputList,fun=function(x){
      # set up network parameters
      netPar=Network.Parameters(
        runID=x$i,simTime=simulationTime,nPatch=numNodes,
        beta=bioParameters$betaK,muAd=bioParameters$muAd,popGrowth=bioParameters$popGrowth,
        tEgg=bioParameters$tEgg,tLarva=bioParameters$tLarva,tPupa=bioParameters$tPupa,
        AdPopEQ=patchPops
      )
      # set seed
      randomSeed=as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
      # pass parameters and run
      MGDrivEv2::stochastic_multiple(seed = randomSeed,
                                     cubeR = driveCube,
                                     parametersR = netPar,
                                     migrationFemaleR = migration,
                                     migrationMaleR = migration,
                                     migrationBatchR = batchMigration,
                                     releasesR = patchReleases,
                                     output = x$folder,
                                     verbose = FALSE)
      gc()
    })
    
    
    
    ###############################################################################################################
    ############################### POST-ANALYSIS #################################################################
    # create vector of just folder names from the outputList above
    FolderVec <- c(vapply(X = OutputList, FUN = "[[", 1, FUN.VALUE = character(repsPerCore)))
    
    # split output and aggregate females. seems to work in one loop, may have to split
    parallel::parLapply(cl=cl,X=FolderVec,fun=function(x){
      MGDrivE::splitOutput(readDir=x, writeDir = NULL, remFile = TRUE, numCores = 1)
      MGDrivE::aggregateFemales(readDir=x, writeDir = NULL, genotypes = driveCube$genotypesID, remFile = TRUE, numCores = 1)
      gc()
    })
    
    # do analysis in parallel, over number of files to use
    parallel::parLapply(cl=cl, X=numFileAnalysis, fun=function(x){
      AnalyzeQuantilesMOD(readDirectory = DataDir,
                          writeDirectory = MGDrivEAnalysisDir,
                          numFiles = x,
                          numPatches = numNodes,
                          FPop = Population)
      gc()
    })
    
    
    
  }# end loop over populations
}# end loop over kernels

# stop cluster
parallel::stopCluster(cl)

# detach packages
# will pull necessary functions directly from packages to keep MGDrive from 
#  covering things in mPlex
detach("package:MGDrivE", unload=TRUE)
detach("package:MGDrivEv2", unload=TRUE)

###############################################################################
###############################################################################
###############################################################################
# mPLEX
###############################################################################
library(mPlexCpp)

# setup cluster
cl=parallel::makePSOCKcluster(names=numCores)
parallel::clusterEvalQ(cl=cl,expr={
  library(mPlexCpp)
})

for(Kernel in Movement){
  for(Population in FPop){
    ###############################################################################
    # Setup Parameters for Network
    ###############################################################################
    # set parameters
    numNodes=NROW(Kernel)
    patchPops=rep(2*Population,numNodes)
    bioParameters=list(betaK=8,tEgg=6,tLarva=11,tPupa=4,popGrowth=1.096,muAd=0.09)
    
    
    #setup alleles to initiate patches
    alleloTypes <- vector(mode = "list", length = 1L) #1 locus
    alleloTypes[[1]]$alleles <- c("W")
    alleloTypes[[1]]$probs <- c(1L)
    
    AllAlleles <- replicate(n = numNodes, expr = alleloTypes, simplify = FALSE)
    
    
    # This sets up a basic CRISPR drive, with perfect homing and no resistance or backgorund mutation
    reproductionReference <- mPlexCpp::MakeReference_Multiplex_mLoci(H = c(0.992),
                                                           R = c(0.003),
                                                           S = c(0.005),
                                                           d = c(0.00),
                                                           omega = c("HH"=0.6, "HR"=0.6, "HS"=0.8, "RR"=0.6, "HW"=0.8, "RW"=0.8, "RS"=0.8))
    
    # movement matrix and batch migration
    migration <- Kernel
    batchMigration <- basicBatchMigration(batchProbs = 0.0, sexProbs = c(0.5,0.5), numPatches = numNodes)

    
    
    ###############################################################################
    # Release Setup
    ###############################################################################
    
    # create Release List
    patchReleases = replicate(n = numNodes,
                              expr = list(maleReleases = NULL,
                                          femaleReleases = NULL,
                                          eggReleases = NULL),
                              simplify = FALSE)
    
    
    # Create release object to pass to patches
    holdRel <- mPlexCpp::Release_basicRepeatedReleases(releaseStart = 100L,
                                             releaseEnd = 109L,
                                             releaseInterval = 1,
                                             genMos = c("HH"),
                                             numMos = c(Population),
                                             minAge = 16L,
                                             maxAge = 24L,
                                             ageDist = rep(x = 1, times = 24-16+1)/9)
    
    patchReleases[[1]]$maleReleases <- holdRel
    
    ###############################################################################
    # Folder Setup
    ###############################################################################
    #numCores and repsPerCore have to divide to whole numbers. Basically, pay attention
    # and balance runs so they make sense?
    repsPerCore <- repetitions/numCores
    repStart <- seq(1, repetitions, repetitions/numCores)
    OutputList <- vector(mode = "list", length = numCores)
    
    # loop over number of cores
    for(i in 1:numCores){
      # vector of folder names for each core to output to
      folderName=paste0(DataDir,"/Rep_", formatC(x = repStart[i]:(repStart[i]+repsPerCore-1), width = 4, format = "d", flag = "0"), sep = "")
      # create the folders if they don't exist, clear them if they do
      for(j in 1:repsPerCore){
        if(!dir.exists(folderName[j])){dir.create(folderName[j])}else{eraseDirectoryMOD(folderName[j])}
      }
      # store in list, the vector of folder names and the number of which rep to start on
      OutputList[[i]]$folder=folderName
      OutputList[[i]]$i=repStart[i]
    } # end loop
    
    
    
    OutputList2=paste0(DataDir2,"/Rep_", formatC(x = 1:repetitions, width = 4, format = "d", flag = "0"), sep = "")
    for(j in 1:repetitions){
      if(!dir.exists(OutputList2[j])){dir.create(OutputList2[j])}else{eraseDirectoryMOD(OutputList2[j])}
    }

    ###############################################################################
    # Run
    ###############################################################################
    # setup parallel cluster and run!
    parallel::clusterExport(
      cl=cl,
      varlist=c("simulationTime","numNodes","bioParameters","patchPops","patchReleases","migration",
                "AllAlleles","batchMigration", "reproductionReference", "AnalyzeQuantilesMOD2", "Population", "DataDir2", "mPlexAnalysisDir")
    )
    
    
    parallel::parLapply(cl=cl,X=OutputList,fun=function(x){
      # set up network parameters
      netPar = NetworkParameters(nPatch = numNodes,
                                 simTime = simulationTime,
                                 alleloTypes = AllAlleles,
                                 AdPopEQ = patchPops,
                                 runID = x$i,
                                 dayGrowthRate = bioParameters$popGrowth,tEgg = bioParameters$tEgg,
                                 tLarva = bioParameters$tLarva,tPupa = bioParameters$tPupa,muAd = bioParameters$muAd,
                                 beta = bioParameters$betaK)
      # set seed
      randomSeed=as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
      # pass parameters and run
      mPlexCpp::mPlex_runRepetitions(seed = randomSeed,networkParameters = netPar,
                                     reproductionReference = reproductionReference,
                                     patchReleases = patchReleases,
                                     migrationMale = migration,
                                     migrationFemale = migration,
                                     migrationBatch = batchMigration,
                                     output_directory = x$folder,
                                     reproductionType = "mPlex_mLoci",
                                     verbose = FALSE)
      
      gc()
    })
    
    
    
    ###############################################################################################################
    ############################### POST-ANALYSIS #################################################################
    # create vector of just folder names from the outputList above
    FolderVec <- c(vapply(X = OutputList, FUN = "[[", 1, FUN.VALUE = character(repsPerCore)))
    
    FolderList <- as.list(data.frame(rbind(FolderVec, OutputList2), stringsAsFactors = FALSE)) # combine with second list, put into nice format. Weird command
    
    # split output and aggregate females. seems to work in one loop, may have to split
    parallel::parLapply(cl=cl,X=FolderList,fun=function(x){
      mPlexCpp::splitOutput(readDirectory = x[1], numCores = 1)
      mPlexCpp::AnalyzeOutput_mLoci_Daisy(readDirectory = x[1],saveDirectory = x[2],
                                          genotypes = list(NULL),collapse = c(FALSE),
                                          numCores = 1)
      gc()
    })
    
    # do analysis in parallel, over number of files to use
    parallel::parLapply(cl=cl, X=numFileAnalysis, fun=function(x){
      AnalyzeQuantilesMOD2(readDirectory = DataDir2,
                           writeDirectory = mPlexAnalysisDir,
                           numFiles = x,
                           numPatches = numNodes,
                           FPop = Population)
      gc()
    })
    
    
  }# end loop over populations
}# end loop over kernels

# stop cluster
parallel::stopCluster(cl)


# detach packages
detach("package:mPlexCpp", unload=TRUE)

###############################################################################
###############################################################################
# ANALYSIS
###############################################################################

#https://stats.stackexchange.com/questions/184101/comparing-two-histograms-using-chi-square-distance
#https://www.researchgate.net/post/What_is_chi-squared_distance_I_need_help_with_the_source_code
# chi-squared distance metric. Compares vectors, will apply over genotypes
X2_distance <- function(x,y){
  hold <- (x!=0) | (y!=0) # both are zero protection
  return(0.5*sum( (x[hold]-y[hold])^2/(x[hold]+y[hold]) ))
}

# plot function
Plot_Data <- function(meanData, varianceData, FStatValue, fileName){
  
  ## set output
  png(filename = fileName, width=3840, height=2160, units="px")
  
  
  # things used often later
  numPatches <- dim(meanData)[2]
  patches <- 1:numPatches
  timeLength <- dim(meanData)[1]
  
  
  #setup plot layout
  lmatrix <- matrix(data = 1:(numPatches*2), nrow = numPatches, ncol = 2, byrow = TRUE)
  layout(mat = lmatrix)
  
  
  
  #plot first patch and the legend
  #mean
  #  par(mar = c(2,3,3,1), las = 1, font.lab = 2, font.axis = 2, font.main = 2, cex.main = 1.75)
  
  par(mar = c(7,5,6,4), mgp = c(5,1.5,0), las = 1, font.lab = 2, font.axis = 2,
      font.main = 2, cex.main = 5, cex.lab = 3, cex.axis = 3)
  
  
  plot(x = 1:timeLength, y = meanData[,1], type = "l", lty = 1,
       main = "Chi-Squared Distance of Mean",
       ylab = "", lwd=4, xlab = "Time (days)",
       ylim = c(0, max(meanData[,1])),
       xlim = c(0, timeLength), yaxs = "i", xaxs = "i",
       col = "black")
  #title(ylab = "Population", line = 2)
  box(lwd = 4)
  grid()
  
  #variance
  par(mar = c(7,5,6,4), las = 1)
  plot(x = 1:timeLength, y = varianceData[,1], type = "l", lty = 1,
       main = "Euclidian Norm of Variance", ylab = "", lwd=4, xlab = "Time (days)",
       ylim = c(0, max(varianceData[,1])),
       xlim = c(0, timeLength), yaxs = "i", xaxs = "i",
       col = "black")
  abline(h = FStatValue, col = "red", lwd = 4)
  mtext(paste0("Patch ", patches[1]), side = 4, line = 2, las = 0, cex = 3, font = 2)
  box(lwd = 4)
  grid()
  
  
  ##rest of the patches
  if(numPatches>1){
    for(patch in patches[-1]){
      par(mar = c(7,5,1,4), las = 1)
      plot(x = 1:timeLength, y = meanData[,patch], type = "l", lty = 1,
           ylab = "", lwd=4, xlab = "Time (days)",
           ylim = c(0, max(meanData[,patch])),
           xlim = c(0, timeLength), yaxs = "i", xaxs = "i",
           col = "black")
      #title(ylab = "Population", line = 2)
      box(lwd = 4)
      grid()
      
      par(mar = c(7,5,1,4))
      plot(x = 1:timeLength, y = varianceData[,patch], type = "l", lty = 1,
           ylab = "", lwd=4, xlab = "Time (days)",
           ylim = c(0, max(varianceData[,patch])),
           xlim = c(0, timeLength), yaxs = "i", xaxs = "i",
           col = "black")
      abline(h = FStatValue, col = "red", lwd = 4)
      mtext(paste0("Patch ", patch), side = 4, line = 2, las = 0, cex = 3, font = 2)
      box(lwd = 4)
      grid()
    }#end patch loop
  }#end if
  
  
  # close output to file
  dev.off()
  
}



# read in all files
MGDriveFiles <- list.files(path = MGDrivEAnalysisDir, full.names = TRUE)
mPlexFiles <- list.files(path = mPlexAnalysisDir, full.names = TRUE)

# holder objects that don't change size
MGDriveHolder <- matrix(data = 0.0, nrow = simulationTime-1, ncol = 10)
mPlexHolder <- matrix(data = 0.0, nrow = simulationTime-1, ncol = 10)
boolHolder <- matrix(data = FALSE, nrow = simulationTime-1, ncol = 10)
FStatHolder <- matrix(data = 0.0, nrow = simulationTime-1, ncol = 10)
counter <- 1
sortOrder <- c(10,4,9,7,1,3,2,8,6,5)


#subset files by patch, then pop, then rep
pNames <- formatC(x = unlist(lapply(X = Movement, FUN = NROW)), width = 2, format = "d", flag = "0")
fNames <- formatC(x = FPop, width = 3, format = "d", flag = "0")
aNames <- formatC(x = numFileAnalysis, width = 3, format = "d", flag = "0")

for(numPatch in pNames){
  
  # holder objects that don't change size
  ChiSquaredDistance <- matrix(data = 0, nrow = simulationTime-1, ncol = as.integer(numPatch))
  FStat <- matrix(data = 0, nrow = simulationTime-1, ncol = as.integer(numPatch))
  
  
  for(femPop in fNames){
    for(AnalysisReps in aNames){
      
      # set Fstatistic value for this number of reps
      FStatValue <- qf(p = 0.975, df1 = as.integer(AnalysisReps), df2 = as.integer(AnalysisReps))
      FStatValue <- sqrt(sum(rep.int(x = FStatValue, times = 10)^2))
      
      
      # set pattern for files to work on
      PATTERN <- paste0("NumPatches_", numPatch, ":[0-9]{2}_FPop_", femPop, "_NumReps_", AnalysisReps)
      
      # get files to work on
      MGDriveCurrentFiles <- grep(pattern = PATTERN, x = MGDriveFiles, value = TRUE)
      mPlexCurrentFiles <- grep(pattern = PATTERN, x = mPlexFiles, value = TRUE)
      
      
      # read and calculate female means/variance analysis
      MGDriveSexFiles <- grep(pattern = "Female", x = MGDriveCurrentFiles, value = TRUE, fixed = TRUE)
      mPlexSexFiles <- grep(pattern = "Female", x = mPlexCurrentFiles, value = TRUE, fixed = TRUE)
      
      counter <- 1
      for(i in 1:(length(MGDriveSexFiles)/2) ){
        
        # read in means
        MGDriveHolder[] <- matrix(data = scan(file = MGDriveSexFiles[counter], what = numeric(), sep = ",", dec = ".", skip = 1, quiet = TRUE),
                                  ncol = 11, byrow = TRUE)[,-1]
        mPlexHolder[] <- matrix(data = scan(file = mPlexSexFiles[counter], what = numeric(), sep = ",", dec = ".", skip = 1, quiet = TRUE),
                                  ncol = 11, byrow = TRUE)[-1,-1][ ,sortOrder]
        
        # calculate X^2 distance
        for(currentRow in 1:(simulationTime-1)){
          ChiSquaredDistance[currentRow,i] <- X2_distance(x = MGDriveHolder[currentRow,], y = mPlexHolder[currentRow,])
        }
        
        # read in variance
        MGDriveHolder[] <- matrix(data = scan(file = MGDriveSexFiles[counter+1], what = numeric(), sep = ",", dec = ".", skip = 1, quiet = TRUE),
                                  ncol = 11, byrow = TRUE)[,-1]
        mPlexHolder[] <- matrix(data = scan(file = mPlexSexFiles[counter+1], what = numeric(), sep = ",", dec = ".", skip = 1, quiet = TRUE),
                                ncol = 11, byrow = TRUE)[-1,-1][ ,sortOrder]
        
        # calculate variance ratios and euclidian norm
          # these calculate variance ratios
        boolHolder[] <- (MGDriveHolder>mPlexHolder) & (mPlexHolder!=0)
        FStatHolder[boolHolder] <- MGDriveHolder[boolHolder]/mPlexHolder[boolHolder]
        
        boolHolder[] <- (mPlexHolder>MGDriveHolder) & (MGDriveHolder!=0)
        FStatHolder[boolHolder] <- MGDriveHolder[boolHolder]/mPlexHolder[boolHolder]
          
          # euclidian norm of variance ratios
          # this may be inaccurate for large/small values of x
        FStat[,i] <- apply(X = FStatHolder, MARGIN = 1, FUN = function(x){sqrt(sum(x^2))})

        # increment counter
        counter <- counter+2
      }
      # clear holders for re-use
      FStatHolder[] <- 0
      
      #build plot name, and place to store it
      plotName <- paste0(AnalysisDir, "/NPatches_", numPatch, "_FPop_", femPop, "_NumReps_", AnalysisReps, "_Female.png")
      
      # plot
      Plot_Data(meanData = ChiSquaredDistance, varianceData = FStat, FStatValue = FStatValue, fileName = plotName)
      
      
      
      
      

      
      
      # read and calculate male means/variance analysis
      MGDriveSexFiles <- grep(pattern = "Male", x = MGDriveCurrentFiles, value = TRUE, fixed = TRUE)
      mPlexSexFiles <- grep(pattern = "Male", x = mPlexCurrentFiles, value = TRUE, fixed = TRUE)
      
      counter <- 1
      for(i in 1:(length(MGDriveSexFiles)/2) ){
        
        # read in means
        MGDriveHolder[] <- matrix(data = scan(file = MGDriveSexFiles[counter], what = numeric(), sep = ",", dec = ".", skip = 1, quiet = TRUE),
                                  ncol = 11, byrow = TRUE)[,-1]
        mPlexHolder[] <- matrix(data = scan(file = mPlexSexFiles[counter], what = numeric(), sep = ",", dec = ".", skip = 1, quiet = TRUE),
                                ncol = 11, byrow = TRUE)[-1,-1][ ,sortOrder]
        
        # calculate X^2 distance
        for(currentRow in 1:(simulationTime-1)){
          ChiSquaredDistance[currentRow,i] <- X2_distance(x = MGDriveHolder[currentRow,], y = mPlexHolder[currentRow,])
        }
        
        # read in variance
        MGDriveHolder[] <- matrix(data = scan(file = MGDriveSexFiles[counter+1], what = numeric(), sep = ",", dec = ".", skip = 1, quiet = TRUE),
                                  ncol = 11, byrow = TRUE)[,-1]
        mPlexHolder[] <- matrix(data = scan(file = mPlexSexFiles[counter+1], what = numeric(), sep = ",", dec = ".", skip = 1, quiet = TRUE),
                                ncol = 11, byrow = TRUE)[-1,-1][ ,sortOrder]
        
        # calculate variance ratios
          # these calculate variance ratios
        boolHolder[] <- (MGDriveHolder>mPlexHolder) & (mPlexHolder!=0)
        FStatHolder[boolHolder] <- MGDriveHolder[boolHolder]/mPlexHolder[boolHolder]
        
        boolHolder[] <- (mPlexHolder>MGDriveHolder) & (MGDriveHolder!=0)
        FStatHolder[boolHolder] <- MGDriveHolder[boolHolder]/mPlexHolder[boolHolder]
        
          # euclidian norm of variance ratios
        FStat[,i] <- apply(X = FStatHolder, MARGIN = 1, FUN = function(x){sqrt(sum(x^2))})
        
        # increment counter
        counter <- counter+2
      }
      # clear holders for re-use
      FStatHolder[] <- 0
      
      
      #build plot name, and place to store it
      plotName <- paste0(AnalysisDir, "/NPatches_", numPatch, "_FPop_", femPop, "_NumReps_", AnalysisReps, "_Male.png")
      
      # plot
      Plot_Data(meanData = ChiSquaredDistance, varianceData = FStat, FStatValue = FStatValue, fileName = plotName)

    }# end loop over number of files used in analysis
  }# end loop over female population
}# end loop of number of patches










































