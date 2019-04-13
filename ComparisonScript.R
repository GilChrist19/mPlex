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

# DOES NEED DATA.TABLE!!!!!!

#library(MGDrivE)
#library(MGDrivEv2)
#library(mPlexCpp)
#library(parallel)

###############################################################################
# Modified quantiles function
###############################################################################
# This sets all data.table functions to a single thread
data.table::setDTthreads(threads = 1)

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
  allMaleFiles <- list.files(path = readDirectory, pattern = "M_", full.names = TRUE)[1:(numFiles*numPatches)] 
  allFemaleFiles <- list.files(path = readDirectory, pattern = "F_", full.names = TRUE)[1:(numFiles*numPatches)]  
    
  
  #generate a list of all patches to run over
  patchList = unique(regmatches(x = allMaleFiles,
                                m = regexpr(pattern = "Patch_[0-9]+",
                                            text = allMaleFiles,
                                            perl = TRUE)))
  
  #read in a file initially to get variables and setup return array
  testFile <- data.table::fread(input = allMaleFiles[1],
                                header = TRUE, verbose = FALSE, showProgress = FALSE,
                                logical01 = FALSE, sep = ",", drop = c("Time", "Other"))
  
  #bunch of constants that get used several times
  columnNames <- c("Time", names(testFile))
  numRow <- dim(testFile)[1]
  numCol <- dim(testFile)[2]+1
  
  #setup input data holder
  popDataMale <- array(data = 0, dim = c(numRow,  numFiles, numCol-1))
  popDataFemale <- array(data = 0, dim = c(numRow, numFiles, numCol-1))
  
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
    maleFiles <- grep(pattern = patch, x = allMaleFiles, value = TRUE,
                      fixed = TRUE, useBytes = TRUE)
    
    femaleFiles <- grep(pattern = patch, x = allFemaleFiles, value = TRUE,
                        fixed = TRUE, useBytes = TRUE)
    
    
    #Read in all repetitions for this patch
    for(repetition in 1:numFiles){
      popDataMale[ ,repetition, ] <- as.matrix(data.table::fread(input = maleFiles[repetition],
                                                                 header = TRUE, verbose = FALSE, showProgress = FALSE,
                                                                 logical01 = FALSE, sep = ",", drop = c("Time", "Other")))
      popDataFemale[ ,repetition, ] <- as.matrix(data.table::fread(input = femaleFiles[repetition],
                                                                   header = TRUE, verbose = FALSE, showProgress = FALSE,
                                                                   logical01 = FALSE, sep = ",", drop = c("Time", "Other")))
    }
    
    
    #do mean
    for(whichCol in 1:(numCol-1)){
      outputDataMale[ ,whichCol+1,1] <- .rowMeans(x = popDataMale[ , ,whichCol],
                                                  m = numRow, n = numFiles)
      outputDataFemale[ ,whichCol+1,1] <- .rowMeans(x = popDataFemale[ , ,whichCol],
                                                    m = numRow, n = numFiles)
    }
    
    
    #write output
    maleFileName <- file.path(writeDirectory,
                              file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                        ":", formatC(x = as.integer(gsub(pattern = "Patch_", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                        "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                        "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                        "_Mean_Male.csv", fsep = "")
    )
    femaleFileName <- file.path(writeDirectory,
                                file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                          ":", formatC(x = as.integer(gsub(pattern = "Patch_", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
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
                                        ":", formatC(x = as.integer(gsub(pattern = "Patch_", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
                                        "_FPop_", formatC(x = FPop, width = 3, format = "d", flag = "0"),
                                        "_NumReps_", formatC(x = numFiles, width = 3, format = "d", flag = "0"),
                                        "_Variance_Male.csv", fsep = "")
    )
    femaleFileName <- file.path(writeDirectory,
                                file.path("NumPatches_", formatC(x = numPatches, width = 2, format = "d", flag = "0"),
                                          ":", formatC(x = as.integer(gsub(pattern = "Patch_", replacement = "", x = patch, fixed = TRUE))+1, width = 2, format = "d", flag = "0"),
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

X2_distance <- function(x,y){
  # chi-squared distance metric. Compares vectors, will apply over genotypes
  
  #https://stats.stackexchange.com/questions/184101/comparing-two-histograms-using-chi-square-distance
  #https://www.researchgate.net/post/What_is_chi-squared_distance_I_need_help_with_the_source_code
  
  hold <- (x!=0) | (y!=0) # both are zero protection
  return(0.5*sum( (x[hold]-y[hold])^2/(x[hold]+y[hold]) ))
}#end function

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
  
}#end function

###############################################################################
# Setup folder structure and shared values
###############################################################################
DataDir = "~/Desktop/OUTPUT/ComparisonData"
DataDir2 = "~/Desktop/OUTPUT/ComparisonData2"
AnalysisDir = "~/Desktop/OUTPUT/ComparisonAnalysis"
MGDrivEAnalysisDir = "~/Desktop/OUTPUT/ComparisonAnalysis/MGDrivE"
mPlexAnalysisDir = "~/Desktop/OUTPUT/ComparisonAnalysis/mPlex"


if(!dir.exists(DataDir)){dir.create(DataDir)}else{eraseDirectoryMOD(DataDir)}
if(!dir.exists(DataDir2)){dir.create(DataDir2)}else{eraseDirectoryMOD(DataDir2)}
if(!dir.exists(AnalysisDir)){dir.create(AnalysisDir)}else{eraseDirectoryMOD(AnalysisDir)}
if(!dir.exists(MGDrivEAnalysisDir)){dir.create(MGDrivEAnalysisDir)}else{eraseDirectoryMOD(MGDrivEAnalysisDir)}
if(!dir.exists(mPlexAnalysisDir)){dir.create(mPlexAnalysisDir)}else{eraseDirectoryMOD(mPlexAnalysisDir)}




# shared values
FPop <- c(10,50,100,500)
Movement <- list(matrix(data = 1,nrow = 1,ncol = 1),matrix(data = 1/9,nrow = 3,ncol = 3),matrix(data = 1/100,nrow = 10,ncol = 10))
numFileAnalysis <- c(10,25,50,100,252)

simulationTime=200 # Number of "days" run in the simulation
repetitions=numFileAnalysis[length(numFileAnalysis)] #252
numCores <- 2#detectCores()/2 #should give number of real cpu cores
    # NUM CORES MUST BE A DIVISOR OF REPETITIONS!!!!!!!!!

# drive parameters, need to be shared by both
cutting <- 0.9
homing <- 0.9
resistance <- 0.9
bioParameters=list(betaK=20,tEgg=6,tLarva=11,tPupa=4,popGrowth=1.096,muAd=0.09)



#######################################
# Run MGDrivE
#######################################


for(Kernel in Movement){
  for(Population in FPop){
    ########################################
    # Setup Parameters for Network
    ########################################
    # set parameters
    numNodes=NROW(Kernel)
    patchPops=rep(2*Population,numNodes)
    
    # set cube with random interesting things
    #eM = 0.984, eF = 0.984, rM = 0.01, bM = 0.006, rF = 0.01, bF = 0.006
    driveCube = MGDrivE::Cube_HomingDrive(cM = cutting, cF = cutting, dF = 0,
                                          chM = homing, crM = resistance,
                                          chF = homing, crF = resistance,
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
    maleReleasesVector1=MGDrivE::generateReleaseVector(driveCube=driveCube,releasesParameters=releasesParameters,sex="M")
    
    patchReleases[[1]]$maleReleases=maleReleasesVector1
    
    
    
    
    
    # vector of folder names for each core to output to
    folderName=paste0(DataDir,"/Rep_",
                      formatC(x = 1:repetitions, width = 4, format = "d", flag = "0"),
                      sep = "")
    
    # create folders
    for(i in folderName){
      if(!dir.exists(i)){dir.create(i)}else{eraseDirectoryMOD(i)}
    }
    
    
    
    
    
    
    # setup network parameters
    netPar=MGDrivE::Network.Parameters(
        runID=1,simTime=simulationTime,nPatch=numNodes,
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
                                   output = folderName,
                                   verbose = FALSE)
      
    
    ###############################################################################################################
    ############################### POST-ANALYSIS #################################################################
    # split data by patch and aggregate females over mate
    MGDrivEv2::SplitAggregateCpp(readDir = DataDir, writeDir = DataDir,
                                 simTime = simulationTime, numPatch = numNodes,
                                 genotypes = driveCube$genotypesID, remFiles = TRUE)  
      
    # run modified analysis function
    #  fast enough, simple without parallel loop
     for(i in numFileAnalysis){
       AnalyzeQuantilesMOD(readDirectory = DataDir,
                          writeDirectory = MGDrivEAnalysisDir,
                          numFiles = i,
                          numPatches = numNodes,
                          FPop = Population)
     }
    
    
  }# end loop over populations
}# end loop over kernels


#######################################
#######################################
#######################################
# mPLEX
#######################################
#library(mPlexCpp)
library(parallel)



#######################################
# Setup mPlex aggregate key
#######################################
# example function for other uses
#  this has to get written to a special place!!!!!
#mPlexCpp::genOI_mLoci_Daisy(outputFile = "~/Desktop/0_aggKey.csv",genotypes = list(NULL),collapse = c(FALSE))
aggKey <- data.frame("Key"=c("HH","HR","HS","HW","RR","RS","RW","SS","SW","WW"),
                     "Group"=c(1,2,3,4,5,6,7,8,9,10))
      
      

# setup cluster
cl=parallel::makePSOCKcluster(names=numCores)
parallel::clusterEvalQ(cl=cl,expr={
  library(mPlexCpp)
})

for(Kernel in Movement){
  for(Population in FPop){
    ########################################
    # Setup Parameters for Network
    ########################################
    # set parameters
    numNodes=NROW(Kernel)
    patchPops=rep(2*Population,numNodes)
    
    
    #setup alleles to initiate patches
    alleloTypes <- vector(mode = "list", length = 1L) #1 locus
    alleloTypes[[1]]$alleles <- c("W")
    alleloTypes[[1]]$probs <- c(1L)
    
    AllAlleles <- replicate(n = numNodes, expr = alleloTypes, simplify = FALSE)
    
    
    # This sets up a basic CRISPR drive, with perfect homing and no resistance or backgorund mutation
    reproductionReference <- mPlexCpp::MakeReference_Multiplex_mLoci(H = c(cutting),
                                                           R = c(cutting*(1-homing)*resistance),
                                                           S = c(cutting*(1-homing)*(1-resistance)),
                                                           d = c(0.00),
                                                           omega = c("HH"=0.6, "HR"=0.6,
                                                                     "HS"=0.8, "RR"=0.6,
                                                                     "HW"=0.8, "RW"=0.8,
                                                                     "RS"=0.8))
    
    # movement matrix and batch migration
    batchMigration <- mPlexCpp::basicBatchMigration(batchProbs = 0.0, sexProbs = c(0.5,0.5), numPatches = numNodes)

    
    
    ########################################
    # Release Setup
    ########################################
    
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
    
    ########################################
    # Folder Setup
    ########################################
    nReps <- repetitions/numCores
    repStarIDt <- seq(1, repetitions, nReps)
    
    # create folders/clear folders
    for(i in c(DataDir,DataDir2)){
      if(!dir.exists(i)){dir.create(i)}else{eraseDirectoryMOD(i)}
    }
    
    # put aggKey in second folder
    write.csv(x = aggKey, file = file.path(DataDir2,"0_AggKey.csv"), row.names = FALSE)
    
    
    

    ########################################
    # Run
    ########################################
    # setup parallel cluster and run!
    parallel::clusterExport(
      cl=cl,
      varlist=c("simulationTime","numNodes","bioParameters","patchPops","patchReleases",
                "Kernel","AllAlleles","batchMigration","reproductionReference",
                "AnalyzeQuantilesMOD2","Population","DataDir2","mPlexAnalysisDir","DataDir", "nReps")
    )
    
    
    parallel::parLapply(cl=cl,X=repStarIDt,fun=function(x){
      # set up network parameters
      netPar = NetworkParameters(nPatch = numNodes,
                                 simTime = simulationTime,
                                 alleloTypes = AllAlleles,
                                 AdPopEQ = patchPops,
                                 runID = x,
                                 dayGrowthRate = bioParameters$popGrowth,tEgg = bioParameters$tEgg,
                                 tLarva = bioParameters$tLarva,tPupa = bioParameters$tPupa,muAd = bioParameters$muAd,
                                 beta = bioParameters$betaK)
      # set seed
      randomSeed=as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
      # pass parameters and run
      mPlexCpp::mPlex_runRepetitions(seed = randomSeed,
                                     numReps = nReps,
                                     networkParameters = netPar,
                                     reproductionReference = reproductionReference,
                                     patchReleases = patchReleases,
                                     migrationMale = Kernel,
                                     migrationFemale = Kernel,
                                     migrationBatch = batchMigration,
                                     outputDirectory = DataDir,
                                     reproductionType = "mPlex_mLoci",
                                     verbose = FALSE)
      
      gc()
    })
    
    
    ########################################
    # Post-Analysis
    ########################################
    # aggregate by key
    mPlexCpp::SimAggregation(readDirectory = DataDir, writeDirectory = DataDir2, simTime = simulationTime)
    
    # do analysis.
    #  serial, it's fast enough and this is simpler/safer
    for(i in numFileAnalysis){
      AnalyzeQuantilesMOD2(readDirectory = DataDir2,
                           writeDirectory = mPlexAnalysisDir,
                           numFiles = i,
                           numPatches = numNodes,
                           FPop = Population)
     }
    
    
  }# end loop over populations
}# end loop over kernels

# stop cluster
parallel::stopCluster(cl)


# detach packages
detach("package:parallel", unload=TRUE)



###############################################################################
###############################################################################
# ANALYSIS
###############################################################################



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



