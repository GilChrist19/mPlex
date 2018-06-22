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

numPatch <- 100
migration <- matrix(data = runif(numPatch*numPatch), nrow = numPatch, ncol = numPatch)
migration <- migration/rowSums(migration)
patchPops = rep(100L,numPatch)
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







detach("package:mPlexCpp", unload=TRUE)




















###############################################################################
# write split function
###############################################################################
library(data.table)




#' Plot mPlex
#'
#' Plots the analyzed output of mPlex.
#'
#' @usage Plot_mPlex(file, totalPop=TRUE)
#'
#' @param file Path to .gzip file from analyze function
#' @param totalPop Boolean, to plot the total population or not.
#'
#' @details This function plots output from the analyze function. Setting totalPop
#' to FALSE keeps it from plotting the total population.
#'
#' @export
Plot_mPlex <- function(file = NULL, totalPop = TRUE){
  
  #keep old plot parameters to reset later
  oldPar <- par(no.readonly = TRUE)
  
  #Get the data to plot
  Data <- readRDS(file = file)
  
  #Get genotypes and number of patches
  genotypes <- dimnames(Data$maleData)[[2]][-1]
  patches <- dimnames(Data$maleData)[[3]]
  numPatches <- length(patches)
  numGen <- length(genotypes)
  
  if(!totalPop){numGen <- numGen-1}
  
  col <- ggCol_utility(n = numGen)
  
  #setup plot layout
  lmatrix <- matrix(data = 1:(numPatches*3), nrow = numPatches, ncol = 3, byrow = TRUE)
  if(numPatches>1){
    #fill in rest of plot labels
    lmatrix[2:numPatches, c(1,2)] <- matrix(data = 4:(3+2*(numPatches-1)),
                                            ncol = 2, byrow = TRUE)
    #legend gets whole right side
    lmatrix[,3] <- 3
  }
  
  layout(lmatrix, widths = c(3,3,1))
  
  #plot first patch and the legend
  #male
  par(mar = c(2,3,3,1), las = 1, font.lab = 2, font.axis = 2, font.main = 2, cex.main = 1.75)
  matplot(Data$femaleData[,1+(1:numGen), patches[1]], type = "l", lty = 1,
          main = "Female Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(Data$femaleData[,1+(1:numGen), patches[1]])),
          xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
          col = col)
  title(ylab = "Population", line = 2)
  box(lwd = 2)
  grid()
  
  #female
  par(mar = c(2,2,3,1), las = 1)
  matplot(Data$maleData[,1+(1:numGen),patches[1]], type = "l", lty = 1,
          main = "Male Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(Data$maleData[,1+(1:numGen), patches[1]])),
          xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
          col = col)
  mtext(patches[1], side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
  box(lwd = 2)
  grid()
  
  #legend
  par(mar = c(0,0,0,0), font=2)
  plot.new()
  legend(x = "left", legend = genotypes[1:numGen] , col = col,
         bty = "n", bg = "transparent",lty = 1, lwd=3,cex = 1)
  
  
  ##rest of the patches
  if(numPatches>1){
    for(patch in patches[-1]){
      par(mar = c(2,3,1,1), las = 1)
      matplot(Data$femaleData[,1+(1:numGen), patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(Data$femaleData[,1+(1:numGen), patch])),
              xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
              col = col)
      title(ylab = "Population", line = 2)
      box(lwd = 2)
      grid()
      
      par(mar = c(2,2,1,1))
      matplot(Data$maleData[,1+(1:numGen),patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(Data$maleData[,1+(1:numGen), patch])),
              xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
              col = col)
      mtext(patch, side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
      box(lwd = 2)
      grid()
    }#end patch loop
  }#end if
  
  #reset par()
  par(oldPar)
}












Plot_mPlex2 <- function(directory = NULL, totalPop = TRUE){
  
  #keep old plot parameters to reset later
  oldPar <- par(no.readonly = TRUE)
  
  #Get the data to plot
  dirFiles = list.files(path = directory, pattern = ".*\\.csv$")
  
  #test file to get sizes
  testFile <- read.csv(file = file.path(directory, dirFiles[1]), header = TRUE, stringsAsFactors = FALSE)
  
  # get genotypes, num patches, etc
  patches = unique(x = regmatches(x = dirFiles, m = regexpr(pattern = "Patch[0-9]+", text = dirFiles)))
  genotypes <- dimnames(testFile)[[2]][-1]
  numPatches <- length(patches)
  numGen <- length(genotypes)
  
  # create data list holder, and then fill it
  Data = list("maleData"=array(data = 0, dim = c(nrow(testFile),numGen+1,numPatches),dimnames = list(NULL, dimnames(testFile)[[2]], patches)),
              "femaleData"=array(data = 0, dim = c(nrow(testFile),numGen+1,numPatches),dimnames = list(NULL, dimnames(testFile)[[2]], patches)))
  
  
  maleData=array(data = 0, dim = c(nrow(testFile),numGen+1,numPatches),dimnames = list(NULL, dimnames(testFile)[[2]], patches))
  femaleData=array(data = 0, dim = c(nrow(testFile),numGen+1,numPatches),dimnames = list(NULL, dimnames(testFile)[[2]], patches)) 
                   
                   
  for(patch in patches){
    names = grep(pattern = patch, x = dirFiles, fixed = TRUE, value = TRUE)
    femaleData[,,patch] = read.csv(file = file.path(directory, names[1]), header = TRUE, stringsAsFactors = FALSE)
    maleData[,,patch] = read.csv(file = file.path(directory, names[2]), header = TRUE, stringsAsFactors = FALSE)
  }
  
  
  
  
  
  
  
  Data <- readRDS(file = file)
  
  
  
  
  
  #Get genotypes and number of patches
  genotypes <- dimnames(Data$maleData)[[2]][-1]
  #  patches <- dimnames(Data$maleData)[[3]]
  numPatches <- length(patches)
  numGen <- length(genotypes)
  
  
  
  
  
  if(!totalPop){numGen <- numGen-1}
  
  col <- ggCol_utility(n = numGen)
  
  #setup plot layout
  lmatrix <- matrix(data = 1:(numPatches*3), nrow = numPatches, ncol = 3, byrow = TRUE)
  if(numPatches>1){
    #fill in rest of plot labels
    lmatrix[2:numPatches, c(1,2)] <- matrix(data = 4:(3+2*(numPatches-1)),
                                            ncol = 2, byrow = TRUE)
    #legend gets whole right side
    lmatrix[,3] <- 3
  }
  
  layout(lmatrix, widths = c(3,3,1))
  
  #plot first patch and the legend
  #male
  par(mar = c(2,3,3,1), las = 1, font.lab = 2, font.axis = 2, font.main = 2, cex.main = 1.75)
  matplot(Data$femaleData[,1+(1:numGen), patches[1]], type = "l", lty = 1,
          main = "Female Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(Data$femaleData[,1+(1:numGen), patches[1]])),
          xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
          col = col)
  title(ylab = "Population", line = 2)
  box(lwd = 2)
  grid()
  
  #female
  par(mar = c(2,2,3,1), las = 1)
  matplot(Data$maleData[,1+(1:numGen),patches[1]], type = "l", lty = 1,
          main = "Male Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(Data$maleData[,1+(1:numGen), patches[1]])),
          xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
          col = col)
  mtext(patches[1], side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
  box(lwd = 2)
  grid()
  
  #legend
  par(mar = c(0,0,0,0), font=2)
  plot.new()
  legend(x = "left", legend = genotypes[1:numGen] , col = col,
         bty = "n", bg = "transparent",lty = 1, lwd=3,cex = 1)
  
  
  ##rest of the patches
  if(numPatches>1){
    for(patch in patches[-1]){
      par(mar = c(2,3,1,1), las = 1)
      matplot(Data$femaleData[,1+(1:numGen), patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(Data$femaleData[,1+(1:numGen), patch])),
              xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
              col = col)
      title(ylab = "Population", line = 2)
      box(lwd = 2)
      grid()
      
      par(mar = c(2,2,1,1))
      matplot(Data$maleData[,1+(1:numGen),patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(Data$maleData[,1+(1:numGen), patch])),
              xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
              col = col)
      mtext(patch, side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
      box(lwd = 2)
      grid()
    }#end patch loop
  }#end if
  
  #reset par()
  par(oldPar)
}








system.time(AnalyzeOutput_oLocus(readDirectory = "~/Desktop/HOLD/MGDrivE/",
                                 saveDirectory = "~/Desktop/HOLD/mPlex/",
                                 alleles = list(list(c(NULL)),list(c(NULL))),
                                 collapse = list(c(F),c(F)),
                                 numCores = 1))





Rprof(filename = "~/Desktop/HOLD/profiling.out", interval = 0.01, line.profiling = TRUE)
AnalyzeOutput_oLocus(readDirectory = "~/Desktop/HOLD/MGDrivEHOLD/",
                     saveDirectory = "~/Desktop/HOLD/mPlex/",
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













