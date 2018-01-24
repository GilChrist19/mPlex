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

migration = diag(5) #migration matrix
N = nrow(migration) #number of patches
patchPops = rep(10L,N) #population of eachpatch
directory <- "~/Desktop/HOLD"

    #setup alleles to initiate patches
alleloTypes <- vector(mode = "list", length = 1) #3 loci
alleloTypes[[1]]$alleles <- c("W")
alleloTypes[[1]]$probs <- c(1)
# alleloTypes[[2]]$alleles <- c("W","H")
# alleloTypes[[2]]$probs <- c(1,0)
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
reproductionReference <- MakeReference_Multiplex_mLoci(H = c(0.98),
                                                       R = c(0.0001),
                                                       S = c(0.00003),
                                                       d = c(0))

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
patchReleases[[1]]$maleReleases <- Release_basicRepeatedReleases(releaseStart = 50,
                                                                 releaseEnd = 70,
                                                                 releaseInterval = 2,
                                                                 genMos = c("HH"),
                                                                 numMos = c(10),
                                                                 minAge = 16,
                                                                 maxAge = 24,
                                                                 ageDist = rep(x = 1, times = 24-16+1)/9)

###############################################################################
# Calculate parameters and initialize network
###############################################################################

    # calculate network parameters, auxiliary function
netPar = Network.Parameters(nPatch = N,simTime = 1000L,
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
                          genotypes = list(NULL),
                          collapse = c(F))

Run1 <- readRDS(file = "~/Desktop/HOLD/20180119_Run1_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")
Run2 <- readRDS(file = "~/Desktop/HOLD/20180119_Run2_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")
Run3 <- readRDS(file = "~/Desktop/HOLD/20180119_Run3_(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW)(HH|HR|HS|HW|RR|RS|RW|SS|SW|WW).rds")
Run4 <- readRDS(file = "~/Desktop/HOLD/20180119_Run4_HH_HR_HS_HW_RR_RS_RW_SS_SW_WW.rds")
test <- readRDS(file = "~/Desktop/HOLD/20180120_Run4_HH_HR_HS_HW_RR_RS_RW_SS_SW_WW.rds")



file <- "~/Desktop/HOLD/20180120_Run4_HH_HR_HS_HW_RR_RS_RW_SS_SW_WW.rds"
library(viridisLite)

Plot_mPlex <- function(file = NULL){

  #keep old plot parameters to reset later
  oldPar <- par(no.readonly = TRUE)

  #Get the data to plot
  Data <- readRDS(file = file)

  #Get genotypes and number of patches
  genotypes <- dimnames(Data$maleData)[[2]][-1]
  patches <- dimnames(Data$maleData)[[3]]

  #setup plot layout
  lmatrix <- matrix(data = 1:(length(patches)*3), nrow = length(patches), ncol = 3, byrow = TRUE)
  if(length(patches)>1){
    #fill in rest of plot labels
    lmatrix[2:length(patches), c(1,2)] <- matrix(data = 4:(3+2*(length(patches)-1)),
                                                 ncol = 2, byrow = TRUE)
    #legend gets whole right side
    lmatrix[,3] <- 3
  }

  layout(lmatrix, widths = c(3,3,1))

  #plot first patch and the legend
  #male
  par(mar = c(2,3,3,1), las = 1, font.lab = 2, font.axis = 2, font.main = 2, cex.main = 1.75)
  matplot(Data$femaleData[,1+(1:length(genotypes)), patches[1]], type = "l", lty = 1,
          main = "Female Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(Data$femaleData[,1+(1:length(genotypes)), patches[1]])),
          xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
          col = 1:length(genotypes))
  title(ylab = "Population", line = 2)
  box(lwd = 2)
  grid()

  #female
  par(mar = c(2,2,3,1), las = 1)
  matplot(Data$maleData[,1+(1:length(genotypes)),patches[1]], type = "l", lty = 1,
          main = "Male Mosquitoes", ylab = "", lwd=2,
          ylim = c(0, max(Data$maleData[,1+(1:length(genotypes)), patches[1]])),
          xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
          col = 1:length(genotypes))
  mtext(patches[1], side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
  box(lwd = 2)
  grid()

  #legend
  par(mar = c(0,0,0,0), font=2)
  plot.new()
  legend(x = "left", legend = genotypes , col = 1:length(genotypes),
         bty = "n", bg = "transparent",lty = 1, lwd=3,cex = 1)


  ##rest of the patches
  if(length(patches)>1){
    for(patch in patches[-1]){
      par(mar = c(2,3,1,1), las = 1)
      matplot(Data$femaleData[,1+(1:length(genotypes)), patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(Data$femaleData[,1+(1:length(genotypes)), patch])),
              xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
              col = 1:length(genotypes))
      title(ylab = "Population", line = 2)
      box(lwd = 2)
      grid()

      par(mar = c(2,2,1,1))
      matplot(Data$maleData[,1+(1:length(genotypes)),patch], type = "l", lty = 1,
              ylab = "", lwd=2,
              ylim = c(0, max(Data$maleData[,1+(1:length(genotypes)), patch])),
              xlim = c(0, dim(Data$maleData)[1]), yaxs = "i", xaxs = "i",
              col = 1:length(genotypes))
      mtext(patch, side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
      box(lwd = 2)
      grid()
    }#end patch loop
  }#end if

  #reset par()
  par(oldPar)
}





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






