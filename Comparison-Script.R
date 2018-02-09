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
#
###############################################################################

#this is to get the split and aggregated files after a run.
MGDrivE=retrieveOutput(directory="~/Desktop/HOLD/MGDrivE/",genotypes=driveCube$genotypesID)
mPlex=

#This gets the names for each run
Runs <- names(MGDrivE)

#separate male and female
Females <- lapply(X = Runs, FUN = function(run){MGDrivE[[run]]$F})
Males <- lapply(X = Runs, FUN = function(run){MGDrivE[[run]]$M})

#get how many patches and names for each
#get time of simulation and number of genotypes
Patches <- names(Females[[1]])
Time_numGen <- dim(Females[[1]]$Patch1)





#setup return lists for female/male mean/standard deviation
# each list hold every patch
Fmean <- vector(mode = "list", length = length(Patches))
Mmean <- vector(mode = "list", length = length(Patches))
Fsd <- vector(mode = "list", length = length(Patches))
Msd <- vector(mode = "list", length = length(Patches))

for(i in 1:length(Patches)){

  #get patch i from each run
  Fhold <- lapply(X = Females, '[[', i)
  Mhold <- lapply(X = Males, '[[', i)

  #sum all runs of patch i, divide by total number of runs
  Fmean[[i]] <- Reduce(f = "+", x = Fhold)/length(Runs)
  Mmean[[i]] <- Reduce(f = "+", x = Mhold)/length(Runs)

  #This converts the lists into an array and works.
  #I found it on the internet.
  Fsd[[i]] <- apply(array(unlist(Fhold), c(Time_numGen, length(Runs))), c(1,2), sd)
  Msd[[i]] <- apply(array(unlist(Fhold), c(Time_numGen, length(Runs))), c(1,2), sd)

}

#metaData
Note <- "This is data for comparing the master vs FastUnsafe branches.
It was run on 1/26/2018. There are 44 nodes, run for 6000 time steps,
for 50 runs. This data is the mean and sd of each patch over all the repititions."

#Save the data as a nice object
saveRDS(object = list(Fmean = Fmean,
                      Mmean = Mmean,
                      Fsd = Fsd,
                      Msd = Msd,
                      Note = Note),
        file = "~/Desktop/HOLD1/AnalyzedData_UDmel_new_1000pop",
        compress = "gzip")

###############################################################################
#Use this for one run over a network
#This assumes nultiple identical patches that you're averaging over.
###############################################################################

#this is to get the split and aggregated files after a run.
output=retrieveOutput(directory="~/Desktop/HOLD1/UDMel_old/",genotypes=driveCube$genotypesID)
mPlexThing <- readRDS(file = "~/Desktop/HOLD/20180205_Run1_HH_HR_HS_HW_RR_RS_RW_SS_SW_WW.rds")

Patches <- dimnames(mPlexThing$femaleData)[[3]]

#This gets the names for each run
#get time of simulation and number of genotypes
Time_numGen <- dim(mPlexThing$maleData[,,"Patch1"])

Fmean <- apply(X = mPlexThing$femaleData, MARGIN = c(1,2), FUN = mean)
Mmean <- apply(X = mPlexThing$maleData, MARGIN = c(1,2), FUN = mean)

Fsd <- apply(X = mPlexThing$femaleData, c(1,2), sd)
Msd <- apply(X = mPlexThing$maleData, c(1,2), sd)

#metaData
Note <- "This is data for comparing MGDrivE to mPlex. It was run on 2/5/2018.
There are 50 nodes run for 2500 time steps, no migration, so each node is a separate
experiment. 5 releases were done, starting at t=1000 and ending at t=1010, of
10 HH males each time. 5 more releases were dont from t=1200-t=1210 of 10 R2R2
males each time."

#Save the data as a nice object
saveRDS(object = list(Fmean = Fmean,
                      Mmean = Mmean,
                      Fsd = Fsd,
                      Msd = Msd,
                      Note = Note),
        file = "~/Desktop/HOLD/mPlex_50Patch_1000Pop",
        compress = "gzip")


###############################################################################
#Read things in a qqplot them.
###############################################################################

meanOld <- readRDS(file = "~/Desktop/HOLD1/AnalyzedData_UDmel_old_1000pop")
meanNew <- readRDS(file = "~/Desktop/HOLD1/AnalyzedData_UDmel_new_1000pop")

#pick any patch in the list, find the difference, remove extra zeros for better fit.
#too many zeros makes in heavy in the middle, which is good, but looks terrible
#  in a qqplot.
patch <- 31
hold <- meanOld$Fmean[[patch]][,1]-meanNew$Fmean[[patch]][,1]
hold <- hold[!(hold==0)]
qqnorm(hold)
qqline(hold)








mplexSmall <- readRDS(file = "~/Desktop/HOLD/mPlex_50Patch_1000Pop")
MGDrivESmall <- readRDS(file = "~/Desktop/HOLD/MGDrivE_50Patch_1000Pop")



par(mfrow=c(2,1))
matplot(MGDrivESmall$Mmean, main = "MGDrivE: 1000 Pop, 50 Reps",
        xlab = "Days", ylab="Population", type = "l",
        lty = 1, lwd = 4, col = 1:dim(MGDrivESmall$Mmean)[2])
segments(x0 = 1:dim(MGDrivESmall$Mmean)[1], y0 = MGDrivESmall$Mmean-MGDrivESmall$Msd,
         x1 = 1:dim(MGDrivESmall$Mmean)[1], y1 = MGDrivESmall$Mmean+MGDrivESmall$Msd,
         col = rgb(red = 192/255, green = 192/255,blue = 192/255, alpha = 0.5), lty = 1, lwd = 1)
box(lwd=2)
grid()
legend(x = "right", legend = colnames(MGDrivESmall$Mmean),
       col = 1:dim(MGDrivESmall$Mmean)[2], bty = "n", lty = 1,
       lwd = 3, cex = 1)


matplot(mplexSmall$Mmean[,2:11], main = "mPlex: 1000 Pop, 50 Reps",
        xlab = "Days", ylab="Population", type = "l",
        lty = 1, lwd = 4, col = 1:dim(mplexSmall$Mmean[,2:11])[2])
segments(x0 = 1:dim(mplexSmall$Mmean[,2:11])[1], y0 = mplexSmall$Mmean[,2:11]-mplexSmall$Msd[,2:11],
         x1 = 1:dim(mplexSmall$Mmean[,2:11])[1], y1 = mplexSmall$Mmean[,2:11]+mplexSmall$Msd[,2:11],
         col = rgb(red = 192/255, green = 192/255,blue = 192/255, alpha = 0.5), lty = 1, lwd = 1)
box(lwd=2)
grid()
legend(x = "right", legend = colnames(mplexSmall$Mmean[,2:11]),
       col = 1:dim(mplexSmall$Mmean[,2:11])[2], bty = "n", lty = 1,
       lwd = 3, cex = 1)


matplot(mplexSmall$Mmean[,2:11])





matplot(something$Mmean[,2:11])
something$Mmean

matplot(out[,c("PREV","PREV_H","PREV_L")], main = paste("Prevalence Test: cL =", sweep[i], sep = ""),
        xlab = "Years", ylab="Fraction of Population", type = "l", xaxt = 'n',
        lty = 1, lwd = 2, col = 1:3)
axis(side=1, at=seq(0,15*12*10,12*10), labels = seq(0,15,1))
box(lwd=2)
grid()
legend(x = "right", legend = c("PREV","PREV_H","PREV_L"), col = 1:3, bty = "n", lty = 1,
       lwd = 3, cex = 1)









