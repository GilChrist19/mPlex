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
# MGDrivE Crispr 2 resistant alleles against mPlex
###############################################################################
###############################################################################
# Clean environment and source files
###############################################################################
rm(list=ls());gc()
library(MGDrivE)

###############################################################################
# Setup Parameters for Network
###############################################################################

MGDrivE.Setup(stochasticityON = TRUE)
migration = diag(50L)
N = nrow(migration)
patchPops = rep(500L*2L,N)
directory = "~/Desktop/HOLD/MGDrivE/"


###############################################################################
# reproduction Setup
###############################################################################

crisprCube = Cube_CRISPRTwoResistantAllele(e = 0.98, p1 = 0.001, p2 = 0.0003)

###############################################################################
# Release Setup
###############################################################################

patchReleases = replicate(n = N,
                          expr = {list(maleReleases = NULL,
                                       femaleReleases = NULL)},
                          simplify = FALSE)

holdRel = Release_basicRepeatedReleases(genotypes = crisprCube$genotypesID,
                                        releaseStart = 1000L,
                                        releaseEnd = 1010L,
                                        releaseInterval = 2L,
                                        releaseVector = c("HH"=10L,"Hh"=0L,"HR1"=0L,
                                                          "HR2"=0L,"hh"=0L,"hR1"=0L,
                                                          "hR2"=0L,"R1R1"=0L,"R1R2"=0L,
                                                          "R2R2"=0L),
                                        sex = "M")

holdRel2 = Release_basicRepeatedReleases(genotypes = crisprCube$genotypesID,
                                        releaseStart = 1200L,
                                        releaseEnd = 1210L,
                                        releaseInterval = 2L,
                                        releaseVector = c("HH"=0L,"Hh"=0L,"HR1"=0L,
                                                          "HR2"=0L,"hh"=0L,"hR1"=0L,
                                                          "hR2"=0L,"R1R1"=0L,"R1R2"=0L,
                                                          "R2R2"=10L),
                                        sex = "M")

for(i in seq(1,N, 1)){
  patchReleases[[i]]$maleReleases <- c(holdRel, holdRel2)
}

###############################################################################
# Calculate parameters and initialize network
###############################################################################

  #network parameters
netPar = Network.Parameters(nPatch = N,runID = 1L,popGrowth = 1.1,simTime = 2500L,beta = 32L, AdPopEQ = patchPops)

  # initialize network
network = Network$new(networkParameters = netPar,
                      driveCube = crisprCube,
                      patchReleases = patchReleases,
                      migrationMale = migration,
                      migrationFemale = migration,
                      directory = directory)




network$oneRun()
network$reset()
# network$oneRun()
# network$reset()

splitOutput(directory = directory)
aggregateFemales(directory = directory,genotypes = crisprCube$genotypesID,remove = FALSE)

output = retrieveOutput(directory = directory,genotypes = crisprCube$genotypesID)

# plot the first run:
library(reshape2)
library(ggplot2)
datF = reshape2::melt(output$Run1$F)
datM = reshape2::melt(output$Run1$M)

datJoint = rbind(data.frame(datF,sex=rep("F",nrow(datF))),data.frame(datM,sex=rep("M",nrow(datM))))
ggplot(data = datJoint) +
  geom_line(aes(x=Var1,y=value,color=Var2)) +
  facet_grid(L1~sex,scales="free") +
  theme_bw()

detach("package:MGDrivE", unload=TRUE)






