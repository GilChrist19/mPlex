###############################################################################
#                            ____  __          ______
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/
#                                               /_/   /_/
# Jared Bennett
# jared_bennett@berkeley.edu
# December, 2019
###############################################################################
###################
# References
###################
# lots of reading on discrete-time markov chains
#  Sean coded the initial stuff, then I modified for matrix notation

###################
# Setup
###################
library(expm) # matrix power

distMat <- "Distance matrix generated when creating the landscape.
            Assume it comes from 2DGrid.R or sphericalGrid.R"
probsMat <- "Matrix of movement probabilities from one node to another.
             Assume it comes from 2DGrid.R or sphericalGrid.R"
muAD <- "Adult death rate, will be dependent on the experiment being performed"

# row and column vectors for sums
nNodes <- dim(distMat)[1]
rVec <- matrix(data = 1, nrow = 1, ncol = nNodes)
cVec <- t(rVec)

# expected adult lifespan
nA <- round(1/muAD)

# initial conditions
# This is a vector of the probability of starting at any given node
# It is based on the population size at each node, so popNode/totalPop is the
# probability of starting at any given node.
# This will be experiment specific, but the default is equal
x0 <- rep.int(x = 1, times = nNodes)/nNodes



###################
# DISPLACEMENT
###################
# P(x0 -> anywhere | go through nA steps)
xn <- diag(x = x0) %*% (probsMat %^% nA)

# E[displacement]
rVec %*% (distMat*xn) %*% cVec # sum(distMat * xn)



###################
# CUMULATIVE MOVEMENT
###################
# using vector of current location, set initial distance to 0
currentLoc <- x0
cumDistance <- 0

for(i in 1:nA){
  # P(place @ i+1 | place @ i)
  step <- diag(currentLoc) %*% probsMat
  # sum of the distance traveled for all of those paths
  cumDistance <- cumDistance + rVec %*% (distMat*step) %*% cVec
  # P(new location | all starting locations at previous step)
  currentLoc <- c(rVec %*% step)
}


# # this is the same as above, just tweaked since distance doesn't change
# currentLoc <- x0
# D <- matrix(data = 0, nrow = n, ncol = n)
# for(i in 1:nA){
#   # P(place @ i+1 | place @ i)
#   step <- diag(currentLoc) %*% probsMat
#   # D is now a weight of going from place A to B
#   D <- D + step
#   # P(new location | all starting locations at previous step)
#   currentLoc <- c(rVec %*% step)
# }
#
# # add distance out here
# cumDistance <- rVec %*% (distMat*D) %*% cVec





