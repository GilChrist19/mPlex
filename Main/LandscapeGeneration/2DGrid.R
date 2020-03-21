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

###################
# Distance Function
###################
# Euclidean distance function, vectorized over the second point
#  ie, it calculates the distance from an initial point to all other points
eDist <- function(point1, points2){
  return(sqrt((points2[ ,1] - point1[1])^2 +
                (points2[ ,2] - point1[2])^2)
         )
}

###################
# Setup
###################
# setup dimensions of landscape
nRow <- 30
nCol <- 30
rate <- MGDrivE::kernels$exp_rate
p0 <- 0.965
myFile <- "~/Desktop/gridMoveMat.csv"

###################
# point locations
###################
latLongs <- expand.grid(1:nRow, 1:nCol)[ ,c(2,1)]

###################
# Distance Calcs
###################
# distance between all points, output is matrix
distMat <- apply(X = latLongs, MARGIN = 1, FUN = eDist, points2 = latLongs)

# faster euclidean distance function
# sqrt(outer(latLongs[ ,1], latLongs[ ,1], FUN = (x - y)^2) +
#   outer(latLongs[ ,2], latLongs[ ,2], FUN = function(x,y){(x-y)^2})
# )

###################
# Distance Spread
###################
# The current distance matrix is done using unit distance between nodes.
# However, the rate from MGDrivE::kernels$exp_rate is calculated in meters per day,
#  with the median value being ~ 55 meters per day. The probability of staying, per 
#  whole life, is (1-p)^(1/mu), where mu is adult lifespan and p is probability of 
#  moving per day. Given this, the probability to leave during lifespan is 
#  1-(1-p)^(1/mu). Assuming movement is exponential, the distance traveled, given 
#  that you move, is exp(-lambda(n/4)x), where lambda is the rate, n is number of 
#  nodes you want to travel, and x is distance between nodes. n is chosen by user, 
#  I chose 15 for the 30x30 grid.

# this is very specific to the experiment chosen, follow the above paragraph and resolve
# there is a 1% chance that you move 4 times as an adult, and a 1% chance that 
# you move 15 nodes in one jump
distMat <- distMat * 16.6062


###################
# move Probs
###################
probsMat <- mPlexCpp::calcHurdleExpKernel(distMat = distMat, rate = rate, p0 = p0)

###################
# write output
###################
write.table(x = probsMat, file = myFile,
            sep = ",", row.names = FALSE, col.names = FALSE)
















