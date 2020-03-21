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
#  https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
#
# spherical coordinates notes
#  r = spherical radius, = 1 here
#  theta = azimuthal angle in x/y plane, 0 <= theta < 2pi from x axis
#  phi = zenith angle, 0 <= phi <= pi ( 0 is pointing up, pi points down)
#

###################
# Distance Function
###################
# great distance function, vectorized over the second point
#  ie, it calculates the distance from an initial point to all other points
#  as warned on wikipedia (https://en.wikipedia.org/wiki/Great-circle_distance)
#  this has some numerical issues. Not bad, but not perfect is all.
#  will use haversine distance instead
gcDist <- function(point1, points2){
  return(acos(sin(point1[1])*sin(points2[ ,1]) +
         cos(point1[1])*cos(points2[ ,1])*cos(abs(points2[ ,2] - point1[2])) )
         )
}

###################
# Setup
###################
# setup size of landscape
nPoints <- 9
rate <- 1
p0 <- 0.9
myFile <- "~/Desktop/sphericalMoveMat.csv"

###################
# point locations
###################
# phi,theta locations of all points
indices = 0:(nPoints-1) + 0.5
phi = acos(1 - 2*indices/nPoints)
theta = pi * (1 + 5^0.5) * indices

# combine lat/longs for all points
#  phi is technically an azimuthal angle, which is not the same as a latitude
#  The pi/2-phi converts azimuthal to latitude
latLongs <- cbind(pi/2 - phi, theta)

# cartesian coordinates if desired
#  would make plotting easier to check point distribution
# x <- cos(theta) * sin(phi)
# y <- sin(theta) * sin(phi)
# z <- cos(phi)

###################
# Distance Calcs
###################
# distance between all points, output is matrix
#  first is the mathematical function, with some small numerical issues
#  second is using the haversine function, more stable
#    It's written for degrees, not radians, so convert input
#  third is the vincenty formula for a sphere, more accurate
#    It's written for degrees, not radians, so convert input
#  There's effectively no difference between the last two approaches.
#distMat <- apply(X = latLongs, MARGIN = 1, FUN = gcDist, points2 = latLongs)
#distMat <- mPlexCpp::calcHaversine(latLongs = latLongs*(180/pi), r = 1)
distMat <- mPlexCpp::calcVinSph(latLongs = latLongs*(180/pi), r = 1)

###################
# move Probs
###################
probsMat <- mPlexCpp::calcHurdleExpKernel(distMat = distMat, rate = rate, p0 = p0)

###################
# write output
###################
write.table(x = probsMat, file = myFile,
            sep = ",", row.names = FALSE, col.names = FALSE)

















