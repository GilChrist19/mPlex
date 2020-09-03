################################################################################
#   _    __          _ __     
#  | |  / /__  _____(_) /_  __
#  | | / / _ \/ ___/ / / / / /
#  | |/ /  __/ /  / / / /_/ / 
#  |___/\___/_/  /_/_/\__, /  
#                    /____/   
#
# Marshall Lab - CKMR Project
# Jared Bennett
# jared_bennett@berkeley.edu
################################################################################
################################################################################
# 20200902
#  create script
#  plot initial landscape, highlight traps
# 20200903
#  generate exponential kernel with 80% daily probability of staying 
#   This will be a "hurdle-exponential", but semantics
#  Uses the Vincenty Ellipsoid method for calculating distance
#  Assumes and adult death rate of 0.09
#   IF THIS CHANGES, NEED TO REGENERATE KERNEL!
#
################################################################################
# plot to show grid and sampling locations

# load data
load("./c2_centroids_info.rds")


# color houses by trap or normal
setColors <- rep.int(x = "grey", times = nrow(c2_centroids_info))
setColors[as.logical(c2_centroids_info$trap)] <- "magenta"

# print png of figure
png(filename = "~/Desktop/OUTPUT/landscape.png", width = 540, height = 540)

plot(x = c2_centroids_info$xcoord, y = c2_centroids_info$ycoord,
     main = "Control Area 2, Fresno",
     xlab = "X",ylab = "Y", type = "p", pch = 19,
     col = setColors)

dev.off()

################################################################################
# Generate synthetic kernel

# load data
load("./c2_centroids_info.rds")

# calculate distance between points
#  Vincenty ellipsoid method
vDist <- CKMR::calcVinEll(latLongs = as.matrix(c2_centroids_info[,c("ycoord","xcoord")]))

# rate
#  This describes the average movement, given exponential distribution.
#  We want an average distance of 112.5 meters, so the movement rate is the inverse 
#  of that.
rate <- 1/112.5

# p0
#  Given an 80% probability to stay in one place over their entire life, 
#  we calculate the chance of staying daily by raising the power to the life-time 
#  of the mosquito, 
#  Using a deathrate of 0.09
p0 <- 0.8^0.09 

# calculate kernel
c2_kernel_exp80 <- CKMR::calcHurdleExpKernel(distMat = vDist, rate = rate, p0 = p0)

# store for loading into the simulation
save(c2_kernel_exp80, file = "./c2_kernel_exp80.rds", compress = "xz", compression_level = 9)












