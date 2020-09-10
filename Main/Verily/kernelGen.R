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
# 20200909
#  New plot, labeling the traps with their number this time.
#
################################################################################
# plot to show grid and sampling locations

# load data
load("./c2_centroids_info.rds")


####################
# Colored Dots Plot
####################
# color houses by trap or normal
setColors <- rep.int(x = "grey", times = nrow(c2_centroids_info))
setColors[as.logical(c2_centroids_info$trap)] <- "magenta"

# print png of figure
png(filename = "~/Desktop/OUTPUT/landscape.png", width = 540, height = 540)

plot(x = c2_centroids_info$xcoord, y = c2_centroids_info$ycoord,
     main = "Control Area 2, Fresno",
     xlab = "X", ylab = "Y", type = "p", pch = 19,
     col = setColors)

dev.off()


####################
# Numbered Dots Plot
####################
# label the traps with their node number

# get indices for plotting
ind0 <- which(c2_centroids_info$trap == 0)
ind1 <- which(c2_centroids_info$trap == 1)

# print png of figure
png(filename = "~/Desktop/OUTPUT/landscapeLabeled.png", width = 540, height = 540)

plot(x = c2_centroids_info$xcoord[ind0], y = c2_centroids_info$ycoord[ind0],
     main = "Control Area 2, Fresno",
     xlab = "X", ylab = "Y", type = "p", pch = 19,
     col = 'grey')
text(x = c2_centroids_info$xcoord[ind1], y = c2_centroids_info$ycoord[ind1],
     label = as.character(ind1), col = 'magenta',
     offset = 0, font = 2, cex = 0.75)

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
#  Given an 80% probability to stay in one place per day, 
#  we calculate the chance of staying daily by raising the power to the life-time 
#  of the mosquito, 
#  Don't need the deathrate, when daily probs is given
p0 <- 0.8

# calculate kernel
c2_kernel_exp80 <- CKMR::calcHurdleExpKernel(distMat = vDist, rate = rate, p0 = p0)

# store for loading into the simulation
save(c2_kernel_exp80, file = "./c2_kernel_exp80.rds", compress = "xz", compression_level = 9)












