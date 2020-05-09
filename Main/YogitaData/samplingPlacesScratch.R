# plot to show grid and sampling locations
# 1: make basic 25x25 square grid.
#    used 2DGrid.R for this
# 2: pick points. I chose 36 below, far away from the boundary

# points chosen for first row
#  idk, chose these to be away from the edge and well spaces
initPts <- seq(from = 5, to = 20, by = 3)

# generate all points to sample from
#  thse are the initial points, put on every row
sampPts <- as.vector(vapply(X = (25*initPts),
                  FUN = function(x){x + initPts},
                  FUN.VALUE = numeric(length = length(initPts)))
                  )

setColors <- rep.int(x = "grey", times = 25*25)
setColors[sampPts] <- "magenta"


png(filename = "~/Desktop/OUTPUT/landscape.png", width = 540, height = 540)

plot(x = latLongs$Var2, y = latLongs$Var1,
     main = "Circles are houses, Magenta are sampling locations",
     xlab = "X",ylab = "Y", type = "p", pch = 19,
     col = setColors)

dev.off()