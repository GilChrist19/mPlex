## This is a run script for my homework
## Jared Bennett
## Stat 243

library(R6)
source("mFunc.R")
source("mClass.R")




# create holder for all my play objects
allObjects <- vector(mode = "list", length = 10)

# fill my holder with all of my play objects
for (x in 1:length(allObjects)){
  hold <- runif(n = 1, min = 0, max = 1)
  sex <- c("M", "F")
  allObjects[[x]] <- Mosquito$new(genotype = "WW", age = 10, stage = "Adult")
  allObjects[[x]]$set_sex(sex = sex[round(hold)+1])
  
}





allObjects[[1]]$grow_up()




