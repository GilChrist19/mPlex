## This is a run script for my homework
## Jared Bennett
## Stat 243

library(R6)
source("mFunc.R")
source("mClass.R")


#how to set up intelligently?
Need to make possible genotypes.
Need to make homing/n list stuff
  Applies to all things and only needs calculated once
Fill adults/ eggs, larva, pupa with probable numbers and 
  appropriate age distribution, and genotypes.
  Base this off of chosen adult size




# create holder for all my play objects
allObjects <- vector(mode = "list", length = 10)

# fill my holder with all of my play objects
for (x in 1:length(allObjects)){
  hold <- runif(n = 1, min = 0, max = 1)
  sex <- c("M", "F")
  allObjects[[x]] <- Mosquito$new(genotype = "WW", age = 10, stage = "Adult")
  allObjects[[x]]$set_sex(sex = sex[round(hold)+1])
  
}



Mate(allObject = allObjects)





Eggs <- new.env(hash = TRUE, parent = parent.frame())

assign(x = "1", value = Mosquito$new(genotype = "WW", age = 0, stage = "Egg"), pos = Eggs, inherits = F)

#x has to be a string. Has to be unique. 
#how to move that object out, use assign? then how to know what names to use/



Larva <- new.env()
Pupa <- new.env()
Adults <- new.env()








