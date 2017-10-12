## Mosquito class functions
## Jared Bennett

# This file defines the functions used in the mosquito class.

###############################################################################

GrowUp <- function(mObject){
  
  #increase age
  mObject$age <- mObject$age+1
  
  #increarse life stage
  cStage <- mObject$stage
  
  if (cStage == "Egg"){
    if (mObject$age >= 7){
      mObject$stage <- "Larva"
    }
  }
  
  if (cStage == "Larva"){
    if (mObject$age >= 14){
      mObject$stage <- "Pupa"
    }
  }
  
  if (cStage == "Pupa"){
    if (mObject$age >= 21){
      mObject$age <- "Adult"
    }
  }
  
}

Death <- function(mObject){
  
  if (mObject$stage == "Adult" && mObject$age >=42){
    return(NULL)
  }
  
}


###############################################################################
## These functions don't get stuff inside function
## Re-organize somehow. 

Mate <- function(allObject){
  #currently remates all females :()
  
  sexList <- vapply(1:length(allObject), function(x){allObject[[x]]$get_sex()}, FUN.VALUE = character(1))

  mIndex <- which(sexList == "M")
  fIndex <- which(sexList == "F")
  
  mGenotype <- vapply(1:length(mIndex), function(x){allObject[[x]]$get_genotype()}, FUN.VALUE = character(1))
  
  mSample <- sample(x = mGenotype, size = length(fIndex), replace = T)
  
  for (i in 1:length(fIndex)){
    allObject[[ fIndex[i] ]]$set_mate(mate = mSample[i])
  }
  
}
