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


functions
  age(critter, age factor){
    This function will age every critter, with a poisson distribution on their
    current life stage to move them to the next.
    Should also kill for too old?
  }
  
  stage move(critter, next environment){
    get life stage of critter, then move it to the next stage.
    This will only be used on ELP stages, not on adult
  }
  
  clean_names(environment){
    get the list of objects in the environment, find whichi ones are  empty, then
    shift things to fill them in.
    
    get list of names, move objects up in list maybe?
      rm()/remove()
      ls()/objects()

    or, just get the names that are gone, then call new objects by those names? 
      max(list of names) then add there for more names

  }
  
  Density effects(factors, environment){
    These are specific to each stage. apply to whole environment, randomly kill.
    Maybe apply this before we clean environment names?
  }
  
  Mate (adult?){
    mating, done by female that isnt already mated. need list of all males currently
    in environment. So maybe apply this to whole environment? need to only mate
    females that arent already mated, so make part of adult move function?
  }
  
  Offspring(){
    This applies to all adult females. They are already mated, so have mate 
    genotype.
    
    if (male)
      break
    if (female)
      do offspring
    
    push offspring into egg environment. 
    
    
  }
  
  


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
