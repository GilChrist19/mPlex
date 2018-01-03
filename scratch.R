#Things for parameters
numMos <- 10
age <- 10
alleloTypes <- vector(mode = "list", length = 3)

alleloTypes[[1]]$alleles <- c("W","H")
alleloTypes[[1]]$probs <- c(1,0)
alleloTypes[[2]]$alleles <- c("W","H")
alleloTypes[[2]]$probs <- c(1,0)
alleloTypes[[3]]$alleles <- c("W","H")
alleloTypes[[3]]$probs <- c(1,0)




CreateMosquitoes_Distribution_Genotype <- function(numMos, minAge, maxAge, ageDist, aTypes){
  
  population <- vector(mode = "list", length = numMos)
  genotypes <- vector(mode = "list", length = length(aTypes))
  
  for(i in 1:numMos){
    
    #generate genotypes from distribution
    for(locus in 1:length(genotypes)){
      hold <- sample(x = aTypes[[locus]]$alleles, 
                                   size = 2, replace = T, 
                                   prob = aTypes[[locus]]$probs)
      genotypes[[locus]] <- sort(x = hold)
    }
    
    #generate age
    holdAge <- sample(x = minAge:maxAge, size = 1, replace = FALSE, prob = ageDist)
    
    #create new mosquito
    population[[i]] <- Mosquito$new(genotype = paste0(unlist(genotypes), collapse = ""),
                                    age = holdAge)
  }
  
  return(population)
}



## This will pair with the releases function below!
CreateMosquitoes_Defined_Genotype <- function(genMos, numMos, minAge, maxAge, ageDist){
  #genMos is a list of genotypes to relaese
  #numMos is a vector of the number of mosquitoes you want to make, corresponding
  #  to the genotypes of genMos
  
  #minAge, maxAge are the min/max age range. To get a single age, must be length
  # 2 with a c(1,0) vector in ageDist
  #ageDist - probabilities to sample from for age range. must be length
  #  minAge:maxAge
  
  #return list
  population <- vector(mode = "list", length = sum(numMos))

  #external counter
  count = 1L
  
  #loop over each genotype
  for(gen in 1:length(genMos)){
    #loop over number of mosquitoes of that genotype
    for(num in 1:numMos[gen]){
      #generate age
      holdAge <- sample(x = minAge:maxAge, size = 1, replace = FALSE, prob = ageDist)
      
      #create new mosquito
      population[[count]] <- Mosquito$new(genotype = genMos[gen], age = holdAge)
      
      count = count + 1L
    }
  }
  
  return(population)
}


##releases
# release vector must be a list of the mosquitoes to release. Can create using 
#  the functions above. 
Release_basicRepeatedReleases <- function(releaseStart, releaseEnd, releaseInterval, releaseVector, sex="M"){

  # check timing of releases
  if(releaseInterval > (releaseEnd - releaseStart)){
    stop("interval between releases cannot be greater than time between start and end of releases")
  }

  # name and check releaseVector. Initialize release times. Initialize return list
  releaseTimes = seq(from=releaseStart,to = releaseEnd,by = floor(releaseInterval))
  releaseList = vector(mode="list",length=length(releaseTimes))

  # check for male/female/larvae. Fill appropriate list.
  if(sex=="M"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$nuM = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
    }
  } else if(sex=="F"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$nuF = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
    }
  } else if(sex=="L"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$larvae = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
    }
  } else {
    stop(paste0("expected character in 'M','F','L' in argument 'sex', got: ",sex))
  }
  return(releaseList)
}




test <- CreateMosquitoes_Defined_Genotype(numMos = c(25000,25000,25000), genMos = c("AABBCC", "DDDDDD", "DSDF"),minAge = 10, maxAge = 11,ageDist = c(1,0))
listedTest <- replicate(n = 10, expr = test, simplify = FALSE)
listedTest2 <- replicate(n = 10, expr = listedTest, simplify = FALSE)


numNodes = 10
hold <- vector(mode = "list", length = numNodes)
for(i in 1:numNodes){
  for(j in 1:numNodes){
    hold[[i]] <- c(hold[[i]], listedTest2[[i]][[j]])
  }
}




hold <- Reduce(f = c, x = listedTest2)

########################################################################3######
# Figure out intialize and equilibrium stuffs
###############################################################################



any(length(alleloTypes[[1]]) != lengths(alleloTypes))

