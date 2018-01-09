###############################################################################
#                                    ____  _
#                          _ __ ___ |  _ \| | _____  __
#                         | '_ ` _ \| |_) | |/ _ \ \/ /
#                         | | | | | |  __/| |  __/>  <
#                         |_| |_| |_|_|   |_|\___/_/\_\
#
###############################################################################
###############################################################################
# Release Functions
###############################################################################

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

###############################################################################
# Network Initialization Functions
###############################################################################
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

CreateMosquitoes_Eggs <- function(genMos, numMos){
  #genMos is a list of genotypes to relaese
  #numMos is a vector of the number of mosquitoes you want to make, corresponding
  #  to the genotypes of genMos


  #return list
  population <- vector(mode = "list", length = sum(numMos))

  #external counter
  count = 1L

  #loop over each genotype
  for(gen in 1:length(genMos)){
    #loop over number of mosquitoes of that genotype
    for(num in 1:numMos[gen]){
      #create new mosquito
      population[[count]] <- Mosquito$new(genotype = genMos[gen], age = 0)

      count = count + 1L
    }
  }

  #list of new mosquitoes
  return(population)

}

###############################################################################
# Post-processing of Output
###############################################################################
#' Split Output by Patch
#'
#' Split output into multiple files by patches.
#'
#' @param directory directory where output was written to; must not end in path seperator
#'
#' @export
splitOutput <- function(directory){
  dirFiles = list.files(path = directory)

  # for each file read it in
  for(file in dirFiles){
    cat("processing ",file,"\n",sep="")
    fileIn = read.csv(file.path(directory, file))
    # for each file, get all the patches and split into multiple files
    for(patch in unique(fileIn$Patch)){
      patchIn = fileIn[fileIn$Patch==patch,]
      patchName = gsub(pattern = ".csv",replacement = paste0("_Patch",patch,".csv"),x = file)
      write.csv(x = patchIn,file = file.path(directory, patchName),row.names = FALSE)
    }
    cat("removing ",file,"\n",sep="")
    file.remove(file.path(directory, file))
  }
}




###############################################################################
# Random Others
###############################################################################

JaredDirichlet <- function(n=1,alpha){
  Gam <- matrix(0,n,length(alpha))
  for(i in 1:length(alpha)) {Gam[,i] <- rgamma(n,shape=alpha[i])}
  Gam/rowSums(Gam)
}

















