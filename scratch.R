###############################################################################
#                                    ____  _           
#                          _ __ ___ |  _ \| | _____  __
#                         | '_ ` _ \| |_) | |/ _ \ \/ /
#                         | | | | | |  __/| |  __/>  < 
#                         |_| |_| |_|_|   |_|\___/_/\_\                           
#                                
###############################################################################

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


for(critter in test){
  
  critter <- NULL
  
}





test <- CreateMosquitoes_Defined_Genotype(numMos = c(2,2,2), genMos = c("AABBCC", "DDDDDD", "DSDF"),minAge = 10, maxAge = 11,ageDist = c(1,0))
listedTest <- replicate(n = 10, expr = test, simplify = FALSE)
listedTest2 <- replicate(n = 10, expr = listedTest, simplify = FALSE)


numNodes = 10
hold <- vector(mode = "list", length = numNodes)
for(i in 1:numNodes){
  for(j in 1:numNodes){
    hold[[i]] <- c(hold[[i]], listedTest2[[i]][[j]])
  }
}



eggsHist <- rlnorm(n = 10000, meanlog = log(x = 1), sdlog = log(x = 2))
larvaHist <- rlnorm(n = 1000000, meanlog = log(x = 1+14+1), sdlog = log(x = 1.3))
larvaHist <- rpois(n = 10000, lambda = 15)
larvaHist <- rgamma(n = 10000, shape = 15, scale = 1)


hold <- hist(x = larvaHist, breaks = 50, plot = FALSE)

hold$counts <- hold$counts/sum(hold$counts)

plot(hold, xlim = c(0, max(larvaHist)))

