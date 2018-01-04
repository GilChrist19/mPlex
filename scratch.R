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





