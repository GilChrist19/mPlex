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

reference <- MakeReferenceMultiplex_oLocus(H=c(0.9, 0.9,0.9),R=c(0,0,0), S=c(0,0,0)/3, d=c(0,0,0))






###############################################################################
# Original one-locus multiplex Function
###############################################################################

#given the alleles, this will create a vector of all unique genotypes, and
#  fill a cube as traditionally used in MGDrivE

#assume the inputs are
#this tells me there are 3 sites, with 2 potential resistant alleles, and background mutation
##

Multiplex_FULL <- function(H=c(0.9, .6), R=0.0001, S=R/3, d=0.0001){
  #     H=homing alleles, given as a vector of homing rate. The length of that
  #         vector is the multiplex number
  #     R=resistant allele, damaging, from NHEJ
  #     S=resistant allele, not damaging, from NHEJ. defined as 1/3 the rate of R
  #     d=background mutations

  gtype <- c("W", "H", "R", "S")

  library("gtools")
  library("stringi")
  library("matrixStats")

  #experimental
  library("slam")
  library(tractor.base)


  ###############################################################################
  #generate all genotypes, set up vectors and matrices

  #gets all allele possibilities.
  # expand.grid is faster, but requires a list.
  openGtype <- gtools::permutations(n = length(gtype),r =length(H) ,v = gtype, set = FALSE, repeats.allowed = TRUE )

  #convert to list, openGtype is not always the same dimension
  # ncol() is always equal to length(H)
  openGtype <- lapply(X = seq_len(length.out = ncol(openGtype)), FUN = function(i){openGtype[,i]})

  #stitch the alleles together
  alleleType <- do.call(what = stringi::stri_join, args = openGtype)

  #get all possible combinations of allelles, always some rows by 2 columns
  openAllele <- gtools::combinations(n = length(alleleType), r = 2, v = alleleType, set = FALSE, repeats.allowed = TRUE)

  #final allele combinations
  genotypes <- do.call(what = stringi::stri_join, args = list(openAllele[,1], openAllele[,2]))

  #create a thing to see if this works
  #these work for H=2
  setNames(object = numeric(length(genotypes)), nm = genotypes)
  matrix(data = 0, nrow = length(genotypes), ncol = length(genotypes), dimnames = list(genotypes, genotypes))
  tMatrix <- array(data = 0, dim = c(length(genotypes),length(genotypes),length(genotypes)), dimnames = list(genotypes,genotypes,genotypes))


  ###############################################################################
  #setup all probability matrices

  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(H), dimnames = list(gtype, NULL))

  homingProbs[1, ] <- 1-d-H #chance to stay W is 1-homing-background mutation
  homingProbs[2, ] <- H*(1-R-S) #chance to become homing is H*1-H*R-H*S
  #homing probs, - NHEJ probs
  homingProbs[3, ] <- H*R   #NHEJ caused resistance, detrimentalt allele
  homingProbs[4, ] <- d+H*S #good resistant allele, from NHEJ and background mutation rate

  #set up lists to hold probabilities
  homProbsList <- vector(mode = "list", length = length(H))
  noHomProbsList <- vector(mode = "list", length = length(H))

  #fill the lists
  for(i in 1:length(H)){
    homProbsList[[i]]$W <- homingProbs[,i]
    homProbsList[[i]]$H <- setNames(object = c(1-d,d), nm = c("H", "S"))
    homProbsList[[i]]$R <- 1
    homProbsList[[i]]$S <- 1

    noHomProbsList[[i]]$W <- setNames(object = c(1-d,d), nm = c("W", "S"))
    noHomProbsList[[i]]$H <- setNames(object = c(1-d,d), nm = c("H", "S"))
    noHomProbsList[[i]]$R <- 1
    noHomProbsList[[i]]$S <- 1
  }

  ###############################################################################
  #fill transition matrix

  #check if homing allele is present
  hCheck <- grepl(pattern = "H", x = genotypes, ignore.case = FALSE, perl = FALSE, fixed = TRUE)
  nHom <- length(H)

  #loop over all matings, male and female
  for (fi in 1:length(genotypes)){
    for (mi in 1:length(genotypes)){

      #female and male genotype for this mating
      fSplit <- strsplit(genotypes[fi], "")[[1]]
      mSplit <- strsplit(genotypes[mi], "")[[1]]

      #split into parental alleles
      fGList <- list(allele1 = fSplit[1:nHom], allele2 = fSplit[(nHom+1):(2*nHom)] )
      mGList <- list(allele1 = mSplit[1:nHom], allele2 = mSplit[(nHom+1):(2*nHom)] )

      #set up allele and probability lists
      fAllele1 <- vector(mode = "list", length = nHom)
      fProbs1 <- vector(mode = "list", length = nHom)

      fAllele2 <- vector(mode = "list", length = nHom)
      fProbs2 <- vector(mode = "list", length = nHom)

      mAllele1 <- vector(mode = "list", length = nHom)
      mProbs1 <- vector(mode = "list", length = nHom)

      mAllele2 <- vector(mode = "list", length = nHom)
      mProbs2 <- vector(mode = "list", length = nHom)

      #fill female lists




      if(hCheck[fi]){
        #homing in females
        #lists of allele letter and probabilities
        for (i in 1:nHom){

          #female allele1
          if(fGList[[1]][i]=="W"){
            fAllele1[[i]] <- c("W", "H", "R", "S")
            fProbs1[[i]] <- homProbsList[[i]]$W
          } else if(fGList[[1]][i]=="H"){
            fAllele1[[i]] <- c("H", "S")
            fProbs1[[i]] <- homProbsList[[i]]$H
          } else if(fGList[[1]][i]=="R"){
            fAllele1[[i]] <- "R"
            fProbs1[[i]] <- homProbsList[[i]]$R
          } else if(fGList[[1]][i]=="S"){
            fAllele1[[i]] <- "S"
            fProbs1[[i]] <- homProbsList[[i]]$S
          }

          #female allele2
          if(fGList[[2]][i]=="W"){
            fAllele2[[i]] <- c("W", "H", "R", "S")
            fProbs2[[i]] <- homProbsList[[i]]$W
          } else if(fGList[[2]][i]=="H"){
            fAllele2[[i]] <- c("H", "S")
            fProbs2[[i]] <- homProbsList[[i]]$H
          } else if(fGList[[2]][i]=="R"){
            fAllele2[[i]] <- "R"
            fProbs2[[i]] <- homProbsList[[i]]$R
          } else if(fGList[[2]][i]=="S"){
            fAllele2[[i]] <- "S"
            fProbs2[[i]] <- homProbsList[[i]]$S
          }

        }

      } else {
        #no homing in females
        #lists of allele letter and probabilities
        for (i in 1:nHom){

          #female allele1
          if(fGList[[1]][i]=="W"){
            fAllele1[[i]] <- c("W", "S")
            fProbs1[[i]] <- noHomProbsList[[i]]$W
          } else if(fGList[[1]][i]=="H"){
            fAllele1[[i]] <- c("H", "S")
            fProbs1[[i]] <- noHomProbsList[[i]]$H
          } else if(fGList[[1]][i]=="R"){
            fAllele1[[i]] <- "R"
            fProbs1[[i]] <- noHomProbsList[[i]]$R
          } else if(fGList[[1]][i]=="S"){
            fAllele1[[i]] <- "S"
            fProbs1[[i]] <- noHomProbsList[[i]]$S
          }

          #female allele2
          if(fGList[[2]][i]=="W"){
            fAllele2[[i]] <- c("W", "S")
            fProbs2[[i]] <- noHomProbsList[[i]]$W
          } else if(fGList[[2]][i]=="H"){
            fAllele2[[i]] <- c("H", "S")
            fProbs2[[i]] <- noHomProbsList[[i]]$H
          } else if(fGList[[2]][i]=="R"){
            fAllele2[[i]] <- "R"
            fProbs2[[i]] <- noHomProbsList[[i]]$R
          } else if(fGList[[2]][i]=="S"){
            fAllele2[[i]] <- "S"
            fProbs2[[i]] <- noHomProbsList[[i]]$S
          }

        }

      }

      #fill male lists
      if(hCheck[mi]){
        #homing in males
        #lists om allele letter and probabilities
        for (i in 1:nHom){

          #male allele1
          if(mGList[[1]][i]=="W"){
            mAllele1[[i]] <- c("W", "H", "R", "S")
            mProbs1[[i]] <- homProbsList[[i]]$W
          } else if(mGList[[1]][i]=="H"){
            mAllele1[[i]] <- c("H", "S")
            mProbs1[[i]] <- homProbsList[[i]]$H
          } else if(mGList[[1]][i]=="R"){
            mAllele1[[i]] <- "R"
            mProbs1[[i]] <- homProbsList[[i]]$R
          } else if(mGList[[1]][i]=="S"){
            mAllele1[[i]] <- "S"
            mProbs1[[i]] <- homProbsList[[i]]$S
          }

          #male allele2
          if(mGList[[2]][i]=="W"){
            mAllele2[[i]] <- c("W", "H", "R", "S")
            mProbs2[[i]] <- homProbsList[[i]]$W
          } else if(mGList[[2]][i]=="H"){
            mAllele2[[i]] <- c("H", "S")
            mProbs2[[i]] <- homProbsList[[i]]$H
          } else if(mGList[[2]][i]=="R"){
            mAllele2[[i]] <- "R"
            mProbs2[[i]] <- homProbsList[[i]]$R
          } else if(mGList[[2]][i]=="S"){
            mAllele2[[i]] <- "S"
            mProbs2[[i]] <- homProbsList[[i]]$S
          }

        }

      } else {
        #no homing in males
        #lists om allele letter and probabilities
        for (i in 1:nHom){

          #male allele1
          if(mGList[[1]][i]=="W"){
            mAllele1[[i]] <- c("W", "S")
            mProbs1[[i]] <- noHomProbsList[[i]]$W
          } else if(mGList[[1]][i]=="H"){
            mAllele1[[i]] <- c("H", "S")
            mProbs1[[i]] <- noHomProbsList[[i]]$H
          } else if(mGList[[1]][i]=="R"){
            mAllele1[[i]] <- "R"
            mProbs1[[i]] <- noHomProbsList[[i]]$R
          } else if(mGList[[1]][i]=="S"){
            mAllele1[[i]] <- "S"
            mProbs1[[i]] <- noHomProbsList[[i]]$S
          }

          #male allele2
          if(mGList[[2]][i]=="W"){
            mAllele2[[i]] <- c("W", "S")
            mProbs2[[i]] <- noHomProbsList[[i]]$W
          } else if(mGList[[2]][i]=="H"){
            mAllele2[[i]] <- c("H", "S")
            mProbs2[[i]] <- noHomProbsList[[i]]$H
          } else if(mGList[[2]][i]=="R"){
            mAllele2[[i]] <- "R"
            mProbs2[[i]] <- noHomProbsList[[i]]$R
          } else if(mGList[[2]][i]=="S"){
            mAllele2[[i]] <- "S"
            mProbs2[[i]] <- noHomProbsList[[i]]$S
          }

        }

      }



      #get all combinations of those lists, separately letters and numbers, but in same order
      fKid1 <- expand.grid(fAllele1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      fPKid1 <- expand.grid(fProbs1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      fKid2 <- expand.grid(fAllele2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      fPKid2 <- expand.grid(fProbs2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      mKid1 <- expand.grid(mAllele1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      mPKid1 <- expand.grid(mProbs1, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      mKid2 <- expand.grid(mAllele2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      mPKid2 <- expand.grid(mProbs2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      #add letters together, multiply probabilities
      fKid1Collapse <- do.call(what = stringi::stri_join, args = fKid1)
      fPKid1Collapse <- matrixStats::rowProds(x = as.matrix(fPKid1), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")
      fKid2Collapse <- do.call(what = stringi::stri_join, args = fKid2)
      fPKid2Collapse <- matrixStats::rowProds(x = as.matrix(fPKid2), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")

      mKid1Collapse <- do.call(what = stringi::stri_join, args = mKid1)
      mPKid1Collapse <- matrixStats::rowProds(x = as.matrix(mPKid1), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")
      mKid2Collapse <- do.call(what = stringi::stri_join, args = mKid2)
      mPKid2Collapse <- matrixStats::rowProds(x = as.matrix(mPKid2), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")

      #since these are 2 alleles in same parent, ofspring can't get both, so combine to 1 long list of all possible alleles
      # reorder allele labels, expand.grid calls them backwords from my definition
      fPoss <- stringi::stri_reverse(str = c(fKid1Collapse, fKid2Collapse))
      fProbsPoss <- c(fPKid1Collapse, fPKid2Collapse)

      mPoss <- stringi::stri_reverse(str = c(mKid1Collapse, mKid2Collapse))
      mProbsPoss <- c(mPKid1Collapse, mPKid2Collapse)

      #all combinations of male and female alleles
      offspringgenotypes <- expand.grid(fPoss, mPoss, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      offspringprobabilities <- expand.grid(fProbsPoss, mProbsPoss, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      #recombine into genotypes and their probabilities, 2 lists in same order.
      # collapse is 2 to 1 because expand.grid orders vectors backwords from combinations and permutations above
      offspringgenocollapse <- do.call(what = stringi::stri_join, args = list(offspringgenotypes[,2], offspringgenotypes[,1]))

      offspringprobcollapse <- matrixStats::rowProds(x = as.matrix(offspringprobabilities), rows = NULL, cols = NULL, na.rm = TRUE, method = "expSumLog")/4



      #normalize probabilities
      #offspringprobcollapse <- offspringprobcollapse/4

      #reverse allele permutations that don't match allele combinations
      index <- which((offspringgenocollapse %in% genotypes)==FALSE)

      offspringgenocollapse[index] <- do.call(what = stringi::stri_join, args = list(offspringgenotypes[index,1], offspringgenotypes[index,2]))

      #set names on vector for better access
      nIndex <- match(x = offspringgenocollapse, table = genotypes)


     # offspringprobcollapse <- setNames(object = offspringprobcollapse, nm = offspringgenocollapse)



      #tMatrix[fi,mi, nIndex ] <- tMatrix[fi,mi, nIndex ] + offspringprobcollapse[nIndex]



      #store in matrix
      for (i in 1:length(offspringgenocollapse)){

        tMatrix[fi,mi, nIndex[i] ] <- tMatrix[fi,mi, nIndex[i] ] + offspringprobcollapse[i]

      }

    }
  }

  return(tMatrix)

}

###############################################################################
# better mosquito
###############################################################################

#original
Mosquito <- R6::R6Class(classname = "mosquito",
                        portable = TRUE,
                        cloneable = FALSE,
                        lock_class = FALSE,
                        lock_objects = FALSE,

                        # public memebers
                        public = list(

                          # constructor
                          initialize = function(genotype=NULL, age=NULL){
                            private$age = age
                            private$mate = NULL
                            private$genotype = genotype

                          }, # end constructor

                          # setters
                          set_age = function(age=NULL){private$age = age},
                          set_mate = function(mate=NULL){private$mate = mate},
                          set_genotype = function(genotype=NULL){private$genotype = genotype},

                          age_one_day = function() {private$age = private$age + 1},
                          print_female = function(){file.path(private$age, private$genotype, private$mate, fsep = ",")},
                          print_male = function(){file.path(private$age, private$genotype, fsep = ",")},
                          #file.path can't handle nulls

                          #getters
                          get_age = function(){private$age},
                          get_mate = function(){private$mate},
                          get_genotype = function(){private$genotype}

                        ), # end public

                        private = list(

                          # fields
                          age = NULL,
                          mate = NULL,
                          genotype = NULL

                        ) # end private

) # end class definition

#new
# exact replica of Mosquito class
# drop in, ready to go
NewMosquito <- function(genotype=NULL, age=NULL){
  #last attribute
  mate = NULL

  #getters/setters
  set_age <- function(newAge=NULL){age <<- newAge}
  set_mate <- function(matGen=NULL){mate <<- matGen}
  set_genotype <- function(newGen=NULL){genotype <<- newGen}

  get_age <- function(){age}
  get_mate <- function(){mate}
  get_genotype <- function(){genotype}

  #functions
  age_one_day <- function(){age <<- age + 1}
  print_female = function(){file.path(age, genotype, mate, fsep = ",")}
  print_male = function(){file.path(age, genotype, fsep = ",")}

  #this is an R environment
  environment()
}

# Smallest, fastest
NewMosquitoSMALLEST <- function(genotype=NULL, age=NULL){
  #last attribute
  mate = NULL

  # no getters/setters required

  #functions
  age_one_day = function(){age <<- age + 1}
  print_female = function(){file.path(age, genotype, mate, fsep = ",")}
  print_male = function(){file.path(age, genotype, fsep = ",")}

  #this is an R environment
  environment()
}


R6Private_Portable_NOCLASSCLONE <- R6::R6Class("R6Private_Portable_NOCLASSCLONE",
                         portable = TRUE,
                         class = FALSE,
                         cloneable = FALSE,

                         public = list(
                           initialize = function(x = 1) private$x <- x,
                           getx = function() x,
                           inc = function(n = 1) x <<- x + n
                         ),

                         private = list(
                           x = NULL
                           )

)

R6Private_NOCLASSCLONE <- R6::R6Class("R6Private_NOCLASSCLONE",
                                               portable = FALSE,
                                               class = FALSE,
                                               cloneable = FALSE,

                                               public = list(
                                                 initialize = function(x = 1) private$x <- x,
                                                 getx = function() x,
                                                 inc = function(n = 1) x <<- x + n
                                               ),

                                               private = list(
                                                 x = NULL
                                               )

)



test <- lapply(X = 1:100, FUN = function(X){Mosquito$new(genotype = "AA", age = X)})
test2 <- lapply(X = 1:100, FUN = function(x){NewMosquitoSMALLEST(genotype = "AA", age = x)})

hold <- lapply(X = test2, FUN = '[[', "age")

ages <- rlnorm(n = length(test2), meanlog = log(15), sdlog = log(1.2))

oldMos <- Mosquito$new(genotype = "AA", age = 10)
newMos <- NewMosquito(genotype = "AA", age = 10)
smallMos <- NewMosquitoSMALLEST(genotype = "AA", age = 10)

matured <- integer(length = 100)

microbenchmark::microbenchmark(for(i in 1:100){
                                matured[i] <- test2[i]["age"]
                               },
                               vapply(X = test2, FUN = '[[', "age", FUN.VALUE = integer(length = 1L)),
                               times = 1000)

###############################################################################
# Rcpp Mating Function
###############################################################################

library(Rcpp)
source("./R/4_DaisyGeneratingFunction.R")


fGen <- "HWHHHWWWWW"
mGen <- "WWWWWWWWWW"
reproductionReference <- MakeReference_DaisyDrive(H = c(0.98, 0.98, 0.7, 0.98, 0.7),
                                                       R = c(0.0001,0.0001,0.0001,0.0001,0.0001),
                                                       S = c(0.00003,0.00003,0.00003,0.00003,0.00003),
                                                       d = c(0,0,0,0,0))


fGen <- "WH"
mGen <- "WW"
reproductionReference <- MakeReference_DaisyDrive(H = c(0.98),
                                                        R = c(0.0001),
                                                        S = c(0.00003),
                                                        d = c(0))




sourceCpp("./src/4_DaisyGeneratingFunction.cpp")
HoldTEST <- DaisyOffspring_C(fGen = fGen, mGen = mGen, reference = reproductionReference)
HoldREF <- DaisyOffspring(fGen = fGen, mGen = mGen, reference = reproductionReference)




testList = list(c("A","B","C"), c("D","E"), c("F","G","H","I"))
ListTest(myList = testList)

expand.grid(testList)


microbenchmark::microbenchmark(DaisyOffspring_C(fGen = fGen, mGen = mGen, reference = reproductionReference),
                               DaisyOffspring(fGen = fGen, mGen = mGen, reference = reproductionReference),
                               times = 100)




