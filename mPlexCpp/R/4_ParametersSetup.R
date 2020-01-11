###############################################################################
#                            ____  __          ______          
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
#                                               /_/   /_/      
###############################################################################
###############################################################################
# Release Functions
###############################################################################

#' Make List of Mosquito Releases
#'
#' Sets up a release schedule for a single patch, returns a list to be used in
#' \code{\link{oneDay_Releases_Patch}}
#'
#' @param releaseStart Day releases start
#' @param releaseEnd Day releases end
#' @param releaseInterval Interval between releases
#' @param genMos List of genotypes for new Mosquitoes
#' @param numMos Integer number of Mosquitoes to create
#' @param minAge Integer specifying the minimum age of Mosquitoes
#' @param maxAge Integer specifying the maximum age of Mosquiotes
#' @param ageDist Distribution for ages of Mosquitoes. Must be length(maxAge-minAge+1)
#'
#' @details This function creates a release schedule and the genotypes and ages
#' of mosquitoes to be released at that time. The release schedule is determined
#' by releaseStart, releaseEnd, and releaseInterval. The mosquitoes to be released
#' are determined by genMos (a list of genotypes), numMos (a vector of how many of
#' each genotype to release), and minAge/maxAge/ageDist. The last three determine
#' the age distribution of the mosquitoes to be released.
#'
#' @return List of release dates and the population to be released on that day.
#'
#' @examples
#' # to setup for 3 patches but only release in the first with a defined release schedule:
#'
#' patchReleases = replicate(n = 3,
#'                           expr = list(maleReleases = NULL,femaleReleases = NULL,larvaeReleases=NULL),
#'                           simplify = FALSE)
#'
#' set.seed(42)
#' patchReleases[[1]]$femaleReleases = Release_basicRepeatedReleases(releaseStart = 5,
#' releaseEnd = 30,
#' releaseInterval = 5,
#' genMos = c("A", "B", "C"),
#' numMos = c(10, 20, 30),
#' minAge = 1,
#' maxAge = 10,
#' ageDist = rep(x = 1, times = 10-1+1)/(10-1+1)) #uniform distribution of ages
#'
#' set.seed(42)
#' patchReleases[[1]]$maleReleases = Release_basicRepeatedReleases(releaseStart = 50,
#' releaseEnd = 60,
#' releaseInterval = 1,
#' genMos = c("A", "B", "C"),
#' numMos = c(10, 20, 30),
#' minAge = 5,
#' maxAge = 10,
#' ageDist = c(0,0,0,0,0,1)) #all mosqutioes will be 10
#'
#' set.seed(42)
#' patchReleases[[1]]$larvaeReleases = Release_basicRepeatedReleases(releaseStart = 1,
#' releaseEnd = 5,
#' releaseInterval = 1,
#' genMos = c("A", "B", "C"),
#' numMos = c(10, 20, 30),
#' minAge = 5,
#' maxAge = 10,
#' ageDist = c(0,1,0,1,0,0)) #half of mosquitoes are 6, half are 8
#'
#' @export
Release_basicRepeatedReleases <- function(releaseStart, releaseEnd, releaseInterval, genMos, numMos, minAge, maxAge, ageDist){
  #genMos is a list of genotypes to relaese
  #numMos is a vector of the number of mosquitoes you want to make, corresponding
  #  to the genotypes of genMos
  
  #minAge, maxAge are the min/max age range. To get a single age, must be length
  # 2 with a c(1,0) vector in ageDist
  #ageDist - probabilities to sample from for age range. must be length
  #  minAge:maxAge
  
  # check timing of releases
  if(releaseInterval > (releaseEnd - releaseStart)){
    stop("interval between releases cannot be greater than time between start and end of releases")
  }
  
  # Initialize release times. Initialize return list
  releaseTimes = seq(from=releaseStart,to = releaseEnd,by = floor(releaseInterval))
  releaseList = vector(mode="list",length=length(releaseTimes))
  
  # Initialize genotypes and ages
  if(length(genMos) != length(numMos)) stop("Please specify how many of each genotype of mosquitoes")
  genotypeVector = rep.int(x = genMos, times = numMos)
  ageVector <- sample(x = minAge:maxAge, size = sum(numMos), replace = TRUE, prob = ageDist)
  
  # check for male/female/larvae. Fill appropriate list.
  for(tx in 1:length(releaseTimes)){
    releaseList[[tx]]$genVec = genotypeVector
    releaseList[[tx]]$ageVec = ageVector
    releaseList[[tx]]$tRelease = releaseTimes[[tx]]
  }
  
  return(releaseList)
}

###############################################################################
# Batch Migration Setup
###############################################################################

#' Make List of Batch Migration Parameters
#'
#' Sets up a list containing the probability of a batch migration, the fractional amount of males/females
#' that migrate, and the weighted probabilities for where to migrate.
#'
#' @param batchProbs Probability of a batch migration, either 1 number or vector of length equal to the number of patches
#' @param sexProbs Population fraction of males and females that migration. Either vector c(M,F) or matrix of 2 columns
#' @param moveProbs Movement matrix. There must be 0 probability of moving to your own patch. Default is equal movement to all other patches.
#' @param numPatches Number of patches in the simulation
#'
#' @examples
#' # to setup for 3 patches
#' batchMigration = basicBatchMigration(batchProbs = 1e-5, sexProbs = c(0.1, 0.01), numPatches = 3)
#'
#' @export
basicBatchMigration <- function(batchProbs = 1e-5, sexProbs = c(0.01, 0.01),
                                moveProbs, numPatches = 1){
  
  ##########
  # batchProbs
  ##########
  # check values of probs
  if(!all(batchProbs<=1, batchProbs>=0)){
    stop("Probability of batch migration must be between 0 and 1.")
  }
  
  # check length of probs
  if(length(batchProbs) == 1){
    batchProbs = rep(x = batchProbs, numPatches)
  } else if(length(batchProbs) != numPatches){
    stop("Length of batchProbs vector must be 1 or equal to the number of patches.")
  }
  
  
  ##########
  # sexProbs
  ##########
  # check values of sexProbs
  if(!all(sexProbs<=1, sexProbs>=0)){
    stop("Sex specific movement fraction must be between 0 and 1")
  }
  
  # check dimensions of sexProbs
  sexProbsError <- c("SexProbs is misspecified.\n",
    "If 1 number is provided, male/female probabilities are the same for all patches.\n",
    "If a length 2 vector is provided, male/female probabilities differ, but all patches are the same.\n",
    "If a matrix is provided, it must have rows equal to the number of patches, and 2 columns.")
  
  if(is.null(dim(sexProbs))){
    # ensure dimensions are 1 or 2
    if(!any(length(sexProbs) == c(1,2))) stop(sexProbsError)
    
    sexProbs = matrix(data = sexProbs, nrow = numPatches, ncol = 2,
                         byrow = TRUE, dimnames = list(NULL, c("M","F")))
    
  } else if(!all(dim(sexProbs) == c(numPatches, 2))) stop(sexProbsError)

  
  ##########
  # moveProbs
  ##########
  moveProbsError <- c("moveProbs is misspecified.\n",
    "The default movement is equal probabilities to all other patches.\n",
    "If it is user specified, moveProbs must have dimensions equal to the number of patches, ",
    "all diagonal values must be zero, and the sum of each row must be 1.")
  
  if(missing(moveProbs)){
    # setup movement matrix
    moveProbs <- matrix(data = 1L, nrow = numPatches, ncol = numPatches)
    diag(moveProbs) <- 0L
    moveProbs <- moveProbs/rowSums(x = moveProbs)
    
  } else {
    # dimensions of moveProbs
    if(!all(dim(as.matrix(moveProbs)) == numPatches) ) stop(moveProbsError)
    
    # diagValues
    if(!all(diag(moveProbs) == 0) ) stop(moveProbsError)
    
    # rowSums
    if(!all(rowSums(as.matrix(moveProbs))) == 1 ) stop(moveProbsError)
  }

  
  # return basic batch migration
  return(list("batchProbs" = batchProbs,
              "sexProbs" = sexProbs,
              "moveProbs" = moveProbs)
  )
}

###############################################################################
# Network Initialization Functions
###############################################################################

#' Network Parameters
#'
#' Generate parameters for simulation on a \code{\link{Network}}.
#'
#' @details Parameters average generation time \eqn{g}, population growth rate \eqn{R_{m}},
#' aquatic mortality \eqn{\mu_{Aq}}, and aquatic survival \eqn{\theta_{Aq}}
#' are shared between patches and calculated by \code{\link{calcAverageGenerationTime}},
#' \code{\link{calcPopulationGrowthRate}}, \code{\link{calcLarvalStageMortalityRate}},
#' and \code{\link{calcAquaticStagesSurvivalProbability}}. \cr
#' Patch-specific parameters \eqn{\alpha} and \eqn{L_{eq}}
#' are calculated for each patch by \code{\link{calcDensityDependentDeathRate}}
#' and \code{\link{calcLarvalPopEquilibrium}}.
#'
#' @param nPatch number of \code{\link{Patch}}
#' @param simTime maximum time to run simulation
#' @param sampTime how often output is written. Default is 1 (i.e. every day)
#' @param moveVar variance of stochastic movement (not used in diffusion model of migration).
#' It affects the concentration of probability in the Dirchlet simplex, small
#' values lead to high variance and large values lead to low variance.
#' @param tEgg length of egg stage
#' @param tLarva length of larval instar stage
#' @param tPupa length of pupal stage
#' @param beta female egg batch size of wild-type
#' @param muAd wild-type daily adult mortality (1/muAd is average wild-type lifespan)
#' @param dayGrowthRate daily population growth rate (used to calculate equilibrium)
#' @param AdPopEQ a single number or vector of adult population size at equilibrium
#' @param runID begin counting runs with this set of parameters from this value
#'
#' @return List(nPatch=int, simTime=int, sampTime=int, moveVar=numeric,
#' runID=int, stageTime=vec int, beta=int, dayGrowthRate=numeric, AdPopEq=int vec,
#' genTime=numeric, genGrowthRate=numeric, mu=vec numeric,
#' thetaAq=vec numeric, alpha=vec numeric, Leq=vec int)
#'
#' @examples
#' This is an example
#' So is this
#'
#' @export
NetworkParameters <- function(
  nPatch,
  simTime,
  sampTime = 1,
  moveVar = 1000,
  tEgg = 1L,
  tLarva = 14L,
  tPupa = 1L,
  beta = 32,
  muAd = 0.123,
  dayGrowthRate = 1.096,
  AdPopEQ,
  runID = 1L
){

  # check required parameters
  if(any(missing(nPatch),missing(simTime),missing(AdPopEQ))){
    stop("nPatch, simTime, and AdPopEQ must be provided by the user.")
  }
  
  # make empty parameter list
  pars = list()

  # fill list
  pars$nPatch = nPatch
  pars$simTime = simTime
  pars$sampTime = sampTime
  pars$moveVar = moveVar
  pars$runID = runID

  # biological parameters
  tAdult = as.integer(x = round(x = 1/muAd))
  pars$stageTime = setNames(object = c(tEgg, tLarva, tPupa, tAdult),
                       nm = c("E","L","P","A"))
  pars$beta = beta

  # initial parameters
  pars$dayGrowthRate = dayGrowthRate

  if(length(AdPopEQ) == 1){
    AdPopEQ = rep.int(x = AdPopEQ, times = nPatch)
  } else if(length(AdPopEQ)!=nPatch){
    stop("length of AdPopEQ vector must be 1 or nPatch (number of patches)")
  }
  pars$AdPopEQ = AdPopEQ


  # derived parameters
  pars$genTime = calcAverageGenerationTime(pars$stageTime[c("E","L","P")],muAd)
  pars$genGrowthRate = calcPopulationGrowthRate(dayGrowthRate,pars$genTime)

  #aquatic daily death rate. Equation is for larval stage, but the assumption
  # is that all are equal
  muAq = calcLarvalStageMortalityRate(pars$genGrowthRate,muAd,beta,pars$stageTime[c("E","L","P")])
  pars$mu = setNames(object = c(muAq, muAq, muAq, muAd),
                     nm = c("E", "L", "P", "A"))

  # Survival probability for each state. Holdover from MGDrivE, useful in the equations
  pars$thetaAq = c("E"=calcAquaticStageSurvivalProbability(muAq,tEgg),
                   "L"=calcAquaticStageSurvivalProbability(muAq,tLarva),
                   "P"=calcAquaticStageSurvivalProbability(muAq,tPupa))

  # patch-specific derived parameters
  pars$alpha = calcDensityDependentDeathRate(beta, pars$thetaAq, pars$stageTime["L"],
                                             AdPopEQ, pars$genGrowthRate)
  pars$Leq = calcLarvalPopEquilibrium(pars$alpha,pars$genGrowthRate)
  
  
  # check the list
  invisible(Map(f = check, pars))


  return(pars)
}

########################################################################
# Equilibrium Parameters
########################################################################

# check for positive parameter values
check <- function(x){
  if(is.numeric(x)||is.integer(x)){
    if(any(x < 0)){
      stop("only nonnegative parameter values allowed")
    }
  }
}

#' Calculate Average Generation Time
#'
#' Calculate \eqn{g}, average generation time, given by: \deqn{g=T_e+T_l+T_p+\frac{1}{\mu_{ad}}}
#'
#' @param stagesDuration vector of lengths of aquatic stages, \eqn{T_{e}, T_{l}, T_{p}}
#' @param adultMortality adult mortality rate, \eqn{\mu_{ad}}
#'
#' @export
calcAverageGenerationTime <- function(stagesDuration, adultMortality){
  return(sum(stagesDuration) + (1.0 / adultMortality))
}

#' Calculate Population Growth Rate
#'
#' Calculate \eqn{R_{m}}, population growth in absence of density-dependent mortality, given by: \deqn{(r_{m})^{g}}
#'
#' @param dailyPopGrowthRate daily population growth rate, \eqn{r_{m}}
#' @param averageGenerationTime see \code{\link{calcAverageGenerationTime}}
#'
#' @export
calcPopulationGrowthRate <- function(dailyPopGrowthRate, averageGenerationTime){
  return(dailyPopGrowthRate^averageGenerationTime)
}

#' Calculate Larval Stage Mortality Rate
#'
#' Calculate \eqn{\mu_{l}}, the larval mortality, given by \deqn{\mu_l=1-\Bigg( \frac{R_m * \mu_{ad}}{1/2 * \beta_k * (1-\mu_m)} \Bigg)^{\frac{1}{T_e+T_l+T_p}}}
#'
#' @param generationPopGrowthRate see \code{\link{calcPopulationGrowthRate}}
#' @param adultMortality adult mortality rate, \eqn{\mu_{ad}}
#' @param fertility number of eggs per oviposition for wild-type females, \eqn{\beta_{k}}
#' @param aquaticStagesDuration vector of lengths of aquatic stages, \eqn{T_{e}, T_{l}, T_{p}}
#'
#' @export
calcLarvalStageMortalityRate <- function(generationPopGrowthRate, adultMortality, fertility, aquaticStagesDuration){
  a = 2*generationPopGrowthRate*adultMortality
  b = fertility*(1-adultMortality)
  c = sum(aquaticStagesDuration)
  return(1-(a/b)^(1/c))
}

#' Calculate Aquatic Stage Surival Probability
#'
#' Calculate \eqn{\theta_{st}}, density-independent survival probability, given by: \deqn{\theta_{st}=(1-\mu_{st})^{T_{st}}}
#'
#' @param mortalityRate daily mortality probability, \eqn{\mu_{st}}
#' @param stageDuration duration of aquatic stage, \eqn{T^{st}}
#'
#' @export
calcAquaticStageSurvivalProbability <- function(mortalityRate, stageDuration){
  return((1-mortalityRate)^stageDuration)
}

#' Calculate Density-dependent Larval Mortality
#'
#' Calculate \eqn{\alpha}, the strength of density-dependent mortality during the larval stage, given by: \deqn{\alpha=\Bigg( \frac{1/2 * \beta_k * \theta_e * Ad_{eq}}{R_m-1} \Bigg) * \Bigg( \frac{1-(\theta_l / R_m)}{1-(\theta_l / R_m)^{1/T_l}} \Bigg)}
#'
#' @param fertility number of eggs per oviposition for wild-type females, \eqn{\beta_{k}}
#' @param thetaAq vector of density-independent survival probabilities of aquatic stages, \eqn{\theta_{e}, \theta_{l}}
#' @param timeAq Just the larval stage time, \eqn{T_l}
#' @param adultPopSizeEquilibrium adult population size at equilbrium, \eqn{Ad_{eq}}
#' @param populationGrowthRate population growth in absence of density-dependent mortality \eqn{R_{m}}
#'
#' @export
calcDensityDependentDeathRate <- function(fertility, thetaAq, timeAq, adultPopSizeEquilibrium, populationGrowthRate){
    prodA = (fertility*thetaAq[["E"]]*adultPopSizeEquilibrium) / (2*(populationGrowthRate-1))
    prodB_numerator = 1-(thetaAq[["L"]] / populationGrowthRate)
    prodB_denominator = 1-((thetaAq[["L"]]/populationGrowthRate)^(1/timeAq))
    return(prodA*(prodB_numerator/prodB_denominator))
}

#' Calculate Equilibrium Larval Population
#'
#' Equilibrium larval population to sustain population.
#'
#' @param alpha see \code{\link{calcDensityDependentDeathRate}}
#' @param Rm see \code{\link{calcPopulationGrowthRate}}
#'
#' @export
calcLarvalPopEquilibrium <- function(alpha,Rm){
  return(as.integer(round(alpha * (Rm-1))))
}
