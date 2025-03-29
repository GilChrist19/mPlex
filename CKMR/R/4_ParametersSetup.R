###############################################################################
#     ________ __ __  _______ 
#    / ____/ //_//  |/  / __ \
#   / /   / ,<  / /|_/ / /_/ /
#  / /___/ /| |/ /  / / _, _/ 
#  \____/_/ |_/_/  /_/_/ |_|  
#    
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
basicRepeatedReleases <- function(releaseStart, releaseEnd, releaseInterval, genMos, numMos, minAge, maxAge, ageDist){
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
#' @param numPatches Number of patches in the simulation
#'
#' @examples
#' # to setup for 3 patches
#' batchMigration = basicBatchMigration(batchProbs = 1e-5, sexProbs = c(0.1, 0.01), numPatches = 3)
#'
#' @export
basicBatchMigration <- function(batchProbs = 1e-5, sexProbs = c(0.01, 0.01),
                                numPatches = 1){
  
  # check length of probs
  if(!all(batchProbs<1)){
    stop("Probability of batch migration must be less than 1")
  }
  if(length(batchProbs) != numPatches){
    batchProbs = rep(x = batchProbs, numPatches)
  }
  
  # check length of sexes, make sure less than 1
  if(!all(sexProbs<1)){
    stop("Sex specific movement fraction must be less than 1")
  }
  if(is.null(dim(sexProbs))){
    sexProbsMat = matrix(data = sexProbs, nrow = numPatches, ncol = 2,
                         byrow = TRUE, dimnames = list(NULL, c("M","F")))
  }

  # setup movement matrix
  moveMat <- matrix(data = 1L, nrow = numPatches, ncol = numPatches)
  diag(moveMat) <- 0L
  moveMat <- moveMat/rowSums(x = moveMat)
  
  # return basic batch migration
  return(list("batchProbs" = batchProbs,
              "sexProbs" = sexProbsMat,
              "moveProbs" = moveMat)
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
#' @param moveVar variance of stochastic movement (not used in diffusion model of migration).
#' It affects the concentration of probability in the Dirchlet simplex, small
#' values lead to high variance and large values lead to low variance.
#' @param tEgg length of egg stage
#' @param tLarva length of larval instar stage
#' @param tPupa length of pupal stage
#' @param beta female egg batch size of wild-type
#' @param beta_const use a constant beta or draw from poisson?
#' @param muAd wild-type daily adult mortality (1/muAd is average wild-type lifespan)
#' @param maleMaxAge maximum possible age of adult male
#' @param femaleMaxAge maximum possible age of adult female
#' @param dayGrowthRate daily population growth rate (used to calculate equilibrium)
#' @param AdPopEQ vector of adult population size at equilibrium
#' @param runID begin counting runs with this set of parameters from this value
#'
#' @return List(nPatch=int, simTime=vec int, beta_const=logical, moveVar=numeric,
#' runID=int, stageTime=vec int, beta=int, dayGrowthRate=numeric, AdPopEq=int vec,
#' alleloTypes=list, genTime=numeric, genGrowthRate=numeric, mu=vec numeric,
#' maleMaxAge=int, femaleMaxAge=int
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
  moveVar = 1000,
  tEgg = 1L,
  tLarva = 14L,
  tPupa = 1L,
  beta = 32,
  beta_const = FALSE,
  muAd = 0.123,
  maleMaxAge = 9999,
  femaleMaxAge = 9999,
  dayGrowthRate = 1.096,
  AdPopEQ,
  runID = 1L
){
  
  # make empty parameter list
  pars = list()

  # fill list
  # nPatch safety check
  if(nPatch > 999999){
    stop("Sim is setup to perform fewer than 1,000,000 patches. \n\tPlease stop being silly.")
  }
  pars$nPatch = nPatch
  pars$simTime = simTime
  pars$moveVar = moveVar
  pars$runID = runID

  # biological parameters
  tAdult = as.integer(x = round(x = 1/muAd))
  pars$stageTime = setNames(object = c(tEgg, tLarva, tPupa, tAdult),
                       nm = c("E","L","P","A"))
  pars$beta = beta
  pars$beta_const = beta_const

  # initial parameters
  pars$dayGrowthRate = dayGrowthRate

  # derived parameters
  pars$genTime = calcAverageGenerationTime(pars$stageTime[c("E","L","P")],muAd)
  pars$genGrowthRate = calcPopulationGrowthRate(dayGrowthRate,pars$genTime)

  #aquatic daily death rate. Equation is for larval stage, but the assumption
  # is that all are equal
  muAq = calcLarvalStageMortalityRate(pars$genGrowthRate,muAd,beta,pars$stageTime[c("E","L","P")])
  pars$mu = setNames(object = c(muAq, muAq, muAq, muAd),
                     nm = c("E", "L", "P", "A"))
  pars$maleMaxAge = maleMaxAge
  pars$femaleMaxAge = femaleMaxAge

  # Survival probability for each state. Holdover from MGDrivE, useful in the equations
  pars$thetaAq = setNames(object = c(calcAquaticStageSurvivalProbability(muAq,tEgg),
                                calcAquaticStageSurvivalProbability(muAq,tLarva),
                                calcAquaticStageSurvivalProbability(muAq,tPupa)),
                     nm = c("E","L","P"))

  # patch-specific derived parameters
  #  this is updated to check the adult population shape, and calculate timve-varying
  #  larval parameters
  
  # larval population vector (Leq) needs to remain a vector
  #  it is used to setup the initial population, as well as reserve space for 
  #  better performance (read, less memory re-allocation) later.
  # However, b/c it's used to setup the initial population, this has to match 
  #  T=0.
  # Maybe, if performance becomes a real issue, we can output extra 
  #  vectors for max adult population size and max larval population size, for 
  #  reserving memory.
  
  if(NROW(AdPopEQ)==1 && NCOL(AdPopEQ)==1){
    # single patch, constant population
    
    # safety check
    if(length(AdPopEQ) != nPatch){
      stop("A single population size was specified by AdPopEQ, but the number of patches (nPatch) is not one.")
    }
    # set adult initial size
    pars$AdPopEQ <- as.integer(AdPopEQ)
    
    # calc patch-specific derived parameters
    pars$alpha <- matrix(data = calcDensityDependentDeathRate(fertility=beta,
                                                               thetaAq=pars$thetaAq,
                                                               timeAq=pars$stageTime["L"],
                                                               adultPopSizeEquilibrium=AdPopEQ,
                                                               populationGrowthRate=pars$genGrowthRate),
                         nrow = simTime, ncol = nPatch)

  } else if(NROW(AdPopEQ)!=1 && NCOL(AdPopEQ)!=1){
    # multiple patch, time-varying population
    
    # safety check
    if(NROW(AdPopEQ) != simTime){
      stop("The number of rows in AdPopEQ is not equal to the simulation time (simTime).")
    }
    if(NCOL(AdPopEQ) != nPatch){
      stop("The number of columns in AdPopEQ is not equal to the number of patches (nPatch).")
    }
    # set adult initial size
    pars$AdPopEQ <- as.integer(AdPopEQ[1, ,drop=TRUE])
    
    # calc patch-specific derived parameters
    pars$alpha <- calcDensityDependentDeathRate(fertility=beta,
                                                thetaAq=pars$thetaAq,
                                                timeAq=pars$stageTime["L"],
                                                adultPopSizeEquilibrium=AdPopEQ,
                                                populationGrowthRate=pars$genGrowthRate)

  } else if(NROW(AdPopEQ)!=1 || NCOL(AdPopEQ)!=1){
    # Either multiple patch, constant population
    #  or single patch, time-varying population
    
    # safety check
    if(nPatch == simTime){
      stop("The number of patches (nPatch) is equal to the simulation length (simTime). Therefore, it is impossible to determine the appropriate shape of AdPopEQ.
           Please change the simulation time or provide AdPopEQ as a matrix.")
    }
    
    if(length(AdPopEQ) == nPatch){
      # multiple patch, constant population
      
      # no safety checks for simTime
      
      # set adult initial size
      #  just in case someone puts in a 1-row/col matrix, drop dimensions
      pars$AdPopEQ <- as.integer(AdPopEQ)
      
      # calc patch-specific derived parameters
      pars$alpha <- matrix(data = calcDensityDependentDeathRate(fertility=beta,
                                                                 thetaAq=pars$thetaAq,
                                                                 timeAq=pars$stageTime["L"],
                                                                 adultPopSizeEquilibrium=AdPopEQ,
                                                                 populationGrowthRate=pars$genGrowthRate),
                           byrow=TRUE, nrow=simTime, ncol=nPatch)

    } else if(length(AdPopEQ) == simTime){
      # single patch, time-varying population
      
      # safety check
      if(nPatch != 1){
        stop("A single population was supplied for AdPopEQ, but the number of patches (nPatch) is not 1.")
      }
      
      # set adult initial size
      pars$AdPopEQ <- as.integer(AdPopEQ[1])
      
      # calc patch-specific derived parameters
      pars$alpha <- matrix(data = calcDensityDependentDeathRate(fertility=beta,
                                                                 thetaAq=pars$thetaAq,
                                                                 timeAq=pars$stageTime["L"],
                                                                 adultPopSizeEquilibrium=AdPopEQ,
                                                                 populationGrowthRate=pars$genGrowthRate),
                           nrow=simTime, ncol=nPatch)
    
    } else {
      stop("AdPopEQ was supplied as a vector, but it is neither the length of the simulation time (simTime) nor the number of patches (nPatch).")
    } # end vector AdPopEQ check
  } # end full AdPopEQ check
  # set larval equilibrium size - see note above
  pars$Leq <- calcLarvalPopEquilibrium(alpha=pars$alpha[1, ,drop=TRUE],
                                       Rm=pars$genGrowthRate)
  
  
  return(pars)
}

########################################################################
# Equilibrium Parameters
########################################################################

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
#' @param timeAq vector of lengths of aquatic stages, \eqn{T_{e}, T_{l}, T_{p}}
#' @param adultPopSizeEquilibrium adult population size at equilibrium, \eqn{Ad_{eq}}
#' @param populationGrowthRate population growth in absence of density-dependent mortality \eqn{R_{m}}
#'
#' @export
calcDensityDependentDeathRate <- function(fertility, thetaAq, timeAq, adultPopSizeEquilibrium, populationGrowthRate){
  # original, tweaked for better matrix performance
  # no mathematical change
  #prodA = (fertility*thetaAq[["E"]]*adultPopSizeEquilibrium) / (2*(populationGrowthRate-1))
  prodA = adultPopSizeEquilibrium * ((fertility*thetaAq[["E"]]) / (2*(populationGrowthRate-1)))
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
