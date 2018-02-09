###############################################################################
#                                    ____  _
#                          _ __ ___ |  _ \| | _____  __
#                         | '_ ` _ \| |_) | |/ _ \ \/ /
#                         | | | | | |  __/| |  __/>  <
#                         |_| |_| |_|_|   |_|\___/_/\_\
#
###############################################################################
#    ____       _       _        ____ _
#   |  _ \ __ _| |_ ___| |__    / ___| | __ _ ___ ___
#   | |_) / _` | __/ __| '_ \  | |   | |/ _` / __/ __|
#   |  __/ (_| | || (__| | | | | |___| | (_| \__ \__ \
#   |_|   \__,_|\__\___|_| |_|  \____|_|\__,_|___/___/
#
###############################################################################
#######################################
# Daily Simulation
#######################################

#' Daily Population Dynamics for a Patch
#'
#' Run population dynamics (not including migration) for this patch.
#'
oneDay_PopDynamics_Patch <- function(){

  ################
  # AGE
  ################
  self$oneDay_Age(Population = private$eggs)
  self$oneDay_Age(Population = private$larva)
  self$oneDay_Age(Population = private$pupa)
  self$oneDay_Age(Population = private$adult_male)
  self$oneDay_Age(Population = private$adult_female)

  ################
  # DEATH
  ################
  self$oneDay_EggDeath()
  self$oneDay_LarvalDeath()
  self$oneDay_PupaDeath()
  self$oneDay_AdultDeath()

  ################
  # MATURE
  ################
  self$oneDay_PupaMaturation()
  self$oneDay_LarvaMaturation()
  self$oneDay_EggMaturation()

  ################
  # MATE
  ################
  self$oneDay_Mate()

  ################
  # LAY EGGS
  ################
  self$oneDay_Reproduction()

  ################
  # Releases
  ################
  self$oneDay_Releases()

}
Patch$set(which = "public",name = "oneDay_PopDynamics",
          value = oneDay_PopDynamics_Patch, overwrite = TRUE
)

###############################################################################
# Functions
###############################################################################
#######################################
# Age Function
#######################################
#' Daily Aging
#'
#' Age all Mosquitoes in the population by 1 day
#'
#' @param Population The population to age
#'
oneDay_Age_Patch <- function(Population = NULL){

  for(critter in Population){
    critter$age_one_day()
  }

}
Patch$set(which = "public",name = "oneDay_Age",
          value = oneDay_Age_Patch, overwrite = TRUE
)

#######################################
# Death Functions
#######################################
#' Daily Egg Death
#'
#' Kill eggs based on a chance of death specified by mu("E")
#'
oneDay_EggDeath_Patch <- function(){

  #calculate other dependent factors

  #calculate death
  private$death <- rbinom(n = length(private$eggs),
                    size = 1,
                    prob = private$NetworkPointer$get_mu("E"))

  private$eggs[as.logical(private$death)] <- NULL
}
Patch$set(which = "public",name = "oneDay_EggDeath",
          value = oneDay_EggDeath_Patch, overwrite = TRUE
)

#' Daily Larvae Death
#'
#' Kill larvae based on a chance of death specified by mu("L") and a density
#' dependend parameter alpha
#'
oneDay_LarvalDeath_Patch <- function(){

  #calculate density-dependent portion
  stageDur <- private$NetworkPointer$get_stageTime(stage = "L")
  alpha = private$NetworkPointer$get_alpha(private$patchID)

  private$DenDep <- (alpha/(alpha+length(private$larva)) )^(1/stageDur)

  #calculate death
  # This is done as 1-probLife
  private$death <- rbinom(n = length(private$larva),
                  size = 1,
                  prob = 1-(private$DenDep*(1-private$NetworkPointer$get_mu("L"))))

  private$larva[as.logical(private$death)] <- NULL
}
Patch$set(which = "public",name = "oneDay_LarvalDeath",
          value = oneDay_LarvalDeath_Patch, overwrite = TRUE
)

#' Daily Pupa Death
#'
#' Density independent pupa death specified by mu("P")
#'
oneDay_PupaDeath_Patch <- function(){

  #calculate other dependent factors

  #calculate death
  private$death <- rbinom(n = length(private$pupa),
                    size = 1,
                    prob = private$NetworkPointer$get_mu("P"))

  private$pupa[as.logical(private$death)] <- NULL
}
Patch$set(which = "public",name = "oneDay_PupaDeath",
          value = oneDay_PupaDeath_Patch, overwrite = TRUE
)

#' Daily Adult Death
#'
#' Density independent adult death specified by mu("A")
#'
oneDay_AdultDeath_Patch <- function(){

  #MALE
  #calculate other dependent factors

  #calculate death
  private$death <- rbinom(n = length(private$adult_male),
                    size = 1,
                    prob = private$NetworkPointer$get_mu("A"))

  private$adult_male[as.logical(private$death)] <- NULL


  #FEMALE
  #calculate other dependent factors

  #calculate death
  private$death <- rbinom(n = length(private$adult_female),
                    size = 1,
                    prob = private$NetworkPointer$get_mu("A"))

  private$adult_female[as.logical(private$death)] <- NULL
}
Patch$set(which = "public",name = "oneDay_AdultDeath",
          value = oneDay_AdultDeath_Patch, overwrite = TRUE
)

#######################################
# Mature Functions
#######################################

#' Daily Egg Maturation
#'
#' Eggs mature into larva. The average time spent as an egg is specified in
#' stageTime("E"). This is used to parametrize a lognormal function with
#' standard deviation log(1.2), slightly right-skewing the distribution.
#'
oneDay_EggMature_Patch <- function(){

  private$genericCounter <- length(private$eggs)
  #only do if there are eggs
  if(private$genericCounter > 0){

    #set mean age and sd. THESE NEED TO BE CHECKED or PROVEN
    private$meanAge <- log(x = private$NetworkPointer$get_stageTime(stage = "E"))
    private$sdAge <- log(x = 1.2)

    #draw all the private$ages since this can probably be done.
    # if not, do it one at a time in the loop.
    # This could maybe become a patch variable, so we don't reallocate space.??
    private$ages <- rlnorm(n = private$genericCounter, meanlog = private$meanAge, sdlog = private$sdAge)

    #Get ages
    private$matured <- vapply(X = private$eggs, FUN = '[[', "age", FUN.VALUE = integer(length = 1L))

    #see who does mature
    private$matured <- private$matured > private$ages

    #put matured eggs into larva
    private$larva <- c(private$larva, private$eggs[private$matured])

    #remove matured eggs from eggs
    private$eggs[private$matured] <- NULL

  }

}
Patch$set(which = "public",name = "oneDay_EggMaturation",
          value = oneDay_EggMature_Patch, overwrite = TRUE
)

#' Daily Larva Maturation
#'
#' Larva mature into pupa. The average time spent as a larva is specified in
#' stageTime("L"). This is used to parametrize a lognormal function with
#' standard deviation log(1.2), slightly right-skewing the distribution.
#'
oneDay_LarvaMature_Patch <- function(){

  private$genericCounter <- length(private$larva)
  #only do if things are there
  if(private$genericCounter > 0){

    #set mean age and sd. THESE NEED TO BE CHECKED or PROVEN
    private$meanAge <- log(x = sum(private$NetworkPointer$get_stageTime(stage = c("E", "L") ) ) )
    private$sdAge <- log(x = 1.2)

    #draw all the private$ages since this can probably be done.
    # if not, do it one at a time in the loop.
    # This could maybe become a patch variable, so we don't reallocate space.??
    private$ages <- rlnorm(n = private$genericCounter, meanlog = private$meanAge, sdlog = private$sdAge)

    #get ages
    private$matured <- vapply(X = private$larva, FUN = '[[', "age", FUN.VALUE = integer(length = 1L))

    #see who does mature
    private$matured <- private$matured > private$ages

    #put matured eggs into larva
    private$pupa <- c(private$pupa, private$larva[private$matured])

    #remove matured eggs from eggs
    private$larva[private$matured] <- NULL

  }

}
Patch$set(which = "public",name = "oneDay_LarvaMaturation",
          value = oneDay_LarvaMature_Patch, overwrite = TRUE
)

#' Daily Pupa Maturation
#'
#' Pupae mature into adults. The average time spent as a pupa is specified in
#' stageTime("P"). This is used to parametrize a lognormal function with
#' standard deviation log(1.2), slightly right-skewing the distribution. A binomial
#' distribtution with p=0.5 specifies the sex of maturing pupae. If male, pupa are
#' immediately added to the adult_male population. If female, pupa are added to
#' the unmated_female population.
#'
oneDay_PupaMature_Patch <- function(){

  private$genericCounter <- length(private$pupa)
  #only do things if there are pupa
  if(private$genericCounter > 0){

    #set mean age and sd. THESE NEED TO BE CHECKED or PROVEN
    private$meanAge <- log(x = sum(private$NetworkPointer$get_stageTime(stage = c("E", "L", "P")) ) )
    private$sdAge <- log(x = 1.2)

    #draw all the private$ages since this can probably be done.
    # if not, do it one at a time in the loop.
    # This could maybe become a patch variable, so we don't reallocate space.??
    private$ages <- rlnorm(n = private$genericCounter, meanlog = private$meanAge, sdlog = private$sdAge)

    #get ages
    private$matured <- vapply(X = private$pupa, FUN = '[[', "age", FUN.VALUE = integer(length = 1L))

    #see who does mature
    private$matured <- private$matured > private$ages

    #binomial over sex. Need to add genotype deviance
    sex <- as.logical(rbinom(n = sum(private$matured), size = 1, prob = 0.5))

    #if true, make an unmated female, otherwise, make male
    private$unmated_female <- c(private$unmated_female, private$pupa[private$matured][sex])
    private$adult_male <- c(private$adult_male, private$pupa[private$matured][!sex])

    #remove matured eggs from eggs
    private$pupa[private$matured] <- NULL

  }#end if

}
Patch$set(which = "public",name = "oneDay_PupaMaturation",
          value = oneDay_PupaMature_Patch, overwrite = TRUE
)

#######################################
# Mate Function
#######################################

#' Daily Mating
#'
#' Freshly matured pupa that become females and female releases exist as unmated
#' females. This function mates unmated_females with any member of the current male
#' population, with the possibility that males can mate multiple times per mating.
#' After mating, unmated_females are put into the general adult_female population.
#'
oneDay_Mate_Patch <- function(){

  #number of unwed females
  private$numUnweds <- length(private$unmated_female)
  private$numMates <- length(private$adult_male)

  if(private$numUnweds != 0 && private$numMates != 0){

    #get males for mates, randomly sampled from population
    private$mates <- vapply(X = private$adult_male, FUN = '[[', "genotype", FUN.VALUE = character(length = 1L))

    private$mates <- sample(x = private$mates, size = private$numUnweds, replace = TRUE)

    #set mates
    for(i in 1:private$numUnweds){
      private$unmated_female[[i]]$set_mate(matGen=private$mates[[i]])
    }

    #add to adult females, clear unmated females
    private$adult_female <- c(private$adult_female,private$unmated_female)
    private$unmated_female = NULL
  }

}
Patch$set(which = "public",name = "oneDay_Mate",
          value = oneDay_Mate_Patch, overwrite = TRUE
)

#######################################
# Lay Eggs Function
#######################################

#' Daily Egg Laying
#'
#' The number of eggs laid per female is a Poisson distribution with a mean of
#' beta, the fertility. A multinomial is then used to distribute the number of
#' eggs laid over the offspring distribution (see \code{\link{DaisyOffspring}},
#' \code{\link{MultiplexOffspring_mLoci}}, or \code{\link{MultiplexOffspring_oLocus}}).
#'
oneDay_Reproduction_Patch <- function(){

  if(length(private$adult_female)!=0L){
    #get number of new mosquitoes in each batch
    # reusing mates from above. It was a variable length character vector.
    private$mates <- vapply(X = private$adult_female, FUN = '[[', "genotype", FUN.VALUE = character(length = 1L))
    private$mates <- vapply(X = private$mates,
                            FUN = function(y){calc_parameter(paramList = private$NetworkPointer$get_s(), genotype = y)},
                            FUN.VALUE = numeric(1))
    private$mates <- rpois(n = length(private$adult_female), lambda = private$mates*private$NetworkPointer$get_beta())

    #Generate new mosquitoe population and set counter
    private$newEggs <- vector(mode = "list", length = sum(private$mates))
    private$genericCounter = 1L

    for(critter in 1:length(private$adult_female)){

      #This fills offpsring list, a list of 2 lists: Alleles, Probabilities
      # This function is set during initialization
      private$offspring <- self$offspringDistribution(fGen = private$adult_female[[critter]][["genotype"]],
                                 mGen = private$adult_female[[critter]][["mate"]],
                                 reference = private$NetworkPointer$get_reference())

      #This generates an egg distribution
      private$eggNumber <- rmultinom(n = 1,
                             size = private$mates[critter],
                             prob = private$offspring$Probabilities)

      #loop over each genotype
      for(gen in 1:length(private$offspring$Alleles)){
        #skip if there are no mosquitoes of this genotype
        if(private$eggNumber[gen]==0L){next}

        #loop over number of mosquitoes of that genotype
        for(num in 1:private$eggNumber[gen]){
          #create new mosquito
          private$newEggs[[private$genericCounter]] <- Mosquito(genotype = private$offspring$Alleles[gen], age = 0L)
          #increment counter
          private$genericCounter = private$genericCounter + 1L

        }# end egg number loop
      }# end allele loop
    }#end critter loop

    # add new eggs to the pile
    private$eggs <- c(private$eggs, private$newEggs)

  }#end if statement
}#end function
Patch$set(which = "public",name = "oneDay_Reproduction",
          value = oneDay_Reproduction_Patch, overwrite = TRUE
)
