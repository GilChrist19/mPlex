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

  cat("\nEggs:", length(private$eggs), "\n")
  cat("Larva:", length(private$larva), "\n")
  cat("Pupa:", length(private$pupa), "\n")
  cat("Male:", length(private$adult_male), "\n")
  cat("Female:", length(private$adult_female), "\n")

  ################
  # DEATH
  ################
  #self$oneDay_EggDeath()
  #self$oneDay_LarvalDeath()
  #self$oneDay_PupaDeath()
  #self$oneDay_AdultDeath()

  ################
  # MATURE
  ################
  #self$oneDay_PupaMaturation()
  #self$oneDay_LarvaMaturation()
  #self$oneDay_EggMaturation()

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
  self$oneDay_larvaeReleases()
  self$oneDay_maleReleases()
  self$oneDay_femaleReleases()
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

oneDay_LarvalDeath_Patch <- function(){

  #calculate density-dependent portion
  stageDur <- private$NetworkPointer$get_timeAq("L")
  alpha = private$NetworkPointer$get_alpha(private$patchID)

  private$DenDep <- (alpha/(alpha+length(private$larva)) )^(1/stageDur)

  #calculate death
  ### CAN I MULTIPLY THIS?????? ########
  private$death <- rbinom(n = length(private$larva),
                  size = 1,
                  prob = private$DenDep*private$NetworkPointer$get_mu("L"))

  private$larva[as.logical(private$death)] <- NULL
}
Patch$set(which = "public",name = "oneDay_LarvalDeath",
          value = oneDay_LarvalDeath_Patch, overwrite = TRUE
)

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
oneDay_EggMature_Patch <- function(){

  num <- length(private$eggs)
  #only do if there are eggs
  if(num > 0){

    #set mean age and sd. THESE NEED TO BE CHECKED or PROVEN
    meanAge <- log(x = private$NetworkPointer$get_timeAq(stage = "E"))
    sdAge <- log(x = 1.2)

    #draw all the ages since this can probably be done.
    # if not, do it one at a time in the loop.
    # This could maybe become a patch variable, so we don't reallocate space.??
    ages <- rlnorm(n = num, meanlog = meanAge, sdlog = sdAge)
    matured <- logical(length = num)

    #loop over all eggs
    #This maybe could be  vectorized in R, getting all ages is the only problem.
    for(i in 1:num){

      #check if they mature, move it up, then delete egg
      if(ages[i] <= private$eggs[[i]]$get_age()){
        private$larva <- c(private$larva, private$eggs[[i]])
        matured[i] <- TRUE
      }

    }

    private$eggs[matured] <- NULL
  }

}
Patch$set(which = "public",name = "oneDay_EggMaturation",
          value = oneDay_EggMature_Patch, overwrite = TRUE
)

oneDay_LarvaMature_Patch <- function(){

  num <- length(private$larva)
  #only do if things are there
  if(num > 0){

    #set mean age and sd. THESE NEED TO BE CHECKED or PROVEN
    meanAge <- log(x = private$NetworkPointer$get_timeAq(stage = "E")+
                     private$NetworkPointer$get_timeAq(stage = "L"))
    sdAge <- log(x = 1.2)

    #draw all the ages since this can probably be done.
    # if not, do it one at a time in the loop.
    # This could maybe become a patch variable, so we don't reallocate space.??
    ages <- rlnorm(n = num, meanlog = meanAge, sdlog = sdAge)
    matured <- logical(length = num)

    #loop over all eggs
    #This maybe could be  vectorized in R, getting all ages is the only problem.
    for(i in 1:num){

      #check if they mature, move it up, then delete egg
      if(ages[i] <= private$larva[[i]]$get_age()){
        private$pupa <- c(private$pupa, private$larva[[i]])
        matured[i] <- TRUE
      }

    }

    #remove matured larva
    private$larva[matured] <- NULL
  }

}
Patch$set(which = "public",name = "oneDay_LarvaMaturation",
          value = oneDay_LarvaMature_Patch, overwrite = TRUE
)

oneDay_PupaMature_Patch <- function(){

  num <- length(private$pupa)
  #only do things if there are pupa
  if(num > 0){

    #set mean age and sd. THESE NEED TO BE CHECKED or PROVEN
    meanAge <- log(x = private$NetworkPointer$get_timeAq())
    sdAge <- log(x = 1.2)

    #draw all the ages since this can probably be done.
    # if not, do it one at a time in the loop.
    # This could maybe become a patch variable, so we don't reallocate space.??
    ages <- rlnorm(n = num, meanlog = meanAge, sdlog = sdAge)
    matured <- logical(length = num)

    #loop over all eggs
    #This maybe could be  vectorized in R, getting all ages is the only problem.
    for(i in 1:num){

      #check if they mature, move it up, then delete egg
      if(ages[i] <= private$pupa[[i]]$get_age()){

        #binomial over sex. Need to add genotype deviance
        sex <- as.logical(x = rbinom(n = 1, size = 1, prob = 0.5))

        if(sex){
          #if true, make an unmated female
          private$unmated_female <- c(private$unmated_female, private$pupa[[i]])
        } else {
          #make adult male
          private$adult_male <- c(private$adult_male, private$pupa[[i]])
        }

        matured[i] <- TRUE
      }

    }

    #remove matured pupa
    private$pupa[matured] <- NULL
  }

}
Patch$set(which = "public",name = "oneDay_PupaMaturation",
          value = oneDay_PupaMature_Patch, overwrite = TRUE
)

#######################################
# Mate Function
#######################################

oneDay_Mate_Patch <- function(){

  #number of unwed females
  numUnweds <- length(private$unmated_female)

  if(numUnweds != 0){

    mates <- vector(mode = "list", length = length(private$adult_male))

    #get males for mates, randomly sampled from population
    for(i in 1:length(private$adult_male)){
      mates[[i]] <- private$adult_male[[i]]$get_genotype()
    }

    mates <- sample(x = mates, size = numUnweds, replace = TRUE)

    #set mates
    for(i in 1:numUnweds){
      private$unmated_female[[i]]$set_mate(mate=mates[[i]])
    }

    #add to adult females
    private$adult_female <- c(private$adult_female,private$unmated_female)
  }

  #migration would null unmated_females. currently, this may just grow.
}
Patch$set(which = "public",name = "oneDay_Mate",
          value = oneDay_Mate_Patch, overwrite = TRUE
)





































#######################################
# Lay Eggs Function
#######################################
oneDay_Reproduction_Patch <- function(){

  for(critter in private$adult_female){

    #This produces a list of 2 lists: Alleles, Probabilities
    # This function is set during initialization
    offspring <- self$offspringDistribution(fGen = critter$get_genotype(),
                                      mGen = critter$get_mate(),
                                      reference = private$NetworkPointer$get_reference())

    #This generates an egg distribution
    eggNumber <- rmultinom(n = 1,
                           size = rpois(n = 1, lambda = private$NetworkPointer$get_beta()),
                           prob = offspring$Probabilities)

    #Generate new mosquitoes, put them in eggs class
    private$eggs <- c(private$eggs,
                      CreateMosquitoes_Eggs(genMos = offspring$Alleles,
                                            numMos = eggNumber)
                      )

  }#end loop
}#end function
Patch$set(which = "public",name = "oneDay_Reproduction",
          value = oneDay_Reproduction_Patch, overwrite = TRUE
)



