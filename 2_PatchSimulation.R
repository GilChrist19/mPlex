

###############################################################################
# Daily Simulation
###############################################################################

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
  
  
  
  ################
  # MATE
  ################
  self$oneDay_Mate()




}
Patch$set(which = "public",name = "oneDay_PopDynamics",
          value = oneDay_PopDynamics_Patch, overwrite = TRUE
)

###############################################################################
# Age Function
###############################################################################
oneDay_Age_Patch <- function(Population = NULL){
  
  for(critter in Population){
    critter$age_one_day()
  }
  
}
Patch$set(which = "public",name = "oneDay_Age",
          value = oneDay_Age_Patch, overwrite = TRUE
)

###############################################################################
# Death Functions
###############################################################################
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










###############################################################################
# Mate Function
###############################################################################

oneDay_Mate_Patch <- function(){

  #number of unwed females
  numUnweds <- length(private$unmated_female)
  
  if(numUnweds != 0){
    
    mates <- vector(mode = "list", length = numUnweds)
    
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

