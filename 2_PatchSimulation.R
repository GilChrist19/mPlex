

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

  
  
  ################
  # MATURE
  ################
  
  
  
  ################
  # 
  ################





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








