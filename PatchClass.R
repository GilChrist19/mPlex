#Patch Class
#Jared Bennett

#heavily borrows from patch class in MGDrive package.

#################################################
# Well-mixed Patch
#################################################

# A patch is defined as a fully mixed population of mosquitoes habitating in a given place.
#   Each patch can be mapped to a node within a network of connected populations in
#   which adult mosquitoes can migrate. This is a fully independent and self-contained unit,
#   so it can be treated as an atomic entity.
Patch <- R6::R6Class(classname = "Patch",
            portable = TRUE,
            cloneable = FALSE,
            lock_class = FALSE,
            lock_objects = FALSE,

            # public memebers
            public = list(

                # Constructor
                initialize = function(patchID, simTime, eggs_t0, larva_t0, pupa_t0,
                                      adult_male_t0, Adult_female_t0, maleReleases = NULL,
                                      femaleReleases = NULL, larvaeReleases = NULL){

                  # ID
                  private$patchID = patchID

                  # initial populations
                  private$eggs_t0 = eggs_t0
                  private$larva_t0 = larva_t0
                  private$pupa_t0 = pupa_t0
                  private$adult_male_t0 = adult_male_t0
                  private$adult_female_t0 = Adult_female_t0
                  private$unmated_female_t0 = NULL

                  # Population Lists
                  private$eggs = eggs_t0
                  private$larva = larva_t0
                  private$pupa = pupa_t0
                  private$adult_male = adult_male_t0
                  private$adult_female = Adult_female_t0
                  private$unmated_female = NULL

                  # Mosquito Releases
                  private$maleReleases = maleReleases
                  private$femaleReleases = femaleReleases
                  private$larvaeReleases = larvaeReleases

                },# end constructor

                # patch-level parameters
                get_patchID = function(){return(private$patchID)},

                # population-level parameters
                get_egg = function(){private$eggs},
                get_larva = function(){private$larva},
                get_pupa = function(){private$pupa},
                get_adult_male = function(){private$adult_male},
                get_Adult_female = function(){private$adult_female},
                get_unmated_female = function(){private$unmated_female},

                # migration
                get_maleMigration =   function(){private$maleMigration},
                get_femaleMigration = function(){private$femaleMigration},

                # pointers
                set_NetworkPointer = function(NetPoint){private$NetworkPointer = NetPoint},
                get_NetworkPointer = function(){private$NetworkPointer}
              ),# end public

            # private members
            private = list(

                # patch-level parameters
                patchID = NULL,

                # initial populations
                eggs_t0 = NULL,
                larva_t0 = NULL,
                pupa_t0 = NULL,
                adult_male_t0 = NULL,
                adult_female_t0 = NULL,
                unmated_female_t0 = NULL,

                # populations
                eggs = NULL,
                larva = NULL,
                pupa = NULL,
                adult_male = NULL,
                adult_female = NULL,
                unmated_female = NULL,

                # migration
                maleMigration = NULL, # nGenotypes X nPatch matrix
                femaleMigration = NULL, # nGenotypes X nGenotypes X nPatch array

                # pointers
                NetworkPointer = NULL

              )# end private
)

###############################################################################
# Functions
###############################################################################

#######################################
# These are from patch-methods.R
#######################################


## This may be done ##

#' Reset Patch to Initial Conditions
#'
#' Resets a patch to its initial configuration so that a new one does not have to be created and
#'     allocated in the network (for Monte Carlo simulation).
#'
reset_Patch <- function(){

  cat("reset patch ",private$patchID,"\n",sep="")

  # reset initial population size
  private$eggs = private$eggs_t0
  private$larva = private$larva_t0
  private$pupa = private$pupa_t0
  private$adult_male = private$adult_male_t0
  private$adult_female = private$adult_female_t0
  private$unmated_female = unmated_female_t0

  # Reset Mosquito Releases
  private$maleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"M")
  private$femaleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"F")

}

Patch$set(which = "public",name = "reset",
          value = reset_Patch, overwrite = TRUE
)
























####### THESE NEED UPDATED!!!!!!!!!!!!!!
####### MAKE IT WORk on lists of critters




#' Initialize Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}
#'
oneDay_initOutput_Patch <- function(){

  ##this maybe
  for(i in 1:length(private$adult_male)){

        ## use maybe adult_male_t0
    writeLines(text = file.path(1, private$patchID, private$adult_male[[i]]$print_info(), fsep = "\t"),
               con = private$NetworkPointer$get_conADM(),
               sep = "\n")
  }

  for(i in 1:length(private$adult_female)){

    writeLines(text = file.path(1, private$patchID, private$adult_female[[i]]$print_info(), fsep = "\t"),
               con = private$NetworkPointer$get_conAF1(),
               sep = "\n")
  }
  ##

  # write males
  ADMout = paste0(c(1,private$patchID,private$ADM[[1]]),collapse = ",")
  writeLines(text = ADMout,con = private$NetworkPointer$get_conADM(),sep = "\n")

  # write females
  AF1out = paste0(c(1,private$patchID,c(t(private$AF1[[1]]))),collapse = ",")
  writeLines(text = AF1out,con = private$NetworkPointer$get_conAF1(),sep = "\n")

}

Patch$set(which = "public",name = "oneDay_initOutput",
          value = oneDay_initOutput_Patch, overwrite = TRUE
)













#' Write Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}
#'
oneDay_writeOutput_Patch <- function(){

  ## test this
  for(i in 1:length(private$adult_male)){

    writeLines(text = file.path(private$NetworkPointer$get_tNow(), private$patchID,
                                private$adult_male[[i]]$print_info(), fsep = "\t"),
               con = private$NetworkPointer$get_conADM(),
               sep = "\n")
  }

  for(i in 1:length(private$adult_female)){

    writeLines(text = file.path(private$NetworkPointer$get_tNow(), private$patchID,
                                private$adult_female[[i]]$print_info(), fsep = "\t"),
               con = private$NetworkPointer$get_conAF1(),
               sep = "\n")
  }
  ##




  tNow = private$NetworkPointer$get_tNow()

  # write males
  ADMout = paste0(c(tNow,private$patchID,private$ADM[[tNow]]),collapse = ",")
  writeLines(text = ADMout,con = private$NetworkPointer$get_conADM(),sep = "\n")

  # write females
  AF1out = paste0(c(tNow,private$patchID,c(t(private$AF1[[tNow]]))),collapse = ",")
  writeLines(text = AF1out,con = private$NetworkPointer$get_conAF1(),sep = "\n")

}

Patch$set(which = "public",name = "oneDay_writeOutput",
          value = oneDay_writeOutput_Patch, overwrite = TRUE
)




#######################################
# These are from patch-releases.R
#######################################



#' Release Male Mosquitoes in a Patch
#'
#' Based on this patch's release schedule, handle daily releases.
#'
oneDay_maleReleases_Patch <- function(){
  # male releases

  ##test this
  if( (length(private$maleReleases) > 0) && (private$maleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()) ){
    private$adult_male = c(private$adult_male, private$maleReleases[[1]])
    private$maleReleases[[1]] = NULL
  }
  ##


  if(length(private$maleReleases) > 0){
    # browser()
    if(private$maleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()){
      private$ADMnew = private$admSurvival + private$admPupating + private$maleReleases[[1]]$nuM
      private$maleReleases[[1]] = NULL
    } else {
      private$ADMnew = private$admSurvival + private$admPupating
    }

  } else {
    private$ADMnew = private$admSurvival + private$admPupating
  }





}

Patch$set(which = "public",name = "oneDay_maleReleases",
          value = oneDay_maleReleases_Patch, overwrite = TRUE
)








#' Release Female Mosquitoes in a Patch
#'
#' Based on this patch's release schedule, handle daily releases.
#'
oneDay_femaleReleases_Patch <- function(){
  # female releases

  ##test this
  if( (length(private$femaleReleases) > 0) && (private$femaleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()) ){
    private$unmated_female = c(private$unmated_female, private$femaleReleases[[1]])
    private$femaleReleases[[1]] = NULL
  }
  ##




  if(length(private$femaleReleases) > 0){
    # browser()
    if(private$femaleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()){
      private$af1Pupation = private$af1Pupation + private$femaleReleases[[1]]$nuF
      private$femaleReleases[[1]] = NULL
    }

  }

}

Patch$set(which = "public",name = "oneDay_femaleReleases",
          value = oneDay_femaleReleases_Patch, overwrite = TRUE
)




#######################################
# These are from patch-migration.R
#######################################



## figure out this one ##

#' Stochastic Oubound Migration
#'
#' Stochastic model of migration of \code{AF1new} females from this patch, fills up the \code{femaleMigration} array.
#' Migration is modeled as a Dirichlet-Multinomial process parameterized by \code{moveVar} multiplied by the row corresponding to this
#' patch from the stochastic matrix. A Dirichlet distributed random variate is sampled from \code{\link[MCMCpack]{rdirichlet}} according to that
#' parameter vector and then movement is sampled from \code{\link[stats]{rmultinom}}.
#'
oneDay_migrationOut_stochastic_Patch <- function(){

  # values that are used more than once
  mVar <- private$NetworkPointer$get_moveVar()
  pNum <- private$patchID

  # Males
  # return matrix
  private$maleMigration <-MaleMigrationC(Population = private$ADMnew,
                                        migrationPoint = private$NetworkPointer$get_migrationMaleRow(pNum)*mVar)


  # Females
  # return array
  private$femaleMigration <- FemaleMigrationC(Population = private$AF1new,
                                             migrationPoint = private$NetworkPointer$get_migrationFemaleRow(pNum)*mVar)

}







##THis one might be done ##

#' Inbound Migration
#'
#' Accumulate all inbound migration to this patch.
#'
#' @param maleIn vector of inbound migration
#' @param femaleIn matrix of inbound migration
#'
oneDay_migrationIn_Patch <- function(maleIn, femaleIn){
  private$adult_male = maleIn
  private$adult_female = femaleIn
}

Patch$set(which = "public",name = "oneDay_migrationIn",
          value = oneDay_migrationIn_Patch, overwrite = TRUE
)







