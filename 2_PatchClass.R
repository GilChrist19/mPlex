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
# Well-mixed Patch
#######################################

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
                                      adult_male_t0, unmated_female_t0, maleReleases = NULL,
                                      femaleReleases = NULL, larvaeReleases = NULL){

                  # ID
                  private$patchID = patchID

                  # initial populations
                  private$eggs_t0 = eggs_t0
                  private$larva_t0 = larva_t0
                  private$pupa_t0 = pupa_t0
                  private$adult_male_t0 = adult_male_t0
                  private$adult_female_t0 = NULL
                  private$unmated_female_t0 = unmated_female_t0

                  # Population Lists
                  private$eggs = eggs_t0
                  private$larva = larva_t0
                  private$pupa = pupa_t0
                  private$adult_male = adult_male_t0
                  private$adult_female = NULL
                  private$unmated_female = unmated_female_t0

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
                get_adult_female = function(){private$adult_female},
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
                genericCounter = integer(length = 1),

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
                maleMigration = NULL,
                femaleMigration = NULL,
                numMigrate = integer(length = 1),
                migrationDist = NULL,

                # Mosquito Releases
                maleReleases = NULL,
                femaleReleases = NULL,
                larvaeReleases = NULL,

                # pointers
                NetworkPointer = NULL,
                DenDep = numeric(length = 1),
                death =NULL

              )# end private
)

###############################################################################
# Functions
###############################################################################

#######################################
# These are from patch-methods.R
#######################################

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
  private$unmated_female = private$unmated_female_t0

  # Reset Mosquito Releases
  private$maleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"M")
  private$femaleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"F")

}

Patch$set(which = "public",name = "reset",
          value = reset_Patch, overwrite = TRUE
)

#' Initialize Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}
#'
oneDay_initOutput_Patch <- function(){

  # Write Males
  for(mosquito in private$adult_male){
    writeLines(text = file.path(1, private$patchID, mosquito$print_male(), fsep = ","),
               con = private$NetworkPointer$get_conADM(),
               sep = "\n")
  }

  # Write Females
  for(mosquito in private$unmated_female){
    writeLines(text = file.path(1, private$patchID, mosquito$print_male(), "NULL", fsep = ","),
               con = private$NetworkPointer$get_conAF1(),
               sep = "\n")
  }

}

Patch$set(which = "public",name = "oneDay_initOutput",
          value = oneDay_initOutput_Patch, overwrite = TRUE
)

#' Write Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}
#'
oneDay_writeOutput_Patch <- function(){

  # Write Males
  for(mosquito in private$adult_male){

    writeLines(text = file.path(private$NetworkPointer$get_tNow(), private$patchID,
                                mosquito$print_male(), fsep = ","),
               con = private$NetworkPointer$get_conADM(),
               sep = "\n")
  }

  # Write Females
  for(mosquito in private$adult_female){

    writeLines(text = file.path(private$NetworkPointer$get_tNow(), private$patchID,
                                mosquito$print_female(), fsep = ","),
               con = private$NetworkPointer$get_conAF1(),
               sep = "\n")
  }

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

  if( (length(private$maleReleases) > 0) && (private$maleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()) ){
    private$adult_male = c(private$adult_male, private$maleReleases[[1]]$nuM)
    private$maleReleases[[1]] = NULL
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

  #clear unmated females. This will run every time.
  private$unmated_female = NULL

  if( (length(private$femaleReleases) > 0) && (private$femaleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()) ){
    private$unmated_female = c(private$unmated_female, private$femaleReleases[[1]]$nuF)
    private$femaleReleases[[1]] = NULL
  }

}
Patch$set(which = "public",name = "oneDay_femaleReleases",
          value = oneDay_femaleReleases_Patch, overwrite = TRUE
)

#' Release larvae in a Patch
#'
#' Based on this patch's release schedule, handle daily releases.
#'
oneDay_larvaeReleases_Patch <- function(){
  # female releases

  if( (length(private$larvaeReleases) > 0) && (private$larvaeReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()) ){
    private$larva = c(private$larva, private$larvaeReleases[[1]]$larvae)
    private$larvaeReleases[[1]] = NULL
  }

}
Patch$set(which = "public",name = "oneDay_larvaeReleases",
          value = oneDay_larvaeReleases_Patch, overwrite = TRUE
)
#######################################
# These are from patch-migration.R
#######################################

oneDay_migrationOut_Patch <- function(){

  #set empty return lists
  private$maleMigration <- private$femaleMigration <- vector(mode = "list", length = private$NetworkPointer$get_nPatch())

  #MALE
  if(length(private$adult_male)>0 && any(private$NetworkPointer$get_migrationMale(private$patchID)[-private$patchID]!=0)){
    #number who migrate: population times migration fraction
    private$numMigrate <- as.integer(x = round(x = length(private$adult_male)*private$NetworkPointer$get_migrationFractionMale(private$patchID)))

    #migration distribtuion, removing the current patch
    private$migrationDist <- JaredDirichlet(n = 1, alpha = private$NetworkPointer$get_migrationMale(private$patchID)[-private$patchID]*private$NetworkPointer$get_moveVar())

    #generate a random sample of the population who will migrate
    male_who_migrate <- hold <- sample(x = 1:length(private$adult_male), size = private$numMigrate, replace = FALSE)

    #place how many migrate where
    male_num_migrate <- rmultinom(n = 1, size = private$numMigrate, prob = private$migrationDist)

    private$genericCounter <- 1
    #loop over other patches, not your own
    for(patch in private$NetworkPointer$get_listPatch()[-private$patchID]){
      #if no mosquitoes go there, skip it
      if(male_num_migrate[private$genericCounter]==0){next}

      #get all males who migrate, then remove those males and that number from lists
      private$maleMigration[[patch]] <- private$adult_male[male_who_migrate[1:male_num_migrate[1]]]
      male_who_migrate <- male_who_migrate[-1:male_num_migrate[1]]
      male_num_migrate <- male_num_migrate[-1]
    }

    #remove all the mosquitoes who  migrated
    private$adult_male[hold] <- NULL

  }#end male migration

  #FEMALE
  if(length(private$adult_female)>0 && any(private$NetworkPointer$get_migrationFemale(private$patchID)[-private$patchID]!=0)){

    #number who migrate: population times migration fraction
    private$numMigrate <- as.integer(x = round(x = length(private$adult_female)*private$NetworkPointer$get_migrationFractionFemale(private$patchID)))

    #migration distribtuion, removing the current patch
    private$migrationDist <- JaredDirichlet(n = 1, alpha = private$NetworkPointer$get_migrationFemale(private$patchID)[-private$patchID]*private$NetworkPointer$get_moveVar())

    #generate a random sample of the population who will migrate
    female_who_migrate <- hold <- sample(x = 1:length(private$adult_female), size = private$numMigrate, replace = FALSE)

    #place how many migrate where
    female_num_migrate <- rmultinom(n = 1, size = private$numMigrate, prob = private$migrationDist)

    private$genericCounter <- 1
    #loop over other patches, not your own
    for(patch in private$NetworkPointer$get_listPatch()[-private$patchID]){

      #if no mosquitoes go there, skip it
      if(female_num_migrate[private$genericCounter]==0) next

      #get all males who migrate, then remove those males and that number from lists
      private$femaleMigration[[patch]] <- private$adult_female[female_who_migrate[1:female_num_migrate[1]]]
      female_who_migrate <- female_who_migrate[-1:female_num_migrate[1]]
      female_num_migrate <- female_num_migrate[-1]
    }

    #remove all the mosquitoes who  migrated
    private$adult_female[hold] <- NULL

  }#end female migration

}
Patch$set(which = "public",name = "oneDay_migrationOut",
          value = oneDay_migrationOut_Patch, overwrite = TRUE
)

#' Inbound Migration
#'
#' Accumulate all inbound migration to this patch.
#'
#' @param maleIn vector of inbound migration
#' @param femaleIn matrix of inbound migration
#'
oneDay_migrationIn_Patch <- function(maleIn, femaleIn){
  private$adult_male = c(private$adult_male, maleIn)
  private$adult_female = c(private$adult_female, femaleIn)
}
Patch$set(which = "public",name = "oneDay_migrationIn",
          value = oneDay_migrationIn_Patch, overwrite = TRUE
)







