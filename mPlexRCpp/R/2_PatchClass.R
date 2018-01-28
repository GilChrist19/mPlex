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

#' Patch Class Definition
#'
#' A Patch represents a single, well-mixed population of \code{\link{Mosquito}}
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#'
#' @section **Constructor**:
#'  * patchID: Integer ID of this patch
#'  * simTime: Maximum time of simulation
#'  * eggs_t0: List of Mosquito objects as the initial egg population, \eqn{L_{eq}}
#'  * larva_t0: List of Mosquito objects as the initial larval population, \eqn{L_{eq}}
#'  * pupa_t0: List of Mosquito objects as the initial pupae population, \eqn{L_{eq}}
#'  * adult_male_t0: List of Mosquito objects as the initial adult male population, \eqn{Ad_{eq}}
#'  * unmated_female_t0: List of Mosquito objects as the initial adult female population, \eqn{Ad_{eq}}
#'  * maleReleases: Male release schedule for this patch, see \code{\link{Release_basicRepeatedReleases}}
#'  * femaleReleases: Female release schedule for this patch, see \code{\link{Release_basicRepeatedReleases}}
#'  * larvaeReleases: Larvae release schedule for this patch, see \code{\link{Release_basicRepeatedReleases}}
#'
#' @section **Methods**:
#'  * get_patchID: Return ID of current patch
#'  * get_egg: Return current egg population
#'  * get_larva: Return current larva population
#'  * get_pupa: Return current pupa population
#'  * get_adult_male: Return current adult male population
#'  * get_adult_female: Return current adult female population
#'  * get_unmated_female: Return current adult, unmated female population
#'  * get_maleMigration: Return migratory males (nPatch length list of lists)
#'  * get_femaleMigration: Return migratory females (nPatch length list of lists)
#'  * set_NetworkPointer: Set a referenced to the enclosing \code{\link{Network}} object
#'  * get_NetworkPointer: Return a reference to the enclosing \code{\link{Network}} object
#'  * reset: see \code{\link{reset_Patch}}
#'  * oneDay_initOutput: see \code{\link{oneDay_initOutput_Patch}}
#'  * oneDay_writeOutput: see \code{\link{oneDay_writeOutput_Patch}}
#'  * oneDay_maleReleases: see \code{\link{oneDay_maleReleases_Patch}}
#'  * oneDay_femaleReleases: see \code{\link{oneDay_femaleReleases_Patch}}
#'  * oneDay_larvaeReleases: see \code{\link{oneDay_larvaeReleases_Patch}}
#'  * oneDay_migrationOut: see \code{\link{oneDay_migrationOut_Patch}}
#'  * oneDay_migrationIn: see \code{\link{oneDay_migrationIn_Patch}}
#'  * oneDay_PopDynamics: see \code{\link{oneDay_PopDynamics_Patch}}
#'  * oneDay_Age: see \code{\link{oneDay_Age_Patch}}
#'  * oneDay_EggDeath: see\code{\link{oneDay_EggDeath_Patch}}
#'  * oneDay_LarvalDeath: see \code{\link{oneDay_LarvalDeath_Patch}}
#'  * oneDay_PupaDeath: see \code{\link{oneDay_PupaDeath_Patch}}
#'  * oneDay_AdultDeath: see \code{\link{oneDay_AdultDeath_Patch}}
#'  * oneDay_EggMaturation: see \code{\link{oneDay_EggMature_Patch}}
#'  * oneDay_LarvaMaturation: see \code{\link{oneDay_LarvaMature_Patch}}
#'  * oneDay_PupaMature: see \code{\link{oneDay_PupaMature_Patch}}
#'  * oneDay_Mate: see \code{\link{oneDay_Mate_Patch}}
#'  * oneDay_Reproduction: See \code{\link{oneDay_Reproduction_Patch}}
#'
#' @section **Fields**:
#'  * patchID: Integer ID of current patch
#'  * genericCounter: Generic counter used for indexing
#'  * eggs_t0: List of Mosquitoes comprising the initial egg population
#'  * larva_t0: List of Mosquitoes comprising the initial larval population
#'  * pupa_t0: List of Mosquitoes comprising the initial pupa population
#'  * adult_male_t0: List of Mosquitoes comprising the initial adult male population
#'  * adult_female_t0: List of Mosquitoes comprising the initial adult female population
#'  * unmated_female_t0: List of Mosquitoes comprising the initial unmated female population
#'  * eggs: List of Mosquitoes comprising the current egg population
#'  * larva: List of Mosquitoes comprising the current larval population
#'  * pupa: List of Mosquitoes comprising the current pupa population
#'  * adult_male: List of Mosquitoes comprising the current adult male population
#'  * adult_female: List of Mosquitoes comprising the current adult female population
#'  * unmated_female: List of Mosquitoes comprising the current unmated female population
#'  * maleMigration: List of outbound male Mosquitoes, length nPatch
#'  * femaleMigration: List of outbound female Mosquitoes, length nPatch
#'  * numMigration: Holder for the number of migratory Mosquitoes, integer
#'  * migrationDist: Holder for the distribution of Mosquitoes over patches, vector length nPatch-1
#'  * whoMigrate: Holder for which mosquitoes migrate, vector of indices
#'  * numMigrateWhere: Holder for where and how many mosquitoes migrate, vector of integers
#'  * maleReleases: Male release schedule and list of Mosquitoes to release
#'  * femaleReleases: female release schedule and list of Mosquitoes to release
#'  * larvaeReleases: Larval release schedule and list of Mosquitoes to release
#'  * NetworkPointer: a reference to enclosing \code{\link{Network}}
#'  * DenDep: Density-dependent parameter for larvae
#'  * death: Vector of T/F death at each stage. Slowly grows with population size
#'  * meanAge: Holder for maturation functions
#'  * sdAge: Holder for maturation functions
#'  * ages: Holder for maturation functions
#'  * matured: Holder for maturation functions
#'  * numUnweds: Number of unmated females. Holder for mating function
#'  * numMates: Number of adult males. Holder for mating function
#'  * mates: Character vector of male genotypes. Holder for mating function
#'  * offspring: List(Alleles, Probabilities) specifying new offspring genotype and distribution.
#'  * eggNumber: Integer vector of how many eggs to lay. Holder for reproduction function
#'  * newEggs: list of new Mosquitoes laid. Holder for reproduction function
#'
#' @export
Patch <- R6::R6Class(classname = "Patch",
            portable = TRUE,
            cloneable = FALSE,
            lock_class = FALSE,
            lock_objects = FALSE,
            class = FALSE,

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
                whoMigrate = NULL,
                numMigrateWhere = NULL,

                # Mosquito Releases
                maleReleases = NULL,
                femaleReleases = NULL,
                larvaeReleases = NULL,

                # pointers
                NetworkPointer = NULL,

                # daily death storage
                DenDep = numeric(length = 1),
                death = NULL,

                # daily mature storage
                meanAge = numeric(length = 1),
                sdAge = numeric(length = 1),
                ages = NULL,
                matured = NULL,

                # daily mating storage
                numUnweds = integer(length = 1),
                numMates = integer(length = 1),
                mates = NULL,

                # daily reproduction and reproduction functions
                offspring = NULL,
                eggNumber = NULL,
                newEggs = NULL

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
#' Resets a patch to its initial configuration so that a new one does not have
#' to be created and allocated in the network (for Monte Carlo simulation).
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
#' Writes initial output to the text connections specified in the enclosing
#' \code{\link{Network}}
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
               con = private$NetworkPointer$get_conADF(),
               sep = "\n")
  }

}
Patch$set(which = "public",name = "oneDay_initOutput",
          value = oneDay_initOutput_Patch, overwrite = TRUE
)

#' Write Output from Focal Patch
#'
#' Writes daily output to the text connections specified in the enclosing
#' \code{\link{Network}}
#'
oneDay_writeOutput_Patch <- function(){

  # Write Males
  if(length(private$adult_male)==0){
    writeLines(text = file.path(private$NetworkPointer$get_tNow(), private$patchID,
                                "NULL", "NULL", fsep = ","),
               con = private$NetworkPointer$get_conADM(),
               sep = "\n")

  } else{
    for(mosquito in private$adult_male){

      writeLines(text = file.path(private$NetworkPointer$get_tNow(), private$patchID,
                                  mosquito$print_male(), fsep = ","),
                 con = private$NetworkPointer$get_conADM(),
                 sep = "\n")
    }
  }

  # Write Females
  if(length(private$adult_female)==0){
    writeLines(text = file.path(private$NetworkPointer$get_tNow(), private$patchID,
                                "NULL", "NULL", "NULL", fsep = ","),
               con = private$NetworkPointer$get_conADF(),
               sep = "\n")
  } else{
    for(mosquito in private$adult_female){

      writeLines(text = file.path(private$NetworkPointer$get_tNow(), private$patchID,
                                  mosquito$print_female(), fsep = ","),
                 con = private$NetworkPointer$get_conADF(),
                 sep = "\n")
    }
  }

}
Patch$set(which = "public",name = "oneDay_writeOutput",
          value = oneDay_writeOutput_Patch, overwrite = TRUE
)

#######################################
# These are from patch-releases.R
#######################################

#' Release Mosquitoes in a Patch
#'
#' Based on this patch's release schedule, handle daily releases.
#'
oneDay_Releases_Patch <- function(){

  ################
  # MALE RELEASES
  ################
  if( (length(private$maleReleases) > 0) && (private$maleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()) ){

    # initialize holder list, then fill it with new mosquitoes
    private$newEggs <- vector(mode = "list", length = length(private$maleReleases[[1]]$genVec))
    for(i in 1:length(private$maleReleases[[1]]$genVec)){
      private$newEggs[[i]] <- Mosquito$new(genotype = private$maleReleases[[1]]$genVec[i],
                                           age = private$maleReleases[[1]]$ageVec[i])
    }

    #combine new mosquitoes with the general population, then remove this release from the plan
    private$adult_male = c(private$adult_male, private$newEggs)
    private$maleReleases[[1]] = NULL
  }

  ################
  # FEMALE RELEASES
  ################
  if( (length(private$femaleReleases) > 0) && (private$femaleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()) ){

    # initialize holder list, then fill it with new mosquitoes
    private$newEggs <- vector(mode = "list", length = length(private$femaleReleases[[1]]$genVec))
    for(i in 1:length(private$femaleReleases[[1]]$genVec)){
      private$newEggs[[i]] <- Mosquito$new(genotype = private$femaleReleases[[1]]$genVec[i],
                                           age = private$femaleReleases[[1]]$ageVec[i])
    }

    #combine new mosquitoes with the general population, then remove this release from the plan
    private$unmated_female = c(private$unmated_female, private$newEggs)
    private$femaleReleases[[1]] = NULL
  }

  ################
  # LARVAE RELEASES
  ################
  if( (length(private$larvaeReleases) > 0) && (private$larvaeReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()) ){

    # initialize holder list, then fill it with new mosquitoes
    private$newEggs <- vector(mode = "list", length = length(private$larvaeReleases[[1]]$genVec))
    for(i in 1:length(private$larvaeReleases[[1]]$genVec)){
      private$newEggs[[i]] <- Mosquito$new(genotype = private$larvaeReleases[[1]]$genVec[i],
                                           age = private$larvaeReleases[[1]]$ageVec[i])
    }

    #combine new mosquitoes with the general population, then remove this release from the plan
    private$larva = c(private$larva, private$newEggs)
    private$larvaeReleases[[1]] = NULL
  }

}
Patch$set(which = "public",name = "oneDay_Releases",
          value = oneDay_Releases_Patch, overwrite = TRUE
)

#######################################
# These are from patch-migration.R
#######################################

#' Oubound Migration
#'
#' Stochastic model of migration of adults from this patch.
#'
#' @details Migration is modeled as a Dirichlet-Multinomial process parameterized
#' by \code{moveVar} multiplied by the row corresponding to this patch from the
#' migration matrix. A Dirichlet distributed random variate is sampled over the
#' other patch frequency and then from \code{JaredDirichlet} according to that
#' parameter vector and then movement is sampled from \code{\link[stats]{rmultinom}}.
#'
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
    private$whoMigrate <- sample(x = 1:length(private$adult_male), size = private$numMigrate, replace = FALSE)

    #place how many migrate where
    private$numMigrateWhere <- cumsum(x = c(1, rmultinom(n = 1, size = private$numMigrate, prob = private$migrationDist)))

    private$genericCounter <- 1
    #loop over other patches, not your own
    for(patch in private$NetworkPointer$get_listPatch()[-private$patchID]){
      #if no mosquitoes go there, skip it
      if(private$numMigrateWhere[private$genericCounter] - private$numMigrateWhere[private$genericCounter+1]==0){next}

      #get all males who migrate
      private$maleMigration[[patch]] <- private$adult_male[private$whoMigrate[private$numMigrateWhere[private$genericCounter]:(private$numMigrateWhere[private$genericCounter+1]-1)]]
    }

    #remove all the mosquitoes who  migrated
    private$adult_male[private$whoMigrate] <- NULL

  }#end male migration

  #FEMALE
  if(length(private$adult_female)>0 && any(private$NetworkPointer$get_migrationFemale(private$patchID)[-private$patchID]!=0)){

    #number who migrate: population times migration fraction
    private$numMigrate <- as.integer(x = round(x = length(private$adult_female)*private$NetworkPointer$get_migrationFractionFemale(private$patchID)))

    #migration distribtuion, removing the current patch
    private$migrationDist <- JaredDirichlet(n = 1, alpha = private$NetworkPointer$get_migrationFemale(private$patchID)[-private$patchID]*private$NetworkPointer$get_moveVar())

    #generate a random sample of the population who will migrate
    private$whoMigrate <- sample(x = 1:length(private$adult_female), size = private$numMigrate, replace = FALSE)

    #place how many migrate where
    private$numMigrateWhere <- cumsum(x = c(1, rmultinom(n = 1, size = private$numMigrate, prob = private$migrationDist)))

    private$genericCounter <- 1
    #loop over other patches, not your own
    for(patch in private$NetworkPointer$get_listPatch()[-private$patchID]){
      #if no mosquitoes go there, skip it
      if(private$numMigrateWhere[private$genericCounter] - private$numMigrateWhere[private$genericCounter+1]==0){next}

      #get all males who migrate, then remove those males and that number from lists
      private$femaleMigration[[patch]] <- private$adult_female[private$whoMigrate[private$numMigrateWhere[private$genericCounter]:(private$numMigrateWhere[private$genericCounter+1]-1)]]

    }

    #remove all the mosquitoes who  migrated
    private$adult_female[private$whoMigrate] <- NULL

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







