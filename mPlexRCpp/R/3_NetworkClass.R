###############################################################################
#                                    ____  _
#                          _ __ ___ |  _ \| | _____  __
#                         | '_ ` _ \| |_) | |/ _ \ \/ /
#                         | | | | | |  __/| |  __/>  <
#                         |_| |_| |_|_|   |_|\___/_/\_\
#
###############################################################################
#    _   _      _                      _       ____ _
#   | \ | | ___| |___      _____  _ __| | __  / ___| | __ _ ___ ___
#   |  \| |/ _ \ __\ \ /\ / / _ \| '__| |/ / | |   | |/ _` / __/ __|
#   | |\  |  __/ |_ \ V  V / (_) | |  |   <  | |___| | (_| \__ \__ \
#   |_| \_|\___|\__| \_/\_/ \___/|_|  |_|\_\  \____|_|\__,_|___/___/
#
###############################################################################

#' Network Class Definition
#'
#' A Network is a collection of \code{\link{Patch}} with migration between them.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#'
#' @section **Constructor**:
#'  * networkParameters: see \code{\link{Network.Parameters}}
#'  * patchReleases: see \code{\link{Release_basicRepeatedReleases}}
#'  * reproductionType: one of c("DaisyDrive", "mPlex_oLocus", "mPlex_mLoci") specifying
#'  the type of reproduction to use. Options are: \code{\link{DaisyOffspring}},
#'  \code{\link{MultiplexOffspring_oLocus}}, or \code{\link{MultiplexOffspring_mLoci}}.
#'  * offspringReference: Inheritance patterns for the offspring. Must match reproductionType.
#'  See \code{\link{MakeReference_DaisyDrive}}, \code{\link{MakeReference_Multiplex_oLocus}},
#'  or \code{\link{MakeReference_Multiplex_mLoci}}.
#'  * migrationMale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationFemale: a stochastic matrix whose dimensions conform to the number of patches
#'  * directory: character string specifying the output directory
#'
#' @section **Methods**:
#'  * get_nPatch: Return the number of patches in the Network
#'  * get_simTime: Return total time for the simulation
#'  * get_parallel: Return boolean if running parallel networks
#'  * get_moveVar: Return numeric variance in Dirchlet-Multinomial movement
#'  * get_runID: Return the run identifier
#'  * get_stageTime: Return duration of life stage
#'  * get_beta: Return wild-type daily fertility
#'  * get_mu: Return stage-specific density-independent death rate
#'  * get_alpha: Return larval density-dependent parameter
#'  * get_listPatch: Integer list of Patches in the Network
#'  * get_patch: Return a list of all patches or a single Patch specified by PatchID
#'  * get_reference: Return reference list for calculating offspring distribution
#'  * get_s: Return genotype-dependent list of fractional fertilities
#'  * get_tNow: Return current simulation time
#'  * get_directory: Return output directory
#'  * get_conADM: Return adult_male file connection
#'  * get_conADF: Return adult_female file connection
#'  * close_connections: Close conADM and conADF
#'  * get_migrationMale: Return the migrationMale matrix or a single row specified by PatchID
#'  * get_migrationFemale: Return the migrationFemale matrix or a single row specified by PatchID
#'  * get_migrationFractionMale: Return fraction of adult_male population that migrates, specified by PatchID
#'  * get_migrationFractionFemale: Return fraction of adult_female population that migrates, specified by PatchID
#'  * get_patchReleases: Return release schedule, specified by PatchID and sex
#'
#' @section **Fields**:
#'  * nPatch: Number of patches
#'  * simTime: Maximum time of simulation
#'  * parallel: Boolean if running several patches as the same time
#'  * moveVar: Variation in Dirichlet sampling for migration
#'  * runID: Run identifier
#'  * stageTime: Integer vector specifying the time spent in each life stage
#'  * beta: Integer wild-type daily fertility
#'  * mu: Integer vector specifying density-independent death rate for each life stage
#'  * alpha: Numeric vector specifying larval density factor for each patch
#'  * listPatch: Integer vector listing all the patches, 1:nPatch
#'  * patches: A list of \code{\link{Patch}} objects
#'  * reference: Inheritance pattern list
#'  * s: Fractional reduction/increase in genotype-dependent fertility
#'  * tNow: Current simulation time
#'  * directory: a character string of where to store output
#'  * conADM: a \code{\link[base]{connection}} to write male population dynamics out to
#'  * conAF1: a \code{\link[base]{connection}} to write female population dynamics out to
#'  * migrationMale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationFractionMale: Numeric vector specifing the fraction of the population that migrates
#'  * migrationFemale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationFractionFemale: Numeric vector specifing the fraction of the population that migrates
#'  * patchReleases: a list of release schedules for each patch
#'  * maleMoveOut: List of outgoing males per patch, exists for memory management
#'  * femaleMofeOut: List of outgoing females per patch, exists for memory management
#'  * maleMoveIn: List of incoming males per patch, exists for memory management
#'  * femaleMoveIn: List of incoming females per patch, exists for memory management
#'  * holdMale: Holder for use getting maleMoveIn. Prevents list issues with setting "NULL"
#'  * holdFemale: Holder for use getting femaleMoveIn. Prevents list issues with setting "NULL"
#'
#' @export
Network <- R6::R6Class(classname = "Network",
            portable = TRUE,
            cloneable = FALSE,
            lock_class = FALSE,
            lock_objects = FALSE,
            class = FALSE,

            # public memebers
            public = list(

                # Constructor
                initialize = function(networkParameters, patchReleases, reproductionType=NULL, offspringReference, migrationMale, migrationFemale, directory){

                  #set basic things from Parameters
                  private$nPatch = networkParameters$nPatch
                  private$simTime = networkParameters$simTime
                  private$parallel = networkParameters$parallel
                  private$moveVar = networkParameters$moveVar
                  private$runID = networkParameters$runID
                  private$stageTime = networkParameters$stageTime
                  private$beta = networkParameters$beta
                  private$mu = networkParameters$mu
                  private$alpha = networkParameters$alpha

                  private$listPatch = 1:networkParameters$nPatch
                  private$patches = vector(mode="list",length=private$nPatch)

                  private$directory = directory

                  #Check releases and set them
                  if(length(patchReleases) != networkParameters$nPatch){
                    stop("length of patchReleases must equal number of patches in networkParameters!")
                  }
                  private$patchReleases = patchReleases

                  #set reproduction type
                  private$reference = offspringReference[1:(length(offspringReference)-1)]
                  private$s = offspringReference$s
                  switch(EXPR = reproductionType,
                         DaisyDrive = {Patch$set(which = "public",
                                                 name = "offspringDistribution",
                                                 value = DaisyOffspring_C,
                                                 overwrite = TRUE)
                           cat("Using DaisyOffspring() for reproduction.\n\n")},
                         mPlex_oLocus = {Patch$set(which = "public",
                                                   name = "offspringDistribution",
                                                   value = MultiplexOffspring_oLocus,
                                                   overwrite = TRUE)
                           cat("Using MultiplexOffspring_oLocus() for reproduction.\n\n")},
                         mPlex_mLoci = {Patch$set(which = "public",
                                                  name = "offspringDistribution",
                                                  value = MultiplexOffspring_mLoci_C,
                                                  overwrite = TRUE)
                           cat("Using MultiplexOffspring_mLoci() for reproduction.\n\n")},
                         stop("\n\nreproductionType must match one of these choices:\n",
                              " DaisyDrive\n mPlex_oLocus\n mPlex_mLoci\n")
                         )

                  #set migration matrices
                  private$migrationMale = migrationMale
                  private$migrationFemale = migrationFemale

                  #set fraction of population that moves. This doesn't change,
                  # so filling at initialization
                  private$migrationFractionMale = numeric(length = networkParameters$nPatch)
                  private$migrationFractionFemale = numeric(length = networkParameters$nPatch)
                  for(i in 1:networkParameters$nPatch){
                    private$migrationFractionMale[i] = sum(migrationMale[i,-i,drop=FALSE])
                    private$migrationFractionFemale[i] = sum(migrationFemale[i,-i,drop=FALSE])
                  }

                  # age ranges became a pain
                  lenAgeEgg <- length(0:networkParameters$stageTime["E"])

                  minAgeLarva <- networkParameters$stageTime["E"] + 1L
                  maxAgeLarva <- minAgeLarva + networkParameters$stageTime["L"]
                  distAgeLarva <- (1-networkParameters$mu["L"])^(minAgeLarva:maxAgeLarva)

                  minAgePupa <- maxAgeLarva + 1L
                  maxAgePupa <- minAgePupa + networkParameters$stageTime["P"]
                  lenAgePupa <- length(minAgePupa:maxAgePupa)

                  minAgeAdult <- maxAgePupa + 1L
                  maxAgeAdult <- minAgeAdult + networkParameters$stageTime["A"]
                  lenAgeAdult <- length(minAgeAdult:maxAgeAdult)

                  # initialize patches
                  for(i in 1:private$nPatch){

                    # initialize patch i
                    cat("initializing patch: ",i," of ",private$nPatch,"\n")

                    private$patches[[i]] = Patch$new(patchID = i,
                                                     simTime = private$simTime,
                                                     eggs_t0 = CreateMosquitoes_Distribution_Genotype(numMos = networkParameters$Leq[i],
                                                                                                      minAge = 0L, maxAge = networkParameters$stageTime["E"],
                                                                                                      ageDist =  rep(x = 1, times = lenAgeEgg)/lenAgeEgg,
                                                                                                      aTypes = networkParameters$alleloTypes[[i]]),


                                                     #approximation to number and age distribution of larva.
                                                     larva_t0 = CreateMosquitoes_Distribution_Genotype(numMos = networkParameters$Leq[i],
                                                                                                       minAge = minAgeLarva,
                                                                                                       maxAge = maxAgeLarva,
                                                                                                       ageDist = distAgeLarva/sum(distAgeLarva),
                                                                                                       aTypes = networkParameters$alleloTypes[[i]]),

                                                     #using AdPopEQ as empirically better
                                                     pupa_t0 = list(),#CreateMosquitoes_Distribution_Genotype(numMos = networkParameters$AdPopEQ[i]/2L,
                                                                #                                      minAge = minAgePupa,
                                                                #                                      maxAge = maxAgePupa,
                                                                #                                      ageDist = rep(x = 1, times = lenAgePupa)/lenAgePupa,
                                                                #                                      aTypes = networkParameters$alleloTypes[[i]]),

                                                     adult_male_t0 = CreateMosquitoes_Distribution_Genotype(numMos = networkParameters$AdPopEQ[i]/2L,
                                                                                                            minAge = minAgeAdult,
                                                                                                            maxAge = maxAgeAdult,
                                                                                                            ageDist = rep(x = 1, times = lenAgeAdult)/lenAgeAdult,
                                                                                                            aTypes = networkParameters$alleloTypes[[i]]),

                                                     unmated_female_t0 = CreateMosquitoes_Distribution_Genotype(numMos = networkParameters$AdPopEQ[i]/2L,
                                                                                                              minAge = minAgeAdult,
                                                                                                              maxAge = maxAgeAdult,
                                                                                                              ageDist = rep(x = 1, times = lenAgeAdult)/lenAgeAdult,
                                                                                                              aTypes = networkParameters$alleloTypes[[i]]),

                                                     maleReleases = patchReleases[[i]]$maleReleases,
                                                     femaleReleases = patchReleases[[i]]$femaleReleases,
                                                     larvaeReleases = patchReleases[[i]]$larvaeReleases
                                                   )

                    # set pointers
                    private$patches[[i]]$set_NetworkPointer(self)
                  }

                  # Output
                  if(!dir.exists(directory)){
                    dir.create(directory)
                  } else {
                    # if running in serial remove files, else do nothing
                    if(!networkParameters$parallel){
                      dirFiles = list.files(path = directory)
                      if(length(dirFiles)>0){
                        selection = utils::menu(c(paste0("delete all (",length(dirFiles),") files found in ",directory), "do not delete any files (files may be overwritten)"))
                        if(selection==1){
                          for(i in dirFiles){
                            cat("removing file: ",file.path(directory, i),"\n", sep = "")
                            file.remove(file.path(directory, i))
                          }
                        }

                      }
                    }

                  }

                }, # end constructor

                #getters and setters
                get_nPatch = function(){private$nPatch},
                get_simTime = function(){private$simTime},
                get_parallel = function(){private$parallel},
                get_moveVar = function(){private$moveVar},
                get_runID = function(){private$runID},
                get_stageTime = function(stage = NULL){private$stageTime[stage]},
                get_beta = function(){private$beta},
                get_mu = function(stage=NULL){private$mu[stage]},
                get_alpha = function(patch=NULL){private$alpha[patch]},
                get_listPatch = function(){private$listPatch},
                get_patch = function(patch = NULL){
                              if(is.null(patch)){private$patches
                              } else {private$patches[[patch]]}
                            },
                get_reference = function(){private$reference},
                get_s = function(){private$s},
                get_tNow  = function(){private$tNow},

                get_directory = function(){private$directory},
                get_conADM = function(){private$conADM},
                get_conADF = function(){private$conADF},
                close_connections = function(){
                                      close(private$conADM)
                                      close(private$conADF)
                                    },

                get_migrationMale = function(patch = NULL){
                                      if(is.null(patch)){private$migrationMale
                                      } else {private$migrationMale[patch, ,drop=FALSE]}
                                      #Does this drop thing matter here?
                                    },
                get_migrationFemale = function(patch = NULL){
                                        if(is.null(patch)){private$migrationFemale
                                        } else {private$migrationFemale[patch, ,drop=FALSE]}
                                      },
                get_migrationFractionMale = function(patch=NULL){private$migrationFractionMale[patch]},
                get_migrationFractionFemale = function(patch=NULL){private$migrationFractionFemale[patch]},

                get_patchReleases = function(patch, sex = "M"){
                                      switch(sex,
                                        M = {private$patchReleases[[patch]]$maleReleases},
                                        F = {private$patchReleases[[patch]]$femaleReleases},
                                        L = {private$patchReleases[[patch]]$larvaeReleases}
                                      )
                                    }
              ),

            # private members
            private = list(
                nPatch = NULL, # number of patches
                simTime = NULL, # max time of sim
                parallel = logical(1),
                moveVar = numeric(1),
                runID = integer(1),
                stageTime = numeric(4),
                beta = integer(1),
                mu = numeric(4),
                alpha = NULL,
                listPatch = NULL, #list of patch numbers
                patches = NULL, # list of patches
                reference = NULL, # reference rates for offspring
                s = NULL, #fractional reduction/increase in fertility
                tNow = 2L, # time starts at 2 because time 1 is the initial condition

                # output
                directory = NULL, # directory to store all patch output
                conADM = NULL,
                conADF = NULL,

                # inter-patch migration
                migrationMale = NULL,
                migrationFractionMale = NULL,
                migrationFemale = NULL,
                migrationFractionFemale = NULL,
                maleMoveOut = NULL, # for memory management
                femaleMoveOut = NULL, # for memory management
                maleMoveIn = NULL, # for memory management
                femaleMoveIn = NULL, # for memory management

                # release schedule
                patchReleases = NULL

              )
)

###############################################################################
# Functions
###############################################################################

#######################################
# These are from network-Simulation.R
#######################################

#' Reset Network
#'
#' Reset a \code{\link{Network}} between runs, useful for Monte Carlo simulation. This calls \code{\link{reset_Patch}} on each patch
#' and resets \code{tNow = 2} and increments the \code{runID}.
#'
reset_Network <- function(){

  cat("reset network\n",sep="")

  for(i in 1:private$nPatch){
    private$patches[[i]]$reset()
  }

  private$tNow = 2L
  private$runID = private$runID + 1L

}

Network$set(which = "public",name = "reset",
          value = reset_Network, overwrite = TRUE
)

#' Run Simulation
#'
#' Run a single simulation on this network.
#'
#' @param conADM an optional \code{\link[base]{connection}} to write male population dynamics to
#' @param conADF an optional \code{\link[base]{connection}} to write female population dynamics to
#'
oneRun_Network <- function(conADM = NULL, conADF = NULL){

  # open connections & write headers
  # parallel
if(private$parallel){

    pid = Sys.getpid()
    if(is.null(conADM)){
      private$conADM = file(description = paste0(private$directory, .Platform$file.sep, "ADM_pid_",pid,"_Run",private$runID,".csv"),open = "wt")
    } else {
      private$conADM = file(description = file.path(private$directory, conADM),open = "wt")
    }

    if(is.null(conADF)){
      private$conADF = file(description = paste0(private$directory, .Platform$file.sep, "ADF_pid_",pid,"_Run",private$runID,".csv"),open = "wt")
    } else {
      private$conADF = file(description = file.path(private$directory, conADF),open = "wt")
    }

  # serial
  } else {
    private$conADM = file(description = paste0(private$directory, .Platform$file.sep, "ADM_Run",private$runID,".csv"),open = "wt")
    private$conADF = file(description = paste0(private$directory, .Platform$file.sep, "ADF_Run",private$runID,".csv"),open = "wt")

  }

  # males
  writeLines(text = file.path("Time", "Patch", "Age", "Genotype", fsep = ","),
             con = private$conADM, sep = "\n")

  # females
  writeLines(text = file.path("Time", "Patch", "Age", "Genotype", "Mate", fsep = ","),
             con = private$conADF, sep = "\n")

  cat("begin run ",private$runID,"\n",sep="")

  # setup output
  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_initOutput()
  }

  # initialize progress bar
  pb = txtProgressBar(min = 0,max = private$simTime,style = 3)

  # run simulation
  while(private$simTime >= private$tNow){

    self$oneDay()

    private$tNow = private$tNow + 1L
    setTxtProgressBar(pb,value = private$tNow)
  }

  #close connections
  close(private$conADM)
  close(private$conADF)

  cat("run ",private$runID," over\n",sep="")
}

Network$set(which = "public",name = "oneRun",
          value = oneRun_Network, overwrite = TRUE
)

#' Run a Single Day on a Network
#'
#' Runs a single day of simulation on a \code{\link{Network}} object,
#' handling population dynamics, migration, population update, and output.
#'
oneDay_Network <- function(){

  # intra-patch population dynamics
  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_PopDynamics()
  }

  # inter-patch migration
  self$oneDay_Migration()

  # write output
  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_writeOutput()
  }

}
Network$set(which = "public",name = "oneDay",
          value = oneDay_Network, overwrite = TRUE
)

#######################################
# These are from network-migration.R
#######################################

#' Inter-Patch Migration
#'
#' Simulate migration between patches.
#'
oneDay_Migration_Network <- function(){

  ######################################
  # migration out
  ######################################

  for(i in private$listPatch){
    private$patches[[i]]$oneDay_migrationOut()
  }

  ######################################
  # migration in
  ######################################

  # grab moving mosquitoes
  private$maleMoveOut = vector(mode="list",length=private$nPatch)
  private$femaleMoveOut = vector(mode="list",length=private$nPatch)

  for(i in private$listPatch){

    private$maleMoveOut[[i]] = private$patches[[i]]$get_maleMigration()
    private$femaleMoveOut[[i]] = private$patches[[i]]$get_femaleMigration()

  }

  # sum over patches
  private$maleMoveIn <- lapply(X = private$listPatch,
                               FUN = function(x){unlist(lapply(X = private$maleMoveOut, FUN = '[[', x))})
  private$femaleMoveIn <- lapply(X = private$listPatch,
                                 FUN = function(x){unlist(lapply(X = private$femaleMoveOut, FUN = '[[', x))})

  ######################################
  # migration in
  ######################################

  for(i in private$listPatch){
    private$patches[[i]]$oneDay_migrationIn(maleIn = private$maleMoveIn[[i]],
                                             femaleIn = private$femaleMoveIn[[i]])
  }

}
Network$set(which = "public",name = "oneDay_Migration",
            value = oneDay_Migration_Network, overwrite = TRUE
)
