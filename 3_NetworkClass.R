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

Network <- R6::R6Class(classname = "Network",
            portable = TRUE,
            cloneable = FALSE,
            lock_class = FALSE,
            lock_objects = FALSE,

            # public memebers
            public = list(

                # Constructor
                initialize = function(networkParameters, patchReleases, migrationMale, migrationFemale, directory){

                  if(length(patchReleases) != networkParameters$nPatch){
                    stop("length of patchReleases must equal number of patches in networkParameters!")
                  }

                  private$parameters = networkParameters
                  private$nPatch = networkParameters$nPatch
                  private$listPatch = 1:networkParameters$nPatch
                  private$patches = vector(mode="list",length=private$nPatch)
                  private$simTime = networkParameters$simTime
                  private$directory = directory
                  private$runID = networkParameters$runID

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

                  private$patchReleases = patchReleases

                  # age ranges became a pain
                  lenAgeEgg <- length(0:networkParameters$timeAq["E"])

                  minAgeLarva <- networkParameters$timeAq["E"] + 1L
                  maxAgeLarva <- minAgeLarva + networkParameters$timeAq["L"]
                  lenAgeLarva <- length(minAgeLarva:maxAgeLarva)

                  minAgePupa <- maxAgeLarva + 1L
                  maxAgePupa <- minAgePupa + networkParameters$timeAq["P"]
                  lenAgePupa <- length(minAgePupa:maxAgePupa)

                  minAgeAdult <- maxAgePupa + 1L
                  maxAgeAdult <- minAgeAdult + networkParameters$tAdult
                  lenAgeAdult <- length(minAgeAdult:maxAgeAdult)

                  # initialize patches
                  for(i in 1:private$nPatch){

                    # initialize patch i
                    cat("initializing patch: ",i," of ",private$nPatch,"\n")

                    private$patches[[i]] = Patch$new(patchID = i,
                                                     simTime = private$simTime,
                                                     eggs_t0 = CreateMosquitoes_Distribution_Genotype(numMos = networkParameters$Leq[i],
                                                                                                      minAge = 0L, maxAge = networkParameters$timeAq["E"],
                                                                                                      ageDist =  rep(x = 1, times = lenAgeEgg)/lenAgeEgg,
                                                                                                      aTypes = networkParameters$alleloTypes[[i]]),

                                                     larva_t0 = CreateMosquitoes_Distribution_Genotype(numMos = networkParameters$Leq[i],
                                                                                                       minAge = minAgeLarva,
                                                                                                       maxAge = maxAgeLarva,
                                                                                                       ageDist = rep(x = 1, times = lenAgeLarva)/lenAgeLarva,
                                                                                                       aTypes = networkParameters$alleloTypes[[i]]),

                                                     pupa_t0 = CreateMosquitoes_Distribution_Genotype(numMos = networkParameters$Leq[i],
                                                                                                      minAge = minAgePupa,
                                                                                                      maxAge = maxAgePupa,
                                                                                                      ageDist = rep(x = 1, times = lenAgePupa)/lenAgePupa,
                                                                                                      aTypes = networkParameters$alleloTypes[[i]]),

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
                get_moveVar = function() {private$parameters$moveVar},
                get_timeAq = function(stage = NULL){
                                if(is.null(stage)){sum(private$parameters$timeAq)
                                } else {private$parameters$timeAq[[stage]]}
                              },
                get_beta = function(){private$parameters$beta},
                get_rm = function(){private$parameters$rm},
                get_AdPopEQ = function(ix){private$parameters$AdPopEQ[ix]},
                get_g = function(){private$parameters$g},
                get_Rm = function(){private$parameters$Rm},
                get_mu = function(stage){private$parameters$mu[[stage]]},
                get_alpha = function(ix){private$parameters$alpha[ix]},
                get_Leq = function(ix){private$parameters$Leq[ix]},
                get_patch = function(patch = NULl){
                              if(is.null(patch)){private$patches
                              } else {private$patches[[patch]]}
                            },
                get_nPatch = function(){private$nPatch},
                get_listPatch = function(){private$listPatch},
                get_directory = function(){private$directory},
                get_simTime = function(){private$simTime},
                get_conADM = function(){private$conADM},
                get_conAF1 = function(){private$conAF1},
                get_tNow  = function(){private$tNow},
                get_patchReleases = function(ix, sex = "M"){
                                      switch(sex,
                                        M = {private$patchReleases[[ix]]$maleReleases},
                                        F = {private$patchReleases[[ix]]$femaleReleases},
                                        L = {private$patchReleases[[ix]]$larvaeReleases}
                                      )
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
                close_connections = function(){
                                      close(private$conADM)
                                      close(private$conAF1)
                                    }
              ),

            # private members
            private = list(

                parameters = NULL,
                patches = NULL, # list of patches
                nPatch = NULL, # number of patches
                listPatch = NULL, #list of patch numbers
                simTime = NULL, # max time of sim
                tNow = 2L, # time starts at 2 because time 1 is the initial condition
                runID = numeric(1),

                # output
                directory = NULL, # directory to store all patch output
                conADM = NULL,
                conAF1 = NULL,

                # inter-patch migration
                migrationMale = NULL,
                migrationFractionMale = NULL,
                migrationFemale = NULL,
                migrationFractionFemale = NULL,

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
#' @param conAF1 an optional \code{\link[base]{connection}} to write female population dynamics to
#'
oneRun_Network <- function(conADM = NULL, conAF1 = NULL){

  # open connections & write headers
  # parallel
if(private$parameters$parallel){

    pid = Sys.getpid()
    if(is.null(conADM)){
      private$conADM = file(description = paste0(private$directory, .Platform$file.sep, "ADM_pid_",pid,"_Run",private$runID,".csv"),open = "wt")
    } else {
      private$conADM = file(description = file.path(private$directory, conADM),open = "wt")
    }

    if(is.null(conAF1)){
      private$conAF1 = file(description = paste0(private$directory, .Platform$file.sep, "AF1_pid_",pid,"_Run",private$runID,".csv"),open = "wt")
    } else {
      private$conAF1 = file(description = file.path(private$directory, conAF1),open = "wt")
    }

  # serial
  } else {
    private$conADM = file(description = paste0(private$directory, .Platform$file.sep, "ADM_Run",private$runID,".csv"),open = "wt")
    private$conAF1 = file(description = paste0(private$directory, .Platform$file.sep, "AF1_Run",private$runID,".csv"),open = "wt")

  }

  # males
  writeLines(text = file.path("Time", "Patch", "Age", "Genotype", fsep = ","),
             con = private$conADM, sep = "\n")

  # females
  writeLines(text = file.path("Time", "Patch", "Age", "Genotype", "Mate", fsep = ","),
             con = private$conAF1, sep = "\n")

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
  close(private$conAF1)

  cat("run ",private$runID," over\n",sep="")
}

Network$set(which = "public",name = "oneRun",
          value = oneRun_Network, overwrite = TRUE
)

#' Run a Single Day on a Network
#'
#' Runs a single day of simulation on a \code{\link{Network}} object, handling population dynamics, migration, population update, and output.
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

  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_migrationOut()
  }

  ######################################
  # migration in
  ######################################

  # grab moving mosquitoes
  maleMoveOut = vector(mode="list",length=private$nPatch)
  femaleMoveOut = vector(mode="list",length=private$nPatch)

  for(ix in 1:private$nPatch){

    maleMoveOut[[ix]] = private$patches[[ix]]$get_maleMigration()
    femaleMoveOut[[ix]] = private$patches[[ix]]$get_femaleMigration()

  }

  # sum over patches
  maleMoveIn = vector(mode="list",length=private$nPatch)
  femaleMoveIn = vector(mode="list",length=private$nPatch)

  for(first in 1:private$nPatch){
    #dummy variables that won't disappear when NULL is appended to them
    holdMale <- NULL
    holdFemale <- NULL
    for(second in 1:private$nPatch){
      holdMale <- c(holdMale, maleMoveOut[[first]][[second]])
      holdFemale <- c(holdFemale, femaleMoveOut[[first]][[second]])
    }

    #only if they aren't null, add them to the list. Otherwise, it destroys
    # elements
    if(!is.null(x = holdMale)){maleMoveIn[[first]] = holdMale}
    if(!is.null(x = holdFemale)){femaleMoveIn[[first]] = holdFemale}
  }

  ######################################
  # migration in
  ######################################

  for(ix in 1:private$nPatch){
    private$patches[[ix]]$oneDay_migrationIn(maleIn = maleMoveIn[[ix]],
                                             femaleIn = femaleMoveIn[[ix]])
  }

}

Network$set(which = "public",name = "oneDay_Migration",
            value = oneDay_Migration_Network, overwrite = TRUE
)







