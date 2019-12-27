###############################################################################
#                            ____  __          ______          
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
#                                               /_/   /_/      
###############################################################################
###############################################################################
# ONE RUN
###############################################################################

#' Run one run of mPlex
#' 
#' R interface to C++ simulation code.
#' 
#' @param seed An integer seed for random number generator
#' @param numThreads An integer specifying the number of threads to parallelize over
#' @param networkParameters A list of simulation parameters
#' @param reproductionReference A list of reproduction and genotype specific parameters
#' @param patchReleases A list of releases
#' @param migrationMale A matrix specifying male migration rates
#' @param migrationFemale A matrix specifying female migration rates
#' @param migrationBatch A list specifing batch migration probabilities, rates, and sex ratios
#' @param outputDirectory String folder name to write output
#' @param reproductionType String specifying type of reproduction model
#' @param verbose Boolean, be chatty?
#' 
#' @export
mPlex_oneRun <- function(seed, numThreads, networkParameters, reproductionReference, patchReleases,
                         migrationMale, migrationFemale, migrationBatch,
                         outputDirectory, reproductionType, verbose){
  
  # expand so c++ can find it
  outputDirectory = path.expand(outputDirectory)
  
  # pass all down for simulation
  mPlexCpp:::run_mPlex_Cpp(seed,
                           numThreads_ = numThreads,
                           networkParameters,
                           reproductionReference,
                           patchReleases,
                           migrationMale,
                           migrationFemale,
                           migrationBatch,
                           outputDirectory,
                           reproductionType,
                           verbose)
}

###############################################################################
# REPETITIONS WRAPPER
###############################################################################

#' Run multiple runs of mPlex
#' 
#' R interface to C++ simulation code.
#' 
#' @param seed An integer seed for random number generator
#' @param numReps An integer specifying the number of repetitions
#' @param numThreads An integer specifying the number of threads to parallelize over
#' @param networkParameters A list of simulation parameters
#' @param reproductionReference A list of reproduction and genotype specific parameters
#' @param patchReleases A list of releases
#' @param migrationMale A matrix specifying male migration rates
#' @param migrationFemale A matrix specifying female migration rates
#' @param migrationBatch A list specifing batch migration probabilities, rates, and sex ratios
#' @param outputDirectory String folder name to write output
#' @param reproductionType String specifying type of reproduction model
#' @param verbose Boolean, be chatty?
#' 
#' @export
mPlex_runRepetitions <- function(seed, numReps, numThreads, networkParameters, reproductionReference,
                                 patchReleases, migrationMale, migrationFemale,
                                 migrationBatch, outputDirectory, reproductionType,
                                 verbose){
  
  # expand so c++ can find it
  outputDirectory = path.expand(outputDirectory)
  
  # pass all down for simulation
  mPlexCpp:::run_mPlex_Cpp_repetitions(seed,
                                       numReps,
                                       numThreads_ = numThreads,
                                       networkParameters,
                                       reproductionReference,
                                       patchReleases,
                                       migrationMale,
                                       migrationFemale,
                                       migrationBatch,
                                       outputDirectory,
                                       reproductionType,
                                       verbose)
}






###############################################################################
# PROFILING WRAPPER
###############################################################################
# need a profiling wrapper