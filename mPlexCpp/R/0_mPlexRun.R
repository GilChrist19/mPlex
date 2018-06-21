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
#' @param networkParameters A list of simulation parameters
#' @param reproductionReference A list of reproduction and genotype specific parameters
#' @param patchReleases A list of releases
#' @param migrationMale A matrix specifying male migration rates
#' @param migrationFemale A matrix specifying female migration rates
#' @param migrationBatch A list specifing batch migration probabilities, rates, and sex ratios
#' @param output_directory String folder name to write output
#' @param reproductionType String specifying type of reproduction model
#' @param verbose Boolean, be chatty?
#' 
#' @export
mPlex_oneRun <- function(seed, networkParameters, reproductionReference, patchReleases,
                         migrationMale, migrationFemale, migrationBatch,
                         output_directory, reproductionType, verbose){
  
  # expand so c++ can find it
  output_directory = path.expand(output_directory)
  
  # pass all down for simulation
  mPlexCpp:::run_mPlex_Cpp(seed,
                           networkParameters,
                           reproductionReference,
                           patchReleases,
                           migrationMale,
                           migrationFemale,
                           migrationBatch,
                           output_directory,
                           reproductionType,
                           verbose)
  
  
}

###############################################################################
# REPETITIONS WRAPPER
###############################################################################
# need a repetitions wrapper






###############################################################################
# PROFILING WRAPPER
###############################################################################
# need a profiling wrapper