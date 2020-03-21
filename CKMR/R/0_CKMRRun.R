###############################################################################
#     ________ __ __  _______ 
#    / ____/ //_//  |/  / __ \
#   / /   / ,<  / /|_/ / /_/ /
#  / /___/ /| |/ /  / / _, _/ 
#  \____/_/ |_/_/  /_/_/ |_|  
#    
###############################################################################
###############################################################################
# RUN WRAPPER
###############################################################################

#' Run mPlex
#' 
#' R interface to C++ simulation code.
#' 
#' @param seed An integer seed for random number generator. Default is 1
#' @param numReps An integer specifying the number of repetitions. Default is 1
#' @param numThreads An integer specifying the number of threads to parallelize over. Default is 1
#' @param networkParameters A list of simulation parameters
#' @param reproductionReference A list of reproduction and genotype specific parameters
#' @param patchReleases A list of releases
#' @param migrationMale A matrix specifying male migration rates
#' @param migrationFemale A matrix specifying female migration rates
#' @param migrationBatch A list specifing batch migration probabilities, rates, and sex ratios
#' @param samplingParameters A list specifying the time between sampling events and sampling coverage
#' @param outputDirectory String folder name to write output
#' @param verbose Boolean, be chatty? Default is FALSE
#' 
#' @export
runMPlex <- function(seed = 1, numReps = 1, numThreads = 1,
                     networkParameters, reproductionReference,
                     patchReleases, migrationMale, migrationFemale,
                     migrationBatch, samplingParameters,
                     outputDirectory, verbose = FALSE){
  
  # expand so c++ can find it
  outputDirectory = path.expand(outputDirectory)
  
  # pass all down for simulation
  run_mPlex(seed,
             numReps,
             numThreads,
             networkParameters,
             reproductionReference,
             patchReleases,
             migrationMale,
             migrationFemale,
             migrationBatch,
             samplingParameters,
             outputDirectory,
             verbose)
}
