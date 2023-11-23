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

#' Run CKMR
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
#' @param samplingParameters A list specifying the time of sampling events and sampling coverage
#' @param outputDirectory String folder name to write output
#' @param verbose Boolean, be chatty? Default is FALSE
#' 
#' @export
runCKMR <- function(seed = 1, numReps = 1, numThreads = 1,
                     networkParameters, reproductionReference,
                     patchReleases, migrationMale, migrationFemale,
                     migrationBatch, samplingParameters,
                     outputDirectory, verbose = FALSE){
  
  # expand so c++ can find it
  outputDirectory = path.expand(outputDirectory)
  
  # numThreads safety check
  if(numThreads > 999){
    stop("Sim is setup to use less than 1000 cores. \n\tPlease stop being silly.")
  }

  ##########
  # Seed Setup
  ##########
  oldSeed <- .GlobalEnv$.Random.seed
  set.seed(seed)
  fourSeed <- sample(x = abs(.GlobalEnv$.Random.seed), size = 4, replace = FALSE)
  on.exit(.GlobalEnv$.Random.seed <- oldSeed)
  
  # pass all down for simulation
  run_CKMR(s1_ = fourSeed[1],
           s2_ = fourSeed[2],
           s3_ = fourSeed[3],
           s4_ = fourSeed[4],
           numReps,
           numThreads,
           networkParameters,
           reproductionReference,
           patchReleases,
           migrationMale,
           migrationFemale,
           migrationBatch,
           # split this so I can specify datatype - avoids casting issue from MGDivE
           samplingParameters$samplingDays,
           samplingParameters$samplingCoverage,
           outputDirectory,
           verbose)
}
