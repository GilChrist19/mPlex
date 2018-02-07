###############################################################################
#    __  __       _ _   _       _                      _
#   |  \/  |_   _| | |_(_)_ __ | | _____  __      ___ | |    ___   ___ _   _ ___
#   | |\/| | | | | | __| | '_ \| |/ _ \ \/ /____ / _ \| |   / _ \ / __| | | / __|
#   | |  | | |_| | | |_| | |_) | |  __/>  <_____| (_) | |__| (_) | (__| |_| \__ \
#   |_|  |_|\__,_|_|\__|_| .__/|_|\___/_/\_\     \___/|_____\___/ \___|\__,_|___/
#                        |_|
###############################################################################

#' mPlex One Locus Offspring Reference
#'
#' Create a list specifying the offspring probability distribution.
#'
#' @usage MakeReference_Multiplex_oLocus(H, R, S, d)
#'
#' @param H Vector of homing rates for each drive piece
#' @param R Vector of deleterious allele generation rates
#' @param S Vector of neutral allele generation rates
#' @param d Vector of background mutation rates at each locus
#'
#' @details This function creates a reference list for \code{\link{MultiplexOffspring_oLocus}}.
#' It assumes multiple targeting gRNAs for 1 locus, all targets segregate together.
#' The length of each allele is specified by the length of H. R, S, and d must
#' be the same length as H, but will generally be the same number replicated that
#' many times. It is assumed that S is R/3, but this can be varied.
#' This function is similar to \code{\link{MakeReference_DaisyDrive}} and
#' \code{\link{MakeReference_Multiplex_mLoci}}
#'
#' @return list of homing and mendelian genotypes and rates
#'
#' @examples
#' H <- c(0.9,0.4,0.7) # This drive has 3 targets at 1 locus
#' R <- c(0.001, 0.002, 0.003)
#' S <- c(0.0003, 0.0006, 0.001)
#' d <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_Multiplex_oLocus(H,R,S,d)
#'
#' @export
MakeReference_Multiplex_oLocus <- function(H=c(0.9),R=c(0.0), S=R/3, d=c(0.0001)){
  #H is homing rate. The length of this vector determines the number of loci
  # in the multiplex drive. Each drive can have the same or different rates.

  #R is the NHEJ rate for deleterious alleles. Must be same length as H,
  # can be the same or different values.

  #S is the NHEJ rate for neutral alleles. Must be the same length as H,
  # can be the same or different values.

  #d is the background mutation rate. Must be the same length as H, can
  # have the same or different values.



  #Safety checks
  if(any( c(length(H),length(R), length(S), length(d)) != length(H))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(H) == length(R) == length(S) == length(d)"))
  }
  if(any(H>1, R>1, S>1, d>1)){
    return(cat("All rates must be less than or equal to 1\n"))
  }
  if(any((d+H) > 1)){
    return(cat("Homing rates plus background mutation rates must be <= 1\n",
               "i.e. H+d <= 1\n"))
  }
  if(any((R+S) > 1)){
    return(cat("Negative and neutral repair rates must sum to <= 1\n",
               "i.s. R+S <= 1"))
  }


  #setup allele letters
  #W = Wild-type
  #H = Homing
  #R = Deleterious resistant
  #S = Neutral resistant
  gtype <- c("W", "H", "R", "S")

  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(H), dimnames = list(gtype, NULL))

  homingProbs[1, ] <- 1-d-H #chance to stay W is 1-homing-background mutation
  homingProbs[2, ] <- H*(1-R-S) #chance to become homing is H*1-H*R-H*S
  homingProbs[3, ] <- H*R   #NHEJ caused resistance, detrimentalt allele
  homingProbs[4, ] <- d+H*S #good resistant allele, from NHEJ and background mutation rate

  #set up lists to hold probabilities
  mendProbsList <- vector(mode = "list", length = length(H))
  homProbsList <- vector(mode = "list", length = length(H))

  mendAlleleList <- vector(mode = "list", length = length(H))
  homAlleleList <- vector(mode = "list", length = length(H))

  #fill the lists
  for(i in 1:length(H)){
    #Mendelian Probabilities
    mendProbsList[[i]]$W <- setNames(object = c(1-d[i], d[i]), nm = c("W", "S"))
    mendProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
    mendProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    mendProbsList[[i]]$S <- setNames(object = 1, nm = "S")

    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = mendProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      mendProbsList[[i]][[j]] <- mendProbsList[[i]][[j]][ logicalHold[[j]] ]
      mendAlleleList[[i]][[j]] <- names(mendProbsList[[i]][[j]])
    }

    #Homing Probabilities
    homProbsList[[i]]$W <- homingProbs[ ,i]
    homProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
    homProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    homProbsList[[i]]$S <- setNames(object = 1, nm = "S")

    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = homProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      homProbsList[[i]][[j]] <- homProbsList[[i]][[j]][ logicalHold[[j]] ]
      homAlleleList[[i]][[j]] <- names(homProbsList[[i]][[j]])
    }
  }

  return(list(
    mendelian = mendProbsList,
    homing = homProbsList,
    mendelianAlleles = mendAlleleList,
    homingAlleles = homAlleleList))

}
