###############################################################################
#    ____        _           ____       _
#   |  _ \  __ _(_)___ _   _|  _ \ _ __(_)_   _____
#   | | | |/ _` | / __| | | | | | | '__| \ \ / / _ \
#   | |_| | (_| | \__ \ |_| | |_| | |  | |\ V /  __/
#   |____/ \__,_|_|___/\__, |____/|_|  |_| \_/ \___|
#                      |___/
###############################################################################

#' Daisy Drive Offspring Reference
#'
#' Create a list specifying the offspring probability distribution.
#'
#' @usage MakeReference_DaisyDrive(H, R, S, d)
#'
#' @param H Vector of homing rates for each drive piece
#' @param R Vector of deleterious allele generation rates
#' @param S Vector of neutral allele generation rates
#' @param d Vector of background mutation rates at each locus
#'
#' @details This function creates a reference list for \code{\link{DaisyOffspring}}.
#' The number of drive elements is specified by the length of H. The final homing
#' rate doesn't matter, as the last piece of the drive has nothing to drive, but
#' it must be there for compatibility. R,S, and d must be the same length as H,
#' but will generally be the same number replicated that many times. It is assumed
#' that S is R/3, but this can be varied. This function is similar to \code{\link{MakeReference_Multiplex_mLoci}}
#' and \code{\link{MakeReference_Multiplex_oLocus}}
#'
#' @return list of homing, cutting, and mendelian genotypes and rates
#'
#' @examples
#' H <- c(0.9,0.4,0) # This drive has 3 pieces
#' R <- c(0.001, 0.002, 0.003)
#' S <- c(0.0003, 0.0006, 0.001)
#' d <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_DaisyDrive(H,R,S,d)
#'
#' @export
MakeReference_DaisyDrive <- function(H=c(0.9, 0.4, 0.7),R=c(0.0, 0.0, 0.0), S=R/3, d=c(0.0001, 0.0001, 0.0001)){

  #H is homing rate. The length of this vector determines the number of pieces
  # in the daisy drive. Each drive can have the same or different rates.

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
  if(any(H[-length(H)]>1, R>1, S>1, d>1)){
    #last H doesn't matter because it may/may not exist and isnt' used
    return(cat("All rates must be less than or equal to 1\n"))
  }
  if(any((d+c(0, H[-length(H)])) > 1)){
    #need the driving piece's haming rate plus the background mutation rate
    #  of the piece being driven into.
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
  Hshift <- c(0, H[-length(H)]) #because each pieces relies on the efficiency of the previous piece

  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(H), dimnames = list(gtype, NULL))
  cuttingProbs <- matrix(data = 0, nrow = 3, ncol = length(H), dimnames = list(gtype[-2], NULL))

  homingProbs[1, ] <- 1-d-Hshift #chance to stay W is 1-homing-background mutation
  homingProbs[2, ] <- Hshift*(1-R-S) #chance to become homing is H*1-H*R-H*S
  homingProbs[3, ] <- Hshift*R   #NHEJ caused resistance, detrimentalt allele
  homingProbs[4, ] <- d+Hshift*S #good resistant allele, from NHEJ and background mutation rate

  cuttingProbs[1, ] <- 1-d-Hshift*(R+S) #chance to stay W, 1-homingrate*mutations-background
  cuttingProbs[2, ] <- Hshift*R #chance to become R, cutting*NHEJ deleterious rate
  cuttingProbs[3, ] <- Hshift*S+d #become S, cutting*NHEJ neutral + background


  #set up lists to hold probabilities
  mendProbsList <- vector(mode = "list", length = length(H))
  cutProbsList <- vector(mode = "list", length = length(H))
  homProbsList <- vector(mode = "list", length = length(H))

  #fill the lists
  for(i in 1:length(H)){
    mendProbsList[[i]]$W <- setNames(object = c(1-d[i], d[i]), nm = c("W", "S"))
    mendProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
    mendProbsList[[i]]$R <- 1
    mendProbsList[[i]]$S <- 1

    cutProbsList[[i]]$W <- cuttingProbs[ ,i]
    cutProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
    cutProbsList[[i]]$R <- 1
    cutProbsList[[i]]$S <- 1

    homProbsList[[i]]$W <- homingProbs[ ,i]
    homProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
    homProbsList[[i]]$R <- 1
    homProbsList[[i]]$S <- 1
  }


  return(list(
    mendelian = mendProbsList,
    cutting = cutProbsList,
    homing = homProbsList))

}
