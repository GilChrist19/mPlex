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
#' @usage MakeReference_Multiplex_oLocus(H, R, S, d, s_frac)
#'
#' @param H Vector of homing rates for each drive piece
#' @param R Vector of deleterious allele generation rates
#' @param S Vector of neutral allele generation rates
#' @param d Vector of background mutation rates at each locus
#' @param s_frac List of lists for genotype-dependent fertility reduction
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
#' s_frac = list(list("HH"=0, "HW"=0.5, "HR"=0, "HS"=0.5))
#'
#' MakeReference_Multiplex_oLocus(H,R,S,d,s_frac)
#'
#' @export
MakeReference_Multiplex_oLocus <- function(H=c(0.9),R=c(0.0), S=R/3, d=c(0.0001), s_frac=NULL){
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

  #set fractional reduction in fertility
  s <- vector(mode = "list", length = length(H))
  if(!is.null(s_frac)){s <- s_frac}

  return(list(
    mendelian = mendProbsList,
    homing = homProbsList,
    mendelianAlleles = mendAlleleList,
    homingAlleles = homAlleleList,
    s = s))

}

#' mPlex One Locus Offspring
#'
#' Create a list of offspring genotypes and their probabilities
#'
#' @usage MultiplexOffspring_oLocus(fGen, mGen, reference)
#'
#' @param fGen Female genotype
#' @param mGen Male genotype
#' @param reference Offspring reference list
#'
#' @details Using the reference generated by \code{\link{MakeReference_Multiplex_oLocus}},
#' this function expands the possible genotypes of the offspring and the ratios
#' that they occur. Similar to \code{\link{DaisyOffspring}} and
#' \code{\link{MultiplexOffspring_mLoci}}.
#'
#' @return List(Alleles, Probabilities)
#'
#' @examples
#' ref <- MakeReferenceDaisy(H = 0.98, R = 0.001, S = 0.0003, d = .00001)
#' fGen <- "WW"
#' mGen <- "WW"
#'
#' MultiplexOffspring_oLocus(fGen, mGen, ref)
#'
#' @export
MultiplexOffspring_oLocus <- function(fGen, mGen, reference){

  if(nchar(fGen) != nchar(mGen)){
    return(cat("Critters have different number of loci.\nCheck their genotypes."))
  }

  #split genotypes
  #This splits all characters.
  fSplit <- strsplit(x = fGen, split = "")[[1]]
  mSplit <- strsplit(x = mGen, split = "")[[1]]

  #get number of alleles. Divide by 2 is because diploid
  numAlleles <- length(fSplit)/2

  #the paste statement combines all targets in the first loci and all targets
  #  in the second loci into 2 complete alleles
  momAlleles <- list(allele1 = fSplit[1:numAlleles], allele2 = fSplit[(numAlleles+1):(2*numAlleles)] )
  dadAlleles <- list(allele1 = mSplit[1:numAlleles], allele2 = mSplit[(numAlleles+1):(2*numAlleles)] )

  #score them, assumes fGen and mGen are 1 character string, not a list of things
  #  If I stop pasting alleles together, we may have to change this!!
  fScore <- grepl(pattern = "H", x = fGen, ignore.case = FALSE, perl = FALSE, fixed = TRUE)
  mScore <- grepl(pattern = "H", x = mGen, ignore.case = FALSE, perl = FALSE, fixed = TRUE)

  #setup offspring allele lists
  # assuem diploid organism with numAlleles homing sites
  fAllele <- rep(x = list(vector(mode = "list", numAlleles)), 2)
  fProbs <- rep(x = list(vector(mode = "list", numAlleles)), 2)
  mAllele <- rep(x = list(vector(mode = "list", numAlleles)), 2)
  mProbs <- rep(x = list(vector(mode = "list", numAlleles)), 2)

  #Females!
  if(fScore) {
    #if homing allele present
    #loop over alleles, assume diploid
    for(i in 1:2){
      #loop over targets within locus
      for(j in 1:numAlleles){
        #Fill target with letter and probs
        if(momAlleles[[i]][[j]]=="W"){
          fAllele[[i]][[j]] <- reference$homingAlleles[[j]][[1]]
          fProbs[[i]][[j]] <- reference$homing[[j]]$W
        } else if(momAlleles[[i]][[j]]=="H"){
          fAllele[[i]][[j]] <- reference$homingAlleles[[j]][[2]]
          fProbs[[i]][[j]] <- reference$homing[[j]]$H
        } else if(momAlleles[[i]][[j]]=="R"){
          fAllele[[i]][[j]] <- reference$homingAlleles[[j]][[3]]
          fProbs[[i]][[j]] <- reference$homing[[j]]$R
        } else if(momAlleles[[i]][[j]]=="S"){
          fAllele[[i]][[j]] <- reference$homingAlleles[[j]][[4]]
          fProbs[[i]][[j]] <- reference$homing[[j]]$S
        }#end if string
      }#end target loop
    }#end allele loop

  } else {
    #if homing allele not present
    #loop over alleles, assume diploid
    for(i in 1:2){
      #loop over targets within locus
      for(j in 1:numAlleles){
        #Fill target with letter and probs
        if(momAlleles[[i]][[j]]=="W"){
          fAllele[[i]][[j]] <- reference$mendelianAlleles[[j]][[1]]
          fProbs[[i]][[j]] <- reference$mendelian[[j]]$W
        } else if(momAlleles[[i]][[j]]=="H"){
          fAllele[[i]][[j]] <- reference$mendelianAlleles[[j]][[2]]
          fProbs[[i]][[j]] <- reference$mendelian[[j]]$H
        } else if(momAlleles[[i]][[j]]=="R"){
          fAllele[[i]][[j]] <- reference$mendelianAlleles[[j]][[3]]
          fProbs[[i]][[j]] <- reference$mendelian[[j]]$R
        } else if(momAlleles[[i]][[j]]=="S"){
          fAllele[[i]][[j]] <- reference$mendelianAlleles[[j]][[4]]
          fProbs[[i]][[j]] <- reference$mendelian[[j]]$S
        }#end if string
      }#end target loop
    }#end allele loop
  }#end Female if statement

  #Males!
  if(mScore) {
    #if homing allele present
    #loop over alleles, assume diploid
    for(i in 1:2){
      #loop over targets within locus
      for(j in 1:numAlleles){
        #Fill target with letter and probs
        if(dadAlleles[[i]][[j]]=="W"){
          mAllele[[i]][[j]] <- reference$homingAlleles[[j]][[1]]
          mProbs[[i]][[j]] <- reference$homing[[j]]$W
        } else if(dadAlleles[[i]][[j]]=="H"){
          mAllele[[i]][[j]] <- reference$homingAlleles[[j]][[2]]
          mProbs[[i]][[j]] <- reference$homing[[j]]$H
        } else if(dadAlleles[[i]][[j]]=="R"){
          mAllele[[i]][[j]] <- reference$homingAlleles[[j]][[3]]
          mProbs[[i]][[j]] <- reference$homing[[j]]$R
        } else if(dadAlleles[[i]][[j]]=="S"){
          mAllele[[i]][[j]] <- reference$homingAlleles[[j]][[4]]
          mProbs[[i]][[j]] <- reference$homing[[j]]$S
        }#end if string
      }#end target loop
    }#end allele loop

  } else {
    #if homing allele not present
    #loop over alleles, assume diploid
    for(i in 1:2){
      #loop over targets within locus
      for(j in 1:numAlleles){
        #Fill target with letter and probs
        if(dadAlleles[[i]][[j]]=="W"){
          mAllele[[i]][[j]] <- reference$mendelianAlleles[[j]][[1]]
          mProbs[[i]][[j]] <- reference$mendelian[[j]]$W
        } else if(dadAlleles[[i]][[j]]=="H"){
          mAllele[[i]][[j]] <- reference$mendelianAlleles[[j]][[2]]
          mProbs[[i]][[j]] <- reference$mendelian[[j]]$H
        } else if(dadAlleles[[i]][[j]]=="R"){
          mAllele[[i]][[j]] <- reference$mendelianAlleles[[j]][[3]]
          mProbs[[i]][[j]] <- reference$mendelian[[j]]$R
        } else if(dadAlleles[[i]][[j]]=="S"){
          mAllele[[i]][[j]] <- reference$mendelianAlleles[[j]][[4]]
          mProbs[[i]][[j]] <- reference$mendelian[[j]]$S
        }#end if string
      }#end target loop
    }#end allele loop
  }#end Female if statement


  fAllLoci <- vector(mode = "list", length = 2)
  fProbsLoci <- vector(mode = "list", length = 2)
  mAllLoci <- vector(mode = "list", length = 2)
  mProbsLoci <- vector(mode = "list", length = 2)

  for(i in 1:2){
    #expand female alleles and probs
    holdAllOne <- expand.grid(fAllele[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    holdProbOne <- expand.grid(fProbs[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

    #expand male alleles and probs
    holdAllTwo <- expand.grid(mAllele[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    holdProbTwo <- expand.grid(mProbs[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

    #paste things in reverse
    fAllLoci[[i]] <- do.call(what = paste0, c(holdAllOne, collapse = NULL))
    fProbsLoci[[i]] <- apply(X = holdProbOne, MARGIN = 1, FUN = prod)

    mAllLoci[[i]] <- do.call(what = paste0, c(holdAllTwo, collapse = NULL))
    mProbsLoci[[i]] <- apply(X = holdProbTwo, MARGIN = 1, FUN = prod)

  }


  #unlist so that alleles within parent don't combine, then expand all
  #  combinations of male/female alleles
  outAList <- expand.grid(unlist(fAllLoci), unlist(mAllLoci), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  outPList <- expand.grid(unlist(fProbsLoci), unlist(mProbsLoci), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  outAList <- do.call(what = paste0, c(outAList, collapse = NULL))
  outPList <- apply(X = outPList, MARGIN = 1, FUN = prod)
  #could do hold[,1]*hold[,2]

  #aggregate and return
  #there is a vapply way to do this in one of the MGDrivE cubes.
  aggregateHold <- aggregate(outPList~outAList, data=data.frame(outAList, outPList), FUN=sum)


  #get proper type and normalize, then return as list.
  return(list(
    Alleles = as.character(aggregateHold$outAList),
    Probabilities = aggregateHold$outPList/sum(aggregateHold$outPList)
  ))

}
