###############################################################################
#    __  __       _ _   _       _                          _               _
#   |  \/  |_   _| | |_(_)_ __ | | _____  __     _ __ ___ | |    ___   ___(_)
#   | |\/| | | | | | __| | '_ \| |/ _ \ \/ /____| '_ ` _ \| |   / _ \ / __| |
#   | |  | | |_| | | |_| | |_) | |  __/>  <_____| | | | | | |__| (_) | (__| |
#   |_|  |_|\__,_|_|\__|_| .__/|_|\___/_/\_\    |_| |_| |_|_____\___/ \___|_|
#                        |_|
###############################################################################

#' mPlex Multiple Loci Offspring Reference
#'
#' Create a list specifying the offspring probability distribution.
#'
#' @usage MakeReference_Multiplex_mLoci(H, R, S, d, s_frac)
#'
#' @param H Vector of homing rates for each drive piece
#' @param R Vector of deleterious allele generation rates
#' @param S Vector of neutral allele generation rates
#' @param d Vector of background mutation rates at each locus
#' @param s_frac List of lists for genotype-dependent fertility reduction
#'
#' @details This function creates a reference list for \code{\link{MultiplexOffspring_mLoci}}.
#' Each drive element targets one locus. Each locus is independent of the rest.
#' The number of targets is specified by the length of H. R,S, and d must
#' be the same length as H, but will generally be the same number replicated that
#' many times. It is assumed that S is R/3, but this can be varied.
#' This function is similar to \code{\link{MakeReference_DaisyDrive}} and
#' \code{\link{MakeReference_Multiplex_oLocus}}
#'
#' @return list of homing and mendelian genotypes and rates
#'
#' @examples
#' H <- c(0.9,0.4,0.7) # This drive targets 3 loci
#' R <- c(0.001, 0.002, 0.003)
#' S <- c(0.0003, 0.0006, 0.001)
#' d <- c(0.0001, 0.0001, 0.0001)
#' s_frac = list(NULL, c("HH"=0, "HW"=0.5, "HR"=0, "HS"=0.5), NULL)
#'
#' MakeReference_Multiplex_mLoci(H,R,S,d,s_frac)
#'
#' @export
MakeReference_Multiplex_mLoci <- function(H=c(0.9, 0.4, 0.7),R=c(0.0, 0.0, 0.0), S=R/3, d=c(0.0001, 0.0001, 0.0001), s_frac=NULL){

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

#' mPlex Multiple Loci Offspring
#'
#' Create a list of offspring genotypes and their probabilities
#'
#' @usage MultiplexOffspring_mLoci(fGen, mGen, reference)
#'
#' @param fGen Female genotype
#' @param mGen Male genotype
#' @param reference Offspring reference list
#'
#' @details Using the reference generated by \code{\link{MakeReference_Multiplex_mLoci}},
#' this function expands the possible genotypes of the offspring and the ratios
#' that they occur. Similar to \code{\link{DaisyOffspring}} and
#' \code{\link{MultiplexOffspring_oLocus}}.
#'
#' @return List(Alleles, Probabilities)
#'
#' @examples
#' ref <- MakeReference_Multiplex_mLoci(H = 0.98, R = 0.001, S = 0.0003, d = .00001)
#' fGen <- "WW"
#' mGen <- "WW"
#'
#' MultiplexOffspring_mLoci(fGen, mGen, ref)
#'
#' @export
MultiplexOffspring_mLoci <- function(fGen, mGen, reference){

  if(nchar(fGen) != nchar(mGen)){
    return(cat("Critters have different number of loci.\nCheck their genotypes."))
  }

  #split mother genotype
  #This splits all characters.
  fSplit <- strsplit(x = fGen, split = "")[[1]]
  mSplit <- strsplit(x = mGen, split = "")[[1]]

  #get number of alleles. Divide by two because diploid
  numAlleles <- length(fSplit)/2

  #make a list of each allele at every locus. This list is length numAlleles, and each
  # sublist has length 2
  momAlleles <- lapply(X = seq.int(from = 1, to = 2*numAlleles, by = 2), FUN = function(X){fSplit[X:(X+1)]})
  dadAlleles <- lapply(X = seq.int(from = 1, to = 2*numAlleles, by = 2), FUN = function(X){mSplit[X:(X+1)]})

  #score them
  fscore <- grepl(pattern = "H", x = momAlleles, ignore.case = FALSE, perl = FALSE, fixed = TRUE)
  mscore <- grepl(pattern = "H", x = dadAlleles, ignore.case = FALSE, perl = FALSE, fixed = TRUE)

  #setup offspring allele lists
  fAllele <- rep(x = list(vector(mode = "list", 2)), numAlleles)
  fProbs <- rep(x = list(vector(mode = "list", 2)), numAlleles)
  mAllele <- rep(x = list(vector(mode = "list", 2)), numAlleles)
  mProbs <- rep(x = list(vector(mode = "list", 2)), numAlleles)

  #loop over all loci
  for(i in 1:numAlleles){

    #Females!
    if(fscore[i]){
      #If there is homing
      #loop over alleles at the locus. Everything is diploid here.
      for(j in 1:2){

        #Fill allele with letter and probs
        if(momAlleles[[i]][[j]]=="W"){
          fAllele[[i]][[j]] <- reference$homingAlleles[[i]][[1]]
          fProbs[[i]][[j]] <- reference$homing[[i]]$W
        } else if(momAlleles[[i]][[j]]=="H"){
          fAllele[[i]][[j]] <- reference$homingAlleles[[i]][[2]]
          fProbs[[i]][[j]] <- reference$homing[[i]]$H
        } else if(momAlleles[[i]][[j]]=="R"){
          fAllele[[i]][[j]] <- reference$homingAlleles[[i]][[3]]
          fProbs[[i]][[j]] <- reference$homing[[i]]$R
        } else if(momAlleles[[i]][[j]]=="S"){
          fAllele[[i]][[j]] <- reference$homingAlleles[[i]][[4]]
          fProbs[[i]][[j]] <- reference$homing[[i]]$S
        }

      }#end allele loop
    } else {
      #If there is not homing
      #loop over alleles at the locus. Everything is diploid here.
      for(j in 1:2){

        #Fill allele with letter and probs
        if(momAlleles[[i]][[j]]=="W"){
          fAllele[[i]][[j]] <- reference$mendelianAlleles[[i]][[1]]
          fProbs[[i]][[j]] <- reference$mendelian[[i]]$W
        } else if(momAlleles[[i]][[j]]=="H"){
          fAllele[[i]][[j]] <- reference$mendelianAlleles[[i]][[2]]
          fProbs[[i]][[j]] <- reference$mendelian[[i]]$H
        } else if(momAlleles[[i]][[j]]=="R"){
          fAllele[[i]][[j]] <- reference$mendelianAlleles[[i]][[3]]
          fProbs[[i]][[j]] <- reference$mendelian[[i]]$R
        } else if(momAlleles[[i]][[j]]=="S"){
          fAllele[[i]][[j]] <- reference$mendelianAlleles[[i]][[4]]
          fProbs[[i]][[j]] <- reference$mendelian[[i]]$S
        }

      }#end allele loop
    }#end female if statement

    #Males!
    if(mscore[i]){
      #If there is homing
      #loop over alleles at the locus. Everything is diploid here.
      for(j in 1:2){

        #Fill allele with letter and probs
        if(dadAlleles[[i]][[j]]=="W"){
          mAllele[[i]][[j]] <- reference$homingAlleles[[i]][[1]]
          mProbs[[i]][[j]] <- reference$homing[[i]]$W
        } else if(dadAlleles[[i]][[j]]=="H"){
          mAllele[[i]][[j]] <- reference$homingAlleles[[i]][[2]]
          mProbs[[i]][[j]] <- reference$homing[[i]]$H
        } else if(dadAlleles[[i]][[j]]=="R"){
          mAllele[[i]][[j]] <- reference$homingAlleles[[i]][[3]]
          mProbs[[i]][[j]] <- reference$homing[[i]]$R
        } else if(dadAlleles[[i]][[j]]=="S"){
          mAllele[[i]][[j]] <- reference$homingAlleles[[i]][[4]]
          mProbs[[i]][[j]] <- reference$homing[[i]]$S
        }

      }#end allele loop
    } else {
      #If there is not homing
      #loop over alleles at the locus. Everything is diploid here.
      for(j in 1:2){

        #Fill allele with letter and probs
        if(dadAlleles[[i]][[j]]=="W"){
          mAllele[[i]][[j]] <- reference$mendelianAlleles[[i]][[1]]
          mProbs[[i]][[j]] <- reference$mendelian[[i]]$W
        } else if(dadAlleles[[i]][[j]]=="H"){
          mAllele[[i]][[j]] <- reference$mendelianAlleles[[i]][[2]]
          mProbs[[i]][[j]] <- reference$mendelian[[i]]$H
        } else if(dadAlleles[[i]][[j]]=="R"){
          mAllele[[i]][[j]] <- reference$mendelianAlleles[[i]][[3]]
          mProbs[[i]][[j]] <- reference$mendelian[[i]]$R
        } else if(dadAlleles[[i]][[j]]=="S"){
          mAllele[[i]][[j]] <- reference$mendelianAlleles[[i]][[4]]
          mProbs[[i]][[j]] <- reference$mendelian[[i]]$S
        }

      }#end allele loop
    }#end male if statement
  }#end loci loop


  #combine each locus into single lists, so that alleles within a locus can't
  # be combined with each other, but do get combined with allelels for the
  # other loci.
  # ie, unlist the sublists in each allele/probs list. This will give single-
  # depth lists the same length as the number of multiplex loci.
  fAllLoci <- lapply(X = fAllele, FUN = unlist, recursive=TRUE)
  fProbsLoci <- lapply(X = fProbs, FUN = unlist, recursive=TRUE)
  mAllLoci <- lapply(X = mAllele, FUN = unlist, recursive=TRUE)
  mProbsLoci <- lapply(X = mProbs, FUN = unlist, recursive=TRUE)

  #combine male and female alleles at each locus.
  # This requires looping through each locus, getting all combinations of
  lociAList <- vector(mode = "list", length = numAlleles)
  lociPList <- vector(mode = "list", length = numAlleles)

  for( i in 1:numAlleles){

    #get all combinationes of male/female for each allele
    holdAllOne <- expand.grid(fAllLoci[[i]], mAllLoci[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    holdProbOne <- expand.grid(fProbsLoci[[i]], mProbsLoci[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

    #sort each combination so they are the same.
    holdAllOne <- apply(X = holdAllOne, MARGIN = 1, FUN = sort)

    #paste alleles togheter
    holdAllTwo <- do.call(what = "paste0", list(holdAllOne[1, ], holdAllOne[2, ]))
    holdProbTwo <- holdProbOne[ ,1]*holdProbOne[ ,2]

    #aggregate and return
    aggregateHold <- aggregate(holdProbTwo~holdAllTwo, data=data.frame(holdAllTwo, holdProbTwo), FUN=sum)

    #fill lists
    lociAList[[i]] <- as.character(aggregateHold$holdAllTwo)
    lociPList[[i]] <- aggregateHold$holdProbTwo

  }


  #get all combinations of each loci. This gives the total genotype
  outAList <- expand.grid(lociAList, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  outPList <- expand.grid(lociPList, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  #combine allele names and probabilities
  outAList <- apply(X = outAList, MARGIN = 1, FUN = paste0, collapse="")
  outPList <- apply(X = outPList, MARGIN = 1, FUN = prod)
  #can use matrixStats::rowProds(x = as.matrix(outPList))


  #normalize probs and return
  return(list(
    Alleles = outAList,
    Probabilities = outPList/sum(outPList)
  ))

}
