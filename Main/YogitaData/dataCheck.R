################################################################################
# 20200507
# Jared Bennett
# Check mPlex
#  Gordana raised concerns about maternal half-sibs
#  Using this script to check basics of the data
#  Can update/expand as necessary
################################################################################
################################################################################

####################
# read female/male data
####################
females <- read.csv(file = "~/Desktop/OUTPUT/mPlex/newFemale.csv",
                    header = TRUE)
males <- read.csv(file = "~/Desktop/OUTPUT/mPlex/newMale.csv",
                  header = TRUE)
totData <- rbind(females,cbind(males,"Mate"=NA))


####################
# check number of sampled patches
####################
length(unique(females$Patch))
length(unique(males$Patch))


####################
# look at age range
####################
#  no max defined, since random death
#  min should be 8 or 9
range(females$Age)
range(males$Age)


####################
# make sure all individual IDs are unique
####################
#  ie, NROW() should be equal to unique(ID)
nrow(females) == length(unique(females$myID))
nrow(males) == length(unique(males$myID))

# check female mates
#  ensure that there are more unique females that unique mates
#  Males can mate more than once, so females could share the mate
nrow(females) >= length(unique(females$Mate))


####################
# check average captures per patch per day.
####################
#  expected ~7 females per patch per day
#  expected ~4.6 males per patch per day
simTime <- unique(females$Time)
patches <- sort(unique(females$Patch))

retArr <- array(data = 0, dim = c(length(simTime),2,length(patches)),
                dimnames = list(simTime,c("F","M"),patches))

# loop over time
for(sT in 1:length(simTime)){
  # subset by time
  #  speeds stuff up later
  fHold <- females$Patch[females$Time == simTime[sT]]
  mHold <- males$Patch[males$Time == simTime[sT]]

  # count how many in each patch.
  retArr[sT,"F", ] <- vapply(X = patches, FUN = function(x){sum(fHold == x)},
                             FUN.VALUE = numeric(length = 1))
  retArr[sT,"M", ] <- vapply(X = patches, FUN = function(x){sum(mHold == x)},
                             FUN.VALUE = numeric(length = 1))

}

# Metrics
#  average over time and patch
mean(retArr[ ,"F", ])
mean(retArr[ ,"M", ])

#  average over time, conditional on patch
colMeans(retArr[ ,"F", ])
colMeans(retArr[ ,"M", ])

#  average over patch, conditional on time
rowMeans(retArr[ ,"F", ])
rowMeans(retArr[ ,"M", ])


####################
# check M/O samples
####################
## The first part counts any instance of Mother/child as one
##  regardless of how many offspring were captured
#  mothers - female ID
#  offspring - momID
#    This goes for female and male data
momDaught <- sum(females$myID %in% females$momID)
momSon <- sum(females$myID %in% males$momID)

#  total M/O
sum(momDaught,momSon)

## The second part counts all Mother/child matches as independent
##  so, a mother with 2 offspring caught would count as 2, etc.
# get female ID who are mothers
mdID <- females$myID[females$myID %in% females$momID]
msID <- females$myID[females$myID %in% males$momID]

# tabulate all momID
dTab <- table(females$momID)
sTab <- table(males$momID)

# get counts of all offspring
mdCount <- sum(dTab[names(dTab) %in% mdID])
msCount <- sum(sTab[names(sTab) %in% msID])

# total
sum(mdCount, msCount)


####################
# check F/O samples
####################
## The first part counts any instance of father/child as one
##  regardless of how many offspring were captured
#  fathers - male ID
#  offspring - dadID
#    This goes for female and male data
dadDaught <- sum(males$myID %in% females$dadID)
dadSon <- sum(males$myID %in% males$dadID)

#  total F/O
sum(dadDaught,dadSon)

## The second part counts all father/child matches as independent
##  so, a father with 2 offspring caught would count as 2, etc.
# get male ID who are fathers
fdID <- males$myID[males$myID %in% females$dadID]
fsID <- males$myID[males$myID %in% males$dadID]

# tabulate all dadID
dTab <- table(females$dadID)
sTab <- table(males$dadID)

# get counts of all offspring
fdCount <- sum(dTab[names(dTab) %in% fdID])
fsCount <- sum(sTab[names(sTab) %in% fsID])

# total
sum(fdCount, fsCount)


####################
# check sibs (full, maternal half)
####################
# To do this, I'm going to need full sibs, so doing all 3
#  First, get momID
#  Then, get all momID where it is used more than once
#  Finally, check dadID
momCount <- table(totData$momID)
momDup <- momCount[momCount > 1] # all instances of siblings, labeled by mother

# to count number of sibling combinations, you need the number of siblings,
#  times the other half, and divide that by two (because A:B is the same as B:A,
#  so all combinations are double counted)
# This provides the formula n*(n-1)/n, where n is the number of siblings per mother
sum(momDup*(momDup-1)/2)


fullSibF <- function(x){
  # subset all sibs father's IDs
  sibs <- totData$dadID[totData$momID == x]
  # see if all daddy's are the same!
  return(all(sibs == sibs[1]))
}


# This takes a shitload of time (~40min on my laptop)
#A <- Sys.time()
hold <- vapply(X = names(momDup), FUN = fullSibF, FUN.VALUE = logical(length = 1))
#difftime(Sys.time(), A)

# all full sibs?
length(hold) == sum(hold)


####################
# check sibs (full, paternal half)
####################
# To do this, I'm going to need full sibs, so doing all 3
#  First, get dadID
#  Then, get all daID where it is used more than once
#  Finally, check momID
dadCount <- table(totData$dadID)
dadDup <- momCount[dadCount > 1] # all instances of siblings, labeled by mother

# to count number of sibling combinations, you need the number of siblings,
#  times the other half, and divide that by two (because A:B is the same as B:A,
#  so all combinations are double counted)
# This provides the formula n*(n-1)/n, where n is the number of siblings per mother
sum(dadDup*(dadDup-1)/2)

fullSibM <- function(x){
  # subset all sibs father's IDs
  sibs <- totData$momID[totData$dadID == x]
  # see if all daddy's are the same!
  return(all(sibs == sibs[1]))
}


# This takes a shitload of time (~40min on my laptop)
A <- Sys.time()
hold <- vapply(X = names(dadDup), FUN = fullSibM, FUN.VALUE = logical(length = 1))
difftime(Sys.time(), A)

# all full sibs?
length(hold) == sum(hold)


expecting 6935 half sibs
ie
length(momDub) - length(dadDup)

