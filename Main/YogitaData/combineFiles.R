###############################################################################
#                            ____  __          ______
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/
#                                               /_/   /_/
###############################################################################
###############################################################################
#                          __    _            ____        __       
#    _________  ____ ___  / /_  (_)___  ___  / __ \____ _/ /_____ _
#   / ___/ __ \/ __ `__ \/ __ \/ / __ \/ _ \/ / / / __ `/ __/ __ `/
#  / /__/ /_/ / / / / / / /_/ / / / / /  __/ /_/ / /_/ / /_/ /_/ / 
#  \___/\____/_/ /_/ /_/_.___/_/_/ /_/\___/_____/\__,_/\__/\__,_/  
#                                                                  
###############################################################################
# 20200220
#  MGDrivE outputs all information into 2 files, one for males and one for females.
#  However, as it is a class-based model, all mosquito classes (aka genotypes) are 
#  printed in one line per patch per time-point, making it easy to split the file into 
#  individual patch information. mPlex is individual based, so there is no simple 
#  way to split a file by patch. I solved this by opening a file handle for each 
#  patch, and writing them separately so no splitting is required. 
#  This script is written to take whichever output you need, combine them into 
#  one file, and sort it based on time. In theory, this script will work on all 
#  output from mPlex.
#  This will breakdown when data gets too large. IDK what "too large" is.
#
# 20200616
#  Turning this into a function, so it can be sourced from the sims script and 
#  run. May be useful if we need to automate repetitions.
#  
###############################################################################
#rm(list = ls());gc()

# set target directory
# mainDir <- "~/Desktop/OUTPUT/mPlex/CKMR/simDir"
# workIndicator <- 25
# fPattern <- list("M","F")

###############################################################################
###############################################################################
# aux function that's way better than anything I was writing.
cLines <- function(file, chunkSize=50e6, ...) {
  # 20200220
  # this is taken from R.utils package by Henrik Bengtsson
  # https://github.com/HenrikBengtsson/R.utils
  
  # Argument 'file':
  if (inherits(file, "connection")) {
    con <- file
  } else {
    file <- as.character(file)
    con <- gzfile(file, open="rb")
    on.exit(close(con))
  }
  
  LF <- as.raw(0x0a)
  CR <- as.raw(0x0d)
  SPC <- as.raw(32L)
  
  isLastCR <- isLastLF <- FALSE
  isEmpty <- TRUE
  nbrOfLines <- 0L
  while(TRUE) {
    bfr <- readBin(con=con, what=raw(), n=chunkSize)
    if (isLastCR) {
      # Don't count LF following a CR in previous chunk.
      if (bfr[1L] == LF)
        bfr[1L] <- SPC
    }
    
    n <- length(bfr)
    if (n == 0L)
      break
    
    isEmpty <- FALSE
    
    # Replace all CRLF:s to become LF:s
    idxsCR <- which(bfr == CR)
    nCR <- length(idxsCR)
    if (nCR > 0L) {
      idxsCRLF <- idxsCR[(bfr[idxsCR + 1L] == LF)]
      if (length(idxsCRLF) > 0L) {
        bfr <- bfr[-idxsCRLF]
        n <- length(bfr)
        idxsCRLF <- NULL; # Not needed anymore
        nCR <- length(which(bfr == CR))
      }
    }
    
    # Count all CR:s and LF:s
    nLF <- length(which(bfr == LF))
    nbrOfLines <- nbrOfLines + (nCR + nLF)
    
    if (n == 0L) {
      isLastCR <- isLastLF <- FALSE
    } else {
      # If last symbol is CR it might be followed by a LF in
      # the next chunk. If so, don't count that next LF.
      bfrN <- bfr[n]
      isLastCR <- (bfrN == CR)
      isLastLF <- (bfrN == LF)
    }
  } # while()
  
  # Count any last line without newline too
  if (!isEmpty) {
    if (!isLastLF) nbrOfLines <- nbrOfLines + 1L
    #attr(nbrOfLines, "lastLineHasNewline") <- isLastLF
  }
  
  return(nbrOfLines)
}

###############################################################################
###############################################################################
# main function, does proper sizing, reads in files outputs one complete one
combineFiles <- function(mainDir, workIndicator=25, fPattern=c("M","F") ){
  
  # grab desired files from directory
  fList <- lapply(X = fPattern, FUN = function(x){
  	list.files(path = mainDir, pattern = paste0(x,"_"), full.names = TRUE)
  })
  
  # protect from incorrect fPattern
  fList <- fList[as.logical(lengths(fList))]
  
  fPatCount <- 1
  
  # loop over sets of files in folder
  for(stage in fList){
  	
  	##########
    # Get number of lines in files
    ##########
    # use cLines to count the number of lines in each file
    #  use sapply because not vectorized
    # the line count includes headers, so subtract one from all
  	nLines <- sapply(X = stage, FUN = cLines, USE.NAMES = FALSE) - 1
  	
  	##########
  	# Get header
  	##########
  	# This is because females are printed with a different header than everyone else. 
  	header <- scan(file = stage[1], what = character(), sep = ",", nlines = 1,quiet = TRUE)
  	nCol <- length(header)
  	
  	##########
  	# Return matrix
  	##########
  	# number of rows is total number of lines in all files
  	# number of columns is the number in the original files, plus the patch label
  	retMat <- matrix(data = 0L, nrow = sum(nLines), ncol = nCol+1,
  									 dimnames = list(NULL,c(header[1], "Patch", header[-1])))
  	
  	# counters or objects reused in loop
  	lStartCount <- 1
  	nFiles <- length(stage)
  	dataCols <- 3:(nCol+1)
  	
  	# loop over all files
  	for(myFile in 1:nFiles){
  		
  	  ##########
  	  # Input matrix
  	  ##########
  	  # Since all things are numbers, scan as integers
  	  #  we can calculate the number of items so setup memory size exactly
  	  # conver to matrix
  	  #  we also know the size of the matrix, so set that
  		# Some of the small populations were dying out as they were sampled too much,
  	  #  to prevent issues, skip all of this if there are no lines
  		if(nLines[myFile] == 0) next
  		dataMat <- matrix(data = scan(file = stage[myFile], what = integer(), sep = ",",
  																	skip = 1, nlines = nLines[myFile], n = nLines[myFile]*nCol,
  																	quiet = TRUE),
  											ncol = nCol, nrow = nLines[myFile], byrow = TRUE)
  
  		##########
  		# Patch number
  		##########
  		# Just to be safe, get patch number from the file name
  		# split on underscore (only in file name)
  		# unlist
  		# take last element
  		# split on period
  		# unlist
  		# take first element
  		# convert to 1-indexing
  		pNum <- as.integer(strsplit(x = tail(x = strsplit(x = stage[myFile],split = "_", fixed = TRUE)[[1]], n = 1),
  																split = ".", fixed = TRUE)[[1]][1]) + 1L
  		
  		##########
  		# Fill return matrix
  		##########
  		# setup row index
  		lineIndex <- lStartCount:(lStartCount + NROW(dataMat)-1L)
  		
  		# fill sampling time
  		retMat[lineIndex,1] <- dataMat[ ,1]
  		
  		# fill patch label
  		retMat[lineIndex,2] <- pNum
  		
  		# fill rest of data into retMat
  		retMat[lineIndex, dataCols] <- dataMat[ ,2:nCol]
  		
  		##########
  		# Final cleanup
  		##########
  		# increment line count
  		lStartCount <- lStartCount + NROW(dataMat)
  		
  		# indicate I'm working
  		if(!(myFile %% workIndicator)) print(myFile)
  		
  	} # end loop over files
  	
  	
  	##########
  	# Output
  	##########
  	# File name
  	fName <- file.path(mainDir,file.path("000_", fPattern[fPatCount],".csv",fsep = ""))
  	
  	# Sort and write
  	write.table(x = retMat[order(retMat[ ,1]), ], file = fName, sep = ",",
  							row.names = FALSE, col.names = TRUE)
  	
  	# update file pattern count
  	fPatCount <- fPatCount+1
  	
  } # end loop over file list

} # end function

