


library(mPlexCpp)

numPatch <- 5
migration <- matrix(data = runif(numPatch*numPatch), nrow = numPatch, ncol = numPatch)
migration <- migration/rowSums(migration)
patchPops = rep(5L,numPatch)



reproductionReference <- MakeReference_DaisyDrive(H = c(0.98, 0.5),
                                                  R = c(0.0001, 0.0001),
                                                  S = c(0.0003, 0.004),
                                                  d = c(0, 0))



netPar = NetworkParameters(nPatch = numPatch,
                           simTime = 1000L,
                           alleloTypes = vector(mode = "list", length = numPatch),
                           AdPopEQ = patchPops,
                           runID = 1L,
                           dayGrowthRate = 1.1,
                           beta = 32L)

migrationBatch <- basicBatchMigration(numPatches = numPatch)


mPlex_oneRun(seed = 10,
             networkParameters = netPar,
             reproductionReference = reproductionReference,
             migrationMale = migration,
             migrationFemale = migration,
             migrationBatch = migrationBatch,
             reproductionType = "DaisyDrive",
             verbose = TRUE)
