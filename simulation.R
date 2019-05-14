library(splatter)

simulation <- function(batchCells, nGenes, nGroups) {
  params <- newSplatParams()
  
  # Set seed to 0
  params <- setParam(params, "seed", 0)
  
  # Set cell=200, genes=5000
  params <- setParams(params, update = list(nGenes = nGenes, batchCells = batchCells))
  
  # Generate random probabilities for each cell type
  probs <- sample(1:10, nGroups)
  probs <- probs/sum(probs)
  
  # Start simulations
  #sim <- splatSimulate(params)
  sim <- splatSimulate(params,
                       group.prob=probs,
                       method = "groups",
                       verbose = TRUE)
  
  groupdist <- as.vector(table(sim@colData@listData$Group))
  truecounts <- as.data.frame(t(counts(sim)))
  cellinfo <- as.data.frame(colData(sim))
  geneinfo <- as.data.frame(rowData(sim))
  
  # Add dropouts
  params <- setParam(params, "dropout.type", "experiment")
  
  sim.drop1 <- splatter:::splatSimDropout(sim, setParam(params, "dropout.mid", 0))
  sim.drop2 <- splatter:::splatSimDropout(sim, setParam(params, "dropout.mid", 2))
  sim.drop3 <- splatter:::splatSimDropout(sim, setParam(params, "dropout.mid", 5))
  
  counts_drop1 <- as.data.frame(t(counts(sim.drop1)))
  counts_drop2 <- as.data.frame(t(counts(sim.drop2)))
  counts_drop3 <- as.data.frame(t(counts(sim.drop3)))
  
  dropout1 <- sum(assays(sim.drop1)$Dropout)/(batchCells*nGenes)
  dropout2 <- sum(assays(sim.drop2)$Dropout)/(batchCells*nGenes)
  dropout3 <- sum(assays(sim.drop3)$Dropout)/(batchCells*nGenes)
  
  # Return
  return(list(counts_drop1 = counts_drop1,
              counts_drop2 = counts_drop2,
              counts_drop3 = counts_drop3,
              dropout1 = dropout1,
              dropout2 = dropout2,
              dropout3 = dropout3,
              truecounts = truecounts,
              cellinfo=cellinfo,
              geneinfo=geneinfo,
              groupdist=groupdist))
}
