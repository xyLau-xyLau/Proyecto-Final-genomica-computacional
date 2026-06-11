library(adegenet)
library(ape)
library(phangorn)
secuencias <- fasta2DNAbin("secuencias/concatenados.fasta")
phydata <- as.phyDat(secuencias)
treeRatchet <- pratchet(phydata, trace = 0)
treeRatchet <- acctran(treeRatchet, phydata)
treeRatchet <- di2multi(treeRatchet)

if (inherits(treeRatchet, "multyPhylo")){
  treeRatchet <- unique(treeRatchet)
}
treeRatchet <-root(treeRatchet, outgroup= "Gallus_gallus")
plotBS(treeRatchet)
