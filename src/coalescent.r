setwd("~/Downloads/edu/fciencias/genomica_20262/proyecto")
library(phangorn)
dna <-c("secuencias/nd2_completos.fasta",
        "coi5p_alineados.fasta")
i <- 0
mtree <- c()
for (k in dna) {
  i <- i + 1
  dna <- fasta2DNAbin(k)
  phydata <- as.phyDat(dna)
  modeltest <- modelTest(phydata)
  pml <- pml_bb(modeltest)
  mtree[[i]] <- root(pml$tree, outgroup= "Gallus_gallus")
}
class(mtree) <- "multiPhylo"
supertree <- superTree(mtree)
supertree <- root(supertree, outgroup = "Gallus_gallus")
plot(supertree)