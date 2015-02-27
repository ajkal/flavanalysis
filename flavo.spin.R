source("scop.R")
library(filehash)
library(bio3d)
library(RCurl)
library(XML)

scop <- load.scop()
flavo.srch <- txt.search("flavoprotein",scop)
sf.ind <- grep("Flavoproteins",flavo.srch[,"des"])
flavo.hits <- list.scop.node(flavo.srch[sf.ind,1],scop)
##scopdbfiles <- get.scop(flavo.hits[,2],path="/Users/ajkal/Work/flavoproteins/scopdbs",full.names=TRUE)

pdbfiles <- get.pdb(unique(flavo.hits[,1]),path="/Users/ajkal/Work/flavoproteins/pdbs")
sepchains <- pdbsplit(pdbfiles,path="/Users/ajkal/Work/flavoproteins/split_chain")

##pdb.anns <- pdb.annotate(unique(flavo.hits[,"pdbid"]))
pdb.anns <- pdb.annotate(substr(basename(sepchains),1,4))

eafe