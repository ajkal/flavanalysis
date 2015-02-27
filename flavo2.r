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
pdb.anns <- pdb.annotate(substr(bqasename(sepchains),1,4))

pdbs <- pdbaln(sepchains)
allatom <- read.all(pdbs)

fitted.files <- fitonlig(sepchains,ligname,pdb.anns)   
fitaln <- pdbaln(fitted.files)
fitalnall <- read.all(fitaln)


getLocalenv <- function (pdbfiles) {
  localenv.vec <- rep(NA,length(pdbfiles))
  mclapply(1:length(pdbfiles), function(i), {
    pdb <- read.pdb(pdbfiles[i])
    
    
    
    pdb.cmap <- cmap(pdb$xyz, grpby=pdb$atom[,"resno"],dcut=6)
    lig.cmapind <- which(pdb$atom[which(!duplicated(pdb$atom[,"resno"])),"resid"]=="FMN")
    prot.cmapind <- grep("(ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL)",
                         pdb$atom[which(!duplicated(pdb$atom[,"resno"])),"resid"],perl=TRUE)
    localenv.vec[i] <- sum(ref.cmap[prot.cmapind,170]==1,na.rm=TRUE)
  }
  return localenv.vec
}

length(subchains)


