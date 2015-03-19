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

pdbs <- pdbaln(sepchains)
allatom <- read.all(pdbs)

ligname <- "FMN|FAD|FNR|RBF"
fitted.files <- fitonlig(sepchains,ligname,pdb.anns)   
fitaln <- pdbaln(fitted.files)
fitalnall <- read.all(fitaln)

getLocalenv <- function (fitalnall) {
  localenv.vec <- rep(NA,length(fitalnall$id))
  mclapply(1:length(pdbfiles), function(i), {
    
    lig <- grep(ligname,unlist(strsplit(pdb.anns$ligandId[i],split=",")),value=TRUE)
    if (length(lig) != 0) {
      print(paste("Found ligand: ",lig,sep=""))
    } else {
      print("No ligand found")
    }
    pdb <- read.pdb(sepchains[i])
    ## Just pick up one of the ligands if there are >1   
    pdb.sel <- atom.select(pdb,resid=lig,elety=atypes)
    pdb.ord <- order(pdb$atom[pdb.sel$atom,"elety"])
    if (!identical(pdb.ord,ref.ord)) {
      print(paste("pdb file: ",basename(pdbfiles[i])," has atoms ordered differently than reference"))
      match.ind <- match(pdb.ord,ref.ord)
      if (fit==TRUE) {
        xyz <- fit.xyz(ref.pdb$xyz, pdb$xyz,
                       ref.sel$xyz, pdb.sel$xyz[atom2xyz(match.ind)])
        write.pdb(pdb, xyz=xyz, file = outfile)
      }
    } else {
      print(paste("pdb file: ",basename(pdbfiles[i])," has atoms ordered the same as reference"))
      if (fit==TRUE) {
        xyz <- fit.xyz(ref.pdb$xyz, pdb$xyz,
                       ref.sel$xyz, pdb.sel$xyz)
        write.pdb(pdb, xyz=xyz, file = outfile)
      } else if (lig.inds)
    }
    
    pdb.cmap <- cmap(pdb$xyz, grpby=pdb$atom[,"resno"],dcut=6)
    lig.cmapind <- which(pdb$atom[which(!duplicated(pdb$atom[,"resno"])),"resid"]=="FMN")
    prot.cmapind <- grep("(ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL)",
                         pdb$atom[which(!duplicated(pdb$atom[,"resno"])),"resid"],perl=TRUE)
    localenv.vec[i] <- sum(ref.cmap[prot.cmapind,170]==1,na.rm=TRUE)
  }
  return localenv.vec
}

length(subchains)

#ftco <- substr(basename(fitted.files),1,4)
#dm.repids <- flavo.hits[,"pdbid"][which(!duplicated(flavo.hits[,"dm"]))]
#dm.ftco <- na.omit(sapply(rep.dm.ids,match,ftco))
#dm.rep <- fitted.files[dm.ftco]
#dm.scop <- which(rownames(flavo.hits)
#                %in% names(dm.ftco))
#dm.txt <- flavo.hits[,"dm.txt"][dm.scop]
#pdb.anns$title[match(ftco,dm.ftco],
#                     tolower(pdb.anns$structureId))]
#unique(flavo.hits[,"fa.txt"][dmscop])

##substr(flavo.hits[,"pdbid"],1,4)[which(!duplicated(unlist(sapply(substr(flavo.hits[,"pdbid"],1,4),grep,unique(substr(basename(fitted.files),1,4))))))]