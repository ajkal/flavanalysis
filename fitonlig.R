fitonlig <- function(pdbfiles, ligname, fit, pdb.anns) {
  atypes <- c("N1","C2","O2","N3","C4","O4","C4A","N5","C5A",
              "C6","C7","C7M","C8","C8M","C9","C9A","N10","C10")
  ## Pick up ligands of interest in the annotations
  lig.list <- lapply(strsplit(pdb.anns$ligandId,split=","),
                     grep,pattern=ligname,perl=TRUE,value=TRUE)
  
  ## Only work with structures that have one ligand for now
  ligmatch <- unlist(lapply(lig.list,length))
  onematch.inds <- which(ligmatch == 1)
  subfiles <- pdbfiles[onematch.inds]
  
  ref.ind <- which.min(pdb.anns$resolution[onematch.inds])
  ref.lig <- lig.list[onematch.inds][[ref.ind]]
  ref.pdb <- read.pdb(subfiles[ref.ind])
  if (length(ref.lig)!=1) { stop("Reference structure has zero or more than 1 ligands...") }
  ref.sel <- atom.select(ref.pdb,resid=ref.lig,elety=atypes)
  ref.ord <- order(ref.pdb$atom[ref.sel$atom,"elety"])
  
  for (i in (1:length(subfiles))[-ref.ind]) {
    outfile <- paste("fitted/",
                     sapply(strsplit(basename(subfiles[i]),split=".pdb"),`[`,1),
                     "_flsq.pdb",
                     sep="")
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
  }
}
    
  
  
  