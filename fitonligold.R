fitonlig <- function(pdbfiles, ligname, pdb.anns) {
  atypes.fmn <- c("N1","C2","O2","N3","C4","O4","C4A","N5","C5A",
                  "C6","C7","C7M","C8","C8M","C9","C9A","N10","C10")
  atypes.fad <- c("N1","C2","O2","N3","C4","O4","C4X","N5","C5X",
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
  ref.sel <- atom.select(ref.pdb,resid=ref.lig,elety=atypes.fmn)
  ref.ord <- order(ref.pdb$atom[ref.sel$atom,"elety"])
  
  for (i in (1:length(subfiles))[-ref.ind]) {
    outfile <- paste("fitted/",
                     sapply(strsplit(basename(subfiles[i]),split=".pdb"),`[`,1),
                     "_flsq.pdb",
                     sep="")
    lig <- lig.list[[onematch.inds[i]]]
    if (length(lig) != 0) {
      print(paste("Found ligand: ",lig,sep=""))
    } else {
      print("No ligand found")
    }
    if (lig == "FAD") {
      atypes.lig <- atypes.fad
    } else {
      atypes.lig <- atypes.fmn
    }
    pdb <- read.pdb(subfiles[i])
    ## Just pick up one of the ligands if there are >1
    pdb.sel <- atom.select(pdb,resid=lig,elety=atypes.lig)
    if (length(pdb.sel$atom) != 18) {
      print("ERROR: LIGAND IS MISSING ATOMS")
      next
    }
    pdb.ord <- order(pdb$atom[pdb.sel$atom,"elety"])
    if (!identical(pdb.ord,ref.ord)) {
      print(paste("pdb file: ",basename(subfiles[i])," has atoms ordered differently than reference"))
    } else {
      print(paste("pdb file: ",basename(subfiles[i])," has atoms ordered the same as reference"))
    }
    xyz <- fit.xyz(fixed=ref.pdb$xyz,
                   mobile=pdb$xyz,
                   fixed.inds=ref.sel$xyz[atom2xyz(ref.ord)],
                   mobile.inds=pdb.sel$xyz[atom2xyz(pdb.ord)])
    write.pdb(pdb, xyz=xyz, file = outfile)
  }
}
  
  
  