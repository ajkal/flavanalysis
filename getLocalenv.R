getLocalenv <- function (fitalnall, ligname) {
  
  atypes <- c("N1","C2","O2","N3","C4","O4","C4A","N5","C5A",
              "C6","C7","C7M","C8","C8M","C9","C9A","N10","C10")
  
  atypes2 <- c("nix","nam","car","amm","cax","azl","ooh","tyr","trp")
  
  localenv <- matrix(NA,nrow=length(fitalnall$id),ncol=length(atypes))
  
  lapply(1:length(fitalnall$id), function(i) {
    pdb <- read.pdb(fitalnall$id[i])
    print(basename(fitalnall$id[i]))
    fulg.ind <- grep(ligname,pdb$atom[,"resid"])
    if (length(fulg.ind) == 0) {
      print(paste("No ligand atoms found",sep=""))
      break()
    } else {
      ioa.ind <- fulg.ind[match(atypes,pdb$atom[fulg.ind,"elety"])]
      prot.sel <- atom.select(pdb,"protein")
      dmat <- dm.xyz(xyz=pdb$xyz)
      for (j in 1:length(ioa.ind)) {
        atmin <- which.min(dmat[prot.sel$atom,ioa.ind[j]])
        localenv[i,j] <<- pdb$atom[prot.sel$atom,"resid"][atmin]
      }
    }
  })
  return(localenv)
}


