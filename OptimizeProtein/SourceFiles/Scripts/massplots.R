rm(list=ls(all=TRUE))

paths = Sys.glob(file.path("Summary*.txt"))

for(path in 1:length(paths)){
  print(paths[path])
  filename <- paths[path]
  stretches <- read.table(filename,header=T)
  name <- as.character(strsplit(filename,"\\.")[[1]][1])
  if(name == "SummaryStretches"){
    r <- 0
  }else{
    MASS <- read.delim(filename, header=T, dec=".")
    Mols <- unique(unlist(MASS$Mol))
    for(mol in Mols){
      dataset = subset(MASS,MASS$Mol==mol)
      pdf(paste("MASSplot",'_',mol,'.pdf',sep=""))
      plot(dataset$dTANGO,dataset$ddG,pch=19,xlab=expression(paste(Delta,"TANGO")),ylab=expression(paste(Delta,Delta,"G")),col="blue",ylim=c(min(dataset$ddG,na.rm=TRUE)-3,max(dataset$ddG,na.rm=TRUE)+3),xlim=c(min(dataset$dTANGO,na.rm=TRUE)-10,max(dataset$dTANGO,na.rm=TRUE)+10))
      title(paste("MASSplot\nMol: ",mol,sep=""))
      text(dataset$dTANGO,dataset$ddG, dataset$Mutation, cex=0.5, pos=4, col="red")
      dev.off()
    }
  }
}
