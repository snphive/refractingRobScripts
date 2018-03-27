rm(list=ls(all=TRUE))
install.packages("plyr")
library(plyr)

Seqdet <- read.delim("SummaryPerResSeqDet.txt",header=T)
Agad <- read.delim("SummaryPerResAgad.txt",header=T)
Seqdet$Res1 <- plyr::mapvalues(Seqdet$amino.acid,from=c("ALA","CYS","ASP","GLU","PHE","GLY","HIS","H2S","H1S","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"),to=c("A","C","D","E","F","G","H","H","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"))
Seqdet$Resnum <- Seqdet$number
Seqdet$Mol <- Seqdet$chain

Summary <- merge(Seqdet,Agad,by=c("Mol","Resnum","Res1"))
APRs <- subset(Summary,APRcount!=0)
APRsums <- aggregate(. ~ Mol + Pdb + Name + APRcount + chain, data=APRs,sum)

APRsums$dG_corrected <- APRsums$total.energy
APRsums$dG_corrected[APRsums$dG_corrected < -5] <- -5
APRsums$dG_corrected[APRsums$dG_corrected > 5] <- 5
APRsums$dG_norm <- ((APRsums$dG_corrected+5)/10)
APRsums$SolubisScore <- APRsums$dG_norm*APRsums$TANGO

SolubisScoreSum<- aggregate(SolubisScore ~ Pdb,data=APRsums,sum)
TangoScoreSum <- aggregate(TANGO ~ Pdb,data=APRsums,sum)

colors <- c("blue","red","black","green3","cyan","orange","magenta","grey","purple")
mols <- unique(unlist(APRsums$Mol))

pdf(file="Stretchplot.pdf",width=10,height=8,onefile=TRUE,bg='white',version="1.7")
par(mar=c(5, 5, 4, 2) + 0.1)
plot(APRsums$dG_corrected,APRsums$TANGO,xlim=c(-5.5,5.5),col=c(colors[1:length(mols)]),pch=19,cex=1.7,cex.lab=1.5,cex.axis=1.2,cex.main=2,ylim=c(0,max(APRsums$TANGO,na.rm = TRUE)),xlab=expression(paste("Summed ",Delta,"G",sep="")),ylab="Summed Tango score",main=APRsums$Name[1])
legend("topleft",legend=mols,col=c(colors[1:length(mols)]),pch=20,pt.cex=1,cex=0.7)
dev.off()

write.table(SolubisScoreSum,file="SolubisScore.txt",quote=FALSE,row.names=FALSE,sep="\t")

# stretches <- read.delim("SummaryStretches.txt", header=T)
# stretches$Name <- as.factor(stretches$Name)

