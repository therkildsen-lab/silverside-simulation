setwd("E:/Research_AJL/Cornell/silverside-simulation/")

require(vcfR)
require(OutFLANK)
require(qqman)
require(RColorBrewer)

ssvcf <- read.vcfR("1614591178478_SilverSide_Inversion.vcf")

ssvcf_clean <- ssvcf[-grep("MT=[2-4]",ssvcf@fix[,8])]

geno <- ssvcf_clean@gt[,-1] # Remove 1st column, which is 'Format'
position <- as.numeric(getPOS(ssvcf_clean)) # Positions in bp

# for(i in 1:length(position)){
#   if(position[i]>34750 & position[i]<104250){
#     position[i]<-104250-position[i]+34750
#   }
#   if(position[i]>173750 & position[i]<243250){
#     position[i]<-243250-position[i]+173750
#       }
#   if(position[i]>312750 & position[i]<382250){
#     position[i]<-382250-position[i]+312750
#       }
# }

G1 <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G1[geno %in% c("0/0", "0|0")] <- 0
G1[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G1[geno %in% c("1/1", "1|1")] <- 2

Gt<-t(G1)
colnames(Gt)<-paste("M",position,sep="")

PopsALL <- NULL
for(j in rep(1:2)){
  for(i in rep(j,100)){
    PopsALL <- c(PopsALL,i)
  }
}

position_inv<-NULL
position_scaled<-NULL
for(i in 1:length(position)){
  if(position[i]>0 & position[i]<139001){
    position_inv<-c(position_inv,1)
    position_scaled<-c(position_scaled,position[i])
  }
  if(position[i]>139000 & position[i]<278001){
    position_inv<-c(position_inv,2)
    position_scaled<-c(position_scaled,position[i]-139000)
  }
  if(position[i]>278000 & position[i]<417001){
    position_inv<-c(position_inv,3)
    position_scaled<-c(position_scaled,position[i]-278000)
  }
  if(position[i]>417000 & position[i]<556001){
    position_inv<-c(position_inv,4)
    position_scaled<-c(position_scaled,position[i]-417000)
  }
  if(position[i]>556000 & position[i]<695001){
    position_inv<-c(position_inv,5)
    position_scaled<-c(position_scaled,position[i]-556000)
  }
  if(position[i]>695000 & position[i]<834001){
    position_inv<-c(position_inv,6)
    position_scaled<-c(position_scaled,position[i]-695000)
  }
}

Gtf <- data.frame(Gt)

getFSTs_diploids = function(popNameList, SNPDataColumn){
  #eliminating the missing data for this locus
  popnames=unlist(as.character(popNameList))
  popNameTemp=popnames[which(SNPDataColumn!=9)]
  snpDataTemp=SNPDataColumn[SNPDataColumn!=9]
  
  HetCounts <- tapply(snpDataTemp, list(popNameTemp,snpDataTemp), length)
  HetCounts[is.na(HetCounts)] = 0
  
  #Case: all individuals are genetically identical at this locus
  if(dim(HetCounts)[2]==1){
    return (list(He=NA,FST=NA, T1=NA, T2=NA,FSTNoCorr=NA, T1NoCorr=NA, T2NoCorr=NA,meanAlleleFreq = NA))
  }
  
  if(dim(HetCounts)[2]==2){
    if(paste(colnames(HetCounts),collapse="")=="01"){HetCounts=cbind(HetCounts,"2"=0)}
    if(paste(colnames(HetCounts),collapse="")=="12"){HetCounts=cbind("0"=0,HetCounts)}
    if(paste(colnames(HetCounts),collapse="")=="02"){HetCounts=cbind(HetCounts[,1],"1"=0, HetCounts[,2])}
  }
  
  out = WC_FST_Diploids_2Alleles(HetCounts)
  return(out)
}

MakeDiploidFSTMat_2<-function(SNPmat,locusNames,popNames){
  locusname <- unlist(locusNames)
  popname <- unlist(popNames)
  snplevs <- levels(as.factor(unlist(SNPmat)))
  if(any(!(snplevs%in%c(0,1,2,9)))==TRUE) {
    print("Error: Your snp matrix has a character other than 0,1,2 or 9")
    break
  }
  if (dim(SNPmat)[1] != length(popname)) {
    print("Error: your population names do not match your SNP matrix")
    break
  }
  if (dim(SNPmat)[2] != length(locusname)) {
    print("Error:  your locus names do not match your SNP matrix")
    break
  }
  writeLines("Calculating FSTs, may take a few minutes...")
  nloci <- length(locusname)
  FSTmat <- matrix(NA, nrow = nloci, ncol = 8)
  for (i in 1:nloci) {
    FSTmat[i, ] = unlist(getFSTs_diploids(popname, SNPmat[,i]))
    if (i%%10000 == 0) {
      print(paste(i, "done of", nloci))
    }
  }
  outTemp = as.data.frame(FSTmat)
  outTemp = cbind(locusname, outTemp)
  colnames(outTemp) = c("LocusName", "He", "FST", "T1", "T2",
                        "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq")
  return(outTemp)
}

Pfst<-MakeDiploidFSTMat_2(SNPmat = Gtf, locusNames = colnames(Gtf), popNames = PopsALL)
(Max <- max(Pfst$FST, na.rm=T))
(Min <- min(Pfst$FST, na.rm = T))
CHR_data <- data.frame (SNP = Pfst$LocusName , CHR=position_inv, BP= position_scaled, P=Pfst$FST)
# Plot !
manhattan(CHR_data , ylim=c(Min-(Min*0.1),Max+(Max*0.1)), suggestiveline = F, genomewideline = F , logp=F, col=brewer.pal(5, "Set2") , main="Silverside simulation\nm=0.001, r=(1e-5)x4, 1e-6, 1e-4", ylab=expression(paste("F"[ST])))
lines(x=c(34750,104250),y=c(Min-(Min*0.1),Min-(Min*0.1)),col="red",lwd=2)
lines(x=c(173750,243250),y=c(Min-(Min*0.1),Min-(Min*0.1)),col="red",lwd=2)
lines(x=c(312750,382250),y=c(Min-(Min*0.1),Min-(Min*0.1)),col="red",lwd=2)
abline(v=c(34750,104250,173750,243250,312750,382250), col="red")
 