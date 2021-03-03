setwd("E:/Research_AJL/Cornell/silverside-simulation/")

require(vcfR)
require(OutFLANK)
require(qqman)
require(RColorBrewer)

ssvcf<-read.vcfR("results/1693883279556_SilverSide_Inversion_filt.recode.vcf")
#ssvcf <- read.vcfR("results/1615143132138_SilverSide_Inversion.vcf")
#ssvcf <- read.vcfR("1826701057597_fullNeut_filt_SilverSide_Inversion.recode.vcf")
grep("MT=[2]",ssvcf@fix[,8])
ssvcf_clean <- ssvcf[-grep("MT=[2]",ssvcf@fix[,8])]

#if neutral (i.e. no MT>1) then the below command will remove everything
#ssvcf_clean <- ssvcf

geno <- ssvcf_clean@gt[,-1] # Remove 1st column, which is 'Format'
position <- as.numeric(getPOS(ssvcf_clean)) # Positions in bp

for(i in 1:length(position)){
  if(duplicated(position)[i]){
    position[i]<-position[i]+0.1
    #print(position1_filt[i])
  }
}

G1 <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G1[geno %in% c("0/0", "0|0")] <- 0
G1[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G1[geno %in% c("1/1", "1|1")] <- 2

Gt<-t(G1)
colnames(Gt)<-paste("M",position,sep="")

PopsALL <- NULL
for(j in rep(1:2)){
  for(i in rep(j,5000)){
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
    FSTmat[i, ] = unlist(getFSTs_diploids2(popname, SNPmat[,i]))
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

getFSTs_diploids2 = function(popNameList, SNPDataColumn){
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
  
  out = WC_FST_Diploids_2Alleles2(HetCounts)
  return(out)
}

WC_FST_Diploids_2Alleles2 <- function (Sample_Mat) 
{
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat)
  r = n_pops
  n_c = (n_pops * n_ave - sum(sample_sizes^2)/(n_pops * n_ave))/(n_pops - 
                                                                   1)
  p_freqs = (Sample_Mat[, 3] + Sample_Mat[, 2]/2)/sample_sizes #fixed this - was Sample_Mat[, 1] which is hom. absent
  p_ave = sum(sample_sizes * p_freqs)/(n_ave * n_pops)
  s2 = sum(sample_sizes * (p_freqs - p_ave)^2)/((n_pops - 1) * 
                                                  n_ave)
  if (s2 == 0) {
    return(0)
    break
  }
  h_freqs = Sample_Mat[, 2]/sample_sizes
  h_ave = sum(sample_sizes * h_freqs)/(n_ave * n_pops)
  a <- n_ave/n_c * (s2 - 1/(n_ave - 1) * (p_ave * (1 - p_ave) - 
                                            ((r - 1)/r) * s2 - (1/4) * h_ave))
  b <- n_ave/(n_ave - 1) * (p_ave * (1 - p_ave) - (r - 1)/r * 
                              s2 - (2 * n_ave - 1)/(4 * n_ave) * h_ave)
  c <- 1/2 * h_ave
  aNoCorr <- n_ave/n_c * (s2)
  bNoCorr <- (p_ave * (1 - p_ave) - (r - 1)/r * s2 - (2 * n_ave)/(4 * 
                                                                    n_ave) * h_ave)
  cNoCorr <- 1/2 * h_ave
  He <- 1 - sum(p_ave^2, (1 - p_ave)^2)
  FST <- a/(a + b + c)
  FSTNoCorr = aNoCorr/(aNoCorr + bNoCorr + cNoCorr)
  return(list(He = He, FST = FST, T1 = a, T2 = (a + b + c), 
              FSTNoCorr = FSTNoCorr, T1NoCorr = aNoCorr, T2NoCorr = (aNoCorr + 
                                                                       bNoCorr + cNoCorr), meanAlleleFreq = p_ave))
}

Pfst<-MakeDiploidFSTMat_2(SNPmat = Gtf, locusNames = colnames(Gtf), popNames = PopsALL)
#write.csv(Pfst,"1615143132138_Pfst.csv",row.names = F)
# Pfst<-read.csv("1826581014997_Pfst.csv")
# Pfst<-read.csv("1826701057597_neut_Pfst.csv")

(Max <- max(Pfst$FST, na.rm=T))
(Min <- min(Pfst$FST, na.rm = T))
CHR_data <- data.frame (SNP = Pfst$LocusName , CHR=position_inv, BP= position_scaled, P=Pfst$FST)
# Plot !
manhattan(CHR_data , ylim=c(Min-(Min*0.1),Max+(Max*0.1)), suggestiveline = F, genomewideline = F , logp=F, col=brewer.pal(5, "Set2") , main="Silverside simulation\nMAF=0.1, m=0.1, r=(1e-5)x4, 1e-6, 1e-4", ylab=expression(paste("F"[ST])))
lines(x=c(34750,104250),y=c(Min-(Min*0.1),Min-(Min*0.1)),col="red",lwd=2)
lines(x=c(173750,243250),y=c(Min-(Min*0.1),Min-(Min*0.1)),col="red",lwd=2)
lines(x=c(312750,382250),y=c(Min-(Min*0.1),Min-(Min*0.1)),col="red",lwd=2)
abline(v=c(34750,104250,173750,243250,312750,382250), col="red")
 
Pop1_afreq<-rowSums(G1[,1:5000])/(2*ncol(G1[,1:5000]))

Pop1_afreq<-data.frame(Pop1_afreq)
#rownames(Pop1_afreq)<-paste("M",position,sep="")

Pop2_afreq<-rowSums(G1[,5001:10000])/(2*ncol(G1[,5001:10000]))

Pop2_afreq<-data.frame(Pop2_afreq)
#rownames(Pop2_afreq)<-paste("M",position,sep="")

Pfst$P1AlleleFreq <- Pop1_afreq
Pfst$P2AlleleFreq <- Pop2_afreq

Pfst$Position <- position
Pfst$Chrom <- position_inv
#write.csv(Pfst, "results/1693883279556_Pfst.csv", row.names = F)




Pfst[Pfst$Position==34750,]


SNPmat=Gtf
locusNames=colnames(Gtf)
popNames=PopsALL
  
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

popNameList=popname
SNPDataColumn=SNPmat[,184]
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


Sample_Mat = HetCounts

sample_sizes = rowSums(Sample_Mat)
n_ave = mean(sample_sizes)
n_pops = nrow(Sample_Mat)
r = n_pops
n_c = (n_pops * n_ave - sum(sample_sizes^2)/(n_pops * n_ave))/(n_pops - 1)
(p_freqs = (Sample_Mat[, 3] + Sample_Mat[, 2]/2)/sample_sizes)
(p_ave = sum(sample_sizes * p_freqs)/(n_ave * n_pops))
s2 = sum(sample_sizes * (p_freqs - p_ave)^2)/((n_pops - 1) * n_ave)
if (s2 == 0) {
  return(0)
  break
}
(h_freqs = Sample_Mat[, 2]/sample_sizes)
(h_ave = sum(sample_sizes * h_freqs)/(n_ave * n_pops))
a <- n_ave/n_c * (s2 - 1/(n_ave - 1) * (p_ave * (1 - p_ave) - 
                                          ((r - 1)/r) * s2 - (1/4) * h_ave))
b <- n_ave/(n_ave - 1) * (p_ave * (1 - p_ave) - (r - 1)/r * 
                            s2 - (2 * n_ave - 1)/(4 * n_ave) * h_ave)
c <- 1/2 * h_ave
aNoCorr <- n_ave/n_c * (s2)
bNoCorr <- (p_ave * (1 - p_ave) - (r - 1)/r * s2 - (2 * n_ave)/(4 * n_ave) * h_ave)
cNoCorr <- 1/2 * h_ave
He <- 1 - sum(p_ave^2, (1 - p_ave)^2)
FST <- a/(a + b + c)
FSTNoCorr = aNoCorr/(aNoCorr + bNoCorr + cNoCorr)
return(list(He = He, FST = FST, T1 = a, T2 = (a + b + c), 
            FSTNoCorr = FSTNoCorr, T1NoCorr = aNoCorr, T2NoCorr = (aNoCorr + 
                                                                     bNoCorr + cNoCorr), meanAlleleFreq = p_ave))
