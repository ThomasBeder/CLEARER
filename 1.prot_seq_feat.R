# attach required packages
library(seqinr)
library(stringr)
library(protr)
Seq <- read.fasta("/media/tbeder/Elements1/Rechner3/Essential_genes_AnG/Olufemi_data/Sc_Prot.fa",forceDNAtolower = F,seqtype = "AA",as.string = T)

# read in protein seqeunces
Seq <- read.fasta("Sc_prot.fa",forceDNAtolower = F,seqtype = "AA",as.string = T)
names <- names(Seq)
Seq <- getSequence(Seq,as.string = T)
Seq <- as.character(unlist(Seq))
names(Seq) <- names

# remove sequences < 60 amino acids if there are any
#len <- c()
#for (i in 1:length(Seq)){
#  len[i] <- length(s2c(as.character(Seq[i])))
#}
#Seq <- Seq[-c(which(len < 60))]
#names <- names[-c(which(len < 60))]

# seqinr features
names <- c()
Tiny <- c()
Small <- c()
Aliphatic <- c()
Aromatic <- c()
Non.polar <- c()
Polar <- c()
Charged <- c()
Basic <- c()
Acidic <- c()
Pi <- c()

for (i in 1:length(Seq)) {
  names[i] <- as.character(names[i])
  Tiny[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Tiny
  Small[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Small
  Aliphatic[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Aliphatic
  Aromatic[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Aromatic
  Non.polar[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Non.polar
  Polar[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Polar
  Charged[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Charged
  Basic[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Basic
  Acidic[i] <- AAstat(s2c(as.character(Seq[i])),plot = F)$Prop$Acidic
  Pi[i] <- AAstat(s2c(as.character(Seq[1])),plot = F)$Pi
}

datalist = list()
for (i in 1:length(Seq)) {
  dat <- AAstat(s2c(as.character(Seq[i])),plot =F)$Compo
  dat$i <- i  
  datalist[[i]] <- dat 
}
big_data = as.table(do.call(rbind, datalist))
dim(big_data)
write.csv(big_data,"AA_numbers_of_genes.csv",row.names = F)
big_data <- read.csv("AA_numbers_of_genes.csv")
sapply(big_data, class)

PA <- c()
PC <- c()
PD <- c()
PE <- c()
PF <- c()
PG <- c()
PH <- c()
PI <- c()
PK <- c()
PL <- c()
PM <- c()
PN <- c()
PP <- c()
PQ <- c()
PR <- c()
PS <- c()
PT <- c()
PV <- c()
PW <- c()
PY <- c()

for (i in 1:length(rownames(big_data))) {
  PA[i] <- big_data[i,2]/sum(big_data[i,c(2:21)])
  PC[i] <- big_data[i,3]/sum(big_data[i,c(2:21)])
  PD[i] <- big_data[i,4]/sum(big_data[i,c(2:21)])
  PE[i] <- big_data[i,5]/sum(big_data[i,c(2:21)])
  PF[i] <- big_data[i,6]/sum(big_data[i,c(2:21)])
  PG[i] <- big_data[i,7]/sum(big_data[i,c(2:21)])
  PH[i] <- big_data[i,8]/sum(big_data[i,c(2:21)])
  PI[i] <- big_data[i,9]/sum(big_data[i,c(2:21)])
  PK[i] <- big_data[i,10]/sum(big_data[i,c(2:21)])
  PL[i] <- big_data[i,11]/sum(big_data[i,c(2:21)])
  PM[i] <- big_data[i,12]/sum(big_data[i,c(2:21)])
  PN[i] <- big_data[i,13]/sum(big_data[i,c(2:21)])
  PP[i] <- big_data[i,14]/sum(big_data[i,c(2:21)])
  PQ[i] <- big_data[i,15]/sum(big_data[i,c(2:21)])
  PR[i] <- big_data[i,16]/sum(big_data[i,c(2:21)])
  PS[i] <- big_data[i,17]/sum(big_data[i,c(2:21)])
  PT[i] <- big_data[i,18]/sum(big_data[i,c(2:21)])
  PV[i] <- big_data[i,19]/sum(big_data[i,c(2:21)])
  PW[i] <- big_data[i,20]/sum(big_data[i,c(2:21)])
  PY[i] <- big_data[i,21]/sum(big_data[i,c(2:21)])
}        

comb <- as.data.frame(cbind(names(Seq),PA,PC,PD,PE,PF,PG,PH,PI,PK,PL,PM,PN,PP,PQ,PR,PS,PT,PV,PW,PY,Tiny,
                            Small,Aliphatic,Non.polar,Polar,Charged,Basic,Acidic,Pi ))
head(comb)
write.csv(comb,"seqinr.csv",row.names = F)                      

# protr features
seqs <- Seq
head(seqs)
## remove * character at the end from each sequence ##
seqs <- str_sub(seqs, 1, str_length(seqs)-1)
names(seqs) <- names
seqs <- seqs[(sapply(seqs, protcheck))]

AAC <- t(sapply(seqs, extractAAC))
head(AAC)
dim(AAC)
DC <- t(sapply(seqs, extractDC))
dim(DC)
TC <- t(sapply(seqs, extractTC))
dim(TC)
MoreauBroto <- t(sapply(seqs, extractMoreauBroto))
dim(MoreauBroto)
Moran <- t(sapply(seqs, extractMoran))
dim(Moran)
Geary <- t(sapply(seqs, extractGeary))
dim(Geary)
CTDC <- t(sapply(seqs, extractCTDC))
dim(CTDC)
CTDT <- t(sapply(seqs, extractCTDT))
dim(CTDT)
CTDD <- t(sapply(seqs, extractCTDD))
dim(CTDD)
CTriad <- t(sapply(seqs, extractCTriad))
dim(CTriad)
SOCN <- t(sapply(seqs, extractSOCN))
dim(SOCN)
QSO <- t(sapply(seqs, extractQSO))
dim(QSO)
PAAC <- t(sapply(seqs, extractPAAC))
dim(PAAC)
APAAC <- t(sapply(seqs, extractAPAAC))
dim(APAAC)

### rownames are all in the same order so there is no need to sort ####

comb <- as.data.frame(cbind(AAC,DC,TC,MoreauBroto,Moran,Geary,CTDC,CTDT,CTDD,CTriad,SOCN,QSO,PAAC,APAAC))


for(i in 1:ncol(comb)){
  comb[is.na(comb[,i]), i] <- mean(comb[,i], na.rm = TRUE)
}
dim(comb)
write.csv(comb,"~/Essential_genes_AnG/protein_features/Ce_protr.csv")



