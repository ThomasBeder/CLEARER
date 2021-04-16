# 1. paralog features ##################################################################################
# run blastn 
system(paste("blastn -subject Sc_DNA.fa -query Sc_DNA.fa",
       "-outfmt", "'6 qseqid sseqid pident length qlen slen evalue bitscore'", 
       "-evalue 1e-30",
       "-out Sc_Sc_blastn",
       "-max_target_seqs 500",
       "-max_hsps 5",
       "-qcov_hsp_perc 10",
       sep = " "))

# generate paralog features

file <- read.csv("Sc_Sc_blastn",header = F, sep = "\t" )
head(file)
comNames <- paste(file$V1,file$V2,sep = ";")
file$comNames <- as.factor(comNames)
file <- file[order(file$comNames),]
PID <- tapply(file[,3], file[,9], median)
algLen <- tapply(file[,4], file[,9], sum)
qlen <- tapply(file[,5], file[,9], median)
slen <- tapply(file[,6], file[,9], median)

df <- as.data.frame(cbind(PID,algLen,qlen,slen))

m <- c()
for (j in 1:nrow(df)) {
  m[j] <- (df[j,3]+df[j,4])/2
}
df$m <- m
score <- c()
for (j in 1:nrow(df)) {
  score[j] <- df[j,2]/df[j,5]*(df[j,1]/100)
}
df$score <- score

names <- rownames(df)
namesShort <- c()
namesShort2 <- c()
for (i in 1:length(names)){
  namesShort[i] <- strsplit(names[i], '[;]')[[1]][1]
  namesShort2[i] <- strsplit(names[i], '[;]')[[1]][2]
}

df$names <- namesShort
df$names2 <- namesShort2

df <- df[-c(which(df$names == df$names2)),]

Pos <- sapply(as.character(unique(df$names)),function(x)which(x == df$names))

H000 <- c()
H005 <- c()
H010 <- c()
H020 <- c()
H030 <- c()
H040 <- c()
H050 <- c()
H060 <- c()
H070 <- c()
H080 <- c()
H090 <- c()
H095 <- c()

for (i in 1:length(Pos)){
  H000[i] <- length(which(df$score[Pos[[i]]]>0))
  H005[i] <- length(which(df$score[Pos[[i]]]>0.05))
  H010[i] <- length(which(df$score[Pos[[i]]]>0.1))
  H020[i] <- length(which(df$score[Pos[[i]]]>0.2))
  H030[i] <- length(which(df$score[Pos[[i]]]>0.3))
  H040[i] <- length(which(df$score[Pos[[i]]]>0.4))
  H050[i] <- length(which(df$score[Pos[[i]]]>0.5))
  H060[i] <- length(which(df$score[Pos[[i]]]>0.6))
  H070[i] <- length(which(df$score[Pos[[i]]]>0.7))
  H080[i] <- length(which(df$score[Pos[[i]]]>0.8))
  H090[i] <- length(which(df$score[Pos[[i]]]>0.9))
  H095[i] <- length(which(df$score[Pos[[i]]]>0.95))
}

para <- as.data.frame(cbind(names(Pos),H000,H005,H010,H020,H030,H040,H050,H060,H070,H080,H090,H095))
colnames(para)[1] <- "Gene"

# 2. homology features ##################################################################################
# attach required packages
library(seqinr)
fasta <- read.fasta("Sc_Prot.fa",forceDNAtolower = F,seqtype = "AA",as.string = T)
head(fasta)

dir.create("Sc_homology/")
dir.create("Sc_homology/Seqs/")

# split individual in fasta sequences
for (i in 1:length(fasta)){
  write.fasta(sequences = getSequence(fasta[i]),
              names = names(fasta[i]),
              file.out = paste("Sc_homology/Seqs/",names(fasta[i]),".fa",sep = ""),
              nbchar = 80
  )
}
write.table(names(fasta),"Sc_homology/Sequence_names", row.names = F)

# download the entire NCBI refseq (protein) database for example as described here:
# https://cran.r-project.org/web/packages/biomartr/vignettes/Database_Retrieval.html

# run psiblast script
system(paste("sh ","Sc_homology/psiblastSc.sh", sep = ""))

# add a line to prevent empty entries
system(paste("sh ","Sc_homology/cat.sh", sep = ""))

# load alignment results and generate features
list_path <-  dir("Sc_homology/Seqs/", pattern ="*.cat", full.names = T)
head(list_path)
list_path_short <-  dir("Sc_homology/Seqs/", pattern ="*.cat", full.names = F)
head(list_path_short)

Idx <- c()
for (i in 1:length(list_path_short)){
  Idx[i] <- strsplit(list_path_short[i], '[.]')[[1]][1]
}

Idx <- as.numeric(Idx)

names <- read.csv("Sc_homology/Sequence_names")

names <- as.character(names$x)
namesShort <- c()
for (i in 1:length(names)){
  namesShort[i] <- strsplit(names[i], '[-]')[[1]][1]
}
namesShort <- namesShort[Idx]
head(namesShort)
tail(namesShort)
##e-value cutoff 10-30,10-20,10-10,10-7,10- 5,10-3 (H30, H20, H10, H7, H5, H3)
gene <- as.character()
H100 <- c()
H080 <- c()
H050 <- c()
H030 <- c()
H020 <- c()
H010 <- c()
H007 <- c()
H005 <- c()
H003 <- c()

for (i in 1:length(list_path)) {
  file <- read.csv(list_path[i],header = F, sep = "\t")
  gene[i] <- namesShort[i]
  H100[i] <- length(which(file$V7 < 1.00e-100))
  H080[i] <- length(which(file$V7 < 1.00e-80))
  H050[i] <- length(which(file$V7 < 1.00e-50))
  H030[i] <- length(which(file$V7 < 1.00e-30))
  H020[i] <- length(which(file$V7 < 1.00e-20))
  H010[i] <- length(which(file$V7 < 1.00e-10))
  H007[i] <- length(which(file$V7 < 1.00e-07))
  H005[i] <-length(which(file$V7 < 1.00e-05))
  H003[i] <- length(which(file$V7 < 1.00e-03))
}

homo1 <- as.data.frame(cbind(H100,H080,H050,H030,H020,H010,H007,H005,H003))
homo1$Gene <- gene

H000 <- c()
H005 <- c()
H010 <- c()
H020 <- c()
H030 <- c()
H040 <- c()
H050 <- c()
H060 <- c()
H070 <- c()
H080 <- c()
H090 <- c()
H095 <- c()

for (i in 1:length(list_path)) {
  file <- read.csv(list_path[i],header = F, sep = "\t")
  comNames <- paste(file$V1,file$V2,sep = ";")
  file$comNames <- as.factor(comNames)
  file <- file[order(file$comNames),]
  PID <- tapply(file[,3], file[,9], median)
  algLen <- tapply(file[,4], file[,9], sum)
  qlen <- tapply(file[,5], file[,9], median)
  slen <- tapply(file[,6], file[,9], median)
  
  df <- as.data.frame(cbind(PID,algLen,qlen,slen))
  
  m <- c()
  for (j in 1:nrow(df)) {
    m[j] <- (df[j,3]+df[j,4])/2
  }
  df$m <- m
  score <- c()
  for (j in 1:nrow(df)) {
    score[j] <- df[j,2]/df[j,5]*(df[j,1]/100)
  }
  
  H000[i] <- length(which(score>0))
  H005[i] <- length(which(score>0.05))
  H010[i] <- length(which(score>0.1))
  H020[i] <- length(which(score>0.2))
  H030[i] <- length(which(score>0.3))
  H040[i] <- length(which(score>0.4))
  H050[i] <- length(which(score>0.5))
  H060[i] <- length(which(score>0.6))
  H070[i] <- length(which(score>0.7))
  H080[i] <- length(which(score>0.8))
  H090[i] <- length(which(score>0.9))
  H095[i] <- length(which(score>0.95))
}

homo2 <- as.data.frame(cbind(H000,H005,H010,H020,H030,H040,H050,H060,H070,H080,H090,H095))

homo2$Gene <-gene

# combine
colnames(homo1) <- paste("homo", colnames(homo1), sep = "_")
colnames(homo1)[10] <- "Gene"
colnames(homo2) <- paste("homoAd", colnames(homo2), sep = "_")
colnames(homo2)[13] <- "Gene"
colnames(para) <- paste("para", colnames(para), sep = "_")
colnames(para)[1] <- "Gene"

df <- merge.data.frame(homo1,homo2,"Gene", all = T)
df <- merge.data.frame(df,para,"Gene", all = T)
summary(df)

head(df)

write.table(df,"conservation_features.csv", row.names = F, sep = "\t")





