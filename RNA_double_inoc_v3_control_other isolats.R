########################################################################################
#' Import transcript-level abundances and estimated counts for gene-level analysis packages
########################################################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")
#install.packages("gplots")
library("gplots")
library("RColorBrewer")
library(Rgraphviz)
library(topGO)
library(edgeR)
library(limma)
#libraries
library('nlme')
#options("scipen"=100,digits=3) 
library('lme4')
library('ggplot2')
library(reshape)
library("multcomp")
#biocLite(c("GO.db", "preprocessCore", "impute"))
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
#install.packages("~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/WGCNA_1.49.tgz",type = "source", repos = NULL, lib=.Library) 
library(WGCNA)
library(pheatmap)
########################################################################################
# ANALYSIS AT GENE LEVEL 
# input data: Kallisto pseudo-alignment. 

tximport <- function(files,
                     type=c("none","kallisto","salmon","sailfish","rsem"),
                     txIn=TRUE,
                     txOut=FALSE,
                     countsFromAbundance=c("no","scaledTPM","lengthScaledTPM"),
                     tx2gene=NULL,
                     reader=read.delim,
                     geneIdCol,
                     txIdCol,
                     abundanceCol,
                     countsCol,
                     lengthCol,
                     importer,
                     collatedFiles,
                     ignoreTxVersion=FALSE) {
  
  type <- match.arg(type, c("none","kallisto","salmon","sailfish","rsem"))
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  stopifnot(all(file.exists(files)))
  
  # kallisto presets
  if (type == "kallisto") {
    geneIdCol="gene_id"
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    importer <- reader
  }
  
  # salmon/sailfish presets
  if (type %in% c("salmon","sailfish")) {
    geneIdCol="gene_id"
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "EffectiveLength"
    importer <- function(x) reader(x, comment='#') 
  }
  
  # rsem presets
  if (type == "rsem") {
    txIn <- FALSE
    geneIdCol <- "gene_id"
    abundanceCol <- "FPKM"
    countsCol <- "expected_count"
    lengthCol <- "effective_length"
    importer <- reader
  }
  
  if (type == "cufflinks") {
    stop("reading from collated files not yet implemented")
  }
  
  # if input is tx-level, need to summarize abundances, counts and lengths to gene-level
  if (txIn) {
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))
      
      #####################################################################
      # some temporary code for detecting older fishes
      if ((i == 1) &
          (type %in% c("salmon","sailfish")) &
          !("EffectiveLength" %in% names(raw))) {
        lengthCol <- "Length" 
        # because the comment lines have the same comment character
        # as the header, need to name the column names
        importer <- function(x) {
          tmp <- reader(x, comment="#")
          names(tmp) <- c("Name","Length","TPM","NumReads")
          tmp
        }
        # re-read the first file
        raw <- as.data.frame(importer(files[i]))
      }
      #####################################################################
      
      # does the table contain gene association or was an external tx2gene table provided?
      if (is.null(tx2gene) & !txOut) {
        # e.g. Cufflinks includes the gene ID in the table
        stopifnot(all(c(geneIdCol, lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          geneId <- raw[[geneIdCol]]
        } else {
          stopifnot(all(geneId == raw[[geneIdCol]]))
        }
      } else {
        # e.g. Salmon and kallisto do not include the gene ID, need an external table
        stopifnot(all(c(lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          txId <- raw[[txIdCol]]
        } else {
          stopifnot(all(txId == raw[[txIdCol]]))
        }
      }
      # create empty matrices
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[txIdCol]]
        colnames(mat) <- names(files)
        abundanceMatTx <- mat
        countsMatTx <- mat
        lengthMatTx <- mat
      }
      abundanceMatTx[,i] <- raw[[abundanceCol]]
      countsMatTx[,i] <- raw[[countsCol]]
      lengthMatTx[,i] <- raw[[lengthCol]]
    }
    message("")
    
    txi <- list(abundance=abundanceMatTx, counts=countsMatTx, length=lengthMatTx,
                countsFromAbundance="no")
    
    # if the user requested just the transcript-level data:
    if (txOut) {
      return(txi)
    }
    
    txi[["countsFromAbundance"]] <- NULL
    txiGene <- summarizeToGene(txi, tx2gene, ignoreTxVersion, countsFromAbundance)
    return(txiGene)  
    
    # e.g. RSEM already has gene-level summaries
    # just combine the gene-level summaries across files
  } else {
    # stating the obvious:
    if (txOut) stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
    
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))
      stopifnot(all(c(geneIdCol, abundanceCol, lengthCol) %in% names(raw)))
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[geneIdCol]]
        colnames(mat) <- names(files)
        abundanceMat <- mat
        countsMat <- mat
        lengthMat <- mat
      }
      abundanceMat[,i] <- raw[[abundanceCol]]
      countsMat[,i] <- raw[[countsCol]]
      lengthMat[,i] <- raw[[lengthCol]]
    }
  } 
  message("")
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance="no"))
}

# summarizeToGene() splits out the summarization functions
# in tximport(), so it can be called by users to summarize
# transcript-level lists of matrices

#' @describeIn tximport Summarize tx-level matrices to gene-level
#' @export
summarizeToGene <- function(txi,
                            tx2gene,
                            ignoreTxVersion=FALSE,
                            countsFromAbundance=c("no","scaledTPM","lengthScaledTPM")
) {
  
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  
  # unpack matrices from list for cleaner code
  abundanceMatTx <- txi$abundance
  countsMatTx <- txi$counts
  lengthMatTx <- txi$length
  
  txId <- rownames(abundanceMatTx)
  stopifnot(all(txId == rownames(countsMatTx)))
  stopifnot(all(txId == rownames(lengthMatTx)))
  
  # need to associate tx to genes
  # potentially remove unassociated transcript rows and warn user
  if (!is.null(tx2gene)) {
    colnames(tx2gene) <- c("tx","gene")
    if (ignoreTxVersion) {
      txId <- sapply(strsplit(as.character(txId), "\\."), "[[", 1)
    }
    tx2gene$gene <- factor(tx2gene$gene)
    tx2gene$tx <- factor(tx2gene$tx)
    # remove transcripts (and genes) not in the abundances
    tx2gene <- tx2gene[tx2gene$tx %in% txId,]
    tx2gene$gene <- droplevels(tx2gene$gene)
    ntxmissing <- sum(!txId %in% tx2gene$tx)
    if (ntxmissing > 0) message("transcripts missing genes: ", ntxmissing)
    sub.idx <- txId %in% tx2gene$tx
    abundanceMatTx <- abundanceMatTx[sub.idx,,drop=FALSE]
    countsMatTx <- countsMatTx[sub.idx,,drop=FALSE]
    lengthMatTx <- lengthMatTx[sub.idx,,drop=FALSE]
    txId <- txId[sub.idx]
    geneId <- tx2gene$gene[match(txId, tx2gene$tx)]
  }
  
  # summarize abundance and counts
  message("summarizing abundance")
  abundanceMat <- fastby(abundanceMatTx, geneId, colSums)
  message("summarizing counts")
  countsMat <- fastby(countsMatTx, geneId, colSums)
  message("summarizing length")
  
  # the next lines calculate a weighted average of transcript length, 
  # weighting by transcript abundance.
  # this can be used as an offset / normalization factor which removes length bias
  # for the differential analysis of estimated counts summarized at the gene level.
  weightedLength <- fastby(abundanceMatTx * lengthMatTx, geneId, colSums)
  lengthMat <- weightedLength / abundanceMat   
  
  # pre-calculate a simple average transcript length
  # for the case the abundances are all zero for all samples.
  # first, average the tx lengths over samples
  aveLengthSamp <- rowMeans(lengthMatTx)
  # then simple average of lengths within genes (not weighted by abundance)
  aveLengthSampGene <- tapply(aveLengthSamp, geneId, mean)
  
  stopifnot(all(names(aveLengthSampGene) == rownames(lengthMat)))
  
  # check for NaN and if possible replace these values with geometric mean of other samples.
  # (the geometic mean here implies an offset of 0 on the log scale)
  # NaN come from samples which have abundance of 0 for all isoforms of a gene, and 
  # so we cannot calculate the weighted average. our best guess is to use the average
  # transcript length from the other samples.
  lengthMat <- replaceMissingLength(lengthMat, aveLengthSampGene)
  
  if (countsFromAbundance != "no") {
    countsSum <- colSums(countsMat)
    if (countsFromAbundance == "lengthScaledTPM") {
      newCounts <- abundanceMat * rowMeans(lengthMat)
    } else {
      newCounts <- abundanceMat
    }
    newSum <- colSums(newCounts)
    countsMat <- t(t(newCounts) * (countsSum/newSum))
  }
  
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance=countsFromAbundance))
}

# this is much faster than by(), a bit slower than dplyr summarize_each()
fastby <- function(m, f, fun) {
  idx <- split(1:nrow(m), f)
  if (ncol(m) > 1) {
    t(sapply(idx, function(i) fun(m[i,,drop=FALSE])))
  } else {
    matrix(sapply(idx, function(i) fun(m[i,,drop=FALSE])),
           dimnames=list(levels(f), colnames(m)))
  }
}

# function for replacing missing average transcript length values
replaceMissingLength <- function(lengthMat, aveLengthSampGene) {
  nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
  if (length(nanRows) > 0) {
    for (i in nanRows) {
      if (all(is.nan(lengthMat[i,]))) {
        # if all samples have 0 abundances for all tx, use the simple average
        lengthMat[i,] <- aveLengthSampGene[i]
      } else {
        # otherwise use the geometric mean of the lengths from the other samples
        idx <- is.nan(lengthMat[i,])
        lengthMat[i,idx] <-  exp(mean(log(lengthMat[i,!idx]), na.rm=TRUE))
      }
    }
  }
  lengthMat
}

#################################################################################
#################################################################################
# transfrom transcripts to genes.
####################################################################
#Mesculenta
setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta') 

filesToProcess <- dir(pattern = "*_abundance.tsv$")  #files to process.
names(filesToProcess)<-gsub('_abundance.tsv$','',filesToProcess,perl=TRUE)
samples<-read.table('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/Mesculenta_V1_10_B1_abundance.tsv',h=T)

tx2gene<-cbind.data.frame(samples$target_id,gsub('.[0-9].v6.1$','',samples$target_id,perl=TRUE))
colnames(tx2gene)<-c('TXNAME','GENEID')
txi_YUCA <- tximport(filesToProcess, type="kallisto", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM") # normalized library size and transcript length
FOUR_VARS_YUCA<-txi_YUCA[[2]]
colnames(FOUR_VARS_YUCA)<-gsub("Mesculenta_","",colnames(FOUR_VARS_YUCA),perl=T)

#Rirregularis
#setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Umap_Rirregularis')   
setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Umap_Rirregularis/Control_other_genomes')
filesToProcess <- dir(pattern = "*_C2ref_abundance.tsv$")  #files to process.
names(filesToProcess)<-gsub('_abundance.tsv$','',filesToProcess,perl=TRUE)
samples<-read.table('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Umap_Rirregularis/Control_other_genomes/m2V513_CAN_C2ref_abundance.tsv',h=T)

tx2gene<-cbind.data.frame(samples$target_id,gsub('.[0-9].v6.1$','',samples$target_id,perl=TRUE))
colnames(tx2gene)<-c('TXNAME','GENEID')

head(tx2gene)
txi_AMF <- tximport(filesToProcess, type="kallisto", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM") # normalized library size and transcript length
FOUR_VARS_AMF<-txi_AMF[[2]]#[, c(-15,-45,-46,-40,-41)] # exclude repetead lib
colnames(FOUR_VARS_AMF)<-gsub("UmapAMF_"," ",colnames(FOUR_VARS_AMF),perl=T)
colnames(FOUR_VARS_AMF)
########################################################################################################################################


VAR<-gsub( "m","",gsub( "m2","",gsub("_\\w+$"," ",colnames(FOUR_VARS_AMF),perl=T)    ) )
TREAT<-c("CAN","CAN","B1","CANB1","B1","CANB1","CAN","B1","CANB1")
#TREAT_AMF<-factor(TREAT,levels=c("CANB1","CAN","B1"))
#TREAT_AMF<-factor(TREAT,levels=c("CAN","B1","CANB1"))

design_AMF <- model.matrix(~TREAT)

DGE_AMF <- DGEList(FOUR_VARS_AMF)

#FILTERING  

keep_A <- rowSums(DGE_AMF$counts>100) >= 3
DGE_AMF<- DGE_AMF[keep_A,]

#NORMALIZATION
DGE_AMF_N <- voom(DGE_AMF, design_AMF,plot=TRUE)


plot(DGE_AMF$counts[grep("g15421.t1",rownames(DGE_AMF$counts)),]~TREAT,las=2)


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/V4_coinoc /REf_C2_g7296t1.pdf", width=10, height=8, useDingbats = F)

genes<-DGE_AMF$counts[grep("g7296.t1",rownames(DGE_AMF$counts)),]
stDevs <-tapply(genes,TREAT,sd)
means<-tapply(genes,TREAT,mean)
mp<-barplot(tapply(genes,TREAT,mean),ylim=c(0,max(means + stDevs)*1.1),
            beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
            )
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
dev.off()


Results_A_vf2[rownames(Results_A_vf2)=="g3179.t1",]




########################################################################################################################################


#annotation infos
########################################################################################################################################
#M esculenta
mesculenta_go<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/Mesculenta_305_v6.1.annotation_info.txt',h=T)
mesculenta_go2<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/annotations_blast2go_mesculenta_FMv2/blast2go_FM_v3',h=T)
colnames(mesculenta_go2)<-c("locusName","transcriptName","GO","DEF")
#R. irregularis
rirregularis_go<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/blast2go_annot_predicted_prot_hint_glomus_nu6_genome_masked.annot',h=F)
colnames(rirregularis_go)<-c('locusName','GO.ID')
head(mesculenta_go)

########################################################################################################################################
####################################################################
# SEQUENCING RESULTS ANALYSIS
seq<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/Sequencing_results_V1.txt",h=T)

seq<-seq[seq$Variety!="V1",]
seq<-seq[seq$Variety!="V3",]
seq<-seq[seq$Variety!="V8",]

seq<-seq[seq$Samples!="V4-16b",];seq<-seq[seq$Samples!="V4-17",]
seq<-seq[seq$Samples!="V48",];seq<-seq[seq$Samples!="m2V518",]
seq<-seq[seq$Samples!="V57",];seq<-seq[seq$Samples!="V56",]
seq<-seq[seq$Samples!="V5-18",];seq<-seq[seq$Samples!="V5-19",]


seq<-seq[seq$Samples!="V5-9",];seq<-seq[seq$Samples!="V5-1",]
seq<-seq[seq$Samples!="V58",];seq<-seq[seq$Samples!="V5-4",]


## NEED TO FILTER THE SAMPLES TO USE
head(seq[,c(1,3,4,5,6,7,9:11,13)])
colnames(seq)[c(1,3,4,5,6,7,9:11,13)]<-c("Samples","Cultivar","AMF_treatment","Total_reads","Reads_cassava","TPM>100_cassava"
                                         ,"Unmapped_reads_cassava", "Reads_R.irregularis","TPM>100_R.irregularis","sites_cov>10_R.irregularis")
write.table(seq[,c(1,3,4,5,6,7,9:11,13)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Supplementary_table1.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)

#ACCUMULATION CURVES

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/0.Seq_effort_Mesculenta.pdf", width=14, height=10)
par(las=1);par(mfrow=c(1,1));par(mar=c(10,10,4,3),mgp=c(8, 2, 0));options(digits=3)
c("dodgerblue3","lightsalmon3","darkmagenta","azure3")

plot(seq$tpm100_M[seq$Treatment=="CTRL"]~seq$pseudoalign_M[seq$Treatment=="CTRL"], ylim=c(0,3000), cex.lab=3,cex.axis=3, pch=17,col="azure3",cex=3,ylab="Nb. genes >100 counts",xlab="Reads")
points(seq$tpm100_M[seq$Treatment=="CAN"]~seq$pseudoalign_M[seq$Treatment=="CAN"], ylim=c(0,3000), cex.lab=3,cex.axis=3, pch=16,col="lightsalmon3",cex=3)
points(seq$tpm100_M[seq$Treatment=="B1"]~seq$pseudoalign_M[seq$Treatment=="B1"], ylim=c(0,3000), cex.lab=3,cex.axis=3, pch=15,col="dodgerblue3",cex=3)
points(seq$tpm100_M[seq$Treatment=="CANB1"]~seq$pseudoalign_M[seq$Treatment=="CANB1"], ylim=c(0,3000), cex.lab=3,cex.axis=3, pch=15,col="darkmagenta",cex=3)
dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/0.Seq_effort_Rirregularis.pdf", width=14, height=10)
par(las=1);par(mfrow=c(1,1));par(mar=c(10,10,4,3),mgp=c(8, 2, 0));options(digits=3)

plot(seq$tpm_100_R[seq$Treatment=="CTRL"] ~seq$pseudoaligned_R[seq$Treatment=="CTRL"], ylim=c(0,2500), xlim=c(0,8000000),cex.lab=3,cex.axis=3, pch=17,col="azure3",cex=3,ylab="Nb. genes >100 counts",xlab="Reads")
points(seq$tpm_100_R[seq$Treatment=="CAN"] ~seq$pseudoaligned_R[seq$Treatment=="CAN"], ylim=c(0,2000), cex.lab=3,cex.axis=3, pch=16,col="lightsalmon3",cex=3)
points(seq$tpm_100_R[seq$Treatment=="B1"] ~seq$pseudoaligned_R[seq$Treatment=="B1"], ylim=c(0,2000), cex.lab=3,cex.axis=3, pch=15,col="dodgerblue3",cex=3)
points(seq$tpm_100_R[seq$Treatment=="CANB1"] ~seq$pseudoaligned_R[seq$Treatment=="CANB1"], ylim=c(0,2000), cex.lab=3,cex.axis=3, pch=15,col="darkmagenta",cex=3)

dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/0.ControlY_TPM100_treat.pdf")
par(las=2);par(mfrow=c(1,1));par(mar=c(10,10,4,3),mgp=c(8, 1, 0));options(digits=3)
boxplot(seq$tpm100_M~seq$Treatment,ylab="Transcripts > 100 counts",ylim=c(0,3000),xlab=NULL,
        col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.lab=3,cex.axis=3)
kruskal.test(seq$tpm100_M~seq$Treatment)
dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/0.ControlA_TPM100_treat.pdf")
par(las=2);par(mfrow=c(1,1));par(mar=c(10,10,4,3),mgp=c(8, 1, 0));options(digits=3)
boxplot(seq$tpm_100_R~seq$Treatment,ylab="Transcripts > 100 counts",ylim=c(0,3000),xlab=NULL,
        col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.lab=3,cex.axis=3)
kruskal.test(seq$tpm_100_R[seq$Treatment!="CTRL"]~seq$Treatment[seq$Treatment!="CTRL"])

dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/0.ControlY_TPM100_var.pdf")
par(las=2);par(mfrow=c(1,1));par(mar=c(10,10,4,3),mgp=c(8, 2, 0));options(digits=3)
seq$Cultivar[seq$Cultivar=="V4_COL2215"]<-"V4"
boxplot(seq$tpm100_M~droplevels(seq$Variety),ylab="Transcripts > 100 counts",ylim=c(0,3000),xlab=NULL,
        col=c("aquamarine4","burlywood","darkorange4"),cex.lab=3,cex.axis=3)
kruskal.test(seq$tpm100_M~droplevels(seq$Variety))
dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/0.ControlA_TPM100_var.pdf")
par(las=2);par(mfrow=c(1,1));par(mar=c(10,10,4,3),mgp=c(8, 2, 0));options(digits=3)
boxplot(seq$tpm_100_R[seq$Treatment!="CTRL"]~droplevels(seq$Variety[seq$Treatment!="CTRL"]),ylab="Transcripts > 100 counts",ylim=c(0,3000),xlab=NULL,
        col=c("aquamarine4","burlywood","darkorange4"),cex.lab=3,cex.axis=3)

slas<-lm(seq$tpm_100_R[seq$Treatment!="CTRL"]~droplevels(seq$Variety[seq$Treatment!="CTRL"]))
summary(slas)
lsmeans(slas, pairwise )

dev.off()

#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/0.Control_interaction_TPM100.pdf")
#boxplot(seq$tpm_100_R~interaction(seq$Cultivar,seq$Treatment),ylab="Transcripts > 100 counts",las=2,ylim=c(0,3000))
#dev.off()
########################################################################################################################################

### CONTROL Symbiosis related genes in dataset.

Fred_YUCA_genes<-c("Manes.06G156700","Manes.12G124500","Manes.01G123300",
                   "Manes.04G002400","Manes.12G124300","Manes.05G053600",
                   "Manes.08G154600","Manes.18G034100","Manes.16G125700",
                   "Manes.06G143100","Manes.06G143100","Manes.14G029600",
                   "Manes.01G071000")
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/0.Control_Symbiosis_genes.pdf',height=8,width = 10)

symb<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Fred_YUCA_genes,]
rownames(symb)<-rownames(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Fred_YUCA_genes,])
pheatmap(symb,cellwidth =8 ,cellheight=8)
dev.off()




########################################################################################################################################
#annotation infos
########################################################################################################################################
#M esculenta
mesculenta_go<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/Mesculenta_305_v6.1.annotation_info.txt',h=T)
mesculenta_go2<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/annotations_blast2go_mesculenta_FMv2/blast2go_FM_v3',h=T)
colnames(mesculenta_go2)<-c("locusName","transcriptName","GO","DEF")
#R. irregularis
rirregularis_go<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/blast2go_annot_predicted_prot_hint_glomus_nu6_genome_masked.annot',h=F)
colnames(rirregularis_go)<-c('locusName','GO.ID')
head(mesculenta_go)

########################################################################################################################################
# CASSAVA
Mercator_Mesculenta<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Mesculenta_database_v4.txt",h=F,quote="\'")
#replace 'NA'	'NA'	(\w+)\"(\b)   (\w+)\'(\b)   (\w+)\'\, (\w+)\'-  (\w+)\'\( (\w+)\'\) \.\'(\w+)   (\w+)\'\] \(\'(\w+)
colnames(Mercator_Mesculenta)<-c("Bin","Function","gene","details","Type")
head(Mercator_Mesculenta)
# AMF
Mercator_Rirregularis<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Rirregularis_database_v4",h=F,quote="\'")
colnames(Mercator_Rirregularis)<-c("Bin","Function","gene","details")

####################################################################
# DATA SELECTION
####################################################################
# Selection all co-inoculation , no controls

#FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("CTRL"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("CTRL"), colnames(FOUR_VARS_AMF),invert =T)]

# Exclude bad coinoculated (Seen in SNP data)
#??V3-3 (CAN), V1-7 (CAN), V8-11 (B1), V3-15 (B1), V3-19A (CANB1)
# Exclude not conlonized ( a lot of NA in SNP DATA)
# V5-18, V5-9A, V3-19A, V8-6, m2V5-18_B1, V8-18_B1
# Exclude not clusterizing in AMF groups
# V6_14_B1, V5_19_CANB1, V4_17_CAN , V8_13
# Exclude BAD CONTROLS
# V4_8, V4_16B, V8_8
# Exclude COINOC LIKE A CONTROL 
# V6_11

dim(FOUR_VARS_YUCA)
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V3_3"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V1_7"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8_11"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V3_15"), colnames(FOUR_VARS_YUCA),invert =T)]

FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8_6"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_9"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_18"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_6"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8_18"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V3_19"), colnames(FOUR_VARS_YUCA),invert =T)]

FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V6_14"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_19"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_7"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V4_17"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8_13"), colnames(FOUR_VARS_YUCA),invert =T)]

FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V4_8"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V4_16b"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8_8"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V6_11"), colnames(FOUR_VARS_YUCA),invert =T)]

FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_4"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_8"), colnames(FOUR_VARS_YUCA),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3_3"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V1_7"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_11"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3_15"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_6"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_9"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_18"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_6"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_18"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3_19"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V6_14"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_19"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_7"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V4_17"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_13"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V4_8"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V4_16b"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_8"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V6_11"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_4"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_8"), colnames(FOUR_VARS_AMF),invert =T)]


colnames(FOUR_VARS_AMF[,grep("V3", colnames(FOUR_VARS_AMF),invert=F)])

# exclude varieties temporal 
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V1"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V1"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V3"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8"), colnames(FOUR_VARS_YUCA),invert =T)]
# for only one cultivar at time
#FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5"), colnames(FOUR_VARS_YUCA),invert =F)]
#FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5"), colnames(FOUR_VARS_AMF),invert =F)]


####################################
# DESIGN DEFINITION
VAR<-gsub( "m","",gsub( "m2","",gsub("_\\w+$"," ",colnames(FOUR_VARS_YUCA),perl=T)    ) )
TREAT<-gsub("^(\\w+)_","",colnames(FOUR_VARS_YUCA),perl=T)
TREAT_YUCA<-factor(TREAT,levels=c("CTRL","CANB1","CAN","B1"))
#TREAT_YUCA<-factor(TREAT,levels=c("CAN","B1","CANB1"))

design_YUCA <- model.matrix(~TREAT_YUCA*VAR)


#design_YUCA <- model.matrix(~TREAT_YUCA)


VAR<-gsub( "m","",gsub( "m2","",gsub("_\\w+$"," ",colnames(FOUR_VARS_AMF),perl=T)    ) )
TREAT<-gsub("^ (\\w+)_","",colnames(FOUR_VARS_AMF),perl=T)
TREAT_AMF<-factor(TREAT,levels=c("CANB1","CAN","B1"))
#TREAT_AMF<-factor(TREAT,levels=c("CAN","B1","CANB1"))

design_AMF <- model.matrix(~TREAT_AMF*VAR)
#design_AMF <- model.matrix(~TREAT_AMF)


####################################################################
#################################################################################
# DATA FILTERING AND NORMALIZATION
DGE_YUCA <- DGEList(FOUR_VARS_YUCA)
DGE_AMF <- DGEList(FOUR_VARS_AMF)

#FILTERING  
keep_Y <- rowSums(DGE_YUCA$counts>100) >= 3
DGE_YUCA<- DGE_YUCA[keep_Y,]

keep_A <- rowSums(DGE_AMF$counts>100) >= 3
DGE_AMF<- DGE_AMF[keep_A,]

#NORMALIZATION
DGE_YUCA <- calcNormFactors(DGE_YUCA)
DGE_AMF <- calcNormFactors(DGE_AMF)
png("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Voom_Mean-var_CASSAVA.png", width=800, height=600, units="px")
DGE_YUCA_N <- voom(DGE_YUCA, design_YUCA,plot=TRUE)
dev.off()
png("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Voom_Mean-var_AMF.png", width=800, height=600, units="px")
DGE_AMF_N <- voom(DGE_AMF, design_AMF,plot=TRUE)
dev.off()


####################################################################
#DATA VISUALIZATION NORMALIZED
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/MDS_2organism.pdf", width=14, height=10)
par(las=1);par(mfrow=c(2,2));par(mar=c(11,11,4,3),mgp=c(6, 2, 0));options(digits=3)
 
plotMDS(DGE_YUCA_N,col=TREAT_YUCA,pch=16,cex=3,cex.lab=3,cex.axis=3) #by TREAT
#legend("topleft",c('CTRL','CAN','B1'),col=c(1,2,3),pch=16)
plotMDS(DGE_YUCA_N,col=as.numeric(gsub("V","",VAR,perl=TRUE)),pch=16,cex=3,cex.lab=3,cex.axis=3) #by VAR
#legend("topleft",c('V1','V3','V4','V5','V6','V8'),col=c(1,3,4,5,6,8),pch=16)

plotMDS(DGE_AMF_N,col=TREAT_AMF,pch=16,cex=3,cex.lab=3,cex.axis=3) #by TREAT
#legend("topleft",c('CAN','B1'),col=c(1,2),pch=16)
plotMDS(DGE_AMF_N,col=as.numeric(gsub("V","",VAR,perl=TRUE)),pch=16,cex=3,cex.lab=3,cex.axis=3) #by VAR
#legend("topleft",c('V1','V3','V4','V5','V6','V8'),col=c(1,3,4,5,6,8),pch=16)
dev.off()
#
#plotMDS(DGE_AMF_N,col=TREAT_COL,pch=16,main="Rirregularis samples differentiation",cex=3) #by TREAT
#legend("topleft",c('CTRL','CAN','B1'),col=c(1,2,3),pch=16)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/heatmap_CASSAVA_norm_samples.pdf",  width=14, height=10)
sampleDistMatrix <- as.matrix( dist(t(DGE_YUCA_N$E) ))#[,grep(c("V1"), colnames(DGE_YUCA_N$E),invert =F)]
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(sampleDistMatrix, fontsize=18, col=colors)
dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/heatmap_AMF_norm_samples.pdf",  width=14, height=10)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sampleDistMatrix <- as.matrix( dist(t(DGE_AMF_N$E) ))
pheatmap(sampleDistMatrix, fontsize=18, col=colors)
dev.off()
####################################################################


####################################################################
## Differential gene expression


# DGE CASSAVA
DGE_YUCA_N <- voom(DGE_YUCA, design_YUCA,plot=TRUE)
#fit_YUCA <- eBayes(lmFit(DGE_YUCA_N,design_YUCA))
#sign_YUCA<-topTable(fit_YUCA,c(2,3,4,7,8,9,10,11,12),number=dim(DGE_YUCA_N$E)[1]) # all genes
#head(sign_YUCA,5)
dim(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])

YUCA_DGE_vf<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]

YUCA_DGE_vf$gene<-rownames(YUCA_DGE_vf)
YUCA_DGE_vf<-merge(YUCA_DGE_vf,Mercator_Mesculenta,by="gene")
YUCA_DGE_vf[,18]<- unlist(lapply(lapply(strsplit(as.character(YUCA_DGE_vf[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
YUCA_DGE_vf[,18]<-as.factor(YUCA_DGE_vf[,18])
colnames(YUCA_DGE_vf)[c(2:10)]<-c("COL2215_B1.CTRL","COL2215_CAN.CTRL","COL2215_CANB1.CTRL",
                                    "BRA337_B1.CTRL","BRA337_CAN.CTRL","BRA337_CANB1.CTRL",
                                   "CM4574-7_B1.CTRL","CM4574-7_CAN.CTRL","CM4574-7_CANB1.CTRL")


write.table( YUCA_DGE_vf[,c(1:14,16,18,17)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Supplementary_table2.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = F)


# REDO GOOD CONTRAST FOR CASSAVA

REPLICA<-unlist(lapply(strsplit(colnames(FOUR_VARS_YUCA),"_"),function (x) x[[2]]))
VAR<-gsub( " ","",gsub( "m","",gsub( "m2","",gsub("_\\w+$","",colnames(FOUR_VARS_YUCA),perl=T)    ) ))
TREAT<-gsub("^(\\w+)_","",colnames(FOUR_VARS_YUCA),perl=T)
TREAT_YUCA<-factor(TREAT,levels=c("CTRL","CAN","B1","CANB1"))

Treatsss <- factor(paste(VAR,TREAT_YUCA,sep="."))
design <- model.matrix(~0+Treatsss)
colnames(design) <- levels(Treatsss)
corfit <- duplicateCorrelation(DGE_YUCA_N,design)
corfit$consensu
fit_YUCA <- lmFit(DGE_YUCA_N,design,correlation=corfit$consensus)
cm <- makeContrasts(V4_B1.CTRL= V4.B1-V4.CTRL,
                    V4_CAN.CTRL = V4.CAN-V4.CTRL,
                    V4_CANB1.CTRL = V4.CANB1-V4.CTRL,
                    V5_B1.CTRL = V5.B1-V5.CTRL,
                    V5_CAN.CTRL = V5.CAN-V5.CTRL,
                    V5_CANB1.CTRL = V5.CANB1-V5.CTRL,
                    V6_B1.CTRL = V6.B1-V6.CTRL,
                    V6_CAN.CTRL = V6.CAN-V6.CTRL,
                    V6_CANB1.CTRL = V6.CANB1-V6.CTRL,
                    levels=design)
fit_YUCA<- eBayes(contrasts.fit(fit_YUCA, cm))
sign_YUCA<-topTable(fit_YUCA, coef=c(1:9),number=dim(DGE_YUCA_N)[1])

dim(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])

vennDiagram(decideTests(fit_YUCA[,c(1,2,3)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)
vennDiagram(decideTests(fit_YUCA[,c(4,5,6)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)
vennDiagram(decideTests(fit_YUCA[,c(7,8,9)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)


vennDiagram(decideTests(fit_YUCA[,c(1,4,7)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)
vennDiagram(decideTests(fit_YUCA[,c(2,5,8)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)
vennDiagram(decideTests(fit_YUCA[,c(3,6,9)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)


V4<-decideTests(fit_YUCA[,c(1,2,3)],adjust.method="BH")
V5<-decideTests(fit_YUCA[,c(4,5,6)],adjust.method="BH")
V6<-decideTests(fit_YUCA[,c(7,8,9)],adjust.method="BH")
only_canb1<-c(rownames(V4@.Data[V4@.Data[,3]!=0&V4@.Data[,1]==0&V4@.Data[,2]==0,])
,rownames(V5@.Data[V5@.Data[,3]!=0&V5@.Data[,1]==0&V5@.Data[,2]==0,])
,rownames(V6@.Data[V6@.Data[,3]!=0&V6@.Data[,1]==0&V6@.Data[,2]==0,]))


vennCounts(decideTests(fit_YUCA[rownames(fit_YUCA$t) %in% only_canb1,c(3,6,9)],adjust.method="BH"))
vennDiagram(decideTests(fit_YUCA[rownames(fit_YUCA$t) %in% only_canb1,c(3,6,9)],adjust.method="BH"))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.CANB1_venn.pdf",width=8, height=6)

vennCounts(decideTests(fit_YUCA[,c(3,6,9)],adjust.method="BH"))
v <- venneuler::venneuler(c("V6"=31,
                            "V5"=82,
                            "V5&V6"=90, 
                            "V4"=16,
                            "V4&V6"=12, 
                            "V4&V5"=48,
                            "V4&V5&V6"=109))

plot(v)
dev.off()

############################################
# to do analysis of CANB1 comparted to CAN and B1 in R. irregularis
to_analyse<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,1:9]

V4_int<-list()
V5_int<-list()
V6_int<-list()
V4_IC_CAN<-list();V4_IC_B1<-list();V5_IC_CAN<-list();V5_IC_B1<-list()
V6_IC_CAN<-list();V6_IC_B1<-list()
com<-1
NORM_EXP_Y<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  rownames(to_analyse),])
NORM_EXP_Y<-NORM_EXP_Y[grep("CTRL",rownames(NORM_EXP_Y),invert=T),]
tratamiento<-unlist(lapply(strsplit(rownames(NORM_EXP_Y),"_"),function (x) x[3]))
cultivar<-gsub("m","",gsub("m2","",unlist(lapply(strsplit(rownames(NORM_EXP_Y),"_"),function (x) x[1]))   ))

for (gene in rownames(to_analyse)) {
  compt<- gene
  boxplot(to_analyse[grep(compt,rownames(to_analyse)),],las=2)
  
  
  Results_A3<-cbind.data.frame(NORM_EXP_Y[,grep(compt,colnames(NORM_EXP_Y))],tratamiento,cultivar)
  #tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)
  V4_int[[compt]]<-findInterval(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[3,1], 
                                sort(c(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[1,1],
                                       tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[2,1])))==1
  
  V5_int[[compt]]<-findInterval(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[3,2], 
                                sort(c(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[1,2],
                                       tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[2,2])))==1
  
  V6_int[[compt]]<-findInterval(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[3,3], 
                                sort(c(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[1,3],
                                       tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[2,3])))==1
  
  
  V4_IC_CAN[[com]]<-findInterval(to_analyse[compt,c(3)], sort(c(to_analyse[compt,c(2)]*1.2,to_analyse[compt,c(2)]*0.8) )) == 1
  V4_IC_B1[[com]]<-findInterval(to_analyse[compt,c(3)], sort(c(to_analyse[compt,c(1)]*1.2,to_analyse[compt,c(1)]*0.8) )) == 1
  V5_IC_CAN[[com]]<-findInterval(to_analyse[compt,c(6)], sort(c(to_analyse[compt,c(5)]*1.2,to_analyse[compt,c(5)]*0.8) )) == 1
  V5_IC_B1[[com]]<-findInterval(to_analyse[compt,c(6)], sort(c(to_analyse[compt,c(4)]*1.2,to_analyse[compt,c(4)]*0.8) )) == 1
  V6_IC_CAN[[com]]<-findInterval(to_analyse[compt,c(9)], sort(c(to_analyse[compt,c(8)]*1.2,to_analyse[compt,c(8)]*0.8) )) == 1
  V6_IC_B1[[com]]<-findInterval(to_analyse[compt,c(9)], sort(c(to_analyse[compt,c(7)]*1.2,to_analyse[compt,c(7)]*0.8) )) == 1
  
  com<-com+1
}

Results_A<-cbind.data.frame(unlist(V4_int),unlist(V5_int),unlist(V6_int))
rownames(Results_A)<-rownames(to_analyse)
Results_A[Results_A==TRUE]<-"Antagonistic"
Results_A[Results_A==FALSE]<-"Synergistic"
colnames(Results_A)<-c("V4","V5","V6")

Results_A2<-cbind.data.frame(unlist(V4_IC_CAN),unlist(V4_IC_B1),
                             unlist(V5_IC_CAN),unlist(V5_IC_B1),
                             unlist(V6_IC_CAN),unlist(V6_IC_B1))
rownames(Results_A2)<-rownames(to_analyse)
colnames(Results_A2)<-c("V4_can","V4_b1","V5_can","V5_b1","V6_can","V6_b1")

Results_A2$V4_can[Results_A2[,1]]<-"Equal 2 CAN"
Results_A2$V4_can[Results_A2[,1]==FALSE]<-"Different 2 CAN"
Results_A2$V4_b1[Results_A2[,2]]<-"Equal 2 B1"
Results_A2$V4_b1[Results_A2[,2]==FALSE]<-"Different 2 B1"
Results_A2$V5_can[Results_A2[,3]]<-"Equal 2 CAN"
Results_A2$V5_can[Results_A2[,3]==FALSE]<-"Different 2 CAN"
Results_A2$V5_b1[Results_A2[,4]]<-"Equal 2 B1"
Results_A2$V5_b1[Results_A2[,4]==FALSE]<-"Different 2 B1"
Results_A2$V6_can[Results_A2[,5]]<-"Equal 2 CAN"
Results_A2$V6_can[Results_A2[,5]==FALSE]<-"Different 2 CAN"
Results_A2$V6_b1[Results_A2[,6]]<-"Equal 2 B1"
Results_A2$V6_b1[Results_A2[,6]==FALSE]<-"Different 2 B1"


Results_A_vf<-cbind.data.frame(Results_A,Results_A2)
Results_A_vf
Results_A_vf[Results_A_vf$V4=="Equal 2 CAN"&Results_A_vf$V4_can=="Equal 2 CAN"]
Results_A_vf[Results_A_vf$V4=="Equal 2 CAN"&Results_A_vf$V4_can=="Equal 2 CAN"]

Results_A_vf$V4_can=="Equal 2 CAN"

head(Results_A_vf)
dim(Results_A_vf)
Results_A_vf2<-cbind.data.frame(rep("V4",dim(Results_A_vf)[1]),rep("V5",dim(Results_A_vf)[1]),rep("V6",dim(Results_A_vf)[1]))
colnames(Results_A_vf2)<-c("V4","V5","V6")
rownames(Results_A_vf2)<-rownames(Results_A_vf)
Results_A_vf2$V4<-as.vector(Results_A_vf2$V4);Results_A_vf2$V5<-as.vector(Results_A_vf2$V5);Results_A_vf2$V6<-as.vector(Results_A_vf2$V6)


Results_A_vf2$V4[Results_A_vf$V4=="Synergistic"]<-"Synergistic"
Results_A_vf2$V4[Results_A_vf$V4=="Antagonistic"]<-"Antagonistic"
Results_A_vf2$V4[Results_A_vf$V4_can=="Equal 2 CAN"]<-"Equal 2 CAN"
Results_A_vf2$V4[Results_A_vf$V4_b1=="Equal 2 B1"]<-"Equal 2 B1"
Results_A_vf2$V4[Results_A_vf$V4_can=="Equal 2 CAN"&Results_A_vf$V4_b1=="Equal 2 B1"]<-"Equal 2 Both"


Results_A_vf2$V5[Results_A_vf$V5=="Synergistic"]<-"Synergistic"
Results_A_vf2$V5[Results_A_vf$V5=="Antagonistic"]<-"Antagonistic"
Results_A_vf2$V5[Results_A_vf$V5_can=="Equal 2 CAN"]<-"Equal 2 CAN"
Results_A_vf2$V5[Results_A_vf$V5_b1=="Equal 2 B1"]<-"Equal 2 B1"
Results_A_vf2$V5[Results_A_vf$V5_can=="Equal 2 CAN"&Results_A_vf$V5_b1=="Equal 2 B1"]<-"Equal 2 Both"

Results_A_vf2$V6[Results_A_vf$V6=="Synergistic"]<-"Synergistic"
Results_A_vf2$V6[Results_A_vf$V6=="Antagonistic"]<-"Antagonistic"
Results_A_vf2$V6[Results_A_vf$V6_can=="Equal 2 CAN"]<-"Equal 2 CAN"
Results_A_vf2$V6[Results_A_vf$V6_b1=="Equal 2 B1"]<-"Equal 2 B1"
Results_A_vf2$V6[Results_A_vf$V6_can=="Equal 2 CAN"&Results_A_vf$V6_b1=="Equal 2 B1"]<-"Equal 2 Both"

Results_A_vf2
head(Results_A_vf2)
Results_A_vf2<-cbind.data.frame(Results_A_vf2,to_analyse)

# See type of interaction

compt<-"Manes.03G123000"
boxplot(NORM_EXP_Y[,grep(compt,colnames(NORM_EXP_Y))]~interaction(tratamiento,cultivar),las=2)
boxplot(to_analyse[grep(compt,rownames(to_analyse)),c(3,2,1,6,5,4,9,8,7)],las=2)
Results_A_vf2[grep(compt,rownames(Results_A2)),]

table(Results_A_vf2[,c(3)])/sum(table(Results_A_vf2[,c(3)]))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.Categories_CANB1.pdf",width=6, height=8)
par(mfrow=c(1,1),mar=c(4,8,1,1),mgp=c(5, 2, 0))
barplot(matrix(c(
  table(Results_A_vf2[,c(1)]),
  table(Results_A_vf2[,c(2)]),
  table(Results_A_vf2[,c(3)])),nrow=5),beside=F,las=1,ylim=c(0,1200),
  col=c("brown","skyblue3","wheat3","azure3","chartreuse4"),names.arg = c("V4","V5","V6"),cex.axis = 3,cex.names = 3)
chisq.test(matrix(c(
  table(Results_A_vf2[,c(1)]),
  table(Results_A_vf2[,c(2)]),
  table(Results_A_vf2[,c(3)])),nrow=5)[c(4),])
dev.off()

# Info about the synergistic/antagonistic effects
Results_A_vf2$gene<- rownames(Results_A_vf2)  

Results_A_vf3<-merge(Results_A_vf2,Mercator_Mesculenta,by="gene")
head(Results_A_vf3)
unlist(lapply(lapply(strsplit(as.character(genes_yuca$Bin),"\\."),function (x) x[1:2]),function (x) paste(x[1],x[2],sep=".")))
Results_A_vf3<-cbind.data.frame(Results_A_vf3[,c(1:4,6)],unlist(lapply(lapply(strsplit(as.character(Results_A_vf3$Bin),"\\."),function (x) x[1:2]),function (x) paste(x[1],x[2],sep=".")))  )
colnames(Results_A_vf3)[dim(Results_A_vf3)[2]]<-"classifier"


transport.metal


sort(table(Results_A_vf3[Results_A_vf3[,2]=="Synergistic"&Results_A_vf3[,3]!="Synergistic"&Results_A_vf3[,4]!="Synergistic",5]))
Results_A_vf3[Results_A_vf3[,2]=="Synergistic"&Results_A_vf3[,3]!="Synergistic"&Results_A_vf3[,4]!="Synergistic"&
                Results_A_vf3[,5]=="misc.short_chain_dehydrogenase/reductase_(SDR)",]

Results_A_vf4<-Results_A_vf3[Results_A_vf3[,2]=="Synergistic"&Results_A_vf3[,3]!="Synergistic"&Results_A_vf3[,4]!="Synergistic",]
dim(Results_A_vf4[!duplicated(Results_A_vf4[,1]),])
sort(table(Results_A_vf4[,5]))


####
### NOt equal to CANB1
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.Interact_VENN_different_CANB1.pdf",width=8, height=6)

vennCounts(
  Results_A_vf2=="Antagonistic")
v <- venneuler::venneuler(c(V4=38,
                            V5=153,
                            V6=522, 
                            "V4&V5"=7,
                            "V4&V6"=19, "V5&V6"=149,
                            "V4&V5&V6"=6))

plot(v)
dev.off()

# send results
Results_A_vf2$gene<-rownames(Results_A_vf2)
Results_A_vf2<-merge(Results_A_vf2,Mercator_Mesculenta,by="gene")
Results_A_vf2[,17]<- unlist(lapply(lapply(strsplit(as.character(Results_A_vf2[,15]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
Results_A_vf2[,17]<-as.factor(Results_A_vf2[,17])
colnames(Results_A_vf2)[c(2,3,4,17)]<-c("COL2215","BRA337","CM4574-7","Type") 
write.table( Results_A_vf2,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Supplementary_table3.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

head(Results_A_vf2)

## Specific genes per variety
V4_effect<-Results_A_vf2[grep("Equal",Results_A_vf2[,2],invert=T),]
V4_effect<-V4_effect[V4_effect$Bin!="35.2"&V4_effect$Bin!="35.1",]
V4_effect$gene
V5_effect<-Results_A_vf2[grep("Equal",Results_A_vf2[,3],invert=T),]
V5_effect<-V5_effect[V5_effect$Bin!="35.2"&V5_effect$Bin!="35.1",]
V5_effect$gene
V6_effect<-Results_A_vf2[grep("Equal",Results_A_vf2[,4],invert=T),]
V6_effect<-V6_effect[V6_effect$Bin!="35.2"&V6_effect$Bin!="35.1",]
V6_effect$gene


V4_top<-table(droplevels(Results_A_vf2[Results_A_vf2$gene %in% setdiff(setdiff(V4_effect$gene,V5_effect$gene),V6_effect$gene),8]))

V5_top<-table(droplevels(Results_A_vf2[Results_A_vf2$gene %in% setdiff(setdiff(V5_effect$gene,V4_effect$gene),V6_effect$gene),8]))

V6_top<-table(droplevels(Results_A_vf2[Results_A_vf2$gene %in% setdiff(setdiff(V6_effect$gene,V4_effect$gene),V5_effect$gene),8]))
length(V4_top)

pheatmap(t(V4_top),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
pheatmap(t(V5_top),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
pheatmap(t(V6_top),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.Exclusive_Syner_ANTA_genes_var.pdf",width=10, height=8)
pheatmap(t(V4_top),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
pheatmap(t(V5_top),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
pheatmap(t(V6_top),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
write.table( V4_effect,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.Exclusive_Syner_ANTA_genes_VAR4.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)
write.table( V5_effect,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.Exclusive_Syner_ANTA_genes_VAR5.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)
write.table( V6_effect,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.Exclusive_Syner_ANTA_genes_VAR6.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

dev.off()
## Common genes among variety
V_common<-Results_A_vf2[grep("Equal",Results_A_vf2[,2],invert=T),]
V_common<-V_common[grep("Equal",V_common[,3],invert=T),]
V_common<-V_common[grep("Equal",V_common[,4],invert=T),]


V_common<-Results_A_vf2[Results_A_vf2[,2]=="Antagonistic",]
V_common<-Results_A_vf2[grep("Equal",Results_A_vf2[,2],invert=T),]
V_common<-V_common[grep("Equal",V_common[,3],invert=T),]

dim(Results_A_vf2[,1:4])
head(V_common)
head(Results_A_vf2[,1:4])


## Antagnoistic , Synergistic and both in 3 cultures

A_S<-subset(Results_A_vf2, COL2215 =="Antagonistic" |COL2215 =="Synergistic" |
         BRA337 =="Antagonistic" |BRA337 =="Synergistic" |
         Results_A_vf2[,4]=="Antagonistic" |Results_A_vf2[,4] =="Synergistic")

A_<-subset(A_S, COL2215 =="Antagonistic"  |
              BRA337 =="Antagonistic"  |
             A_S[,4]=="Antagonistic" )
S_<-subset(A_S, COL2215 =="Synergistic"  |
             BRA337 =="Synergistic"  |
             A_S[,4]=="Synergistic" )

AxS<-subset(A_S, COL2215 =="Synergistic"  & BRA337 =="Antagonistic" |
              COL2215 =="Synergistic"  & A_S[,4] =="Antagonistic" |
              BRA337 =="Synergistic"  & A_S[,4] =="Antagonistic" |
              
              COL2215 =="Antagonistic"  & BRA337 =="Synergistic" |
              COL2215 =="Antagonistic"  & A_S[,4] =="Synergistic" |
              COL2215 =="Antagonistic"  & A_S[,4] =="Synergistic" )


dim(A_pur)
dim(A_S)

A_pur<- A_[!A_$gene %in% AxS$gene,]

S_pur<- S_[!S_$gene %in% AxS$gene,]


length(unique(A_S$gene))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.AntagonistXSynergistics_genes.pdf",width=12, height=8)
pheatmap(t(table(droplevels(AxS[,17]))),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
write.table( AxS,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.AntagonistXSynergistics_genes.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)
length(unique(AxS$gene))
table(droplevels(AxS[,17][!AxS[,17] %in% c(as.vector(S_pur[,17]),as.vector(A_pur[,17]))]))

dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.ONLY_Antagonistic_genes.pdf",width=12, height=8)
pheatmap(t(table(droplevels(A_pur[,17]))),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
write.table( A_pur,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.ONLY_Antagonistic_genes.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)
length(unique(A_pur$gene))
table(droplevels(A_pur[,17][!A_pur[,17] %in% c(as.vector(S_pur[,17]),as.vector(AxS[,17]))]))
A_pur[A_pur$Type=="transport.metabolite_transporters_at_the_envelope_membrane",]

dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.ONLY_Synergistic_genes.pdf",width=12, height=8)
pheatmap(t(table(droplevels(S_pur[,17]))),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
write.table( S_pur,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.ONLY_Synergistic_genes.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)
length(unique(S_pur$gene))
table(droplevels(S_[,17][!S_[,17] %in% c(as.vector(A_[,17]),as.vector(AxS[,17]))]))
head(S_pur,2)
g7584.t1
boxplot(S_pur[S_pur$Type=="cell_wall.precursor_synthesis",5:13],las=2)
S_pur[S_pur$Type=="transport.unspecified_cations",]

dev.off()

########################

########################
########################
# DGE AMF
DGE_AMF_N <- voom(DGE_AMF, design_AMF,plot=TRUE)

#fit_AMF <- eBayes(lmFit(DGE_AMF_N,design_AMF))
#sign_AMF<-topTable(fit_AMF,c(2,3,6,7,8,9),number=dim(DGE_AMF_N$E)[1]) # all genes
head(sign_AMF,50)
dim(sign_AMF[sign_AMF$adj.P.Val<0.05,])

sign_AMF2<- decideTests(fit_AMF[,c(2,3,6,7,8,9)],adjust.method="BH")

REPLICA<-unlist(lapply(strsplit(colnames(FOUR_VARS_AMF),"_"),function (x) x[[2]]))
VAR<-gsub( " ","",gsub( "m","",gsub( "m2","",gsub("_\\w+$","",colnames(FOUR_VARS_AMF),perl=T)    ) ))
TREAT<-gsub("^ (\\w+)_","",colnames(FOUR_VARS_AMF),perl=T)
TREAT_AMF<-factor(TREAT,levels=c("CANB1","CAN","B1"))

Treatsss <- factor(paste(VAR,TREAT_AMF,sep="."))
design <- model.matrix(~0+Treatsss)
colnames(design) <- levels(Treatsss)
corfit <- duplicateCorrelation(DGE_AMF_N,design)
corfit$consensu
fit_AMF <- lmFit(DGE_AMF_N,design,correlation=corfit$consensus)
cm <- makeContrasts(V4_B1.CANB1= V4.B1-V4.CANB1,
            V4_CAN.CANB1 = V4.CAN-V4.CANB1,
            V5_B1.CANB1 = V5.B1-V5.CANB1,
            V5_CAN.CANB1 = V5.CAN-V5.CANB1,
            V6_B1.CANB1 = V6.B1-V6.CANB1,
            V6_CAN.CANB1 = V6.CAN-V6.CANB1,
            levels=design)
fit_AMF<- eBayes(contrasts.fit(fit_AMF, cm))
sign_AMF<-topTable(fit_AMF, coef=c(1:6),number=dim(DGE_AMF_N)[1])
tail(sign_AMF)
dim(sign_AMF[sign_AMF$adj.P.Val<0.05,])


AMF_DGE_vf<-sign_AMF[sign_AMF$adj.P.Val<0.05,]

AMF_DGE_vf$gene<-rownames(AMF_DGE_vf)
AMF_DGE_vf<-merge(AMF_DGE_vf,Mercator_Rirregularis,by="gene")
AMF_DGE_vf[,15]<- unlist(lapply(lapply(strsplit(as.character(AMF_DGE_vf[,13]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
AMF_DGE_vf[,15]<-as.factor(AMF_DGE_vf[,15])
colnames(AMF_DGE_vf)[c(2:7,15)]<-c("COL2215_B1.CANB1","COL2215_CAN.CANB1","BRA337_B1.CANB1","BRA337_CAN.CANB1",
                                   "CM4574-7_B1.CANB1","CM4574-7_CAN.CANB1","Type")
write.table( AMF_DGE_vf[,c(1:11,13,15,14)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Supplementary_table4.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = F)




vennDiagram(decideTests(fit_AMF[,c(1,2)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)
vennDiagram(decideTests(fit_AMF[,c(3,4)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)
vennDiagram(decideTests(fit_AMF[,c(5,6)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)


vennDiagram(decideTests(fit_AMF[,c(1,3,5)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)
vennDiagram(decideTests(fit_AMF[,c(2,4,6)],adjust.method="BH"),include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE Rirregularis",circle.col=1)

common_var<-decideTests(fit_AMF[,c(2,4,6)],adjust.method="BH")
common_var<-rownames(common_var@.Data[common_var@.Data[,1]!=0&common_var@.Data[,2]!=0&common_var@.Data[,3]!=0,])
common_var2<-decideTests(fit_AMF[,c(1,3,5)],adjust.method="BH")
common_var2<-rownames(common_var2@.Data[common_var2@.Data[,1]!=0&common_var2@.Data[,2]!=0,])
unique(c(common_var,common_var2))

V4<-decideTests(fit_AMF[,c(1,2)],adjust.method="BH")
V5<-decideTests(fit_AMF[,c(3,4)],adjust.method="BH")
V6<-decideTests(fit_AMF[,c(5,6)],adjust.method="BH")
rownames(V4@.Data[V4@.Data[,1]!=0|V4@.Data[,2]!=0,])
rownames(V5@.Data[V5@.Data[,1]!=0|V5@.Data[,2]!=0,])
rownames(V6@.Data[V6@.Data[,1]!=0|V6@.Data[,2]!=0,])

v <- venneuler::venneuler(c("CAN"=581,
                            "B1"=1,
                            "B1&CAN"=0))


v <- venneuler::venneuler(c("V6"=414,
                            "V5"=23,
                            "V5&V6"=90, 
                            "V4"=139,
                            "V4&V6"=53, 
                            "V4&V5"=71,
                            "V4&V5&V6"=24))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.CANB1-onlyB1_venn.pdf",width=8, height=6)

plot(v)
dev.off()


write.table(sign_AMF,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.DEG_AMF_ALL.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = F)





toPHEAT<-DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% rownames(sign_AMF[sign_AMF$adj.P.Val<0.001,]),]
pheatmap(toPHEAT)

# to do analysis of CANB1 comparted to CAN and B1 in R. irregularis
to_analyse<-sign_AMF[sign_AMF$adj.P.Val<0.05,1:6]

V4_int<-list()
V5_int<-list()
V6_int<-list()
V4_IC_CAN<-list();V4_IC_B1<-list();V5_IC_CAN<-list();V5_IC_B1<-list()
V6_IC_CAN<-list();V6_IC_B1<-list()
V4_high<-list();V5_high<-list();V6_high<-list()
com<-1
NORM_EXP_A<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  rownames(to_analyse),])
NORM_EXP_A<-NORM_EXP_A[grep("CTRL",rownames(NORM_EXP_A),invert=T),]
tratamiento<-unlist(lapply(strsplit(rownames(NORM_EXP_A),"_"),function (x) x[3]))
cultivar<-gsub("m","",gsub("m2","",unlist(lapply(strsplit(rownames(NORM_EXP_A),"_"),function (x) x[1]))   ))

for (gene in rownames(to_analyse)) {
  compt<- gene
  pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/test.pdf",width=6, height=8)
  
  boxplot(to_analyse[grep(compt,rownames(to_analyse)),],las=2)
  dev.off()
  
    Results_A3<-cbind.data.frame(NORM_EXP_A[,grep(compt,colnames(NORM_EXP_A))],tratamiento,cultivar)
  #tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)
  V4_int[[compt]]<-findInterval(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[3,1], 
               sort(c(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[1,1],
                      tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[2,1])))==1
  
  V4_high[[compt]]<-tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[3,1] > 
    max(c(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[1,1],
           tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[2,1]))*1.1
  
  
  V5_int[[compt]]<-findInterval(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[3,2], 
               sort(c(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[1,2],
                      tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[2,2])))==1
  
  V5_high[[compt]]<-tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[3,2] > 
    max(c(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[1,2],
          tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[2,2]))*1.1
  
  
  V6_int[[compt]]<-findInterval(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[3,3], 
               sort(c(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[1,3],
                      tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[2,3])))==1
  
  V6_high[[compt]]<-tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[3,3] > 
    max(c(tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[1,3],
          tapply(Results_A3[,1],list(Results_A3[,2],Results_A3[,3]),mean)[2,3]))*1.1
  
  
  
  V4_IC_CAN[[com]]<- to_analyse[compt,c(1)] > 0.2 | to_analyse[compt,c(1)] < -0.2 
  V4_IC_B1[[com]]<-to_analyse[compt,c(2)] > 0.2 | to_analyse[compt,c(2)] < -0.2 
  V5_IC_CAN[[com]]<- to_analyse[compt,c(3)] > 0.2 | to_analyse[compt,c(3)] < -0.2 
  V5_IC_B1[[com]]<-to_analyse[compt,c(4)] > 0.2 | to_analyse[compt,c(4)] < -0.2 
  V6_IC_CAN[[com]]<-to_analyse[compt,c(5)] > 0.2 | to_analyse[compt,c(5)] < -0.2 
  V6_IC_B1[[com]]<-to_analyse[compt,c(6)] > 0.2 | to_analyse[compt,c(6)] < -0.2 
  
      com<-com+1
}

Results_A<-cbind.data.frame(unlist(V4_int),unlist(V5_int),unlist(V6_int))
rownames(Results_A)<-rownames(to_analyse)
Results_A[Results_A==TRUE]<-"Antagonistic"
Results_A[Results_A==FALSE]<-"Synergistic"

colnames(Results_A)<-c("V4","V5","V6")

Results_A[Results_A$V4=="Synergistic"&unlist(V4_high)==TRUE,]<-"HIGHER"
Results_A[Results_A$V4=="Synergistic"&unlist(V4_high)==FALSE,]<-"LOWER"
Results_A[Results_A$V5=="Synergistic"&unlist(V5_high)==TRUE,]<-"HIGHER"
Results_A[Results_A$V5=="Synergistic"&unlist(V5_high)==FALSE,]<-"LOWER"
Results_A[Results_A$V6=="Synergistic"&unlist(V6_high)==TRUE,]<-"HIGHER"
Results_A[Results_A$V6=="Synergistic"&unlist(V6_high)==FALSE,]<-"LOWER"

Results_A2<-cbind.data.frame(unlist(V4_IC_CAN),unlist(V4_IC_B1),
                             unlist(V5_IC_CAN),unlist(V5_IC_B1),
                             unlist(V6_IC_CAN),unlist(V6_IC_B1))
rownames(Results_A2)<-rownames(to_analyse)
colnames(Results_A2)<-c("V4_b1","V4_can","V5_b1","V5_can","V6_b1","V6_can")

Results_A2$V4_b1[Results_A2[,1]]<-"Different 2 B1"
Results_A2$V4_b1[Results_A2[,1]==FALSE]<-"Equal 2 B1"
Results_A2$V4_can[Results_A2[,2]]<-"Different 2 CAN"
Results_A2$V4_can[Results_A2[,2]==FALSE]<-"Equal 2 CAN"
Results_A2$V5_b1[Results_A2[,3]]<-"Different 2 B1"
Results_A2$V5_b1[Results_A2[,3]==FALSE]<-"Equal 2 B1"
Results_A2$V5_can[Results_A2[,4]]<-"Different 2 CAN"
Results_A2$V5_can[Results_A2[,4]==FALSE]<-"Equal 2 CAN"
Results_A2$V6_b1[Results_A2[,5]]<-"Different 2 B1"
Results_A2$V6_b1[Results_A2[,5]==FALSE]<-"Equal 2 B1"
Results_A2$V6_can[Results_A2[,6]]<-"Different 2 CAN"
Results_A2$V6_can[Results_A2[,6]==FALSE]<-"Equal 2 CAN"


Results_A_vf<-cbind.data.frame(Results_A,Results_A2)



head(Results_A_vf2)
dim(Results_A_vf)
Results_A_vf2<-cbind.data.frame(rep("V4",dim(Results_A_vf)[1]),rep("V5",dim(Results_A_vf)[1]),rep("V6",dim(Results_A_vf)[1]))
colnames(Results_A_vf2)<-c("V4","V5","V6")
rownames(Results_A_vf2)<-rownames(Results_A_vf)
Results_A_vf2$V4<-as.vector(Results_A_vf2$V4);Results_A_vf2$V5<-as.vector(Results_A_vf2$V5);Results_A_vf2$V6<-as.vector(Results_A_vf2$V6)


Results_A_vf2$V4[Results_A_vf$V4=="Synergistic"]<-"Synergistic"
Results_A_vf2$V4[Results_A_vf$V4=="HIGHER"]<-"HIGHER"
Results_A_vf2$V4[Results_A_vf$V4=="LOWER"]<-"LOWER"
Results_A_vf2$V4[Results_A_vf$V4=="Antagonistic"]<-"Antagonistic"
Results_A_vf2$V4[Results_A_vf$V4_can=="Equal 2 CAN"]<-"Equal 2 CAN"
Results_A_vf2$V4[Results_A_vf$V4_b1=="Equal 2 B1"]<-"Equal 2 B1"
Results_A_vf2$V4[Results_A_vf$V4_can=="Equal 2 CAN"&Results_A_vf$V4_b1=="Equal 2 B1"]<-"Equal 2 BOTH"




Results_A_vf2$V5[Results_A_vf$V5=="Synergistic"]<-"Synergistic"
Results_A_vf2$V5[Results_A_vf$V5=="HIGHER"]<-"HIGHER"
Results_A_vf2$V5[Results_A_vf$V5=="LOWER"]<-"LOWER"
Results_A_vf2$V5[Results_A_vf$V5=="Antagonistic"]<-"Antagonistic"
Results_A_vf2$V5[Results_A_vf$V5_can=="Equal 2 CAN"]<-"Equal 2 CAN"
Results_A_vf2$V5[Results_A_vf$V5_b1=="Equal 2 B1"]<-"Equal 2 B1"
Results_A_vf2$V5[Results_A_vf$V5_can=="Equal 2 CAN"&Results_A_vf$V5_b1=="Equal 2 B1"]<-"Equal 2 BOTH"



Results_A_vf2$V6[Results_A_vf$V6=="Synergistic"]<-"Synergistic"
Results_A_vf2$V6[Results_A_vf$V6=="HIGHER"]<-"HIGHER"
Results_A_vf2$V6[Results_A_vf$V6=="LOWER"]<-"LOWER"
Results_A_vf2$V6[Results_A_vf$V6=="Antagonistic"]<-"Antagonistic"
Results_A_vf2$V6[Results_A_vf$V6_can=="Equal 2 CAN"]<-"Equal 2 CAN"
Results_A_vf2$V6[Results_A_vf$V6_b1=="Equal 2 B1"]<-"Equal 2 B1"
Results_A_vf2$V6[Results_A_vf$V6_can=="Equal 2 CAN"&Results_A_vf$V6_b1=="Equal 2 B1"]<-"Equal 2 BOTH"



Results_A_vf2
table(Results_A_vf2[,3])
Results_A_vf2<-cbind.data.frame(Results_A_vf2,to_analyse)


tratamiento<-unlist(lapply(strsplit(colnames(DGE_AMF$counts),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(DGE_AMF$counts),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(7:9,1:3,4:6)]

boxplot(DGE_AMF_N$E[grep("g3179.t1",rownames(DGE_AMF_N$E)),]~interacts,las=2)
boxplot(DGE_AMF$counts[grep("g3179.t1",rownames(DGE_AMF$counts)),]~interacts,las=2)
Results_A_vf2[rownames(Results_A_vf2)=="g3179.t1",]

breaks<-c(0,10,50,100,200,300,400,500,700)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(9)

Results_A[rownames(Results_A)=="g6073.t1",]
boxplot(to_analyse[rownames(to_analyse)=="g6073.t1",],las=2)

Results_A

table(Results_A_vf2[,c(3)])/sum(table(Results_A_vf2[,c(3)]))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.Categories_CANB1.pdf",width=6, height=8)
par(mfrow=c(1,1),mar=c(4,8,1,1),mgp=c(5, 2, 0))
barplot(matrix(c(
table(Results_A_vf2[,c(1)]),
table(Results_A_vf2[,c(2)]),
table(Results_A_vf2[,c(3)])),nrow=6),beside=F,las=1,ylim=c(0,1000),
col=c("khaki4","skyblue3","wheat3","azure3","chartreuse4","indianred4"),names.arg = c("COL2215","BRA337","CM4574-7"),cex.axis = 2,cex.names = 1.5)
chisq.test(matrix(c(
  table(Results_A_vf2[,c(1)]),
  table(Results_A_vf2[,c(2)]),
  table(Results_A_vf2[,c(3)])),nrow=6)[c(4),])
dev.off()


### NOt equal to CANB1
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.Interact_VENN_different_CANB1.pdf",width=8, height=6)

vennCounts(
Results_A_vf2=="Synergistic")
v <- venneuler::venneuler(c(V4=35,
                            V5=259,
                            V6=241, 
                            "V4&V5"=35,
                            "V4&V6"=22, "V5&V6"=157,
                            "V4&V5&V6"=15))

plot(v)
dev.off()

# send results
Results_A_vf2$gene<-rownames(Results_A_vf2)
Results_A_vf2<-merge(Results_A_vf2,Mercator_Rirregularis,by="gene")
Results_A_vf2[,14]<- unlist(lapply(lapply(strsplit(as.character(Results_A_vf2[,12]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
Results_A_vf2[,14]<-as.factor(Results_A_vf2[,14])

colnames(Results_A_vf2)[c(2,3,4,14)]<-c("COL2215","BRA337","CM4574-7","Type") 
write.table( Results_A_vf2,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Supplementary Table5.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = F)

head(Results_A_vf2,2)


V4_effect<-Results_A_vf2[grep("Equal",Results_A_vf2[,2],invert=T),]
V4_effect<-V4_effect[V4_effect$Bin!="35.2"&V4_effect$Bin!="35.1",]
V4_effect$gene
V5_effect<-Results_A_vf2[grep("Equal",Results_A_vf2[,3],invert=T),]
V5_effect<-V5_effect[V5_effect$Bin!="35.2"&V5_effect$Bin!="35.1",]
V5_effect$gene
V6_effect<-Results_A_vf2[grep("Equal",Results_A_vf2[,4],invert=T),]
V6_effect<-V6_effect[V6_effect$Bin!="35.2"&V6_effect$Bin!="35.1",]
V6_effect$gene


V4_top<-table(droplevels(Results_A_vf2[Results_A_vf2$gene %in% setdiff(setdiff(V4_effect$gene,V5_effect$gene),V6_effect$gene),8]))

V5_top<-table(droplevels(Results_A_vf2[Results_A_vf2$gene %in% setdiff(setdiff(V5_effect$gene,V4_effect$gene),V6_effect$gene),8]))

V6_top<-table(droplevels(Results_A_vf2[Results_A_vf2$gene %in% setdiff(setdiff(V6_effect$gene,V4_effect$gene),V5_effect$gene),8]))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.Exclusive_Syner_ANTA_genes_var.pdf",width=10, height=8)
pheatmap(t(V4_top),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
pheatmap(t(V5_top),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
pheatmap(t(V6_top),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
dev.off()


## Antagnoistic , Synergistic and both in 3 cultures

A_S<-subset(Results_A_vf2, COL2215 =="Antagonistic" |COL2215 =="Synergistic" |
              BRA337 =="Antagonistic" |BRA337 =="Synergistic" |
              Results_A_vf2[,4]=="Antagonistic" |Results_A_vf2[,4] =="Synergistic")

A_<-subset(A_S, COL2215 =="Antagonistic"  |
             BRA337 =="Antagonistic"  |
             A_S[,4]=="Antagonistic" )
S_<-subset(A_S, COL2215 =="Synergistic"  |
             BRA337 =="Synergistic"  |
             A_S[,4]=="Synergistic" )

AxS<-subset(A_S, COL2215 =="Synergistic"  & BRA337 =="Antagonistic" |
              COL2215 =="Synergistic"  & A_S[,4] =="Antagonistic" |
              BRA337 =="Synergistic"  & A_S[,4] =="Antagonistic" |
              
              COL2215 =="Antagonistic"  & BRA337 =="Synergistic" |
              COL2215 =="Antagonistic"  & A_S[,4] =="Synergistic" |
              COL2215 =="Antagonistic"  & A_S[,4] =="Synergistic" )
dim(  A_S)
dim(S_pur)
dim(A_pur)
dim(AxS)

A_pur<- A_[!A_$gene %in% AxS$gene,]

S_pur<- S_[!S_$gene %in% AxS$gene,]

t(table(droplevels(A_pur[,8])))

length(unique(A_S$gene))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.AntagonistXSynergistics_genes.pdf",width=12, height=8)
pheatmap(t(table(droplevels(AxS[,14]))),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
write.table( AxS,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.AntagonistXSynergistics_genes.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)
length(unique(AxS$gene))
table(droplevels(AxS[,14][!AxS[,14] %in% c(as.vector(S_pur[,14]),as.vector(A_pur[,14]))]))

dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.ONLY_Antagonistic_genes.pdf",width=12, height=8)
pheatmap(t(table(droplevels(A_pur[,14]))),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
write.table( A_pur,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.ONLY_Antagonistic_genes.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)
length(unique(A_pur$gene))
table(droplevels(A_pur[,14][!A_pur[,14] %in% c(as.vector(S_pur[,14]),as.vector(AxS[,14]))]))

dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.ONLY_Synergistic_genes.pdf",width=12, height=8)
pheatmap(t(table(droplevels(S_pur[,14]))),cluster_rows = F,cluster_cols = F,cellwidth =8 ,cellheight = 8)
write.table( S_pur,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.ONLY_Synergistic_genes.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)
length(unique(S_pur$gene))
table(droplevels(S_pur[,14][!S_pur[,14] %in% c(as.vector(A_pur[,14]),as.vector(AxS[,14]))]))
head(S_pur,10)
g7584.t1
boxplot(S_pur[S_pur$Function=="transport.sugars",5:10],las=2)
S_pur[grep("acid_and_other_phosphatases",S_pur$Function),]
dev.off()

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################




########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

################################################################
########################################################################
## almost not expression in single-inoculation (DIRECT INTERACTION)
genes_amf<-sign_AMF[sign_AMF$adj.P.Val<0.05,]
genes_amf_compe<-genes_amf[genes_amf$V4_B1.CANB1<(-2)&genes_amf$V4_CAN.CANB1<(-2) |
            genes_amf$V5_B1.CANB1<(-2)&genes_amf$V5_CAN.CANB1<(-2) |
               genes_amf$V6_B1.CANB1<(-2)&genes_amf$V6_CAN.CANB1<(-2),]

genes_amf_compe$gene<-rownames(genes_amf_compe)
genes_amf_compe<-merge(genes_amf_compe,Mercator_Rirregularis,by="gene")
genes_amf_compe[,15]<- unlist(lapply(lapply(strsplit(as.character(genes_amf_compe[,13]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_amf_compe[,15]<-as.factor(genes_amf_compe[,15])

genes_amf_compe2<-DGE_AMF$counts[rownames(DGE_AMF_N$E) %in% genes_amf_compe$gene, ] 
#DGE_AMF$counts

tratamiento<-unlist(lapply(strsplit(colnames(genes_amf_compe2),"_"),function (x) x[3]))
cultivar<-gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_amf_compe2),"_"),function (x) x[1]))   ))


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/FigureX.pdf",width=3, height=10,useDingbats = F)
par(mfrow=c(3,1),mar=c(15,8,1,1),mgp=c(5, 1, 0))

stDevs <-tapply(genes_amf_compe2[1,],interaction(tratamiento,cultivar),sd)
means<-tapply(genes_amf_compe2[1,],interaction(tratamiento,cultivar),mean)
mp<-barplot(tapply(genes_amf_compe2[1,],interaction(tratamiento,cultivar),mean),
            beside=T,las =2,ylim=c(0,250),col=c("dodgerblue3","lightsalmon3","darkmagenta"),cex.axis = 2,cex.names = 2,cex.lab=2,
            ylab="Normlized log gene-expression")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)

stDevs <-tapply(genes_amf_compe2[2,],interaction(tratamiento,cultivar),sd)
means<-tapply(genes_amf_compe2[2,],interaction(tratamiento,cultivar),mean)
mp<-barplot(tapply(genes_amf_compe2[2,],interaction(tratamiento,cultivar),mean),
            beside=T,las =2,ylim=c(0,300),col=c("dodgerblue3","lightsalmon3","darkmagenta"),cex.axis = 2,cex.names = 2,cex.lab=2,
            ylab="Normlized log gene-expression")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)

stDevs <-tapply(genes_amf_compe2[3,],interaction(tratamiento,cultivar),sd)
means<-tapply(genes_amf_compe2[3,],interaction(tratamiento,cultivar),mean)
mp<-barplot(tapply(genes_amf_compe2[3,],interaction(tratamiento,cultivar),mean),
            beside=T,las =2,ylim=c(0,1000),col=c("dodgerblue3","lightsalmon3","darkmagenta"),cex.axis = 2,cex.names = 2,cex.lab=2,
            ylab="Normlized log gene-expression")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
dev.off()
######################################################
######################################################

# phosphate transporters YUCA
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
genes_YUCA_compe<-genes_YUCA_compe[genes_YUCA_compe$Type=="transport.phosphate",]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
#DGE_AMF$counts

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
               gsub("V6","CM4574-7",
               gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Phoshphate_transporters_YUCA.pdf",width=10, height=4,useDingbats = F)
par(mfrow=c(1,3),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in c(1,3,5)) {
  stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
  means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()

######################################################
######################################################


# Nitrate transporters YUCA
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
genes_YUCA_compe<-genes_YUCA_compe[genes_YUCA_compe$Type=="transport.nitrate",]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
#DGE_AMF$counts

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Nitrate_transporters_YUCA.pdf",width=10, height=4,useDingbats = F)
par(mfrow=c(1,3),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:2) {
  stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
  means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()

######################################################
######################################################

# AMMonium transporters YUCA
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
genes_YUCA_compe<-genes_YUCA_compe[grep("transport.ammonium",genes_YUCA_compe$Type),]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
#DGE_AMF$counts

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Ammonium_transporters_YUCA.pdf",width=10, height=8,useDingbats = F)
par(mfrow=c(2,3),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:4) {
  stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
  means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()

######################################################
######################################################
# Sugars transporters YUCA
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
genes_YUCA_compe<-genes_YUCA_compe[grep("transport.sugar",genes_YUCA_compe$Type),]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
#DGE_AMF$counts

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Sugar_transporters_YUCA.pdf",width=13, height=8,useDingbats = F)
par(mfrow=c(2,4),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:7) {
stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
            beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
            ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i]))
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()

######################################################
######################################################
# Potasium transporters YUCA
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
genes_YUCA_compe<-genes_YUCA_compe[grep("transport.potassium",genes_YUCA_compe$Type),]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
#DGE_AMF$counts

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Potassium_transporters_YUCA.pdf",width=13, height=8,useDingbats = F)
par(mfrow=c(2,4),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:2) {
  stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
  means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()

######################################################
######################################################
# MEtal transporters YUCA
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
genes_YUCA_compe<-genes_YUCA_compe[grep("transport.metal",genes_YUCA_compe$Type),]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
#DGE_AMF$counts

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Metal_transporters_YUCA.pdf",width=13, height=12,useDingbats = F)
par(mfrow=c(3,4),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:12) {
  stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
  means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()

######################################################
######################################################
# transporters YUCA
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
genes_YUCA_compe<-genes_YUCA_compe[grep("transport",genes_YUCA_compe$Type),]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
table(genes_YUCA_compe$Type)

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Transport_YUCA.pdf",width=13, height=12,useDingbats = F)
par(mfrow=c(3,4),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:dim(genes_YUCA_compe2)[1]) {
  stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
  means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i],"\n",genes_YUCA_compe$Type[genes_YUCA_compe$gene %in% rownames(genes_YUCA_compe2)][i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()


#################################################### 

# only expressed in controls
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
#genes_YUCA_compe<-genes_YUCA_compe[grep("transport.metal",genes_YUCA_compe$Type),]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
#genes_YUCA_compe2<-genes_YUCA_compe2[genes_YUCA_compe2[,8]<100&genes_YUCA_compe2[,12]<100,]
head(genes_YUCA_compe2)
genes_YUCA_compe$gene

genes_YUCA_compe[genes_YUCA_compe$gene %in% "Manes.01G038100" ,]
# NO expression on ctrl
genes_YUCA_compe2<-genes_YUCA_compe2[genes_YUCA_compe2[,8]<20&genes_YUCA_compe2[,12]<20|
                                       genes_YUCA_compe2[,17]<20|genes_YUCA_compe2[,20]<20& 
                                       genes_YUCA_compe2[,23]<20&genes_YUCA_compe2[,26]<20&
                                       genes_YUCA_compe2[,28]<20&genes_YUCA_compe2[,31]<20
                                       ,]
# Expression only in ctrl

only_CTRL<-c("Manes.01G038100","Manes.01G170500","Manes.01G269100","Manes.01G152100",
  "Manes.02G042400","Manes.02G125400","Manes.02G111300","Manes.03G019200",
  "Manes.03G043500","Manes.04G089800","Manes.05G009800","Manes.04G100100",
  "Manes.05G049200","Manes.05G081700","Manes.05G144100","Manes.06G071500",
  "Manes.06G134900","Manes.06G152800","Manes.07G089100","Manes.07G111900",
  "Manes.08G114700","Manes.10G034100","Manes.11G028100","Manes12G121200",
  "Manes.13G029800","Manes.13G033400","Manes.14G035900",
  "Manes.14G036000","Manes.14G118100","Manes.14G130600","Manes.14G132600",
  "Manes.14G173500","Manes.15G138100","Manes.16G011000","Manes.16G016100",
  "Manes.16G065800","Manes.16G086700","Manes.16G095000","Manes.17G058500",
  "Manes.18G019000","Manes.18G052600","Manes.S053600","Manes.S073300")

genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% only_CTRL, ] 
genes_YUCA_compe$adj.P.Val[genes_YUCA_compe$gene %in% only_CTRL]
table(droplevels(genes_YUCA_compe$details[genes_YUCA_compe$gene %in% only_CTRL]))

table(droplevels(genes_YUCA_compe$gene[genes_YUCA_compe$gene %in% rownames(genes_YUCA_compe2)]))
head(genes_YUCA_compe[genes_YUCA_compe$gene %in% rownames(genes_YUCA_compe2),])
# Expression no funcitonal on vars

no_funky<-c("Manes.01G050800","Manes.02G003200","Manes.03G000100","Manes.05G13970","Manes.07G027400","Manes.08G014200","Manes.08G086500","Manes.12G004600","Manes.16G056100","Manes.17G034700","Manes.17G060300","Manes.17G093900","Manes.S036400","Manes.S111200","Manes.13G05080")
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% no_funky, ] 
genes_YUCA_compe$gene[genes_YUCA_compe$gene %in% no_funky]
table(droplevels(genes_YUCA_compe$Type[genes_YUCA_compe$gene %in% no_funky]))


# Coinoc higher than single

genes_YUCA_compe2<-genes_YUCA_compe2[genes_YUCA_compe2[,2]>genes_YUCA_compe2[,1]&genes_YUCA_compe2[,3]>genes_YUCA_compe2[,10]|
                                       genes_YUCA_compe2[,16]>genes_YUCA_compe2[,15]&genes_YUCA_compe2[,19]>genes_YUCA_compe2[,18]& 
                                       genes_YUCA_compe2[,5]>genes_YUCA_compe2[,29]&genes_YUCA_compe2[,6]>genes_YUCA_compe2[,20],]

table(droplevels(genes_YUCA_compe$Type[genes_YUCA_compe$gene %in% rownames(genes_YUCA_compe2)]))

testss<-genes_YUCA_compe[genes_YUCA_compe$gene %in% rownames(genes_YUCA_compe2),]
testss[grep("transport.phosphate",testss$Type),]

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/METAL_TRANSPORT_YUCA.pdf",width=13, height=12,useDingbats = F)
par(mfrow=c(3,4),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:dim(genes_YUCA_compe2)[1]) {
  stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
  means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
  
}
dev.off()
rownames(genes_YUCA_compe2)
######################################################
######################################################
# ABC transporters YUCA
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
genes_YUCA_compe<-genes_YUCA_compe[grep("transport.ABC_transporters_and_multidrug_resistance_systems",genes_YUCA_compe$Type),]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
table(genes_YUCA_compe$Type)

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/ABC_transporters_YUCA.pdf",width=13, height=12,useDingbats = F)
par(mfrow=c(3,4),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:13) {
  stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
  means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()

######################################################
######################################################
# CELL wall YUCA
genes_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
genes_YUCA_compe<-genes_YUCA
genes_YUCA_compe$gene<-rownames(genes_YUCA_compe)
genes_YUCA_compe<-merge(genes_YUCA_compe,Mercator_Mesculenta,by="gene")
genes_YUCA_compe[,18]<- unlist(lapply(lapply(strsplit(as.character(genes_YUCA_compe[,16]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_YUCA_compe[,18]<-as.factor(genes_YUCA_compe[,18])
genes_YUCA_compe<-genes_YUCA_compe[grep("cell_wall",genes_YUCA_compe$Type),]
genes_YUCA_compe2<-DGE_YUCA$counts[rownames(DGE_YUCA$counts) %in% genes_YUCA_compe$gene, ] 
table(genes_YUCA_compe$Type)

tratamiento<-unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_YUCA_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(9:12,1:4,5:8)]

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/CELL_WALL_YUCA.pdf",width=13, height=12,useDingbats = F)
par(mfrow=c(3,4),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:dim(genes_YUCA_compe2)[1]) {
  stDevs <-tapply(genes_YUCA_compe2[i,],interacts,sd)
  means<-tapply(genes_YUCA_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_YUCA_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_YUCA_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()



######################################################
# GENES AMF
genes_AMF<-sign_AMF[sign_AMF$adj.P.Val<0.05,]
genes_AMF_compe<-genes_AMF
genes_AMF_compe$gene<-rownames(genes_AMF_compe)
genes_AMF_compe<-merge(genes_AMF_compe,Mercator_Rirregularis,by="gene")
genes_AMF_compe[,15]<- unlist(lapply(lapply(strsplit(as.character(genes_AMF_compe[,13]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_AMF_compe[,15]<-as.factor(genes_AMF_compe[,15])
colnames(genes_AMF_compe)[15]<-"Type"
#genes_AMF_compe<-genes_AMF_compe[grep("secondary_metabolism",genes_AMF_compe$Type),]
genes_AMF_compe2<-DGE_AMF$counts[rownames(DGE_AMF$counts) %in% genes_AMF_compe$gene, ] 
table(genes_AMF_compe$Type)

tratamiento<-unlist(lapply(strsplit(colnames(genes_AMF_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_AMF_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)
levels(interacts)<-levels(interacts)[c(7:9,1:3,4:6)]

# NO functional in variety
genes_AMF_compe2<-genes_AMF_compe2[genes_AMF_compe2[,1]<20&genes_AMF_compe2[,12]<20&
                                       genes_AMF_compe2[,13]<20&genes_AMF_compe2[,14]<20& 
                                       genes_AMF_compe2[,15]<20&genes_AMF_compe2[,16]<20&
                                       genes_AMF_compe2[,17]<20&genes_AMF_compe2[,19]<20
                                     ,]


# No functional in one isolate
genes_AMF_compe2<-genes_AMF_compe2[genes_AMF_compe2[,2]<20&genes_AMF_compe2[,3]<20&
                                     genes_AMF_compe2[,4]<20|genes_AMF_compe2[,7]<20& 
                                     genes_AMF_compe2[,8]<20|genes_AMF_compe2[,9]<20|
                                     genes_AMF_compe2[,21]<20&genes_AMF_compe2[,22]<20
                                   ,]
colnames(genes_AMF_compe2)
genes_AMF_compe$details[genes_AMF_compe$gene %in% rownames(genes_AMF_compe2) ]
table(genes_AMF_compe2$Type)
genes_AMF_compe[genes_AMF_compe$gene %in% "g7952.t1" ,]

as.vector(genes_AMF_compe2[,1]<20&genes_AMF_compe2[,12]<20|genes_AMF_compe2[,8]<20&genes_AMF_compe2[,12]<20)

B1_only<-c("g10037.t1","g10150.t1","g10184.t1","g10516.t1","g10543.t1","g10544.t1","g10706.t1","g11014.t1","g11149.t1","g11302.t1","g11367.t1","g11796.t1","g11848.t1","g11910.t1","g11918.t1","g11933.t1","g12393.t1","g12735.t1","g12737.t1","g12819.t1","g12825.t1","g12826.t1","g13061.t1","g13114.t1","g13123.t1","g13141.t1","g13172.t1","g13287.t1","g13327.t1","g13415.t1","g13541.t1","g13557.t1","g13584.t1","g13652.t1","g13685.t1","g13742.t1","g13783.t1","g13969.t1","g13991.t1","g14012.t1","g14030.t1","g14088.t1","g14253.t1","g14254.t1","g14404.t1","g14583.t1","g14593.t1","g14609.t1","g14724.t1","g14727.t1","g14788.t1","g14871.t1","g15124.t1","g15140.t1","g15153.t1","g15389.t1","g15610.t1","g15802.t1","g15889.t1","g15896.t1","g1735.t1","g1843.t1","g1904.t1","g2234.t1","g2484.t1","g2674.t1","g2981.t1","g3780.t1","g4743.t1","g4852.t1","g5904.t1","g6106.t1","g6212.t1","g6314.t1","g6619.t1","g7021.t1","g7098.t1","g7175.t1","g7674.t1","g7810.t1","g782.t1","g8310.t1","g9045.t1","g9111.t1","g9125.t1","g9339.t1","g9614.t1","g9921.t1","g9960.t1")
genes_AMF_compe2<-genes_AMF_compe2[rownames(genes_AMF_compe2) %in% B1_only ,]

CAN_only<-c("g10176.t1","g10244.t1","g10710.t1","g10743.t1","g10846.t1","g10863.t1","g10904.t1","g11074.t1","g11165.t1","g11203.t1","g11209.t1","g11233.t1","g11340.t1","g11370.t1","g11382.t1","g11383.t1","g11384.t1","g11486.t1","g11593.t1","g11730.t1","g11901.t1","g12035.t1","g12035.t1","g12164.t1","g12246.t1","g12414.t1","g12473.t1","g12500.t1","g12680.t1","g12763.t1","g12764.t1","g12873.t1","g12877.t1","g12913.t1","g12952.t1","g13002.t1","g13046.t1","g13119.t1","g13150.t1","g13247.t1","g13248.t1","g13267.t1","g13268.t1","g13306.t1","g13345.t1","g13352.t1","g13449.t1","g13523.t1","g13524.t1","g13534.t1","g13563.t1","g13693.t1","g13805.t1","g13885.t1","g13947.t1","g14198.t1","g14290.t1","g14303.t1","g14344.t1","g14378.t1","g14390.t1","g14403.t1","g1444.t1","g1446.t1","g14514.t1","g14550.t1","g14556.t1","g14581.t1","g14646.t1","g14717.t1","g14785.t1","g1480.t1","g14985.t1","g15034.t1","g15104.t1","g15110.t1","g15208.t1","g15266.t1","g15274.t1","g15397.t1","g15599.t1","g15701.t1","g15733.t1","g15798.t1","g2705.t1","g2939.t1","g3084.t1","g3298.t1","g3959.t1","g3960.t1","g4665.t1","g5313.t1","g5381.t1","g5644.t1","g5790.t1","g6262.t1","g6329.t1","g6330.t1","g6630.t1","g7541.t1","g7644.t1","g7645.t1","g7646.t1","g7952.t1","g8356.t1","g8357.t1","g8358.t1","g8581.t1","g8624.t1","g8625.t1","g8880.t1","g9052.t1","g9162.t1","g9542.t1","g9543.t1","g9562.t1","g9890.t1","g9955.t1","g9970.t1")
genes_AMF_compe2<-genes_AMF_compe2[rownames(genes_AMF_compe2) %in% CAN_only ,]

table(droplevels(genes_AMF_compe$Type[genes_AMF_compe$gene %in% rownames(genes_AMF_compe2)]))
# co inoculation higer than treatments
genes_AMF_compe2<-genes_AMF_compe2[genes_AMF_compe2[,2]>genes_AMF_compe2[,7]&genes_AMF_compe2[,3]>genes_AMF_compe2[,8]&
                                     genes_AMF_compe2[,14]>genes_AMF_compe2[,13]&genes_AMF_compe2[,16]>genes_AMF_compe2[,17]&
                                     genes_AMF_compe2[,5]>genes_AMF_compe2[,20]&genes_AMF_compe2[,6]>genes_AMF_compe2[,21]                  
                                   ,]

# BRA 337 CANB1 higher than CAN, than B1
genes_AMF_compe2<-genes_AMF_compe2[genes_AMF_compe2[,2]>genes_AMF_compe2[,7]&genes_AMF_compe2[,3]>genes_AMF_compe2[,8]&
                                     genes_AMF_compe2[,14]>genes_AMF_compe2[,13]&genes_AMF_compe2[,16]>genes_AMF_compe2[,17]&
                                     genes_AMF_compe2[,5]>genes_AMF_compe2[,20]&genes_AMF_compe2[,6]>genes_AMF_compe2[,21]                  
                                   ,]
CANdiff<-list();B1diff<-list()
for (i in 1:dim(genes_AMF_compe2)[1]) {
CANdiff[[i]]<-mean(genes_AMF_compe2[i,14],genes_AMF_compe2[i,16],genes_AMF_compe2[i,19]) > mean(genes_AMF_compe2[i,1],genes_AMF_compe2[i,12],genes_AMF_compe2[i,17])*1.2
B1diff[[i]]<-mean(genes_AMF_compe2[i,14],genes_AMF_compe2[i,16],genes_AMF_compe2[i,19]) > mean(genes_AMF_compe2[i,13],genes_AMF_compe2[i,15],genes_AMF_compe2[i,18])*1.2
}
table(unlist(CANdiff))
table(unlist(B1diff))
# cell wall
genes_AMF_compe<-genes_AMF_compe[grep("transport",genes_AMF_compe$Type),]
genes_AMF_compe2<-DGE_AMF$counts[rownames(DGE_AMF$counts) %in% genes_AMF_compe$gene, ] 


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Secondary_Metabolism_AMF.pdf",width=13, height=12,useDingbats = F)
par(mfrow=c(3,4),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:dim(genes_AMF_compe2)[1]) {
  stDevs <-tapply(genes_AMF_compe2[i,],interacts,sd)
  means<-tapply(genes_AMF_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_AMF_compe2[i,],interacts,mean),ylim=c(0,max(means + stDevs)*1.1),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","darkmagenta"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_AMF_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()





########################################################################################################################################




########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


################################################################################################################################
# GO ANALYSIS FOR DIFFERENTIALLY EXPRESSED
###########################


#annotation infos
########################################################################################################################################
#M esculenta
mesculenta_go2<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/annotations_blast2go_mesculenta_FMv2/blast2go_FM_v3',h=T)
colnames(mesculenta_go2)<-c("locusName","transcriptName","GO","DEF")
#R. irregularis
rirregularis_go<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/blast2go_annot_predicted_prot_hint_glomus_nu6_genome_masked.annot',h=F)
colnames(rirregularis_go)<-c('locusName','GO.ID')
head(mesculenta_go)

# CASSAVA
######################################################
# 2) ANOTATION DATABASE
######################################################
## CASSAVA
locusName<-rownames(sign_YUCA)
sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
sign_feat_go<-merge(sign_feat2,mesculenta_go2,by='locusName')
sign_feat_go<-na.omit(sign_feat_go)
gene_2_go<-sign_feat_go[,c(1,16)]
head(sign_feat_go)
head(gene_2_go)
write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
            col.names = F, row.names = F)
geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")

######################################################
# 3) SELECT UNIVERSE
######################################################
# TOTAL UNIVERSE CASSAVA
locusName<-rownames(sign_YUCA)
sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
sign_feat_go<-merge(sign_feat2,mesculenta_go2[,c(1,3)],by='locusName') # do not select transcriptname. working at gene level
sign_feat_go<-na.omit(sign_feat_go)
sign_feat_go<-unique(sign_feat_go)
length(unique(sign_feat_go$locusName))
######################################################
# 4) SELECT TOP GENES
######################################################
## DEG SIGNFICANT 
Significant_GENES_Y<-rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])
Significant_GENES_Y<-rownames(V4@.Data[V4@.Data[,3]!=0&V4@.Data[,1]==0&V4@.Data[,2]==0,])

# modules
Significant_GENES_Y<-rownames(gene_names_YUCA[gene_names_YUCA$moduleColors_Y=="black",])

######################################################
# 5) ANNOTATION
######################################################
#AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
geneList <- factor(as.integer(sign_feat_go$locusName %in% Significant_GENES_Y))
names(geneList) <- sign_feat_go$locusName
sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5,ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)

sampleGOdata
resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)

allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                           classicKS = resultKS, elimKS = resultKS.elim,
                           orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)

RES_Y<-allRes
RES_Y
#write.table(RES_Y,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/2.GO_CASSAVA_SPE_CANB1.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = F)
write.table(RES_Y[RES_Y[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/TEMP_module_yuca.txt',sep="_")
            ,quote=FALSE,sep='\t',col.names = F, row.names = F)
#png("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_MODULE_TURQUOISE_CASSAVA_GO_ENRICH.png", width=1400, height=1000, units="px")
#par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
barplot(RES_Y$Significant[RES_Y$classicFisher<0.05], names = paste(RES_Y$GO.ID[RES_Y$classicFisher<0.05],RES_Y$Term[RES_Y$classicFisher<0.05],sep=" "),
        xlab = "Nb. of significant genes",horiz=T,cex.lab=2,las=1,cex.axis=2)
#dev.off()
to_pheat_Y<-RES_Y[RES_Y$classicFisher<0.05,c(6)]
to_pheat_Y<-as.numeric(to_pheat_Y)
names(to_pheat_Y)<-paste(RES_Y[RES_Y$classicFisher<0.05,c(1)],RES_Y[RES_Y$classicFisher<0.05,c(2)],sep=" ")

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.GO_CASSAVA_CANB1-to-single-inoculations.pdf",width=12, height=6)
pheatmap(t(to_pheat_Y),cluster_cols = T, cluster_rows=F, cellwidth=12, cellheight=12,fontsize = 13,
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)) )
dev.off()

#################################################################################
# AMF
######################################################
# 2) ANOTATION DATABASE
######################################################
#AMF
locusName<-rownames(sign_AMF)
sign_feat2<-cbind.data.frame(sign_AMF,locusName)
sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName')
sign_feat_go<-na.omit(sign_feat_go)
gene_2_go<-sign_feat_go[,c(1,12)]

write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
            col.names = F, row.names = F)
geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")


######################################################
# 3) SELECT UNIVERSE
######################################################

#TOTAL UNIVERSE ALL LOCI
#TOTAL UNIVERSE AMF
locusName<-rownames(sign_AMF)
sign_feat2<-cbind.data.frame(sign_AMF,locusName)
sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName') # do not select transcriptname. working at gene level
sign_feat_go<-na.omit(sign_feat_go)
sign_feat_go<-unique(sign_feat_go)
length(unique(sign_feat_go$locusName))
######################################################
# 4) SELECT TOP GENES
######################################################
## DEG SIGNFICANT 
#Significant_GENES_A<-rownames(Results_A[Results_A[,1]=="Different 2 CANB1"|Results_A[,2]=="Different 2 CANB1",])
Significant_GENES_A<-sign_AMF[sign_AMF$adj.P.Val<0.05,]
Significant_GENES_A<-V6@.Data[V6@.Data[,1]!=0|V6@.Data[,2]!=0,]
######################################################
# 5) ANNOTATION
######################################################
#AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
geneList <- factor(as.integer(sign_feat_go$locusName %in% rownames(Significant_GENES_A) ))
names(geneList) <- sign_feat_go$locusName
sampleGOdata <- new("topGOdata", nodeSize = 5,ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO)
sampleGOdata
resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)

allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                           classicKS = resultKS, elimKS = resultKS.elim,
                           orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
RES_A<-allRes
RES_A
############## PRINT RESULTS
#write.table(RES_A,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/1.GO_AMF.txt',quote=FALSE,sep='\t',
        #    col.names = T, row.names = F)
write.table(RES_A[RES_A[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.GO_AMF.txt',sep="_")
            ,quote=FALSE,sep='\t',col.names = F, row.names = F)

#png("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_MODULE_TURQUOISE_CASSAVA_GO_ENRICH.png", width=1400, height=1000, units="px")
#par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
barplot(RES_A$Significant[RES_A$classicFisher<0.05], names = paste(RES_A$GO.ID[RES_A$classicFisher<0.05],RES_A$Term[RES_A$classicFisher<0.05],sep=" "),
        xlab = "Nb. of significant genes",horiz=T,cex.lab=2,las=1,cex.axis=2)
#dev.off()
to_pheat_A<-RES_A[RES_A$classicFisher<0.05,c(6)]
to_pheat_A<-as.numeric(to_pheat_A)
names(to_pheat_A)<-paste(RES_A[RES_A$classicFisher<0.05,c(1)],RES_A[RES_A$classicFisher<0.05,c(2)],sep=" ")

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.GO_AMF_CANB1-to-single-inoculations.pdf",width=8, height=6)
pheatmap(t(to_pheat_A),cluster_cols = T, cluster_rows=F, cellwidth=12, cellheight=12,,fontsize = 13,
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)) )
dev.off()



################################################################################################

################################################################################################
##### Function Mercator tools

sign_YUCA2<-sign_YUCA[rownames(sign_YUCA) %in% Significant_GENES_Y,]

sign_YUCA2$gene<-rownames(sign_YUCA[rownames(sign_YUCA) %in% Significant_GENES_Y,])
sign_YUCA2<-merge(sign_YUCA2,Mercator_Mesculenta,by="gene")
head(sign_YUCA2)
write.table( sign_YUCA2[,c(1:10,14,16:17)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.DEG_Cassava_CANB1-single-inoculation_v1.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)




Significant_GENES_A$gene<-rownames(Significant_GENES_A)
Significant_GENES_A<-merge(Significant_GENES_A,Mercator_Rirregularis,by="gene")
head(Significant_GENES_A)
write.table( Significant_GENES_A[,c(1:7,11,13:14)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/3.DEG_AMF_CANB1-single-inoculation_v1.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)
################################################################################################

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################### PHENOTYPIC MEASUREMENTS
################################
####### 1. COLONIZATION MEASUREMENTS
col<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/Phenodata_vf/colonization_vf.txt",h=T) #tablero
summary(col)
#take out non used data in RNASeq
#col<-col[col$treat!="A",]
col<-col[col$treat!="A",]
col<-col[col$var!="V1_CM6438-14",]
col<-col[col$var!="V2_COL1468",]
col<-col[col$var!="V3_TAI-8",]
col<-col[col$var!="V7_VEN329A",]
col<-col[col$var!="V8_CM523-7",]
col<-col[col$var!="V9_NGA-16",]

col<-droplevels(col)
gsub("//+","",col$treat)
summary(col)
head(col)
col$treat<-droplevels(col$treat)
glm_col<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var*col$treat,family=binomial)
glm_col2<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var+col$treat,family=binomial)
glm_col3<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var,family=binomial)

anova(glm_col,glm_col2,test='Chisq')
anova(glm_col2,glm_col3,test='Chisq')
anova(glm_col)
summary(glm_col)



######################################################################################################
########### FINAL COLONIZATION ANALYSIS
### ANALYSIS WITHOUT CONTROL

col<-col[col$treat!="A",]
#col<-col[col$var!="V2_COL1468",];col<-col[col$var!="V7_VEN329A",];col<-col[col$var!="V9_NGA-16",];col<-col[col$var!="V3_TAI-8",]

summary(col)
col$treat<-droplevels(col$treat)
col$var<-droplevels(col$var)

# using normal anova, not good test
#glm_col<-glm(cbind(Positive,(Obs-Positive))~var*treat,data= col,family=binomial)
#glm_col2<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var+col$treat,family=binomial)
#glm_col3<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var,family=binomial)
#anova(glm_col2,glm_col,test='Chisq')

#anova(glm_col,test="Chisq")

##GOOD STATISTICS
library(lsmeans)
glm_col<-glm(cbind(Positive,(Obs-Positive))~var*treat,data= col,family=binomial)
lsmeans(glm_col, pairwise ~ treat | var)



# multiple comparisons
#glm_col <- glm(cbind(Positive,(Obs-Positive)) ~ var * treat, data= col,family=binomial)
glm_col<-aov(col$Positive/col$Obs~col$treat*col$var)
#TukeyHSD(aov(col$Positive/col$Obs~col$treat*col$var))

summary(glm_col)


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/1.PHENO_colonization_v1.pdf",width=12, height=10,useDingbats = F)
par(mfrow=c(1,1),mar=c(15,8,1,1),mgp=c(5, 1, 0))
boxplot(col$perc~interaction(col$treat,paste("V",col$cultivar,sep="")),las=2,ylim=c(0,1),ylab="% of colonization (n=10)",
     col=c("dodgerblue3","lightsalmon3","darkmagenta"),cex.lab=3,cex.axis=3)
dev.off()

tapply(col$perc, interaction(col$treat,col$var), mean)
###############################
######## 3. FINAL WEIGHT and HEIGHT MEASUREMENTs
###??TO MANUSCRIPT
################################

dry<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/Phenodata_vf/Harvest_VF.txt",h=T)
summary(dry)
dry$BLOCK<-as.factor(dry$BLOCK)
dry$VARIETY<-as.factor(dry$VARIETY)
dry$aerial<-as.numeric(dry$aerial)
dry$roots<-as.numeric(dry$roots)
dry$tuber<-as.numeric(dry$tuber)
dry$total<-as.numeric(dry$total)
dry$ratio<-as.numeric(dry$ratio)
dry$X125dai<-as.numeric(dry$X125dai)

dry$TREATMENT<-relevel(dry$TREATMENT, ref = "C") #order factors to put control as reference
dry$AMF<-relevel(dry$AMF, ref = "CTRL") #order factors to put control as reference

all_root<-dry$roots+dry$tuber
dry<-cbind(dry,all_root)
#dry<-dry[dry$TREATMENT!="A+B",]

#dry<-dry[dry$TREATMENT!="C",]
#dry<-dry[dry$CULTIVARIETY!=2,];dry<-dry[dry$CULTIVARIETY!=7,];dry<-dry[dry$CULTIVARIETY!=9,];dry<-dry[dry$CULTIVARIETY!=3,]
dry<-dry[dry$CULTIVAR!="V1_CM6438-14",]
dry<-dry[dry$CULTIVAR!="V2_COL1468",]
dry<-dry[dry$CULTIVAR!="V3_TAI-8",]
dry<-dry[dry$CULTIVAR!="V7_VEN329A",]
dry<-dry[dry$CULTIVAR!="V8_CM523-7",]
dry<-dry[dry$CULTIVAR!="V9_NGA-16",]

dry<-droplevels(dry)

######################
##GOOD STATISTICS
library(lsmeans)
ALLV_lme<-lme(total~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
lsmeans(ALLV_lme, pairwise ~ AMF | CULTIVAR)
######################
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/1.PHENO_totalDW_v1.pdf",width=12, height=10,useDingbats = F)
par(mfrow=c(1,1),mar=c(15,8,1,1),mgp=c(5, 1, 0))
boxplot(dry$total~interaction(dry$AMF,paste("V",dry$VARIETY,sep="")),las=2,cex.axis=3,cex.lab=3,
        ylab="Total dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/1.PHENO_harvest_v1.pdf",width=18, height=16)
par(mfrow=c(2,2),mar=c(20,8,6,2),mgp=c(6, 2, 0))

#boxplot(dry$total~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
 #       ylab="Total dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
boxplot(dry$aerial~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Aboveground dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
#boxplot(dry$roots~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
 #       ylab="Fine roots dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
#boxplot(dry$tuber~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
#        ylab="Tuberized roots dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
boxplot(dry$all_root~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Underground dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
boxplot(dry$X125dai~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
  ylab="Height (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))

dev.off()
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Supplementary_Figure6.pdf",width=12, height=10,useDingbats = F)
par(mfrow=c(1,1),mar=c(15,8,1,1),mgp=c(5, 1, 0))
  boxplot(dry$ratio~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Ratio root/shoots dry weigth (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
dev.off()
### ALL VARIETYieties

######### total dry weight
ALLV_lme<-lme(total~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
interaction.plot(dry$AMF,dry$CULTIVAR,dry$total)
anova(ALLV_lme)
summary(ALLV_lme)
summary(glht(ALLV_lme, linfct=mcp(CULTIVAR="Tukey")) )
summary(glht(ALLV_lme, linfct=mcp(AMF="Tukey")) )

## differences between CAN B1
#dry$sxt<- interaction(dry$AMF,dry$CULTIVAR)
#ALLV_lme<-lme(aerial~-1+sxt,random=~1|BLOCK,data=dry)

#tmp <- expand.grid(treat = unique(dry$AMF), var = unique(dry$CULTIVAR))
#X <- model.matrix(~ var * treat, data = tmp)

#Tukey <- contrMat(table(dry$AMF), "Tukey")
#K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*1))
#rownames(K1) <- paste(levels(dry$CULTIVAR)[1], rownames(K1), sep = ":")
#K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
#rownames(K2) <- paste(levels(dry$CULTIVAR)[2], rownames(K2), sep = ":")
#K <- rbind(K1, K2)
#colnames(K) <- rep(colnames(Tukey), 2)
#summary(glht(ALLV_lme, linfct = K %*% X,alternative="two.sided"),test=adjusted("Shaffer"))


########################
ALLV_lme<-lme(aerial~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)
summary(glht(ALLV_lme, linfct=mcp(CULTIVAR="Tukey")) )
summary(glht(ALLV_lme, linfct=mcp(AMF="Tukey")) )



ALLV_lme<-lme(roots~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)
summary(glht(ALLV_lme, linfct=mcp(CULTIVAR="Tukey")) )
summary(glht(ALLV_lme, linfct=mcp(AMF="Tukey")) )

ALLV_lme<-lme(tuber~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)
summary(glht(ALLV_lme, linfct=mcp(CULTIVAR="Tukey")) )
summary(glht(ALLV_lme, linfct=mcp(AMF="Tukey")) )


ALLV_lme<-lme(X125dai~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)
summary(glht(ALLV_lme, linfct=mcp(CULTIVAR="Tukey")) )
summary(glht(ALLV_lme, linfct=mcp(AMF="Tukey")) )


#############################
# only samples on RNAseq
onlyRNA<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/Phenodata_vf/Harvest_COINOC.txt",h=T)
onlyRNA<-onlyRNA[onlyRNA$CULTIVAR!="V6_CM4574-7",]
onlyRNA<-onlyRNA[onlyRNA$CULTIVAR!="V8_CM523-7",]

onlyRNA$AMF<-relevel(onlyRNA$AMF, ref = "CTRL") #order factors to put control as reference
onlyRNA<-droplevels(onlyRNA)
all_root<-onlyRNA$roots+onlyRNA$tuber
onlyRNA<-cbind(onlyRNA,all_root)

#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/1.PHENO_harvest_B_v1.pdf",width=18, height=16)
par(mfrow=c(2,3),mar=c(20,8,6,2),mgp=c(6, 2, 0))

boxplot(onlyRNA$total~interaction(onlyRNA$AMF,onlyRNA$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Total dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
boxplot(onlyRNA$aerial~interaction(onlyRNA$AMF,onlyRNA$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Aboveground dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
boxplot(onlyRNA$roots~interaction(onlyRNA$AMF,onlyRNA$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Fine roots dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
boxplot(onlyRNA$tuber~interaction(onlyRNA$AMF,onlyRNA$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Tuberized roots dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
boxplot(onlyRNA$all_root~interaction(onlyRNA$AMF,onlyRNA$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Underground dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))
boxplot(onlyRNA$X125dai~interaction(onlyRNA$AMF,onlyRNA$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Height (n=14)",col=c("azure3","dodgerblue3","lightsalmon3","darkmagenta"))

#dev.off()

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################




################################################################################################
##SUPPLEMENTARY INFO ABOUT RNA-SEQ.
# Correlation replicate 1 vs 2. for each treatment.
#accumulaticon curve

for (vari in c("V1","V3","V4","V5","V6","V8")){
CORR_1<-DGE_YUCA_N$E[,grep(c(vari), colnames(DGE_YUCA_N$E),invert =F)]
CORR_2<-CORR_1[,grep(c("CAN"),colnames(CORR_1))]
CORR_3<-CORR_1[,grep(c("B1"),colnames(CORR_1))]
CORR_4<-CORR_1[,grep(c("CTRL"),colnames(CORR_1))]
png(paste("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor_CASSAVA_",vari,".png",sep=""), width=800, height=600, units="px")
par(mar=c(5,5,3,1),las=1,mgp=c(3, 1, 0))
par(mfrow=c(3,3))
plot(CORR_2[,c(1,2)],cex.lab=2,cex.axis=2);plot(CORR_2[,c(1,3)],cex.lab=2,cex.axis=2);plot(CORR_2[,c(2,3)],cex.lab=2,cex.axis=2)
plot(CORR_3[,c(1,2)],cex.lab=2,cex.axis=2);plot(CORR_3[,c(1,3)],cex.lab=2,cex.axis=2);plot(CORR_3[,c(2,3)],cex.lab=2,cex.axis=2)
plot(CORR_4[,c(1,2)],cex.lab=2,cex.axis=2);plot(CORR_4[,c(1,3)],cex.lab=2,cex.axis=2);plot(CORR_4[,c(2,3)],cex.lab=2,cex.axis=2)
dev.off()
}

par(mfrow=c(3,3))
for (vari in c("V1","V3","V4","V5","V6","V8")){
  CORR_1<-DGE_AMF_N$E[,grep(c(vari), colnames(DGE_AMF_N$E),invert =F)]
  CORR_2<-CORR_1[,grep(c("CAN"),colnames(CORR_1))]
  CORR_3<-CORR_1[,grep(c("B1"),colnames(CORR_1))]
  png(paste("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor_AMF_",vari,".png",sep=""), width=800, height=600, units="px")
  par(mar=c(5,5,3,1),las=1,mgp=c(3, 1, 0))
  par(mfrow=c(2,3))
  plot(CORR_2[,c(1,2)],cex.lab=2,cex.axis=2);plot(CORR_2[,c(1,3)],cex.lab=2,cex.axis=2);plot(CORR_2[,c(2,3)],cex.lab=2,cex.axis=2)
  plot(CORR_3[,c(1,2)],cex.lab=2,cex.axis=2);plot(CORR_3[,c(1,3)],cex.lab=2,cex.axis=2);plot(CORR_3[,c(2,3)],cex.lab=2,cex.axis=2)
  dev.off()
}


###########################################################################################################################
################################################################################
# summary transcripts TPM -> GENE
DTE_G_YUCA<-c('Total','More_than_0','more_than10','more_than100')
for (x in 1:dim(FOUR_VARS_YUCA)[2]){
samp<-c(length(FOUR_VARS_YUCA[,x]),dim(FOUR_VARS_YUCA[FOUR_VARS_YUCA[,x]>0,])[1],dim(FOUR_VARS_YUCA[FOUR_VARS_YUCA[,x]>10,])[1],dim(FOUR_VARS_YUCA[FOUR_VARS_YUCA[,x]>100,])[1])
DTE_G_YUCA<-cbind.data.frame(DTE_G_YUCA,samp)}

colnames(DTE_G_YUCA)[2:(dim(FOUR_VARS_YUCA)[2]+1)]<-colnames(FOUR_VARS_YUCA)
rownames(DTE_G_YUCA)<-c('Total','More_than_0','more_than10','more_than100')
DTE_G_YUCA
write.table(DTE_G_YUCA,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Genes_summary_CASSAVA_sample.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)



DTE_G_AMF<-c('Total','More_than_0','more_than10','more_than100')
for (x in 1:dim(FOUR_VARS_AMF)[2]){
  samp<-c(length(FOUR_VARS_AMF[,x]),dim(FOUR_VARS_AMF[FOUR_VARS_AMF[,x]>0,])[1],dim(FOUR_VARS_AMF[FOUR_VARS_AMF[,x]>10,])[1],dim(FOUR_VARS_AMF[FOUR_VARS_AMF[,x]>100,])[1])
  DTE_G_AMF<-cbind.data.frame(DTE_G_AMF,samp)}
colnames(DTE_G_AMF)[2:(dim(FOUR_VARS_AMF)[2]+1)]<-colnames(FOUR_VARS_AMF)
DTE_G_AMF
rownames(DTE_G_AMF)<-c('Total','More_than_0','more_than10','more_than100')
write.table(DTE_G_AMF,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Genes_summary_CASSAVA_sample.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)



####################################################################
####################################################################
####################################################################
####################################################################
####################################################################



###########################
# select significant compared to reference
YUCA_sign_TREAT<-rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])
FOUR_VARS_YUCA_S<-DGE_YUCA_N$E[rownames(DGE_YUCA_N) %in% YUCA_sign_TREAT,]
FOUR_VARS_YUCA_S<-FOUR_VARS_YUCA_S[,grep("CTRL",colnames(FOUR_VARS_YUCA_S),invert=T)]
head(FOUR_VARS_YUCA_S)
AMF_sign_TREAT<-rownames(sign_AMF[sign_AMF$adj.P.Val<0.05,])
FOUR_VARS_AMF_S<-DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_sign_TREAT,]


####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

## create eigene modules per sample. ##
###################################
#### M esculenta gene counts ######
###################################
colnames(t(FOUR_VARS_YUCA_S))
#create set of powers
powers=c(c(1:10),seq(from=12, to=30,by=2))
#network topology analysis fonction, show power of networks
sft_Y=pickSoftThreshold(t(FOUR_VARS_YUCA_S),powerVector = powers,verbose = 5)
par(mfrow=c(1,2))
plot(sft_Y$fitIndices[,1], -sign(sft_Y$fitIndices[,3])*sft_Y$fitIndices[,2],xlab="SoftThresholdpower)",ylab="ScaleFreeTopologyModelFit,signed^2",main =paste("Scaleindependence"),type="n")
text(sft_Y$fitIndices[,1], -sign(sft_Y$fitIndices[,3])*sft_Y$fitIndices[,2],labels=powers)

plot(sft_Y$fitIndices[,1], sft_Y$fitIndices[,5],type="n")
text(sft_Y$fitIndices[,1], sft_Y$fitIndices[,5],labels=powers)

####################################################
## clustering methods**
####################################################

# create adjacency,
softPower = 12
adjacency_Y = adjacency(t(FOUR_VARS_YUCA_S),power= softPower);
#creat topological overlap matrix
TOM_Y= TOMsimilarity(adjacency_Y)
dissTOM_Y=1-TOM_Y
#Clustering using TOM
geneTree_Y= hclust(as.dist(dissTOM_Y),method="average")
#sizeGrWindow(12,9)
#plot(geneTree_Y,main="Gene clustering based on TOM dissimilarity",labels=F,hang=0.04)
####################################################
#definition of module size.
####################################################
minModuleSize=30
#module identification using dynamic tree cut. module
dynmicMods_Y= cutreeDynamic(dendro=geneTree_Y, distM=dissTOM_Y,deepSplit=2,pamRespectsDendro=F,minClusterSize=minModuleSize)
table(dynmicMods_Y)

dynamicColors_Y=labels2colors(dynmicMods_Y)
table(dynamicColors_Y)
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree_Y,dynamicColors_Y,"Dynamic Tree Cut",dendroLabels= F,hang=0.03,addGuide = T,guideHang = 0.05,main="Gene dendrogram and module colors")
####################################################
#merging modules whose expression profiles are similar
####################################################
#calculate eignegens
MEList_Y=moduleEigengenes(t(FOUR_VARS_YUCA_S),colors=dynamicColors_Y)
MEs_Y=MEList_Y$eigengenes
#calculate dissimilarity od module eigengenes
MEDiss_Y= 1 -cor(MEs_Y)
#cluster module eigengenes
METree_Y=hclust(as.dist(MEDiss_Y),method="average")
sizeGrWindow(7,6)
plot(METree_Y, main="Clustering of module eigengenes") # choose cut for a correlation of 0.75

MEDissThres=0.05
abline(h=MEDissThres,col="red")
merge= mergeCloseModules(t(FOUR_VARS_YUCA_S),dynamicColors_Y,cutHeight = MEDissThres,verbose=3)
mergedColors=merge$colors
mergedMEs_Y=merge$newMEs
# to see effect of merging
sizeGrWindow(12,9)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Module_definition_CASSAVA.pdf", width=10, height=8)
par(las=1)
plotDendroAndColors(geneTree_Y,mergedColors,"Merged modules",dendroLabels = F,hang=0.03,addGuide = T,guideHang = 0.05,
                    cex.colorLabels = 2, cex.dendroLabels = 2, 
                    cex.rowText = 2)
dev.off()
#change names (to work better)
moduleColors_Y=mergedColors
colorOrder=c("grey",standardColors(50))
moduleLabels_Y=match(moduleColors_Y,colorOrder)-1
MEs_Y=mergedMEs_Y
table(moduleColors_Y)
##############################################################################
#other method of clustering: adjacency # do not perform better than DOM
##############################################################################
#hierADJ=hclust(as.dist(1-adjacency_Y), method="average" )
#dynmicMods_Y= cutreeDynamic(dendro=hierADJ, distM=1-adjacency_Y,deepSplit=2,pamRespectsDendro=F,minClusterSize=minModuleSize)
#table(dynmicMods_Y)
#dynamicColors_Y=labels2colors(dynmicMods_Y)
#table(dynamicColors_Y)
#plotDendroAndColors(hierADJ, colors = dynamicColors_Y, dendroLabels = FALSE, hang = 0.03,
#                    main = "Gene hierarchical clustering dendrogram and simulated module colors" ) 


###########################################################################################################################
# PREPARE DATA YUCA FOR quantifying module-module associations 

#define nb genes and samples
nGenes=ncol(t(FOUR_VARS_YUCA_S))
nSamples=nrow(t(FOUR_VARS_YUCA_S))
#recalculate MEs with color labels
MEs0_Y=moduleEigengenes(t(FOUR_VARS_YUCA_S),moduleColors_Y)$eigengenes
MEs_Y=orderMEs(MEs0_Y)
names(MEs_Y)<-paste("Y",names(MEs_Y),sep="_")
rownames(MEs_Y)<-rownames(t(FOUR_VARS_YUCA_S))



###########################################################################################################################
###########################################################################################################################
###################################
#### R irregularis gene modules####
###################################

## create eigene modules per sample. 
#create set of powers
powers=c(c(1:10),seq(from=12, to=40,by=2))
#network topology analysis fx
sft_A=pickSoftThreshold(t(FOUR_VARS_AMF_S),powerVector = powers,verbose = 5)
par(mfrow=c(1,2))
plot(sft_A$fitIndices[,1], -sign(sft_A$fitIndices[,3])*sft_A$fitIndices[,2],xlab="SoftThresholdpower)",ylab="ScaleFreeTopologyModelFit,signed^2",main =paste("Scaleindependence"),type="n")
text(sft_A$fitIndices[,1], -sign(sft_A$fitIndices[,3])*sft_A$fitIndices[,2],labels=powers)

plot(sft_A$fitIndices[,1], sft_A$fitIndices[,5],type="n")
text(sft_A$fitIndices[,1], sft_A$fitIndices[,5],labels=powers)

####################################################
## clustering methods**
####################################################
#adjacency
adjacency_A = adjacency(t(FOUR_VARS_AMF_S),power= 6); # change power
#topological overlap dissimilarity matrix 
dissTOM_A=1-TOMsimilarity(adjacency_A)

#Clustering using TOM
geneTree_A= hclust(as.dist(dissTOM_A),method="average")
#sizeGrWindow(12,9)
#plot(geneTree_A,main="Gene clustering based on TOM dissimilarity",labels=F,hang=0.04)

####################################################
#definition of module size.
####################################################
minModuleSize=30
#module identification using dynamic tree cut
dynmicMods_A= cutreeDynamic(dendro=geneTree_A, distM=dissTOM_A,deepSplit=2,pamRespectsDendro=F,minClusterSize=minModuleSize)
dynamicColors_A=labels2colors(dynmicMods_A)
table(dynamicColors_A)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree_A,dynamicColors_A,"Dynamic Tree Cut",dendroLabels= F,hang=0.03,addGuide = T,guideHang = 0.05,main="Gene dendrogram and module colors")

#merging modules whose expression profiles are similar
#calculate eignegens
MEList_AMF=moduleEigengenes(t(FOUR_VARS_AMF_S),colors=dynamicColors_A)
MEs_AMF=MEList_AMF$eigengenes
#calculate dissimilarity od module eigengenes
MEDiss_AMF= 1 -cor(MEs_AMF)
#cluster module eigengenes
METree_AMF=hclust(as.dist(MEDiss_AMF),method="average")
sizeGrWindow(7,6)
plot(METree_AMF, main="Clustering of module eigengenes") # choose cut for a correlation of 0.75

MEDissThres=0.05
abline(h=MEDissThres,col="red")
merge_AMF= mergeCloseModules(t(FOUR_VARS_AMF_S),dynamicColors_A,cutHeight = MEDissThres,verbose=3)
mergedColors_AMF=merge_AMF$colors
mergedMEs_AMF=merge_AMF$newMEs
# to see effect of merging
sizeGrWindow(12,9)
plotDendroAndColors(geneTree_A,cbind(dynamicColors_A,mergedColors_AMF),c("Dynamic Tree cut","Merged Dynamic"),dendroLabels = F,hang=0.03,addGuide = T,guideHang = 0.05)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Module_definition_AMF.pdf", width=14, height=10)
par(las=1)
plotDendroAndColors(geneTree_A,mergedColors_AMF,"Merged modules",dendroLabels = F,hang=0.03,addGuide = T,guideHang = 0.05,
                   cex.dendroLabels = 2, cex.colorLabels = 2,
                    cex.rowText = 2)
dev.off()


#change names (to work better)
moduleColors_AMF=mergedColors_AMF
colorOrder=c("grey",standardColors(50))
moduleLabels_AMF=match(moduleColors_AMF,colorOrder)-1
MEs_AMF=mergedMEs_AMF

###########################################################################################################################
# PREPARE DATA AMF FOR quantifying module-module associations 

#define nb genes and samples
nGenes=ncol(t(FOUR_VARS_AMF_S))
nSamples=nrow(t(FOUR_VARS_AMF_S))
#recalculate MEs with color labels
MEs0_AMF=moduleEigengenes(t(FOUR_VARS_AMF_S),moduleColors_AMF)$eigengenes

MEs_AMF=orderMEs(MEs0_AMF)
names(MEs_AMF)<-paste("A",names(MEs_AMF),sep="_")
rownames(MEs_AMF)<-rownames(t(FOUR_VARS_AMF_S))
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

#Modules expression in samples
##########################################################################################################################################################
# EIGENGENES
# MOdules expression in different treatments.
# Mesculenta
########??RELATE MODULE EIGENGENES APPARTENSHIP TO TREATMENTS.
modules_YUCA<-MEs_Y

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/module_eigen_sample_heatmap_CASSAVA.pdf", width=20, height=20)
#rownames(modules_YUCA)<-rownames(t(FOUR_VARS_YUCA_topheno))
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
heatmap.2( as.matrix(modules_YUCA), col=colors,scale="none",dendrogram = "both",cexRow =3,cexCol =3, key = F,
           margins=c(23,23), tracecol=NA )
dev.off()


modules_AMF<-MEs_AMF

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/module_eigen_sample_heatmap_AMF.pdf", width=20, height=20)
#rownames(modules_AMF)<-rownames(t(FOUR_VARS_AMF_topheno))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( as.matrix(modules_AMF[,-7]), col=colors,scale="none",dendrogram = "both",cexRow =3,cexCol =3, key = F,
           margins=c(23,23) , tracecol=NA)
dev.off()


################################################################################################

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
################ CORRELATION of MODULES TO PLANT GROWTH
###############################################################################################
# only samples on RNAseq
onlyRNA<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/Phenodata_vf/Harvest_COINOC.txt",h=T)
#onlyRNA<-onlyRNA[onlyRNA$CULTIVAR!="V6_CM4574-7",]
onlyRNA<-onlyRNA[onlyRNA$CULTIVAR!="V8_CM523-7",]

onlyRNA$AMF<-relevel(onlyRNA$AMF, ref = "CTRL") #order factors to put control as reference
onlyRNA<-droplevels(onlyRNA)
all_root<-onlyRNA$roots+onlyRNA$tuber
onlyRNA<-cbind(onlyRNA,all_root)
onlyRNA
onlyRNA2pheno_AMF<-onlyRNA[c(12,4,5,3,   29,28,  6,2,7:8,1,11,18,15,19,16,13,17,14, 31,26 , 25, 30 , 27 ),]

# PLOT CLUSTERING PHENOTYPE REALTIONSHIP
#YUCA
sampleTree2 = hclust(dist(t(FOUR_VARS_YUCA_S)))
plot(sampleTree2)            
traitColors = numbers2colors(onlyRNA2pheno_AMF[,c(11:17)], signed = FALSE);
par(mfrow=c(1,2))
plotDendroAndColors(sampleTree2, traitColors, groupLabels= c("Aboveground DW","Fine roots DW","Tuberized roots DW","Total DW","Underground DW","Undergorund/aboveground ratio","Height"),
                    main = "Sample dendrogram and trait heatmap\n Mesculenta")
#AMF
sampleTree2 = hclust(dist(t(FOUR_VARS_AMF_S)))
plot(sampleTree2)            
traitColors = numbers2colors(onlyRNA2pheno_AMF[,c(11:17)], signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors, groupLabels= c("Aboveground DW","Fine roots DW","Tuberized roots DW","Total DW","Underground DW","Undergorund/aboveground ratio","Height"),
                    main = "Sample dendrogram and trait heatmap\n Rirregularis")

######################################################################################




###########################################################################################################################
#####correlation Modules eigengenes YUCA, phenotype
########################################################
onlyRNA2pheno_AMF<-onlyRNA[c(12,4,5,3,   29,28,  6,2,7:8,1,11,18,15,19,16,13,17,14, 31,26 , 25, 30 , 27 ),]
table(dynamicColors_Y)
### Correlation
moduleTraitCor= cor(MEs_Y,onlyRNA2pheno_AMF[,c(11:17)],use="p")

moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

#Representation
sizeGrWindow(10,8)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",
                 signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))
#disply correlation values within heatmap plot
#X-axis co-expresion modules on AMF, Y-axis co-expresion modules CASSAVA
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(onlyRNA2pheno_AMF[,c(11:17)]),
               yLabels = names(MEs_Y),
               ySymbols = names(MEs_Y),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,cex.lab=1,
               zlim = c(-1,1),
               main = paste("Module CASSAVA-phenotype relationships"))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Modules_CASSAVA_pheno_COMP.pdf", width=20, height=18)
par(mar=c(10,16,3,1))
col_HM<-c(colorRampPalette( rev(brewer.pal(9, "Greens")) )(20),rep("white",10),colorRampPalette( brewer.pal(9, "Blues") )(20)  )
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels =  c("Height","Aboveground DW","Fine roots DW","Tuberized roots DW","Total DW","Underground DW","Undergorund/aboveground ratio"),
               yLabels = names(MEs_Y),
               ySymbols = names(MEs_Y),
               colorLabels = FALSE,
               colors = col_HM,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 2,cex.lab=2,
               zlim = c(-1,1),
               main = paste("Module CASSAVA-module AMF relationships"))
dev.off()


table(dynamicColors_A)


Cor_MOD_CASSAVA_2_CASSAVA<-moduleTraitPvalue[moduleTraitPvalue<0.02]
names(Cor_MOD_CASSAVA_2_CASSAVA)<-c("X125dai Y_MEpurple","X125dai Y_MEgreenyellow","X125dai Y_MEturquoise","X125dai Y_MEmagenta","X125dai Y_MEred",
                                    "roots Y_MEpurple","roots Y_MEgreenyellow","roots Y_MEred",
                                "tuber Y_MEpurple","tuber Y_MEgreen","tuber Y_MEred","tuber Y_MEyellow",
                                "total Y_MEpurple","total Y_MEgreen","total Y_MEred","total Y_MEyellow",
                                "undergrown Y_MEpurple","undergrown Y_MEgreen","undergrown Y_MEred","undergrown Y_MEyellow",
                                "ratio Y_MEpurple","ratio Y_MEgreen","ratio Y_MEred","ratio Y_MEyellow")


Cassava_sign_TREAT<-rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])

sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Cassava_sign_TREAT,]
#take out controls of CASSAVA
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]

phen<-onlyRNA2pheno_AMF[,c(11:17)]
top_gene_sign_intramod_conec_Y<-list()
Import_gene_M_Y<-list()

for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_CASSAVA))  {
  
  gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  phen[,names(phen)==sapply(strsplit(correlated_Modula," "), "[[", 1)]
  
  
  geneModuleMembership_CASSAVA=as.data.frame(cor(t(sign_norm_genes),MEs_Y,use="p"))
  geneTraitSignificance_CASSAVA= as.data.frame(cor(t(sign_norm_genes),
                                               phen[,names(phen)==sapply(strsplit(correlated_Modula," "), "[[", 1)],use="p")) # define correlation with external variable
  
  modNames=substring(names(MEs_Y),5)
  module=gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))
  column= match(module,modNames)
  moduleGenes=moduleColors_Y==module
  ## intramodular conectivity
  ADJ1=abs(cor(t(sign_norm_genes),use="p"))^6
  Alldegrees1=intramodularConnectivity(ADJ1, moduleColors_Y)
  head(Alldegrees1)
  ncol(ADJ1)
  GS1=as.numeric(cor(phen[,names(phen)==sapply(strsplit(correlated_Modula," "), "[[", 1)],t(sign_norm_genes), use="p"))
  GeneSignificance=abs(GS1)
  #GS.6_Y<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0]
  # quantile(abs(geneTraitSignificance_YUCA[moduleGenes,1]),.85)
  table( abs(GS1[moduleGenes])> quantile(abs(GS1[moduleGenes]),.9))
  colorlevels=unique(dynamicColors_Y)
  sizeGrWindow(9,6)
  par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
  par(mar = c(4,5,3,1))
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor2pheno/4.ALL_MODULES_SIGN',correlated_Modula,'.pdf',sep="_"),width=14,height=10,useDingbats = FALSE)
  par(mar=c(8,9,5,2),las=1,mgp=c(6, 2,0))
  
  whichmodule=gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2));
  restrict1 = (dynamicColors_Y==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=dynamicColors_Y[restrict1],
                     main=whichmodule,pch=16,cex=3, cex.axis=3,cex.lab = 3,cex.main = 3,
                     xlab = paste("Connectivity",whichmodule), ylab = paste(paste("Y",whichmodule,sep="_"),"correlation","to",sapply(strsplit(correlated_Modula," "), "[[", 1),sep=" "), abline = TRUE)
  
  dev.off()
  
  
  ####### # 
  #redo filter but with quantile and not absolute value
  #module black
  ##########
  datKME=signedKME(t(sign_norm_genes), MEs_Y, outputColumnName="MM.")
  # FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
  FilterGenes=  abs(GS1)[moduleGenes]> quantile(abs(GS1)[moduleGenes],.8) & abs(datKME[moduleGenes,grep(gsub("Y_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 2)),colnames(datKME))]) >.8
  
  top_sig_connect_high<-cbind.data.frame(rownames(datKME[moduleGenes,][FilterGenes,]),datKME[moduleGenes,grep(gsub("Y_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 2)),colnames(datKME))][FilterGenes])
  as.vector(top_sig_connect_high[,1])
  
  table(FilterGenes)
  
  top_gene_sign_intramod_conec_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-as.vector(top_sig_connect_high[,1]) #dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
  
  Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-sign_YUCA[rownames(sign_YUCA) %in%   as.vector(unlist(top_gene_sign_intramod_conec_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])),]
  Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-cbind.data.frame(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],rownames(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]))
  colnames(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[dim(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[2]]<-"gene"
  Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-merge(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],Mercator_Mesculenta,by="gene")
  
  toprint<-Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]
  write.table( toprint[,c(1:10,14,16:17)],
               paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor2pheno/4.TOP_GENES_FUNCTION',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_"),
               ,quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  summary( toprint)
  
  paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")
  
  ################
  
  MEList_CASSAVA=moduleEigengenes(t(sign_norm_genes),colors=dynamicColors_Y)
  exp_and_modules<-cbind.data.frame(sign_norm_genes,MEList_CASSAVA$validColors)
  colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
  # extract top 20 genes
  exp_and_modules[exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2)),]
  cor_to_module<-as.data.frame(cor(t(sign_norm_genes),
                                   phen[,names(phen)==sapply(strsplit(correlated_Modula," "), "[[", 1)]
                                   ,use="p"))
  top_correlated<-cbind.data.frame(exp_and_modules,cor_to_module)
  colnames(top_correlated)[dim(top_correlated)[2]]<-"CORRELADO"
  
  top_correlated2<-top_correlated[ top_correlated$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
  
  
  #### all markers + p-value for MapMan
  COR_info_module<-cbind.data.frame(rownames(top_correlated2),top_correlated2[,dim(top_correlated2)[2]])
  colnames(COR_info_module)<-c("gene",paste("cor2",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="_"))
  FC_CASSAVA<-cbind.data.frame(rownames(sign_YUCA),sign_YUCA)
  colnames(FC_CASSAVA)[1]<-'gene'
  COR_info_module<-merge(merge(COR_info_module,Mercator_Mesculenta,by="gene"),FC_CASSAVA,by='gene')
  
  write.table(  COR_info_module,
                paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor2pheno/4.GENES_2mapman',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  
  #extract all genes in module
  exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
  
  Genes_module<-t(sign_norm_genes)[,colnames(t(sign_norm_genes)) %in%  rownames(exp_mod2)]
  
  write.table(t(Genes_module), paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor2pheno/4.ALL_GENES_CROSSTALKS',correlated_Modula ,'.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = T, row.names = T)
}

names(Import_gene_M_Y)

################################################################################################
# visualising important genes
All_TOP_genes<-rbind.data.frame(Import_gene_M_Y[[1]],Import_gene_M_Y[[2]],Import_gene_M_Y[[3]],
                                Import_gene_M_Y[[4]],Import_gene_M_Y[[5]],Import_gene_M_Y[[6]],
                                Import_gene_M_Y[[7]],Import_gene_M_Y[[8]],Import_gene_M_Y[[9]],
                                Import_gene_M_Y[[10]],Import_gene_M_Y[[11]],Import_gene_M_Y[[12]],
                                Import_gene_M_Y[[13]],Import_gene_M_Y[[14]],Import_gene_M_Y[[15]],
                                Import_gene_M_Y[[16]],Import_gene_M_Y[[17]],Import_gene_M_Y[[18]],
                                Import_gene_M_Y[[19]],Import_gene_M_Y[[20]],Import_gene_M_Y[[21]],
                                Import_gene_M_Y[[22]],Import_gene_M_Y[[23]],Import_gene_M_Y[[24]])
sort(levels(All_TOP_genes$gene))
All_TOP_genes2<-All_TOP_genes[!duplicated(All_TOP_genes$gene),]
All_TOP_genes3<-All_TOP_genes2
#All_TOP_genes3<-All_TOP_genes2[All_TOP_genes2$Bin!="35.2",]
# individual genes
write.table(All_TOP_genes3[,c(1:10,14,15,16,17)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Genes_summary_CASSAVA_sample.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)

NORM_EXP_TOPCASSAVA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  All_TOP_genes3$gene,])
NORM_EXP_TOPCASSAVA<-NORM_EXP_TOPCASSAVA[grep("CTRL",rownames(NORM_EXP_TOPCASSAVA),invert=T),]
c("g212.t1","g3155.t1","g409.t1","g4691.t1","g532.t1","g6073.t1","g9335.t1")
c("g1884.t1", "g1139.t1", "g2391.t1","g10403.t1","g7542.t1")

correlacion<-cor(NORM_EXP_TOPCASSAVA,phen[,c(1,4:7)])

colsc=c(rgb(241, 54, 23, maxColorValue=255),  "white" , rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space= 'Lab' )
colors = colramp(100)
labo<-data.frame(NORM_EXP_TOPCASSAVA,phen[,c(1,4:7)])
labo
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Correlogram_genes1_CASSAVA.pdf',width=16,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
#labo<-labo[,c(2,1,3:5,7,10,12,6,8:9,11)]
#my.plotcorr(cex =2,cex.lab = 2,cor(labo)[(dim(NORM_EXP_TOPCASSAVA)[2]+1):(dim(NORM_EXP_TOPCASSAVA)[2]+5),1:dim(NORM_EXP_TOPCASSAVA)[2]],
 #           col=colors[((cor(labo)[(dim(NORM_EXP_TOPCASSAVA)[2]+1):(dim(NORM_EXP_TOPCASSAVA)[2]+5),1:dim(NORM_EXP_TOPCASSAVA)[2]] + 1)/2) * 100], diag= 'ellipse' ,  mar=c(0,2,0,0) )
pheatmap(cor(labo)[(dim(NORM_EXP_TOPCASSAVA)[2]+1):(dim(NORM_EXP_TOPCASSAVA)[2]+5),1:dim(NORM_EXP_TOPCASSAVA)[2]],
         cellwidth = 10,cellheight = 10)
dev.off()

tratamiento<-unlist(lapply(strsplit(rownames(NORM_EXP_TOPCASSAVA),"_"),function (x) x[3]))
cultivar<-gsub("m","",gsub("m2","",unlist(lapply(strsplit(rownames(NORM_EXP_TOPCASSAVA),"_"),function (x) x[1]))   ))

#test gene cAssava wiht control
boxplot(NORM_EXP_TOPCASSAVA[,10]~interaction(tratamiento,cultivar),las=2)
plot(NORM_EXP_TOPCASSAVA[,10],onlyRNA[c(12,4,5,3,   29,28,  6,2,7:8,1,11,18,15,19,16,13,17,14, 31,26 , 25, 30 , 27 ),15],las=2)

#plot(NORM_EXP_TOPCASSAVA[,10],t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.12G031900",]),las=2)
pt1<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.12G031900",])
#boxplot(as.vector(pt1)~interaction(tratamiento,cultivar),las=2)
#plot(as.vector(pt1),onlyRNA[c(12,4,5,3,   29,28,  6,2,7:8,1,11,18,15,19,16,13,17,14, 31,26 , 25, 30 , 27 ),15],las=2)
#


gene_impo<-cbind.data.frame(NORM_EXP_TOPCASSAVA,tratamiento,as.factor(cultivar),phen$total)
colnames(gene_impo)[c((dim(gene_impo)[2]-2),(dim(gene_impo)[2]-1),dim(gene_impo)[2])]<-c("tratamiento","cultivar","totalDW")
sort(tapply(gene_impo[,1],gene_impo[,8],mean))
sort(tapply(gene_impo[,1],gene_impo[,8],sd))


pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/all_topgenes_CASSAVA.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0),mfrow=c(3,4))
for (ge in 1:(dim(gene_impo)[2]-3)){
  plot(gene_impo[,ge],phen$total,pch=16,cex=3,ylab="Total Dry weight",xlab=colnames(gene_impo)[ge],cex.lab=3,cex.axis=3)
  abline(lm(phen$total~gene_impo[,ge]))
  #dev.off()
  
  #pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/g212-bySamples_CASSAVA.pdf',width=14,height=10,useDingbats = FALSE)
  #par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
  boxplot(gene_impo[,ge]~interaction(gene_impo$tratamiento,gene_impo$cultivar),las=2,
          ylab="log Gene expression",cex.lab=3,cex.axis=3,col=c("dodgerblue3","lightsalmon3","darkmagenta"))
  
  
}
dev.off()
# single genes
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/topgenes_ubiquitin_AMF.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(15,11,5,2),las=1,mgp=c(5.5, 2,0),mfrow=c(1,1))
ge<-3
colnames(gene_impo)
  plot(gene_impo[,ge],phen$X125dai,pch=16,cex=3,ylab="Height",xlab=colnames(gene_impo)[ge],cex.lab=3,cex.axis=3)
  abline(lm(phen$X125dai~gene_impo[,ge]))
  summary(lm(phen$X125dai~gene_impo[,ge]))
  #dev.off()
  
  #pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/g212-bySamples_CASSAVA.pdf',width=14,height=10,useDingbats = FALSE)
  #par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
  #boxplot(gene_impo[,ge]~interaction(gene_impo$tratamiento,gene_impo$cultivar),las=2,
   #       ylab="log Gene expression",cex.lab=1,cex.axis=1,col=c("dodgerblue3","lightsalmon3","darkmagenta"))

dev.off()

##############################
# extract genes names per module
Cassava_sign_TREAT<-rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])
sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Cassava_sign_TREAT,]
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]

gene_names_YUCA<-cbind.data.frame(sign_norm_genes,moduleColors_Y)
levels(gene_names_YUCA$moduleColors_Y)
rownames(gene_names_YUCA[gene_names_YUCA$moduleColors_Y=="black",])
######################################################################################

######################################################################################
###########################################################################################################################
#####correlation Modules eigengenes AMF, phenotype
########################################################
onlyRNA2pheno_AMF<-onlyRNA[c(12,4,5,3,   29,28,  6,2,7:8,1,11,18,15,19,16,13,17,14, 31,26 , 25, 30 , 27 ),]
table(dynamicColors_A)
### Correlation

moduleTraitCor= cor(MEs_AMF,onlyRNA2pheno_AMF[,c(11:17)],use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

#Representation
sizeGrWindow(10,8)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",
                 signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))
#disply correlation values within heatmap plot
#X-axis co-expresion modules on AMF, Y-axis co-expresion modules CASSAVA
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(onlyRNA2pheno_AMF[,c(11:17)]),
               yLabels = names(MEs_AMF),
               ySymbols = names(MEs_AMF),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,cex.lab=1,
               zlim = c(-1,1),
               main = paste("Module AMF-phenotype relationships"))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Modules_AMF_pheno_COMP.pdf", width=14, height=10)
par(mar=c(10,16,3,1))
col_HM<-c(colorRampPalette( rev(brewer.pal(9, "Greens")) )(20),rep("white",10),colorRampPalette( brewer.pal(9, "Blues") )(20)  )
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels =  c("Height","Aboveground DW","Fine roots DW","Tuberized roots DW","Total DW","Underground DW","Undergorund/aboveground ratio"),
               yLabels = names(MEs_AMF),
               ySymbols = names(MEs_AMF),
               colorLabels = FALSE,
               colors = col_HM,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 2,cex.lab=2,
               zlim = c(-1,1),
               main = paste("Module AMF-module AMF relationships"))
dev.off()


table(dynamicColors_A)


Cor_MOD_CASSAVA_2_AMF<-moduleTraitPvalue[moduleTraitPvalue<0.02]
names(Cor_MOD_CASSAVA_2_AMF)<-c("X125dai A_MEyellow","X125dai A_MEred","X125dai A_MEbrown",
                                "aerial A_MEred",
                                "roots A_MEyellow","roots A_MEred")

#names(Cor_MOD_CASSAVA_2_AMF)<-c("X125dai A_MEgreen","X125dai A_MEbrown","X125dai A_MEyellow","X125dai A_MEred",
 #                               "aerial A_MEred","roots A_MEgreen","roots A_MEbrown",
  #                              "tuber A_MEgrey","total A_MEgrey","undergrown A_MEgrey",
   #                             "ratio A_MEgrey")


AMF_sign_TREAT<-rownames(sign_AMF[sign_AMF$adj.P.Val<0.05,])

sign_norm_genes<-DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_sign_TREAT,]
#take out controls of CASSAVA
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]

phen<-onlyRNA2pheno_AMF[,c(11:17)]
top_gene_sign_intramod_conec_A<-list()
Import_gene_M_A<-list()

for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_AMF))  {
  
  gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  phen[,names(phen)==sapply(strsplit(correlated_Modula," "), "[[", 1)]
  
  
    geneModuleMembership_AMF=as.data.frame(cor(t(sign_norm_genes),MEs_AMF,use="p"))
  geneTraitSignificance_AMF= as.data.frame(cor(t(sign_norm_genes),
                                               phen[,names(phen)==sapply(strsplit(correlated_Modula," "), "[[", 1)],use="p")) # define correlation with external variable
 
    modNames=substring(names(MEs_AMF),5)
  module=gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))
  column= match(module,modNames)
  moduleGenes=moduleColors_AMF==module
    ## intramodular conectivity
  ADJ1=abs(cor(t(sign_norm_genes),use="p"))^6
  Alldegrees1=intramodularConnectivity(ADJ1, moduleColors_AMF)
  head(Alldegrees1)
  ncol(ADJ1)
  GS1=as.numeric(cor(phen[,names(phen)==sapply(strsplit(correlated_Modula," "), "[[", 1)],t(sign_norm_genes), use="p"))
  GeneSignificance=abs(GS1)
  #GS.6_Y<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0]
  # quantile(abs(geneTraitSignificance_YUCA[moduleGenes,1]),.85)
  table( abs(GS1[moduleGenes])> quantile(abs(GS1[moduleGenes]),.9))
  colorlevels=unique(dynamicColors_A)
  sizeGrWindow(9,6)
  par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
  par(mar = c(4,5,3,1))
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor2pheno/4.ALL_MODULES_SIGN',correlated_Modula,'.pdf',sep="_"),width=14,height=10,useDingbats = FALSE)
  par(mar=c(8,9,5,2),las=1,mgp=c(6, 2,0))
  
  whichmodule=gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2));
  restrict1 = (dynamicColors_A==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=dynamicColors_A[restrict1],
                     main=whichmodule,pch=16,cex=3, cex.axis=3,cex.lab = 3,cex.main = 3,
                     xlab = paste("Connectivity",whichmodule), ylab = paste(paste("A",whichmodule,sep="_"),"correlation","to",sapply(strsplit(correlated_Modula," "), "[[", 1),sep=" "), abline = TRUE)
  
  dev.off()
  
  
  ####### # 
  #redo filter but with quantile and not absolute value
  #module black
  ##########
  datKME=signedKME(t(sign_norm_genes), MEs_AMF, outputColumnName="MM.")
  # FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
  FilterGenes=  abs(GS1)[moduleGenes]> quantile(abs(GS1)[moduleGenes],.8) & abs(datKME[moduleGenes,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 2)),colnames(datKME))]) >.8
  
  top_sig_connect_high<-cbind.data.frame(rownames(datKME[moduleGenes,][FilterGenes,]),datKME[moduleGenes,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 2)),colnames(datKME))][FilterGenes])
  as.vector(top_sig_connect_high[,1])
  
  table(FilterGenes)
  
  top_gene_sign_intramod_conec_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-as.vector(top_sig_connect_high[,1]) #dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
  
  Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-sign_AMF[rownames(sign_AMF) %in%   as.vector(unlist(top_gene_sign_intramod_conec_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])),]
  Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-cbind.data.frame(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],rownames(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]))
  colnames(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[dim(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[2]]<-"gene"
  Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-merge(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],Mercator_Rirregularis,by="gene")
  
  toprint<-Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]
  write.table( toprint[,c(1:7,11,13:14)],
               paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor2pheno/4.TOP_GENES_FUNCTION',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_"),
             ,quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  summary( toprint)
  
  paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")
  
  ################
  
  MEList_AMF=moduleEigengenes(t(sign_norm_genes),colors=dynamicColors_A)
  exp_and_modules<-cbind.data.frame(sign_norm_genes,MEList_AMF$validColors)
  colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
  # extract top 20 genes
  exp_and_modules[exp_and_modules$MODULE==gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2)),]
  cor_to_module<-as.data.frame(cor(t(sign_norm_genes),
                                   phen[,names(phen)==sapply(strsplit(correlated_Modula," "), "[[", 1)]
                                   ,use="p"))
  top_correlated<-cbind.data.frame(exp_and_modules,cor_to_module)
  colnames(top_correlated)[dim(top_correlated)[2]]<-"CORRELADO"
  
  top_correlated2<-top_correlated[ top_correlated$MODULE==gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
  
  
  #### all markers + p-value for MapMan
  COR_info_module<-cbind.data.frame(rownames(top_correlated2),top_correlated2[,dim(top_correlated2)[2]])
  colnames(COR_info_module)<-c("gene",paste("cor2",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="_"))
  FC_AMF<-cbind.data.frame(rownames(sign_AMF),sign_AMF)
  colnames(FC_AMF)[1]<-'gene'
  COR_info_module<-merge(merge(COR_info_module,Mercator_Rirregularis,by="gene"),FC_AMF,by='gene')
  
  write.table(  COR_info_module,
                paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor2pheno/4.GENES_2mapman',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  
  #extract all genes in module
  exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
  
  Genes_module<-t(sign_norm_genes)[,colnames(t(sign_norm_genes)) %in%  rownames(exp_mod2)]
  
  write.table(t(Genes_module), paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Cor2pheno/4.ALL_GENES_CROSSTALKS',correlated_Modula ,'.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = T, row.names = T)
}

names(Import_gene_M_A)

################################################################################################
# visualising important genes
All_TOP_genes<-rbind.data.frame(Import_gene_M_A[[1]],Import_gene_M_A[[2]],Import_gene_M_A[[3]],
Import_gene_M_A[[4]],Import_gene_M_A[[5]],Import_gene_M_A[[6]])
levels(All_TOP_genes$gene)
All_TOP_genes2<-All_TOP_genes[!duplicated(All_TOP_genes$gene),]
All_TOP_genes3<-All_TOP_genes2[All_TOP_genes2$Bin!="35.2",]
# individual genes
All_TOP_genes3$gene

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  All_TOP_genes3$gene,])
c("g212.t1","g3155.t1","g409.t1","g4691.t1","g532.t1","g6073.t1","g9335.t1")
c("g1884.t1", "g1139.t1", "g2391.t1","g10403.t1","g7542.t1")

correlacion<-cor(NORM_EXP_TOPAMF,phen[,c(1,4:7)])

colsc=c(rgb(241, 54, 23, maxColorValue=255),  "white" , rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space= 'Lab' )
colors = colramp(100)
labo<-data.frame(NORM_EXP_TOPAMF,phen[,c(1,4:7)])
labo
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/Correlogram_genes1_AMF.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
#labo<-labo[,c(2,1,3:5,7,10,12,6,8:9,11)]
#my.plotcorr(cex =2,cex.lab = 2,cor(labo)[(dim(NORM_EXP_TOPAMF)[2]+1):(dim(NORM_EXP_TOPAMF)[2]+5),1:dim(NORM_EXP_TOPAMF)[2]],
 #           col=colors[((cor(labo)[(dim(NORM_EXP_TOPAMF)[2]+1):(dim(NORM_EXP_TOPAMF)[2]+5),1:dim(NORM_EXP_TOPAMF)[2]] + 1)/2) * 100], diag= 'ellipse' ,  mar=c(0,2,0,0) )
pheatmap(cor(labo)[(dim(NORM_EXP_TOPAMF)[2]+1):(dim(NORM_EXP_TOPAMF)[2]+5),1:dim(NORM_EXP_TOPAMF)[2]],
         cellwidth = 10,cellheight = 10)

dev.off()


tratamiento<-unlist(lapply(strsplit(rownames(NORM_EXP_TOPAMF),"_"),function (x) x[3]))
cultivar<-gsub("m","",gsub("m2","",unlist(lapply(strsplit(rownames(NORM_EXP_TOPAMF),"_"),function (x) x[1]))   ))

gene_impo<-cbind.data.frame(NORM_EXP_TOPAMF,tratamiento,as.factor(cultivar),phen$total)
colnames(gene_impo)[c((dim(gene_impo)[2]-2),(dim(gene_impo)[2]-1),dim(gene_impo)[2])]<-c("tratamiento","cultivar","totalDW")
sort(tapply(gene_impo[,1],gene_impo[,8],mean))
sort(tapply(gene_impo[,1],gene_impo[,8],sd))


pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/all_topgenes_AMF.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0),mfrow=c(3,4))
for (ge in 1:(dim(gene_impo)[2]-3)){
    plot(gene_impo[,ge],phen$total,pch=16,cex=3,ylab="Total Dry weight",xlab=colnames(gene_impo)[ge],cex.lab=3,cex.axis=3)
  abline(lm(phen$total~gene_impo[,ge]))
#dev.off()

#pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/g212-bySamples_AMF.pdf',width=14,height=10,useDingbats = FALSE)
#par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
  boxplot(gene_impo[,ge]~interaction(gene_impo$tratamiento,gene_impo$cultivar),las=2,
        ylab="log Gene expression",cex.lab=3,cex.axis=3,col=c("dodgerblue3","lightsalmon3","darkmagenta"))


}
dev.off()
################################################################################################
# figures WGCNA
table(dynamicColors_A)
table(dynamicColors_Y)



################################################################################################
#### CLASSIFY DEF into mercator classes and explain pattern
genes_yuca<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,]) ,]
genes_yuca<-cbind.data.frame(genes_yuca,rownames(genes_yuca))
colnames(genes_yuca)[dim(genes_yuca)[2]]<-"gene"
genes_yuca<-merge(genes_yuca,Mercator_Mesculenta,by="gene")
head(genes_yuca2)
unlist(lapply(lapply(strsplit(as.character(genes_yuca$Bin),"\\."),function (x) x[1:2]),function (x) paste(x[1],x[2],sep=".")))
genes_yuca2<-cbind.data.frame(genes_yuca[,c(1:33,35)],unlist(lapply(lapply(strsplit(as.character(genes_yuca$Bin),"\\."),function (x) x[1:2]),function (x) paste(x[1],x[2],sep=".")))  )
colnames(genes_yuca2)[dim(genes_yuca2)[2]]<-"classifier"
genes_yuca2<-genes_yuca2[,grep("CTRL",colnames(genes_yuca2),invert=T)]
genes_yuca2<-genes_yuca2[unique(genes_yuca2$gene),]
head(genes_yuca2)
sort(table(genes_yuca2[,27]))
genes_yuca2[genes_yuca2$classifier=="17.5",]
pheatmap(genes_yuca2[genes_yuca2$classifier=="20.1",2:25],labels_row =as.vector(genes_yuca2[genes_yuca2$classifier=="20.1",1]) )
write.table( Coinoc_effect_v5[,c(1:3,5:6)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/2.FUNCTION_CASSAVA_DIFF_CANB1-single-inoculation_v1.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)



################################################################################################
# GO TERMS IN MODULES


## CASSAVA
# extract genes names per module
Cassava_sign_TREAT<-rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])
sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Cassava_sign_TREAT,]
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]

gene_names_YUCA<-cbind.data.frame(sign_norm_genes,moduleColors_Y)
levels(gene_names_YUCA$moduleColors_Y)
rownames(gene_names_YUCA[gene_names_YUCA$moduleColors_Y=="black",])


for (labss in levels(gene_names_YUCA$moduleColors_Y) ) {
locusName<-rownames(sign_YUCA)
sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
sign_feat_go<-merge(sign_feat2,mesculenta_go2,by='locusName')
sign_feat_go<-na.omit(sign_feat_go)
gene_2_go<-sign_feat_go[,c(1,16)]
head(sign_feat_go)
head(gene_2_go)
write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
            col.names = F, row.names = F)
geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")

######################################################
# 3) SELECT UNIVERSE
######################################################
# TOTAL UNIVERSE CASSAVA
locusName<-rownames(sign_YUCA)
sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
sign_feat_go<-merge(sign_feat2,mesculenta_go2[,c(1,3)],by='locusName') # do not select transcriptname. working at gene level
sign_feat_go<-na.omit(sign_feat_go)
sign_feat_go<-unique(sign_feat_go)
length(unique(sign_feat_go$locusName))
######################################################
# 4) SELECT TOP GENES
######################################################

# modules
Significant_GENES_Y<-rownames(gene_names_YUCA[gene_names_YUCA$moduleColors_Y==labss,])

######################################################
# 5) ANNOTATION
######################################################
#AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
geneList <- factor(as.integer(sign_feat_go$locusName %in% Significant_GENES_Y))
names(geneList) <- sign_feat_go$locusName
sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5,ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)

sampleGOdata
resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)

allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                           classicKS = resultKS, elimKS = resultKS.elim,
                           orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)

RES_Y<-allRes
RES_Y

write.table(RES_Y[RES_Y[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/TEMP_module',labss,'yuca.txt',sep="_")
            ,quote=FALSE,sep='\t',col.names = F, row.names = F)

}

## AMF
# extract genes names per module
AMF_sign_TREAT<-rownames(sign_AMF[sign_AMF$adj.P.Val<0.05,])
sign_norm_genes<-DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_sign_TREAT,]
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]

gene_names_AMF<-cbind.data.frame(sign_norm_genes,moduleColors_AMF)
levels(gene_names_AMF$moduleColors_AMF)
rownames(gene_names_AMF[gene_names_AMF$moduleColors_AMF=="black",])


for (labss in levels(gene_names_AMF$moduleColors_AMF)[5:7] ) {
  #AMF
  locusName<-rownames(sign_AMF)
  sign_feat2<-cbind.data.frame(sign_AMF,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName')
  sign_feat_go<-na.omit(sign_feat_go)
  gene_2_go<-sign_feat_go[,c(1,12)]
  
  write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
              col.names = F, row.names = F)
  geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")
  
  
  ######################################################
  # 3) SELECT UNIVERSE
  ######################################################
  
  #TOTAL UNIVERSE ALL LOCI
  #TOTAL UNIVERSE AMF
  locusName<-rownames(sign_AMF)
  sign_feat2<-cbind.data.frame(sign_AMF,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName') # do not select transcriptname. working at gene level
  sign_feat_go<-na.omit(sign_feat_go)
  sign_feat_go<-unique(sign_feat_go)
  length(unique(sign_feat_go$locusName))
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  
  # modules
  Significant_GENES_Y<-rownames(gene_names_AMF[gene_names_AMF$moduleColors_AMF==labss,])
  
  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  geneList <- factor(as.integer(sign_feat_go$locusName %in% Significant_GENES_Y))
  names(geneList) <- sign_feat_go$locusName
  sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5,ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
  
  sampleGOdata
  resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
  resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
  resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
  
  allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                             classicKS = resultKS, elimKS = resultKS.elim,
                             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
  
  RES_Y<-allRes
  RES_Y
  
  write.table(RES_Y[RES_Y[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/TEMP_module',labss,'yuca.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = F, row.names = F)
  
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################################################################################################################
################################################################################################
################################################################################################
