################################################################
################################################################
# Descriptions and pre-processing of each NGS dataset type for the U937-PR9 experiments
# Each section describes the basic pre-processing steps used to create 
# the R objects used for all downstream analysis
# This script can be invoked in R using:
# source("./code/U937_PR9_Processed_Data_Descriptions.R") assuming all libraries are installed
################################################################
################################################################

#Load Packages required for the creation of all R objects for downstream analysis
suppressPackageStartupMessages({
	library(GenomicRanges)
	library(GenomicInteractions)
	library(edgeR)
	library(limma)
	library(GenomicRanges)
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	library(Organism.dplyr)
	library(ChIPpeakAnno)
	data(TSS.human.GRCh37)
	library(limma)
	library(org.Hs.eg.db)
	library(data.table)

})

# load custom functions
#source("./code/Functions_Paper_4_2_20.R")
graphs.dir <- "/graphs/"


################################################################
################################################################
# ATAC-seq
# 2 uninduced U937-PR9 datasets and 2 induced U937-PR9 datasets
# Bam files (trimmomatic-bowtie2(hg19)-blacklist-filtering) had peaks called using macs2 FDR <= 0.01:
# All 4 peak sets called from Macs2 were merged into a consesnsus set for diffbind differential calling using bedtools:
# consensus peak set: $ cat Control_8_FDR_0.01.bed Control_9_FDR_0.01.bed Zinc_8_FDR_0.01.bed Zinc_9_FDR_0.01.bed | bedtools merge -i stdin > ATAC_Merged_conesnsus_bedtools.bed
# ATAC_Merged_conesnsus_bedtools.bed was then used as the peak set to call differential regions 
# using R diffbind: with settings fragmentSize = 300, summits = FALSE, cutoff = 1, contrast = 'Condition'
# Processed bam files for each experiment were used for diffbind read counting
# The final diffbind output created Final_Combined_ATAC.txt 
# Final_Differential_ATAC.txt is the filtered Final_Combined_ATAC.txt based on FDR adjusted p-value
################################################################
################################################################

#read in processed ATAC-seq
Final_Differential_ATAC <- read.table("Final_Differential_ATAC.txt", sep = "\t", head =F)
colnames(Final_Differential_ATAC)[c(1,2,3)] <- c("chromosome","start","end")
Final_Differential_ATAC <- GRanges(Final_Differential_ATAC)

Final_Combined_ATAC <- read.csv("GSE173751_Final_Combined_ATAC-seq.csv.gz")
Final_Combined_ATAC[,1] <- NULL
Final_Combined_ATAC <- GRanges(Final_Combined_ATAC)

seqlevels(Final_Differential_ATAC,pruning.mode="coarse") <- c("chr1","chr2","chr3","chr4",
	"chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
	"chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX"
)

seqlevels(Final_Combined_ATAC,pruning.mode="coarse") <- c("chr1","chr2","chr3","chr4",
	"chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
	"chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX"
)

#Final_Combined_ATAC is the ATAC-seq object used for downstream analysis
cat("#########################","\n")
cat("loaded ATAC-seq","\n")
cat("#########################","\n")


################################################################
################################################################
# Cut&Run
# For peak calling, SEACR was used: https://github.com/FredHutch/SEACR
# Bam files (trimmomatic-bowtie2(hg19)-blacklist-filtering) were converted to bed file:
# $ for p in *bam ;do bamToBed -i $bam > $bam'.bed'
# convert to bedgraph format
# genomeCoverageBed -i $bam'.bed' -bg -g $genometextfile > $bam'.bed.bedgraph' done
# Run SEACR using downloaded SEACR_scripts from https://github.com/FredHutch/SEACR
# bash SEACR.sh experiment.bedgraph control.bedgraph norm union experiment.out
################################################################
################################################################

#load the merged induced U937-PR9 Cut&Run dataset called using the above SEACR code
Data <- read.table('Induced_Fusion_Merged.NU.out.auc.threshold.merge.bed', sep = "\t", head = F)

#convert to GRanges object
SEACR_Peaks <- GRanges(seqnames = Data$V1,
		ranges = IRanges(start = Data$V2, 
		end = Data$V3), q_value = Data$V6,
		summit = Data$V4, score = Data$V5		
)

seqlevels(SEACR_Peaks,pruning.mode="coarse") <- c("chr1","chr2","chr3","chr4",
	"chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
	"chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
	"chr21","chr22","chrX"
)

GetENTREZ <- function(LIST){
load('FINAL_U937_RNA-seq_PAIRED_END.Rdata')

	library(org.Hs.eg.db)
	hs <- org.Hs.eg.db
	my.symbols <- as.character(LIST)
	genes_ID <- select(hs, 
	keys = my.symbols,
    	columns = c("ENTREZID", "SYMBOL"),
       	keytype = "SYMBOL"
    )
    
    genes_ID <- genes_ID[which(genes_ID$ENTREZID %in% y$genes$GeneID),]
	return(genes_ID)
}

#assign each peak a gene by overlapping with 2500bp promoters generated from TxDb.Hsapiens.UCSC.hg19.knownGene
src <- src_organism("TxDb.Hsapiens.UCSC.hg19.knownGene")
keytypes(src)
columns(src)
keytypes(TxDb.Hsapiens.UCSC.hg19.knownGene)
columns(TxDb.Hsapiens.UCSC.hg19.knownGene)
txid = keys(TxDb.Hsapiens.UCSC.hg19.knownGene, "TXID")
df = select(TxDb.Hsapiens.UCSC.hg19.knownGene, txid, "GENEID", "TXID")
PR <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream=2000, downstream=500)

CONVERT <- select(org.Hs.eg.db, df$GENEID, c("SYMBOL", "GENENAME"))
df$Symbol <- CONVERT$SYMBOL[match(df$GENEID, CONVERT$ENTREZID)]

PR$SYMBOL <- df$Symbol[match(PR$tx_id, df$TXID)]
PR$ENTREZ <- df$GENEID[match(PR$tx_id, df$TXID)]

SEACR_Peaks$gene <- 'NA'
SEACR_Peaks[subjectHits(findOverlaps(PR,SEACR_Peaks))]$gene <- PR[queryHits(findOverlaps(PR,SEACR_Peaks))]$SYMBOL
SEACR_Peaks$ENTREZ <- 'NA'
SEACR_Peaks[subjectHits(findOverlaps(PR,SEACR_Peaks))]$ENTREZ <- PR[queryHits(findOverlaps(PR,SEACR_Peaks))]$ENTREZ

#SEACR_Peaks is the final PML-RARA induced set of 15,412 peaks
cat("#########################","\n")
cat("loaded Cut&Run","\n")
cat("#########################","\n")


################################################################
################################################################
# RNA-seq
# Fastq files were trimmed of adaptors and low quality reads with Trimmomatic
# Files were then aligned to the hg19 genome using STAR 2 Pass aligner: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
# Gene level counts per Aligned BAM files were generated using featureCounts() from the Rsubread package
# DGElist object was created using the counts generated and the hg19 genome annotation file: 
# > y <- DGEList(counts = Counts, genes = Ann)
# counts normalised using:
# reads normalised by log2 counts per million were extracted and saved as: GSE173754_Log2_Counts_Per_Million_All_RNA-seq.csv.gz
########
########
# Differential analysis 2 Induced vs 2 uninduced - removing genes with low numbers of reads
# > isexpr <- rowSums(cpm(y) > 1) >= (length(y$samples$group)/100*45)
# > y <- y[isexpr, , keep.lib.sizes = FALSE];
# > y <- calcNormFactors(y,method='TMM')
# > Test <- factor(c('control', 'control', 'zinc', 'zinc'))
# > design <- model.matrix(~Test)
# > v <- voom(y,design,plot=TRUE)
# > fit <- lmFit(v, design)
# > fit.de <- eBayes(fit, robust=TRUE)
# > topTable(fit.de, coef = 'Testzinc', n = Inf, sort = "p")[1:10,]
# > significant <- topTable(fit.de, coef = 'Testzinc', n = Inf, sort = "p")
# > significantDelta <- significant[which(significant$adj.P.Val <= 0.05),]
# > save(y, file = 'FINAL_U937_RNA-seq_PAIRED_END.Rdata')
# > write.table(significant,"FINAL_U937_RNA-seq_PAIRED_END_all_genes.txt", sep = "\t", row.names = F, quote = F)
# > write.table(significantDelta,"GSE173754_Control_vs_Induced_RNA-seq_Limma_DEG.txt.gz", sep = "\t", row.names = F, quote = F)
################################################################
################################################################

significant <- read.table("FINAL_U937_RNA-seq_PAIRED_END_all_genes.txt", sep = "\t", head = T)	
significantDelta <- significant[which(significant$adj.P.Val < 0.05),]
load('FINAL_U937_RNA-seq_PAIRED_END.Rdata')
y <- calcNormFactors(y,method='TMM')
v <- voom (y)	

Up <- significantDelta$GeneID[which(significantDelta$logFC > 0)] 
Down <- significantDelta$GeneID[which(significantDelta$logFC < 0)] 
UP <- significantDelta[which(significantDelta$logFC > 0),] 
DOWN <- significantDelta[which(significantDelta$logFC < 0),] 

#y is the RNA-seq object used for downstream analysis
#significantDelta is the DEG set used for downstream analysis
#Up,up,Down,down are the filtered DEG sets used for downstream analysis

cat("#########################","\n")
cat("loaded RNA-seq","\n")
cat("#########################","\n")

################################################################
################################################################
# Capture Hi-C significant interaction data
# fastq files were processed using the HiCUP pipeline: https://www.bioinformatics.babraham.ac.uk/projects/hicup/
# BAM files for all 3 replicates per condition were merged prior to interaction calling
# Merged BAM files had significant interactions called using 
# CHiCAGO using default settings with a significant interaction with a score >= 5
################################################################
################################################################

#Read in single significance Data 
cDIR <- list.files(pattern = "_Chicago_Score_5", recursive=T, full.names = T)
cDIR <- c(cDIR[grep("Control",cDIR)],cDIR[grep("Zinc",cDIR)])

chicago <- list()
chicago_sig <- list()

for(i in 1:length(cDIR)){
	chicago_sig[[i]] <- read.csv(cDIR[i])
	cat('read',i,'of',length(cDIR), "\n")
}

#load GRanges objects containing the bait coordinates and hindIII fragment coordinates
load("./Annotations/promoterBaitsHuman")
load("./Annotations/hindGR")
files <- list.files(pattern = '.chinput')
testDesignDir <- file.path("./Annotations/")
#get gen names for baitID
baitmap <- read.table(paste0(testDesignDir,"hindIII.baitmap"), sep = "\t", head = F)
hindmap <- read.table(paste0(testDesignDir,"hindIII.rmap"), sep = "\t", head = F)
colnames(hindmap)[c(1:3)] <- c('chromosome','start','end')
hindmap.GR <- GRanges(hindmap)

#create GenomicInteraction objects for downstream analysis
makeGIfromChicago <- function(data){
	
	#again the below code made mistakes
	#utilise how chicago was made orginally, using the unique code for each bait
	one <- hindmap[match(data$baitID, hindmap$V4),]
	two <- hindmap[match(data$otherEndID, hindmap$V4),]
	one$uniqueID <- paste0(one$V4,":",two$V4)
	
	colnames(one) <- c('chromosome','start','end','ID','uniqueID')
	colnames(two) <- c('chromosome','start','end','ID')
	one.GR <- GRanges(one)
	two.GR <- GRanges(two)
	
	GI <- GenomicInteractions(
		one.GR, two.GR)
	GI <- GI[order(GI$anchor1.uniqueID),]
	data$uniqueID <- paste0(data$baitID,":",data$otherEndID)
	data <- data[order(data$uniqueID),]
	
		GI$counts = data$N[match(GI$anchor1.uniqueID,data$uniqueID)]
		GI$score = data$score[match(GI$anchor1.uniqueID,data$uniqueID)]
		GI$bait1 = data$bait1[match(GI$anchor1.uniqueID,data$uniqueID)]
		GI$bait2 = data$bait2[match(GI$anchor1.uniqueID,data$uniqueID)]
		GI$oneID = data$baitID[match(GI$anchor1.uniqueID,data$uniqueID)]
		GI$twoID = data$otherEndID[match(GI$anchor1.uniqueID,data$uniqueID)]
	
	#sanity: GI stupidly changes the order of the data
	#SAMD coordinates 850619
	#chicago[[4]][which(chicago[[4]]$start2 == 850619)]
	#chicago_sig[[4]][which(start(anchorTwo(chicago_sig[[4]])) == 850619)]

	cat('taken and interaction bedpe and converted to GI', "\n")
	
	GI <- unique(GI)
	
	return(GI)
}

names(chicago_sig) <- c("Control_all","Zinc_all")
Control.GI <- makeGIfromChicago(chicago_sig[[1]])
Zinc.GI <- makeGIfromChicago(chicago_sig[[2]])


# get entrez gene IDs for baits and add to each GI object
Eq <- GetENTREZ(Control.GI$bait1)
Control.GI$ENTREZ_1 <- Eq$ENTREZ[match(Control.GI$bait1,Eq$SYMBOL)]
Eq <- GetENTREZ(Control.GI$bait2)
Control.GI$ENTREZ_2 <- Eq$ENTREZ[match(Control.GI$bait2,Eq$SYMBOL)]
Eq <- GetENTREZ(Zinc.GI$bait1)
Zinc.GI$ENTREZ_1 <- Eq$ENTREZ[match(Zinc.GI$bait1,Eq$SYMBOL)]
Eq <- GetENTREZ(Zinc.GI$bait2)
Zinc.GI$ENTREZ_2 <- Eq$ENTREZ[match(Zinc.GI$bait2,Eq$SYMBOL)]

#create vchicago containing characters for each interaction, used for overlap analysis
vchicago <- list()
for(i in 1:length(cDIR)){
	vchicago[[i]] <- unique(paste0(chicago_sig[[i]]$chr1,":",chicago_sig[[i]]$start1,"-", chicago_sig[[i]]$chr2,":",chicago_sig[[i]]$start2))
}

# Control.GI is the uninduced CHiC object used for downstream analysis
# Zinc.GI is the induced CHiC object used for downstream analysis
cat("#########################","\n")
cat("loaded CHiCAGO interactions","\n")
cat("#########################","\n")


################################################################
################################################################
# Capture Hi-C Differential interaction data
# Hicup BAM files had differential  interactions called using GOTHiC differential:
# https://github.com/biscuit13161/GOTHiC 
# Comparative GOTHiC was used to compare each replicate pair (induced rep1 vs uninduced rep1 etc)
# Settings using an ihw cut off of 0.01 were used
# Differential interactions were kept if an interaction was present in at least one other
# comparison and in the same direction, to create a consensus set for downstream analysis
################################################################
################################################################

#read in conesnsus sets
CONSENSUS <- read.csv('GSE173752_Consensus_Gained_Interactions_Zinc_vs_Control.csv.gz')
CONSENSUSDEC <- read.csv('GSE173752_Consensus_Lost_Interactions_Zinc_vs_Control.csv.gz')

#create GI object from consensus sets
makeGIfromGOTHiC_Differential <- function(data){
	
	data$ID <- 1:length(data[,1])
	one.GR <- GRanges(paste0(data$seqnames1,':',data$start1,'-',data$end1),ID=data$ID)
	one.GR.hind <- hindmap.GR[queryHits(findOverlaps(hindmap.GR, one.GR)),]
	one.GR.hind$ID <- one.GR[subjectHits(findOverlaps(hindmap.GR, one.GR)),]$ID
	two.GR <- GRanges(paste0(data$seqnames2,':',data$start2,'-',data$start2),ID=data$ID)
	two.GR.hind <- hindmap.GR[queryHits(findOverlaps(hindmap.GR, two.GR)),]
	two.GR.hind$ID <- two.GR[subjectHits(findOverlaps(hindmap.GR, two.GR)),]$ID
	
	two.GR.hind <- two.GR.hind[order(two.GR.hind$ID)]
	one.GR.hind <- one.GR.hind[order(one.GR.hind$ID)]

	GI <- GenomicInteractions(
		one.GR.hind, two.GR.hind)	
	
		GI$counts = data$counts[match(GI$anchor1.ID,data$ID)]
		GI$bait1 = data$bait1[match(GI$anchor1.ID,data$ID)]
		GI$bait2 = data$bait2[match(GI$anchor1.ID,data$ID)]
		GI$dir = data$dir[match(GI$anchor1.ID,data$ID)]
		GI$ihw = data$qvalue[match(GI$anchor1.ID,data$ID)]
		

	cat('taken and GOTHiC and converted to GI', "\n")
	
	GI <- unique(GI)
	
	return(GI)
}

CONSENSUS_GI <- makeGIfromGOTHiC_Differential(CONSENSUS)
CONSENSUSDEC_GI <- makeGIfromGOTHiC_Differential(CONSENSUSDEC)
ALL.GI <- c(CONSENSUS_GI,CONSENSUSDEC_GI)
Eq <- GetENTREZ(CONSENSUS_GI$bait1)
CONSENSUS_GI$ENTREZ_1 <- Eq$ENTREZ[match(CONSENSUS_GI$bait1,Eq$SYMBOL)]
Eq <- GetENTREZ(CONSENSUS_GI$bait2)
CONSENSUS_GI$ENTREZ_2 <- Eq$ENTREZ[match(CONSENSUS_GI$bait2,Eq$SYMBOL)]
Eq <- GetENTREZ(CONSENSUSDEC_GI$bait1)
CONSENSUSDEC_GI$ENTREZ_1 <- Eq$ENTREZ[match(CONSENSUSDEC_GI$bait1,Eq$SYMBOL)]
Eq <- GetENTREZ(CONSENSUSDEC_GI$bait2)
CONSENSUSDEC_GI$ENTREZ_2 <- Eq$ENTREZ[match(CONSENSUSDEC_GI$bait2,Eq$SYMBOL)]

# CONSENSUSDEC_GI is the lost interaction GI object used for downstream analysis
# CONSENSUS_GI is the gained interaction GI object used for downstream analysis


cat("#########################","\n")
cat("loaded GOTHiC Differential","\n")
cat("#########################","\n")

cat("#########################","\n")
cat("#########################","\n")
cat("loaded all U937-PR9 data required for analysis","\n")
cat("#########################","\n")
cat("#########################","\n")




