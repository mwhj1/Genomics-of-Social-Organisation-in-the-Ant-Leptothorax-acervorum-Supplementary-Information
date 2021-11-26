### BRAKER2 Predicted Genes - BLAST searches for ID ###

# The BRAKER2 gene prediction pipeline was used to predict genes
# in P and FM Supernova assemblies of L. acervorum.  Amino acid 
# and nucleotide sequences of predicted genes were quereied against
# the SwissProt and nt databases (max E-value 1e-10) to aid in 
# identification.  This script analyses those BLAST results.


#1. Setwd(), read in data and format

#Setwd()
setwd("/Users/u1693640/Documents/PhD/2019/RNAseq_for_Annotation/Results/BLAST_searches_of_predicted_AA_sequences")
#Read in data
raw_FM_aa <- read.table(file="BLASTp.augustus.hints.OT3b.aa.txt", quote="", header = FALSE, 
                        stringsAsFactors = FALSE, sep = '\t')
raw_FM_nt <- read.table(file="BLASTn.augustus.hints.OT3b.codingseq.txt", quote="", header = FALSE, 
                        stringsAsFactors = FALSE, sep = '\t')
raw_P_aa <- read.table(file="BLASTp.augustus.hints.PF.aa.txt", quote="", header = FALSE, 
                       stringsAsFactors = FALSE, sep = '\t')
raw_P_nt <- read.table(file="BLASTn.augustus.hints.PF.codingseq.txt", quote="", header = FALSE, 
                       stringsAsFactors = FALSE, sep = '\t')
#Set column names
column_names <- c("qseqid", "sseqid", "evalue", "length", "qlen", "slen", 
                  "qstart", "qend", "sstart", "send", "pident", "salltitles")
colnames(raw_FM_aa) <- column_names
colnames(raw_FM_nt) <- column_names
colnames(raw_P_aa) <- column_names
colnames(raw_P_nt) <- column_names
#Sort by query sequence ID
raw_FM_aa <- as.data.frame(raw_FM_aa[order(raw_FM_aa$qseqid), ])
raw_FM_nt <- as.data.frame(raw_FM_nt[order(raw_FM_nt$qseqid), ])
raw_P_aa <- as.data.frame(raw_P_aa[order(raw_P_aa$qseqid), ])
raw_P_nt <- as.data.frame(raw_P_nt[order(raw_P_nt$qseqid), ])


#2. Filter to high confidence hits

#First filter to E-values of 1e-40
e40_FM_aa <- raw_FM_aa[raw_FM_aa$evalue <= 1e-40, ]
e40_FM_nt <- raw_FM_nt[raw_FM_nt$evalue <= 1e-40, ]
e40_P_aa <- raw_P_aa[raw_P_aa$evalue <= 1e-40, ]
e40_P_nt <- raw_P_nt[raw_P_nt$evalue <= 1e-40, ]
#Then filter to records with a hit to 75% or more of the query
#First create a new column containing the % of the query
#covered by the alignment
percent_qry_aln <- function(BLAST_results){
  BLAST_results$perc_qry_aln <- BLAST_results$length / (BLAST_results$qlen / 100)
  return(BLAST_results)
}
filter_75pc_qry <- function(BLAST_results2){
  result <- percent_qry_aln(BLAST_results2)
  result <- result[result$perc_qry_aln >= 75, ]
  return(result)
}
qry75_e40_FM_aa <- filter_75pc_qry(e40_FM_aa)
qry75_e40_FM_nt <- filter_75pc_qry(e40_FM_nt)
qry75_e40_P_aa <- filter_75pc_qry(e40_P_aa)
qry75_e40_P_nt <- filter_75pc_qry(e40_P_nt)


#Filter on qseqid and sseqid
qsequniq <- function(BLAST_results){
  result <- BLAST_results[!duplicated(BLAST_results[,c('qseqid')]),]
  return(result)
}
ssequniq <- function(BLAST_results){
  result <- BLAST_results[!duplicated(BLAST_results[,c('sseqid')]),]
  return(result)
}
S.Q.sequniq <- function(BLAST_results){
  result <- ssequniq(qsequniq(BLAST_results))
  return(result)
}
qsequniq_FM_aa <- qsequniq(qry75_e40_FM_aa)
qsequniq_FM_nt <- qsequniq(qry75_e40_FM_nt)
qsequniq_P_aa <- qsequniq(qry75_e40_P_aa)
qsequniq_P_nt <- qsequniq(qry75_e40_P_nt)
unique_FM_aa <- S.Q.sequniq(qry75_e40_FM_aa)
unique_FM_nt <- S.Q.sequniq(qry75_e40_FM_nt)
unique_P_aa <- S.Q.sequniq(qry75_e40_P_aa)
unique_P_nt <- S.Q.sequniq(qry75_e40_P_nt)


#3. Get information on genes 
# (BLAST results are for all transcripts)

#Add a column of gene IDs 
#This just entails stripping off the ".*" characters from
#the end of qseqid
qgene <- function(BLAST_results){
  qgenevector <- c()
  for(i in 1:nrow(BLAST_results)){
    qry <- BLAST_results$qseqid[i]
    len <- nchar(qry)
    stoppos <- len-3
    gene <- substr(qry, 1, stoppos)
    qgenevector <- c(qgenevector, gene)
  }
  BLAST_results$qgeneid <- qgenevector
  return(BLAST_results)
}
qsequniq_genes_FM_aa <- qgene(qsequniq_FM_aa)
qsequniq_genes_FM_nt <- qgene(qsequniq_FM_nt)
qsequniq_genes_P_aa <- qgene(qsequniq_P_aa)
qsequniq_genes_P_nt <- qgene(qsequniq_P_nt)
genes_FM_aa <- qgene(unique_FM_aa)
genes_FM_nt <- qgene(unique_FM_nt)
genes_P_aa <- qgene(unique_P_aa)
genes_P_nt <- qgene(unique_P_nt)


#4. Report some results

#Get raw counts for all transcripts
nrow(raw_P_aa)
nrow(raw_FM_aa)
nrow(raw_P_nt)
nrow(raw_FM_nt)

#Get filterd counts for all transcripts
#1st round of filtering (qsequniq)
nrow(qsequniq_genes_FM_aa)
nrow(qsequniq_genes_FM_nt)
nrow(qsequniq_genes_P_aa)
nrow(qsequniq_genes_P_nt)
#2nd round of filtering (ssequniq)
nrow(unique_P_aa)
nrow(unique_FM_aa)
nrow(unique_P_nt)
nrow(unique_FM_nt)

#Get unique query gene (as opposed to transcript) counts
#From 1st round of filtering
length(unique(qsequniq_genes_FM_aa$qgeneid))
length(unique(qsequniq_genes_FM_nt$qgeneid))
length(unique(qsequniq_genes_P_aa$qgeneid))
length(unique(qsequniq_genes_P_nt$qgeneid))
#From 2nd round of filtering
length(unique(genes_FM_aa$qgeneid))
length(unique(genes_FM_nt$qgeneid))
length(unique(genes_P_aa$qgeneid))
length(unique(genes_P_nt$qgeneid))

#Get number of transcript hits (sseqids) not seen in cognate hits
compare_sseq <- function(df1, df2){
  sseq1 <- unique(df1$sseqid)
  sseq2 <- unique(df2$sseqid)
  oneINtwo <- length(which(sseq1 %in% sseq2))
  twoINone <- length(which(sseq2 %in% sseq1))
  perc_oneIntwo <- oneINtwo / (length(sseq1)/100)
  perc_twoInone <- twoINone / (length(sseq2)/100)
  perc_missing_1in2 <- 100 - perc_oneIntwo
  perc_missing_2in1 <- 100 - perc_twoInone
  result <- paste("N from 1 in 2 = ", oneINtwo, " (", perc_oneIntwo, "%) (", perc_missing_2in1, "% missing). ", 
                  "N from 2 in 1 = ", twoINone, " (", perc_twoInone, "%) (", perc_missing_2in1, "% missing).", 
                  sep = "")
  result
}
#From 1st round of filtering
compare_sseq(qsequniq_genes_FM_aa, qsequniq_genes_P_aa)
compare_sseq(qsequniq_genes_P_aa, qsequniq_genes_FM_aa)
compare_sseq(qsequniq_genes_FM_nt, qsequniq_genes_P_nt)
compare_sseq(qsequniq_genes_P_nt, qsequniq_genes_FM_nt)
#From 2nd round of filtering
compare_sseq(genes_FM_aa, genes_P_aa)
compare_sseq(genes_P_aa, genes_FM_aa)
compare_sseq(genes_FM_nt, genes_P_nt)
compare_sseq(genes_P_nt, genes_FM_nt)

#Get number of gene hits (sseqids) not seen in cognate hits
qgeneuniq <- function(BLAST_results){
  result <- BLAST_results[!duplicated(BLAST_results[,c('qgeneid')]),]
  return(result)
}
#From 1st round of filtering
compare_sseq(qgeneuniq(qsequniq_genes_FM_aa), qgeneuniq(qsequniq_genes_P_aa))
compare_sseq(qgeneuniq(qsequniq_genes_P_aa), qgeneuniq(qsequniq_genes_FM_aa))
compare_sseq(qgeneuniq(qsequniq_genes_FM_nt), qgeneuniq(qsequniq_genes_P_nt))
compare_sseq(qgeneuniq(qsequniq_genes_P_nt), qgeneuniq(qsequniq_genes_FM_nt))
#From 2nd round of filtering
compare_sseq(qgeneuniq(genes_FM_aa), qgeneuniq(genes_P_aa))
compare_sseq(qgeneuniq(genes_P_aa), qgeneuniq(genes_FM_aa))
compare_sseq(qgeneuniq(genes_FM_nt), qgeneuniq(genes_P_nt))
compare_sseq(qgeneuniq(genes_P_nt), qgeneuniq(genes_FM_nt))


#Compare the numbers of unique gene IDs from the first
#and second rounds of filtering
dupFMaa <- length(unique(qsequniq_genes_FM_aa$qgeneid)) - length(unique(genes_FM_aa$qgeneid))
dupPaa <- length(unique(qsequniq_genes_P_aa$qgeneid)) - length(unique(genes_P_aa$qgeneid))
dupFMnt <- length(unique(qsequniq_genes_FM_nt$qgeneid)) - length(unique(genes_FM_nt$qgeneid))
dupPnt <- length(unique(qsequniq_genes_P_nt$qgeneid)) - length(unique(genes_P_nt$qgeneid))
se <- function(x) sd(x)/sqrt(length(x))
mean(c(dupFMaa, dupPaa, dupFMnt, dupPnt))
se(c(dupFMaa, dupPaa, dupFMnt, dupPnt))
