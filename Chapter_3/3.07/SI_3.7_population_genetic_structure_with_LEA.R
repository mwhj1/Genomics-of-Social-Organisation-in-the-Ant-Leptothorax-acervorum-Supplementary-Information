#####################################################################################
########################### Population Structure with LEA ###########################
#####################################################################################
### Adapted from http://membres-timc.imag.fr/Olivier.Francois/tutoRstructure.pdf ####
#####################################################################################

#Packages
#
#install.packages(c("fields","RColorBrewer","mapplots"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")
#install.packages("mapplots")
#install.packages("rworldmap")
#install.packages("gridExtra")
library(gridExtra)
library(rworldmap)
library(LEA)
library(parallel)
library(dplyr)
library(mapplots)
library(ggplot2)
library(tidyverse)
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
#
#####################################################################################
#
#1. Set working directory and prepare data
#
top_level_dir = "/Users/u1693640/Documents/PhD/2020/bcftools_Fst"
#
master_VCF_dir = "/Users/u1693640/Documents/PhD/2020/bcftools_Fst/VCFs"
setwd(master_VCF_dir)
#
#Read sample-population-coordinates file
poptable = read.table("./all_sample_names_coordinates.txt", 
                      stringsAsFactors = FALSE, header = TRUE)
#
#Main VCFs - these contain all SNPs - one file for each reference
#(These VCFs were used to produce per-SNP Fst estimates, in turn used 
#for Fst outlier analysis and social marker SNP identification).
#Convert VCF to geno format:
OT_VCF = "OTref.biallelic.DP3g1maf05.recode.vcf"
PF_VCF = "PFref.biallelic.DP3g1maf05.recode.vcf"
output = vcf2geno(OT_VCF)
output = vcf2geno(PF_VCF)
OT_geno = "./OTref.biallelic.DP3g1maf05.recode.geno"
PF_geno = "./PFref.biallelic.DP3g1maf05.recode.geno"
#
#In addition to analysing all SNPs, separate analyses will
#be conducted on SNPs in two subsets of scaffolds from each reference:
#Subset 1: all SNPs on scaffolds containing marker SNPs removed.
#Subset 2: all SNPs on scaffolds not containing marker SNPs removed.
#The rationale of subset 1 is to examine structure using only loci which 
#are completely un-linked to any SNPs identified as social markers, 
#which produces a conservative estimate of structure unrelated to
#selection on social phenotype (assuming no false negatives in the marker SNPs)
#The rationale of subset 2 is to complement subset 1 as well as the full SNP set.
#
#To prepare the subsets, the Fst outlier results are read
setwd("/Users/u1693640/Documents/PhD/2020/bcftools_Fst/Results")
OT_master <- read.table(file = "OTref_MASTER_Fst_results.tsv", 
                        fill = TRUE, stringsAsFactors = FALSE, 
                        header = TRUE, sep = "\t", quote = "")
PF_master <- read.table(file = "PFref_MASTER_Fst_results.tsv", 
                        fill = TRUE, stringsAsFactors = FALSE, 
                        header = TRUE, sep = "\t", quote = "")
#Select only outlier SNPs
OT_markers <- OT_master %>% 
  select(c(CHROM, POS, outlier)) %>% 
  filter(outlier == "outlier") %>% 
  select(-c(outlier))
PF_markers <- PF_master %>% 
  select(c(CHROM, POS, outlier)) %>% 
  filter(outlier == "outlier") %>% 
  select(-c(outlier))
#Filter all rows from *_master dataframes
#which relate to scaffolds in the *_markers subsets
OT_marker_scaffolds <- OT_master[which(OT_master$CHROM %in% OT_markers$CHROM), ] %>% select(c(CHROM, POS))
PF_marker_scaffolds <- PF_master[which(PF_master$CHROM %in% PF_markers$CHROM), ] %>% select(c(CHROM, POS))
#Write to file
write.table(OT_marker_scaffolds, file = "OT_marker_scaffold_SNPs.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(PF_marker_scaffolds, file = "PF_marker_scaffold_SNPs.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
############################
#Externally pass to VCFtools to filter the main VCFs

#Subset 1 (e.g.):
#vcftools --vcf in.vcf \
#--out out_prefix \
#--exclude-positions SNP_positions.txt \
#--recode
#...
#OT: After filtering, kept 1160192 out of a possible 2062621 Sites
#PF: After filtering, kept 1294948 out of a possible 2054974 Sites

#Subset 2 (e.g.):
#vcftools --vcf in.vcf \
#--out out_prefix \
#--positions SNP_positions.txt \
#--recode
#...
#OT: After filtering, kept 902429 out of a possible 2062621 Sites
#PF: After filtering, kept 760026 out of a possible 2054974 Sites
############################
#
#Convert filtered VCFs to geno
#Subset 1
subset1_dir = "/Users/u1693640/Documents/PhD/2020/bcftools_Fst/VCFs/No_marker_scaffolds"
setwd(subset1_dir)
OT_VCF_subset1 = "OTref.biallelic.DP3g1maf05.recode.no_marker_scaffolds.recode.vcf"
PF_VCF_subset1 = "PFref.biallelic.DP3g1maf05.recode.no_marker_scaffolds.recode.vcf"
output = vcf2geno(OT_VCF_subset1)
output = vcf2geno(PF_VCF_subset1)
OT_geno_subset1 = "./OTref.biallelic.DP3g1maf05.recode.no_marker_scaffolds.recode.geno"
PF_geno_subset1 = "./PFref.biallelic.DP3g1maf05.recode.no_marker_scaffolds.recode.geno"
#Subset 2
subset2_dir = "/Users/u1693640/Documents/PhD/2020/bcftools_Fst/VCFs/Only_marker_scaffolds"
setwd(subset2_dir)
OT_VCF_subset2 = "OTref.biallelic.DP3g1maf05.recode.only_marker_scaffolds.recode.vcf"
PF_VCF_subset2 = "PFref.biallelic.DP3g1maf05.recode.only_marker_scaffolds.recode.vcf"
output = vcf2geno(OT_VCF_subset2)
output = vcf2geno(PF_VCF_subset2)
OT_geno_subset2 = "./OTref.biallelic.DP3g1maf05.recode.only_marker_scaffolds.recode.geno"
PF_geno_subset2 = "./PFref.biallelic.DP3g1maf05.recode.only_marker_scaffolds.recode.geno"
#
#####################################################################################
#
#2. Determine correct K value (number of ancestral populations)
# for all three data sets (all SNPs, and the two subsets) using the 
# snmf() function on a range of K from 1:10.  
#
#All SNPs
setwd(master_VCF_dir)
OT_snmf_K_1_10 = snmf(OT_geno, K = 1:10, ploidy = 2, entropy = T, alpha = 0, 
                      CPU = 3, project = "new")
PF_snmf_K_1_10 = snmf(PF_geno, K = 1:10, ploidy = 2, entropy = T, alpha = 0, 
                      CPU = 3, project = "new")
#Subset 1
setwd(subset1_dir)
OT_subset1_snmf_K_1_10 = snmf(OT_geno_subset1, K = 1:10, ploidy = 2, entropy = T, alpha = 0, 
                              CPU = 3, project = "new")
PF_subset1_snmf_K_1_10 = snmf(PF_geno_subset1, K = 1:10, ploidy = 2, entropy = T, alpha = 0, 
                              CPU = 3, project = "new")
#Subset 2
setwd(subset2_dir)
OT_subset2_snmf_K_1_10 = snmf(OT_geno_subset2, K = 1:10, ploidy = 2, entropy = T, alpha = 0, 
                              CPU = 3, project = "new")
PF_subset2_snmf_K_1_10 = snmf(PF_geno_subset2, K = 1:10, ploidy = 2, entropy = T, alpha = 0, 
                              CPU = 3, project = "new")
#
#Plot - do the points platteau (or "elbow")?
#The K value where they do should be the K value to use
setwd(master_VCF_dir)
plot(OT_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "OTref all SNPs")
plot(PF_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "PFref all SNPs")
setwd(subset1_dir)
plot(OT_subset1_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "OTref_subset1 SNPs on scaffolds with no markers")
plot(PF_subset1_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "PFref_subset1 SNPs on scaffolds with no markers")
setwd(subset2_dir)
plot(OT_subset2_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "OTref_subset2 only SNPs on scaffolds with markers")
plot(PF_subset2_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "PFref_subset2 only SNPs on scaffolds with markers")
#
#Arrange plots on one page (these need to be manually exported in the right-hand pane)
par(mfrow = c(3,2))
setwd(master_VCF_dir)
plot(OT_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "OTref all SNPs")
plot(PF_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "PFref all SNPs")
setwd(subset1_dir)
plot(OT_subset1_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "OTref_subset1 SNPs on scaffolds with no markers")
plot(PF_subset1_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "PFref_subset1 SNPs on scaffolds with no markers")
setwd(subset2_dir)
plot(OT_subset2_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "OTref_subset2 only SNPs on scaffolds with markers")
plot(PF_subset2_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "PFref_subset2 only SNPs on scaffolds with markers")
par(mfrow=c(1,1))
#
#####################################################################################
#
#3. Determine correct alpha value
#Take the lowest value of K for each dataset, and re-run snfm() while
#varying alpha (the "regularization parameter")
#
#All SNPs
setwd(master_VCF_dir)
OT_snmf_K_2_a_100 = snmf(OT_geno, K = 2, ploidy = 2, entropy = T, alpha = 100, 
                         CPU = 3, project = "new")
PF_snmf_K_2_a_100 = snmf(PF_geno, K = 2, ploidy = 2, entropy = T, alpha = 100, 
                         CPU = 3, project = "new")
#Subset 1
setwd(subset1_dir)
OT_subset1_snmf_K_2_a_100 = snmf(OT_geno_subset1, K = 2, ploidy = 2, entropy = T, alpha = 100, 
                                 CPU = 3, project = "new")
PF_subset1_snmf_K_2_a_100 = snmf(PF_geno_subset1, K = 2, ploidy = 2, entropy = T, alpha = 100, 
                                 CPU = 3, project = "new")
#Subset 2
setwd(subset2_dir)
OT_subset2_snmf_K_2_a_100 = snmf(OT_geno_subset2, K = 2, ploidy = 2, entropy = T, alpha = 100, 
                                 CPU = 3, project = "new")
PF_subset2_snmf_K_2_a_100 = snmf(PF_geno_subset2, K = 2, ploidy = 2, entropy = T, alpha = 100, 
                                 CPU = 3, project = "new")
#
#All SNPs
setwd(master_VCF_dir)
OT_snmf_K_2_a_500 = snmf(OT_geno, K = 2, ploidy = 2, entropy = T, alpha = 500, 
                         CPU = 3, project = "new")
PF_snmf_K_2_a_500 = snmf(PF_geno, K = 2, ploidy = 2, entropy = T, alpha = 500, 
                         CPU = 3, project = "new")
#Subset 1
setwd(subset1_dir)
OT_subset1_snmf_K_2_a_500 = snmf(OT_geno_subset1, K = 2, ploidy = 2, entropy = T, alpha = 500, 
                                 CPU = 3, project = "new")
PF_subset1_snmf_K_2_a_500 = snmf(PF_geno_subset1, K = 2, ploidy = 2, entropy = T, alpha = 500, 
                                 CPU = 3, project = "new")
#Subset 2
setwd(subset2_dir)
OT_subset2_snmf_K_2_a_500 = snmf(OT_geno_subset2, K = 2, ploidy = 2, entropy = T, alpha = 500, 
                                 CPU = 3, project = "new")
PF_subset2_snmf_K_2_a_500 = snmf(PF_geno_subset2, K = 2, ploidy = 2, entropy = T, alpha = 500, 
                                 CPU = 3, project = "new")
#
#All SNPs
setwd(master_VCF_dir)
OT_snmf_K_2_a_1000 = snmf(OT_geno, K = 2, ploidy = 2, entropy = T, alpha = 1000, 
                          CPU = 3, project = "new")
PF_snmf_K_2_a_1000 = snmf(PF_geno, K = 2, ploidy = 2, entropy = T, alpha = 1000, 
                          CPU = 3, project = "new")
#Subset 1
setwd(subset1_dir)
OT_subset1_snmf_K_2_a_1000 = snmf(OT_geno_subset1, K = 2, ploidy = 2, entropy = T, alpha = 1000, 
                                  CPU = 3, project = "new")
PF_subset1_snmf_K_2_a_1000 = snmf(PF_geno_subset1, K = 2, ploidy = 2, entropy = T, alpha = 1000, 
                                  CPU = 3, project = "new")
#Subset 2
setwd(subset2_dir)
OT_subset2_snmf_K_2_a_1000 = snmf(OT_geno_subset2, K = 2, ploidy = 2, entropy = T, alpha = 1000, 
                                  CPU = 3, project = "new")
PF_subset2_snmf_K_2_a_1000 = snmf(PF_geno_subset2, K = 2, ploidy = 2, entropy = T, alpha = 1000, 
                                  CPU = 3, project = "new")
#
#Place all cross-entropy results from one set of SNPs / reference pair into numeric vectors
OT_all_alpha_vect = c(as.data.frame(summary(OT_snmf_K_1_10)[2][[1]])[2,2], 
                      as.data.frame(summary(OT_snmf_K_2_a_100)[2][[1]])[2,1], 
                      as.data.frame(summary(OT_snmf_K_2_a_500)[2][[1]])[2,1], 
                      as.data.frame(summary(OT_snmf_K_2_a_1000)[2][[1]])[2,1])
PF_all_alpha_vect = c(as.data.frame(summary(PF_snmf_K_1_10)[2][[1]])[2,2], 
                      as.data.frame(summary(PF_snmf_K_2_a_100)[2][[1]])[2,1], 
                      as.data.frame(summary(PF_snmf_K_2_a_500)[2][[1]])[2,1], 
                      as.data.frame(summary(PF_snmf_K_2_a_1000)[2][[1]])[2,1])
OT_subset1_alpha_vect = c(as.data.frame(summary(OT_subset1_snmf_K_1_10)[2][[1]])[2,2], 
                          as.data.frame(summary(OT_subset1_snmf_K_2_a_100)[2][[1]])[2,1], 
                          as.data.frame(summary(OT_subset1_snmf_K_2_a_500)[2][[1]])[2,1], 
                          as.data.frame(summary(OT_subset1_snmf_K_2_a_1000)[2][[1]])[2,1])
PF_subset1_alpha_vect = c(as.data.frame(summary(PF_subset1_snmf_K_1_10)[2][[1]])[2,2], 
                          as.data.frame(summary(PF_subset1_snmf_K_2_a_100)[2][[1]])[2,1], 
                          as.data.frame(summary(PF_subset1_snmf_K_2_a_500)[2][[1]])[2,1], 
                          as.data.frame(summary(PF_subset1_snmf_K_2_a_1000)[2][[1]])[2,1])
OT_subset2_alpha_vect = c(as.data.frame(summary(OT_subset2_snmf_K_1_10)[2][[1]])[2,2], 
                          as.data.frame(summary(OT_subset2_snmf_K_2_a_100)[2][[1]])[2,1], 
                          as.data.frame(summary(OT_subset2_snmf_K_2_a_500)[2][[1]])[2,1], 
                          as.data.frame(summary(OT_subset2_snmf_K_2_a_1000)[2][[1]])[2,1])
PF_subset2_alpha_vect = c(as.data.frame(summary(PF_subset2_snmf_K_1_10)[2][[1]])[2,2], 
                          as.data.frame(summary(PF_subset2_snmf_K_2_a_100)[2][[1]])[2,1], 
                          as.data.frame(summary(PF_subset2_snmf_K_2_a_500)[2][[1]])[2,1], 
                          as.data.frame(summary(PF_subset2_snmf_K_2_a_1000)[2][[1]])[2,1])
#Place in one dataframe
alpha_df = as.data.frame(c(0, 100, 500, 1000), stringsAsFactors = FALSE)
alpha_df$OT_all = OT_all_alpha_vect
alpha_df$PF_all = PF_all_alpha_vect
alpha_df$OT_subset1 = OT_subset1_alpha_vect
alpha_df$PF_subset1 = PF_subset1_alpha_vect
alpha_df$OT_subset2 = OT_subset2_alpha_vect
alpha_df$PF_subset2 = PF_subset2_alpha_vect
colnames(alpha_df)[1] <- "alpha"
alpha_df = alpha_df %>% pivot_longer(!alpha, names_to = "dataset", values_to = "CE")
#Visualise with a plot
ggplot(alpha_df) + 
  geom_line(aes(x = alpha, y = CE, color = dataset))
#Clearly this parameter does not make a lot of difference, I will leave it 
#on the default value of 10
#
#####################################################################################
#
#4. Using K=2 and default alpha (10) each dataset will be run through
#snmf() with 10 repititions to get reliable results
#
#All SNPs
setwd(master_VCF_dir)
OT_snmf_K_2_a_10_10reps = snmf(OT_geno, K = 2, ploidy = 2, entropy = T,  
                               CPU = 3, project = "new", repetitions = 10)
PF_snmf_K_2_a_10_10reps = snmf(PF_geno, K = 2, ploidy = 2, entropy = T, 
                               CPU = 3, project = "new", repetitions = 10)
#Subset 1
setwd(subset1_dir)
OT_snmf_subset1_K_2_a_10_10reps = snmf(OT_geno_subset1, K = 2, ploidy = 2, entropy = T,  
                                       CPU = 3, project = "new", repetitions = 10)
PF_snmf_subset1_K_2_a_10_10reps = snmf(PF_geno_subset1, K = 2, ploidy = 2, entropy = T, 
                                       CPU = 3, project = "new", repetitions = 10)
#Subset 2
setwd(subset2_dir)
OT_snmf_subset2_K_2_a_10_10reps = snmf(OT_geno_subset2, K = 2, ploidy = 2, entropy = T,  
                                       CPU = 3, project = "new", repetitions = 10)
PF_snmf_subset2_K_2_a_10_10reps = snmf(PF_geno_subset2, K = 2, ploidy = 2, entropy = T, 
                                       CPU = 3, project = "new", repetitions = 10)
#
#Select the repetition from each which has the lowest cross-entropy
OT_snmf_best = which.min(cross.entropy(OT_snmf_K_2_a_10_10reps))
PF_snmf_best = which.min(cross.entropy(PF_snmf_K_2_a_10_10reps))
OT_snmf_subset1_best = which.min(cross.entropy(OT_snmf_subset1_K_2_a_10_10reps))
PF_snmf_subset1_best = which.min(cross.entropy(PF_snmf_subset1_K_2_a_10_10reps))
OT_snmf_subset2_best = which.min(cross.entropy(OT_snmf_subset2_K_2_a_10_10reps))
PF_snmf_subset2_best = which.min(cross.entropy(PF_snmf_subset2_K_2_a_10_10reps))
#
#Create Q matrices
OT_qmatrix = Q(OT_snmf_K_2_a_10_10reps, K = 2, run = OT_snmf_best)
PF_qmatrix = Q(PF_snmf_K_2_a_10_10reps, K = 2, run = PF_snmf_best)
OT_subset1_qmatrix = Q(OT_snmf_subset1_K_2_a_10_10reps, K = 2, run = OT_snmf_subset1_best)
PF_subset1_qmatrix = Q(PF_snmf_subset1_K_2_a_10_10reps, K = 2, run = PF_snmf_subset1_best)
OT_subset2_qmatrix = Q(OT_snmf_subset2_K_2_a_10_10reps, K = 2, run = OT_snmf_subset2_best)
PF_subset2_qmatrix = Q(PF_snmf_subset2_K_2_a_10_10reps, K = 2, run = PF_snmf_subset2_best)
#
#####################################################################################
#
#5. Visualise results
#
#Get sample names
setwd(top_level_dir)
names_table = read.table(file = "sample_populations.txt", 
                         header = FALSE, 
                         stringsAsFactors = FALSE)
colnames(names_table) = c("sample", "population")
sample_names = names_table$sample
#Strip out all the unnecessary characters
for(i in 1:length(sample_names)){
  sample_names[i] = substr(sample_names[i], 1, 11)
}
#There are some extra characters which should be removed, this will be done by hand
write.table(sample_names, file = "sample_names.txt", quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)
#[manual editing]
sample_names = read.table("sample_names.txt", 
                          header = FALSE, 
                          stringsAsFactors = FALSE)
sample_names = sample_names$V1
#
#Set plotting directory
plotdir = "/Users/u1693640/Documents/PhD/2020/bcftools_Fst/Structure_results/Plots"
setwd(plotdir)
#
#Simple bar-plots
pdf(file = "barplots.pdf", paper = "a4r", title = "Individual ancestry coefficients", pointsize = 6, width = 14)
par(mfrow = c(6,2))
barplot(t(OT_qmatrix), col = c("orange", "lightgreen"), border = NA, space = 0, 
        xlab = "", ylab = "Admixture coefficients", main = "OTref", 
        names.arg = sample_names, las = 2, cex.names = 0.7)
#Note that for this plot the order of colours is switched - it is clear from the figure
#that the same broad patterns are seen with both references, but the ordering of ancestral
#populations appears to be different.  As the sequences aligned to each reference 
#(before alignment cleaning) were the same, this is likely an artefact of the references
#(although the fact that the other two PFref barplots don't behave require this colour
#switching suggests it might be a result of the analysis).
#Also note when viewing the figures that while switching the colours has made the area
#of each individual corresponding to each of the 2 ancestral populations (i.e. the
#area of each colour in each bar of the plot) comparable between references, the ordering 
#is still inverted - the patterns seen are the same but for the full SNP set plots corresponding
#to alternate references appear inverted with respect to each other even after swapping
#the colour order.
barplot(t(PF_qmatrix), col = c("lightgreen", "orange"), border = NA, space = 0, 
        xlab = "", ylab = "Admixture coefficients", main = "PFref", 
        names.arg = sample_names, las = 2, cex.names = 0.7)
barplot(t(OT_subset1_qmatrix), col = c("orange", "lightgreen"), border = NA, space = 0, 
        xlab = "", ylab = "Admixture coefficients", main = "OTref scaffolds without markers", 
        names.arg = sample_names, las = 2, cex.names = 0.7)
barplot(t(PF_subset1_qmatrix), col = c("orange", "lightgreen"), border = NA, space = 0, 
        xlab = "", ylab = "Admixture coefficients", main = "PFref scaffolds without markers", 
        names.arg = sample_names, las = 2, cex.names = 0.7)
barplot(t(OT_subset2_qmatrix), col = c("orange", "lightgreen"), border = NA, space = 0, 
        xlab = "", ylab = "Admixture coefficients", main = "OTref only marker scaffolds", 
        names.arg = sample_names, las = 2, cex.names = 0.7)
barplot(t(PF_subset2_qmatrix), col = c("orange", "lightgreen"), border = NA, space = 0, 
        xlab = "", ylab = "Admixture coefficients", main = "PFref only marker scaffolds", 
        names.arg = sample_names, las = 2, cex.names = 0.7)
dev.off()
#
#A function to take the mean of each population's admixture coefficients
get_mean_coefficients <- function(qmatrix, population_table){
  #Coerce to dataframes
  qmatrix = as.data.frame(qmatrix)
  population_table = as.data.frame(population_table)
  #Create result as an empty dataframe
  result = as.data.frame(matrix(NA, ncol = ncol(qmatrix), 
                                nrow = length(unique(population_table$Pop))))
  #Get population names in result column
  result$Pop = unique(population_table$Pop)
  #And in qmatrix
  qmatrix$Pop = population_table$Pop
  #Loop over populations
  for(i in 1:length(result$Pop)){
    pop = result$Pop[i]
    #Subset qmatrix (dropping Pop column leaving ncol = n clusters)
    qmatrix_subset = qmatrix[which(qmatrix$Pop == pop), ] %>% select(-Pop)
    #Loop over clusters
    for(j in 1:ncol(qmatrix_subset)){
      result[i,j] = mean(qmatrix_subset[,j])
    }
  }
  #Rename result columns
  for(i in 1:(ncol(result)-1)){
    colnames(result)[i] = paste("cluster", i, sep = "_")
  }
  return(result)
}
#
#Call the function
OT_mean_coefficients = get_mean_coefficients(OT_qmatrix, poptable)
PF_mean_coefficients = get_mean_coefficients(PF_qmatrix, poptable)
OT_subset1_mean_coefficients = get_mean_coefficients(OT_subset1_qmatrix, poptable)
PF_subset1_mean_coefficients = get_mean_coefficients(PF_subset1_qmatrix, poptable)
OT_subset2_mean_coefficients = get_mean_coefficients(OT_subset2_qmatrix, poptable)
PF_subset2_mean_coefficients = get_mean_coefficients(PF_subset2_qmatrix, poptable)
#Stick in coordinates
OT_mean_coefficients$Lat = unique(poptable$Latt)
OT_mean_coefficients$Long = unique(poptable$Long)
PF_mean_coefficients$Lat = unique(poptable$Latt)
PF_mean_coefficients$Long = unique(poptable$Long)
OT_subset1_mean_coefficients$Lat = unique(poptable$Latt)
OT_subset1_mean_coefficients$Long = unique(poptable$Long)
PF_subset1_mean_coefficients$Lat = unique(poptable$Latt)
PF_subset1_mean_coefficients$Long = unique(poptable$Long)
OT_subset2_mean_coefficients$Lat = unique(poptable$Latt)
OT_subset2_mean_coefficients$Long = unique(poptable$Long)
PF_subset2_mean_coefficients$Lat = unique(poptable$Latt)
PF_subset2_mean_coefficients$Long = unique(poptable$Long)
#
#Map plots (these are saved individually because faceting with par(mfrow) breaks
#the mapRegion command)
pdf(file = "OT_all_map.pdf", paper = "a4r", title = "Population mean ancestry coefficients", pointsize = 6, width = 14)
par(mfrow = c(1,1))
#OT
mapPies(OT_mean_coefficients, 
        nameX = "Long", 
        nameY = "Lat", 
        nameZs = c(names(OT_mean_coefficients)[1], names(OT_mean_coefficients)[2]), 
        zColours = c("orange", "lightgreen"), 
        mapRegion = "Spain", 
        landCol = "gray92", 
        symbolSize = 1.8)
title(main=paste("OTref all SNPs"),
      cex=3)
dev.off()
#PF
#(Note that the colours are switched around for this plot - each run was calculated separately
#so the colouring is arbitrary - it is clear though that the ordering of ancestral populations
#has been reversed for this run, and switching the order of colours assigned to cluseters makes
#comparing figures easier.)
pdf(file = "PF_all_map.pdf", paper = "a4r", title = "Population mean ancestry coefficients", pointsize = 6, width = 14)
par(mfrow = c(1,1))
mapPies(PF_mean_coefficients, 
        nameX = "Long", 
        nameY = "Lat", 
        nameZs = c(names(PF_mean_coefficients)[1], names(PF_mean_coefficients)[2]), 
        zColours = c("lightgreen", "orange"), 
        mapRegion = "Spain", 
        landCol = "gray92", 
        symbolSize = 1.8)
title(main=paste("PFref all SNPs"),
      cex=3)
dev.off()
#OT subset1
pdf(file = "OT_no_marker_scaffolds_map.pdf", paper = "a4r", title = "Population mean ancestry coefficients", pointsize = 6, width = 14)
par(mfrow = c(1,1))
mapPies(OT_subset1_mean_coefficients, 
        nameX = "Long", 
        nameY = "Lat", 
        nameZs = c(names(OT_subset1_mean_coefficients)[1], names(OT_subset1_mean_coefficients)[2]), 
        zColours = c("orange", "lightgreen"), 
        mapRegion = "Spain", 
        landCol = "gray92", 
        symbolSize = 1.8)
title(main=paste("OTref no marker scaffolds"),
      cex=3)
dev.off()
#PF subset1
pdf(file = "PF_no_marker_scaffolds_map.pdf", paper = "a4r", title = "Population mean ancestry coefficients", pointsize = 6, width = 14)
par(mfrow = c(1,1))
mapPies(PF_subset1_mean_coefficients, 
        nameX = "Long", 
        nameY = "Lat", 
        nameZs = c(names(PF_subset1_mean_coefficients)[1], names(PF_subset1_mean_coefficients)[2]), 
        zColours = c("orange", "lightgreen"), 
        mapRegion = "Spain", 
        landCol = "gray92", 
        symbolSize = 1.8)
title(main=paste("PFref no markers scaffolds"),
      cex=3)
dev.off()
#OT subset2
pdf(file = "OT_only_marker_scaffolds_map.pdf", paper = "a4r", title = "Population mean ancestry coefficients", pointsize = 6, width = 14)
par(mfrow = c(1,1))
mapPies(OT_subset2_mean_coefficients, 
        nameX = "Long", 
        nameY = "Lat", 
        nameZs = c(names(OT_subset2_mean_coefficients)[1], names(OT_subset2_mean_coefficients)[2]), 
        zColours = c("orange", "lightgreen"), 
        mapRegion = "Spain", 
        landCol = "gray92", 
        symbolSize = 1.8)
title(main=paste("OTref only marker scaffolds"),
      cex=3)
dev.off()
#PF subset2
pdf(file = "PF_only_marker_scaffolds_map.pdf", paper = "a4r", title = "Population mean ancestry coefficients", pointsize = 6, width = 14)
par(mfrow = c(1,1))
mapPies(PF_subset2_mean_coefficients, 
        nameX = "Long", 
        nameY = "Lat", 
        nameZs = c(names(PF_subset2_mean_coefficients)[1], names(PF_subset2_mean_coefficients)[2]), 
        zColours = c("orange", "lightgreen"), 
        mapRegion = "Spain", 
        landCol = "gray92", 
        symbolSize = 1.8)
title(main=paste("PFref only markers scaffolds"),
      cex=3)
dev.off()
#
#Save copies of the first graphs (used to select K) too
pdf(file = "K_selection.pdf", paper = "a4", title = "Cross-entropy as a function of K", pointsize = 6, width = 14)
par(mfrow = c(6,2))
plot(OT_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "OTref all SNPs")
plot(PF_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "PFref all SNPs")
plot(OT_subset1_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "OTref_subset1 no SNPs on scaffolds with markers")
plot(PF_subset1_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "PFref_subset1 no SNPs on scaffolds with markers")
plot(OT_subset2_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "OTref_subset2 only SNPs on scaffolds with markers")
plot(PF_subset2_snmf_K_1_10, col = "blue4", cex = 1.4, pch = 19, main = "PFref_subset2 only SNPs on scaffolds with markers")
dev.off()
#
#####################################################################################
######################################## PCA ########################################
#####################################################################################
#
#PCA
#The LEA package can also do PCA analyses using geno objects as input
#
#Create PCA objects
setwd(master_VCF_dir)
OT_all_PCA = pca(OT_geno, center = FALSE)
PF_all_PCA = pca(PF_geno, center = FALSE)
setwd(subset1_dir)
OT_subset1_PCA = pca(OT_geno_subset1, center = FALSE)
PF_subset1_PCA = pca(PF_geno_subset1, center = FALSE)
setwd(subset2_dir)
OT_subset2_PCA = pca(OT_geno_subset2, center = FALSE)
PF_subset2_PCA = pca(PF_geno_subset2, center = FALSE)
#
#Perform Tracy-Widom tests on PCA objects to identify important components
OT_all_tw = tracy.widom(OT_all_PCA)
PF_all_tw = tracy.widom(PF_all_PCA)
OT_subset1_all_tw = tracy.widom(OT_subset1_PCA)
PF_subset1_all_tw = tracy.widom(PF_subset1_PCA)
OT_subset2_all_tw = tracy.widom(OT_subset2_PCA)
PF_subset2_all_tw = tracy.widom(PF_subset2_PCA)
#
#Visualise with plots
pcadir = "/Users/u1693640/Documents/PhD/2020/bcftools_Fst/PCA_results"
setwd(pcadir)
#Example plots
#PCA projections (PC1:PC2)
ggplot(as.data.frame(OT_all_PCA$projections)) + 
  geom_point(aes(x=V1, y=V2)) + 
  labs(x="PC1", y="PC2", title="OTref all SNPs")
#PCA standard deviations
ggplot(as.data.frame(OT_all_PCA$sdev)) + 
  geom_point(aes(x=1:47, y=V1)) + 
  labs(x="Principal components", y="Standard deviation", title="OTref all SNPs")
#% variance explained by each component
ggplot(as.data.frame(OT_all_tw$percentage)) + 
  geom_point(aes(x=1:47, y=OT_all_tw$percentage)) + 
  labs(x="Principal components", y="Proportion of variance explained", title="OTref all SNPs")
#Tracy-Widom p-values
ggplot(as.data.frame(OT_all_tw$pvalues)) + 
  geom_col(aes(x=1:47, y=OT_all_tw$pvalues)) + 
  labs(x="Principal components", y="p-values", title="OTref all SNPs")
#
#Save plots
#
#Main PCA plots - PCs 1:9 are significant according to the Tracy-Widom results,
#but the proportion of variance explained by each PC indicates it is not necessary
#to go beyond PC4.
#
#All SNPs
#
#OTref
pdf(file = "PCA_OTref.pdf", paper = "a4", title = "PCA_OTref", pointsize = 6, width = 14)
grid.arrange(
  ggplot(as.data.frame(OT_all_PCA$projections)) + 
    geom_point(aes(x=V1, y=V2)) + 
    labs(x="PC1", y="PC2", title="OTref all SNPs - PC1 PC2") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_all_PCA$projections)) + 
    geom_point(aes(x=V1, y=V3)) + 
    labs(x="PC1", y="PC3", title="OTref all SNPs - PC1 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_all_PCA$projections)) + 
    geom_point(aes(x=V1, y=V4)) + 
    labs(x="PC1", y="PC4", title="OTref all SNPs - PC1 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_all_PCA$projections)) + 
    geom_point(aes(x=V2, y=V3)) + 
    labs(x="PC2", y="PC3", title="OTref all SNPs - PC2 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_all_PCA$projections)) + 
    geom_point(aes(x=V2, y=V4)) + 
    labs(x="PC2", y="PC4", title="OTref all SNPs - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_all_PCA$projections)) + 
    geom_point(aes(x=V3, y=V4)) + 
    labs(x="PC3", y="PC4", title="OTref all SNPs - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ncol=2
)
dev.off()
#PFref
pdf(file = "PCA_PFref.pdf", paper = "a4", title = "PCA_PFref", pointsize = 6, width = 14)
grid.arrange(
  ggplot(as.data.frame(PF_all_PCA$projections)) + 
    geom_point(aes(x=V1, y=V2)) + 
    labs(x="PC1", y="PC2", title="PFref all SNPs - PC1 PC2") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_all_PCA$projections)) + 
    geom_point(aes(x=V1, y=V3)) + 
    labs(x="PC1", y="PC3", title="PFref all SNPs - PC1 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_all_PCA$projections)) + 
    geom_point(aes(x=V1, y=V4)) + 
    labs(x="PC1", y="PC4", title="PFref all SNPs - PC1 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_all_PCA$projections)) + 
    geom_point(aes(x=V2, y=V3)) + 
    labs(x="PC2", y="PC3", title="PFref all SNPs - PC2 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_all_PCA$projections)) + 
    geom_point(aes(x=V2, y=V4)) + 
    labs(x="PC2", y="PC4", title="PFref all SNPs - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_all_PCA$projections)) + 
    geom_point(aes(x=V3, y=V4)) + 
    labs(x="PC3", y="PC4", title="PFref all SNPs - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ncol=2
)
dev.off()
#
#Subset1 (no SNPs on scaffolds with marker SNPs)
#
#OTref
pdf(file = "PCA_OTref_no_SNPs_on_marker_scaffolds.pdf", paper = "a4", title = "PCA_OTref_no_SNPs_on_marker_scaffolds", pointsize = 6, width = 14)
grid.arrange(
  ggplot(as.data.frame(OT_subset1_PCA$projections)) + 
    geom_point(aes(x=V1, y=V2)) + 
    labs(x="PC1", y="PC2", title="OTref no SNPs on marker scaffolds - PC1 PC2") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset1_PCA$projections)) + 
    geom_point(aes(x=V1, y=V3)) + 
    labs(x="PC1", y="PC3", title="OTref no SNPs on marker scaffolds - PC1 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset1_PCA$projections)) + 
    geom_point(aes(x=V1, y=V4)) + 
    labs(x="PC1", y="PC4", title="OTref no SNPs on marker scaffolds - PC1 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset1_PCA$projections)) + 
    geom_point(aes(x=V2, y=V3)) + 
    labs(x="PC2", y="PC3", title="OTref no SNPs on marker scaffolds - PC2 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset1_PCA$projections)) + 
    geom_point(aes(x=V2, y=V4)) + 
    labs(x="PC2", y="PC4", title="OTref no SNPs on marker scaffolds - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset1_PCA$projections)) + 
    geom_point(aes(x=V3, y=V4)) + 
    labs(x="PC3", y="PC4", title="OTref no SNPs on marker scaffolds - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ncol=2
)
dev.off()
#PFref
pdf(file = "PCA_PFref_no_SNPs_on_marker_scaffolds.pdf", paper = "a4", title = "PCA_PFref_no_SNPs_on_marker_scaffolds", pointsize = 6, width = 14)
grid.arrange(
  ggplot(as.data.frame(PF_subset1_PCA$projections)) + 
    geom_point(aes(x=V1, y=V2)) + 
    labs(x="PC1", y="PC2", title="PFref no SNPs on marker scaffolds - PC1 PC2") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset1_PCA$projections)) + 
    geom_point(aes(x=V1, y=V3)) + 
    labs(x="PC1", y="PC3", title="PFref no SNPs on marker scaffolds - PC1 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset1_PCA$projections)) + 
    geom_point(aes(x=V1, y=V4)) + 
    labs(x="PC1", y="PC4", title="PFref no SNPs on marker scaffolds - PC1 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset1_PCA$projections)) + 
    geom_point(aes(x=V2, y=V3)) + 
    labs(x="PC2", y="PC3", title="PFref no SNPs on marker scaffolds - PC2 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset1_PCA$projections)) + 
    geom_point(aes(x=V2, y=V4)) + 
    labs(x="PC2", y="PC4", title="PFref no SNPs on marker scaffolds - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset1_PCA$projections)) + 
    geom_point(aes(x=V3, y=V4)) + 
    labs(x="PC3", y="PC4", title="PFref no SNPs on marker scaffolds - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ncol=2
)
dev.off()
#
#Subset2 (only SNPs on scaffolds with marker SNPs)
#
#OTref
pdf(file = "PCA_OTref_only_SNPs_on_marker_scaffolds.pdf", paper = "a4", title = "PCA_OTref_only_SNPs_on_marker_scaffolds", pointsize = 6, width = 14)
grid.arrange(
  ggplot(as.data.frame(OT_subset2_PCA$projections)) + 
    geom_point(aes(x=V1, y=V2)) + 
    labs(x="PC1", y="PC2", title="OTref only SNPs on marker scaffolds - PC1 PC2") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset2_PCA$projections)) + 
    geom_point(aes(x=V1, y=V3)) + 
    labs(x="PC1", y="PC3", title="OTref only SNPs on marker scaffolds - PC1 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset2_PCA$projections)) + 
    geom_point(aes(x=V1, y=V4)) + 
    labs(x="PC1", y="PC4", title="OTref only SNPs on marker scaffolds - PC1 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset2_PCA$projections)) + 
    geom_point(aes(x=V2, y=V3)) + 
    labs(x="PC2", y="PC3", title="OTref only SNPs on marker scaffolds - PC2 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset2_PCA$projections)) + 
    geom_point(aes(x=V2, y=V4)) + 
    labs(x="PC2", y="PC4", title="OTref only SNPs on marker scaffolds - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset2_PCA$projections)) + 
    geom_point(aes(x=V3, y=V4)) + 
    labs(x="PC3", y="PC4", title="OTref only SNPs on marker scaffolds - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ncol=2
)
dev.off()
#PFref
pdf(file = "PCA_PFref_only_SNPs_on_marker_scaffolds.pdf", paper = "a4", title = "PCA_PFref_only_SNPs_on_marker_scaffolds", pointsize = 6, width = 14)
grid.arrange(
  ggplot(as.data.frame(PF_subset2_PCA$projections)) + 
    geom_point(aes(x=V1, y=V2)) + 
    labs(x="PC1", y="PC2", title="PFref only SNPs on marker scaffolds - PC1 PC2") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset2_PCA$projections)) + 
    geom_point(aes(x=V1, y=V3)) + 
    labs(x="PC1", y="PC3", title="PFref only SNPs on marker scaffolds - PC1 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset2_PCA$projections)) + 
    geom_point(aes(x=V1, y=V4)) + 
    labs(x="PC1", y="PC4", title="PFref only SNPs on marker scaffolds - PC1 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset2_PCA$projections)) + 
    geom_point(aes(x=V2, y=V3)) + 
    labs(x="PC2", y="PC3", title="PFref only SNPs on marker scaffolds - PC2 PC3") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset2_PCA$projections)) + 
    geom_point(aes(x=V2, y=V4)) + 
    labs(x="PC2", y="PC4", title="PFref only SNPs on marker scaffolds - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset2_PCA$projections)) + 
    geom_point(aes(x=V3, y=V4)) + 
    labs(x="PC3", y="PC4", title="PFref only SNPs on marker scaffolds - PC2 PC4") + 
    theme(text = element_text(size=8)), 
  ncol=2
)
dev.off()
#
#Supplementary plots
#
#OTref
pdf(file = "PCA_OTref_additional_plots.pdf", paper = "a4r", title = "PCA_OTref_additional_plots", pointsize = 6, width = 14)
grid.arrange(
  #Standard deviation
  ggplot(as.data.frame(OT_all_PCA$sdev)) + 
    geom_point(aes(x=1:47, y=V1)) + 
    labs(x="Principal components", y="Standard deviation", title="OTref all SNPs") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset1_PCA$sdev)) + 
    geom_point(aes(x=1:47, y=V1)) + 
    labs(x="Principal components", y="Standard deviation", title="OTref no SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset2_PCA$sdev)) + 
    geom_point(aes(x=1:47, y=V1)) + 
    labs(x="Principal components", y="Standard deviation", title="OTref only SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  #Variance explained
  ggplot(as.data.frame(OT_all_tw$percentage)) + 
    geom_point(aes(x=1:47, y=OT_all_tw$percentage)) + 
    labs(x="Principal components", y="Proportion of variance explained", title="OTref all SNPs") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset1_all_tw$percentage)) + 
    geom_point(aes(x=1:47, y=OT_subset1_all_tw$percentage)) + 
    labs(x="Principal components", y="Proportion of variance explained", title="OTref no SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset2_all_tw$percentage)) + 
    geom_point(aes(x=1:47, y=OT_subset2_all_tw$percentage)) + 
    labs(x="Principal components", y="Proportion of variance explained", title="OTref only SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  #P-values
  ggplot(as.data.frame(OT_all_tw$pvalues)) + 
    geom_col(aes(x=1:47, y=OT_all_tw$pvalues)) + 
    labs(x="Principal components", y="p-values", title="OTref all SNPs") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset1_all_tw$pvalues)) + 
    geom_col(aes(x=1:47, y=OT_subset1_all_tw$pvalues)) + 
    labs(x="Principal components", y="p-values", title="OTref no SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(OT_subset2_all_tw$pvalues)) + 
    geom_col(aes(x=1:47, y=OT_subset2_all_tw$pvalues)) + 
    labs(x="Principal components", y="p-values", title="OTref only SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  ncol=3
)
dev.off()
#PFref
pdf(file = "PCA_PFref_additional_plots.pdf", paper = "a4r", title = "PCA_PFref_additional_plots", pointsize = 6, width = 14)
grid.arrange(
  #Standard deviation
  ggplot(as.data.frame(PF_all_PCA$sdev)) + 
    geom_point(aes(x=1:47, y=V1)) + 
    labs(x="Principal components", y="Standard deviation", title="PFref all SNPs") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset1_PCA$sdev)) + 
    geom_point(aes(x=1:47, y=V1)) + 
    labs(x="Principal components", y="Standard deviation", title="PFref no SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset2_PCA$sdev)) + 
    geom_point(aes(x=1:47, y=V1)) + 
    labs(x="Principal components", y="Standard deviation", title="PFref only SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  #Variance explained
  ggplot(as.data.frame(PF_all_tw$percentage)) + 
    geom_point(aes(x=1:47, y=PF_all_tw$percentage)) + 
    labs(x="Principal components", y="Proportion of variance explained", title="PFref all SNPs") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset1_all_tw$percentage)) + 
    geom_point(aes(x=1:47, y=PF_subset1_all_tw$percentage)) + 
    labs(x="Principal components", y="Proportion of variance explained", title="PFref no SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset2_all_tw$percentage)) + 
    geom_point(aes(x=1:47, y=PF_subset2_all_tw$percentage)) + 
    labs(x="Principal components", y="Proportion of variance explained", title="PFref only SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  #P-values
  ggplot(as.data.frame(PF_all_tw$pvalues)) + 
    geom_col(aes(x=1:47, y=PF_all_tw$pvalues)) + 
    labs(x="Principal components", y="p-values", title="PFref all SNPs") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset1_all_tw$pvalues)) + 
    geom_col(aes(x=1:47, y=PF_subset1_all_tw$pvalues)) + 
    labs(x="Principal components", y="p-values", title="PFref no SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  ggplot(as.data.frame(PF_subset2_all_tw$pvalues)) + 
    geom_col(aes(x=1:47, y=PF_subset2_all_tw$pvalues)) + 
    labs(x="Principal components", y="p-values", title="PFref only SNPs on marker scaffolds") + 
    theme(text = element_text(size=8)), 
  ncol=3
)
dev.off()