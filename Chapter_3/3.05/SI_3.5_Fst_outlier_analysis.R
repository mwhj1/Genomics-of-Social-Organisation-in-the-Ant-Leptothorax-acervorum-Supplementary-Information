### Fst from bcftools / vcftools ###
###### Social Marker Outliers ######

#install.packages("cumstats")
#install.packages("quantreg")
#install.packages("ggridges")
#library(cumstats)
library(tidyverse)
library(dplyr)
library(magrittr)
library(data.table)
library(scales)
library(quantreg)
library(ggridges)
library(purrr)

#######################
#Part 1 - prepare data#
#######################

#1. Read in data
setwd("/Users/u1693640/Documents/PhD/2020/bcftools_Fst")
#OTref Fst
OTref_NR_L_raw <- read_tsv("./OTref/Raw_Fst/OTrefNR_L_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
OTref_NR_PF_raw <- read_tsv("./OTref/Raw_Fst/OTrefNR_PF_fst.weir.fst", 
                            col_types = "cnn", na = "-nan")
OTref_NR_V_raw <- read_tsv("./OTref/Raw_Fst/OTrefNR_V_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
OTref_OT_L_raw <- read_tsv("./OTref/Raw_Fst/OTrefOT_L_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
OTref_OT_NR_raw <- read_tsv("./OTref/Raw_Fst/OTrefOT_NR_fst.weir.fst", 
                            col_types = "cnn", na = "-nan")
OTref_OT_PF_raw <- read_tsv("./OTref/Raw_Fst/OTrefOT_PF_fst.weir.fst", 
                            col_types = "cnn", na = "-nan")
OTref_OT_V_raw <- read_tsv("./OTref/Raw_Fst/OTrefOT_V_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
OTref_PF_L_raw <- read_tsv("./OTref/Raw_Fst/OTrefPF_L_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
OTref_PF_V_raw <- read_tsv("./OTref/Raw_Fst/OTrefPF_V_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
OTref_V_L_raw <- read_tsv("./OTref/Raw_Fst/OTrefV_L_fst.weir.fst", 
                          col_types = "cnn", na = "-nan")
#PFref Fst
PFref_NR_L_raw <- read_tsv("./PFref/Raw_Fst/PFrefNR_L_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
PFref_NR_PF_raw <- read_tsv("./PFref/Raw_Fst/PFrefNR_PF_fst.weir.fst", 
                            col_types = "cnn", na = "-nan")
PFref_NR_V_raw <- read_tsv("./PFref/Raw_Fst/PFrefNR_V_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
PFref_OT_L_raw <- read_tsv("./PFref/Raw_Fst/PFrefOT_L_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
PFref_OT_NR_raw <- read_tsv("./PFref/Raw_Fst/PFrefOT_NR_fst.weir.fst", 
                            col_types = "cnn", na = "-nan")
PFref_OT_PF_raw <- read_tsv("./PFref/Raw_Fst/PFrefOT_PF_fst.weir.fst", 
                            col_types = "cnn", na = "-nan")
PFref_OT_V_raw <- read_tsv("./PFref/Raw_Fst/PFrefOT_V_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
PFref_PF_L_raw <- read_tsv("./PFref/Raw_Fst/PFrefPF_L_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
PFref_PF_V_raw <- read_tsv("./PFref/Raw_Fst/PFrefPF_V_fst.weir.fst", 
                           col_types = "cnn", na = "-nan")
PFref_V_L_raw <- read_tsv("./PFref/Raw_Fst/PFrefV_L_fst.weir.fst", 
                          col_types = "cnn", na = "-nan")

#2. Create one dataframe per reference, containing all Fst data
#   from all pairwise comparisons

#First check that all raw Fst files contain identical positional information
#(CHROM and POS columns) - if they are all identical then this will be considerably easier
identical(OTref_NR_L_raw$CHROM, OTref_NR_V_raw$CHROM)
identical(OTref_NR_L_raw$POS, OTref_NR_V_raw$POS)
identical(OTref_PF_V_raw$CHROM, OTref_OT_PF_raw$CHROM)
identical(OTref_PF_V_raw$POS, OTref_OT_PF_raw$POS)
identical(PFref_NR_L_raw$CHROM, PFref_NR_V_raw$CHROM)
identical(PFref_NR_L_raw$POS, PFref_NR_V_raw$POS)
identical(PFref_PF_V_raw$CHROM, PFref_OT_PF_raw$CHROM)
identical(PFref_PF_V_raw$POS, PFref_OT_PF_raw$POS)
#These all came out TRUE and I think that is a reasonable check - moreover all
#raw Fst files have the same number of lines (within reference-specific sets).
#If the assumption is wrong I think things will break downstream.

#OTref:
#Start with the first set of raw Fst
OTref_MASTER_raw <- OTref_NR_L_raw
#Re-name Fst column to indicate which population pair it refers to
OTref_MASTER_raw <- rename(OTref_MASTER_raw, NR_L = WEIR_AND_COCKERHAM_FST)
#Then add in all remaining Fst data with renaming in one step
OTref_MASTER_raw <- OTref_MASTER_raw %>% add_column(NR_PF = OTref_NR_PF_raw$WEIR_AND_COCKERHAM_FST, 
                                                    NR_V = OTref_NR_V_raw$WEIR_AND_COCKERHAM_FST, 
                                                    OT_L = OTref_OT_L_raw$WEIR_AND_COCKERHAM_FST, 
                                                    OT_NR = OTref_OT_NR_raw$WEIR_AND_COCKERHAM_FST, 
                                                    OT_PF = OTref_OT_PF_raw$WEIR_AND_COCKERHAM_FST, 
                                                    OT_V = OTref_OT_V_raw$WEIR_AND_COCKERHAM_FST, 
                                                    PF_L = OTref_PF_L_raw$WEIR_AND_COCKERHAM_FST, 
                                                    PF_V = OTref_PF_V_raw$WEIR_AND_COCKERHAM_FST, 
                                                    V_L = OTref_V_L_raw$WEIR_AND_COCKERHAM_FST)
#Finally remove all rows with missing data
OTref_MASTER_raw <- OTref_MASTER_raw %>% drop_na()

#PFref (exactly the same):
#Start with the first set of raw Fst
PFref_MASTER_raw <- PFref_NR_L_raw
#Re-name Fst column to indicate which population pair it refers to
PFref_MASTER_raw <- rename(PFref_MASTER_raw, NR_L = WEIR_AND_COCKERHAM_FST)
#Then add in all remaining Fst data with renaming in one step
PFref_MASTER_raw <- PFref_MASTER_raw %>% add_column(NR_PF = PFref_NR_PF_raw$WEIR_AND_COCKERHAM_FST, 
                                                    NR_V = PFref_NR_V_raw$WEIR_AND_COCKERHAM_FST, 
                                                    OT_L = PFref_OT_L_raw$WEIR_AND_COCKERHAM_FST, 
                                                    OT_NR = PFref_OT_NR_raw$WEIR_AND_COCKERHAM_FST, 
                                                    OT_PF = PFref_OT_PF_raw$WEIR_AND_COCKERHAM_FST, 
                                                    OT_V = PFref_OT_V_raw$WEIR_AND_COCKERHAM_FST, 
                                                    PF_L = PFref_PF_L_raw$WEIR_AND_COCKERHAM_FST, 
                                                    PF_V = PFref_PF_V_raw$WEIR_AND_COCKERHAM_FST, 
                                                    V_L = PFref_V_L_raw$WEIR_AND_COCKERHAM_FST)
#Finally remove all rows with missing data
PFref_MASTER_raw <- PFref_MASTER_raw %>% drop_na()

#2.b It is useful to get some extra information into these
# dataframes now - scaffold length in particular.
# First write out a file of all scaffold IDs
# present in each dataset.  These will be used
# to query the assembly fasta files (externally)

OT_included_scaffolds <- unique(OTref_MASTER_raw$CHROM)
PF_included_scaffolds <- unique(PFref_MASTER_raw$CHROM)
#write.table(OT_included_scaffolds, file = "OT_included_scaffolds.txt", 
#            quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(PF_included_scaffolds, file = "PF_included_scaffolds.txt", 
#            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ... external python scripts ...

#Read back in the scaffold length files
OT_scaffold_lengths <- read.table("OT_included_scaffolds_lengths.txt", 
                                  sep = " ", stringsAsFactors = FALSE)
PF_scaffold_lengths <- read.table("PF_included_scaffolds_lengths.txt", 
                                  sep = " ", stringsAsFactors = FALSE)
colnames(OT_scaffold_lengths) <- c("CHROM", "LENGTH")
colnames(PF_scaffold_lengths) <- c("CHROM", "LENGTH")
#Add in Length as new column to existing master Fst dataframes
OTref_MASTER_raw <- OTref_MASTER_raw %>% dplyr::left_join(OT_scaffold_lengths)
PFref_MASTER_raw <- PFref_MASTER_raw %>% dplyr::left_join(PF_scaffold_lengths)

#2.c Alignments were done against unmasked references in 
# order to not bias where a read best aligns.  BED files
# of repeat-masked ranges in each reference will allow
# identification of any SNPs in low-complexity regions

OTref_repeat_ranges <- read.table("Masked_regions/clear_OT18_12_MP4_re-do_S2_assembly_3b.fasta.ranges.bed", 
                                  stringsAsFactors = FALSE)
PFref_repeat_ranges <- read.table("Masked_regions/clear_PF_15_M1_S3.fasta.ranges.bed", 
                                  stringsAsFactors = FALSE)
colnames(OTref_repeat_ranges) <- c("CHROM", "start", "stop")
colnames(PFref_repeat_ranges) <- c("CHROM", "start", "stop")

#A function to add a column indicating if a SNP is within
#a range defined in a bed file
#(warning - slow and horrible)
low.complex.index <- function(df.Fst, df.ranges){
  df.Fst$low_complexity <- 0
  df.ranges <- filter(df.ranges, CHROM %in% df.Fst$CHROM)
  range <- 1:nrow(df.Fst)
  for(i in range){
    scaff <- df.Fst$CHROM[i]
    print(c(i, scaff))
    set <- filter(df.ranges, CHROM == scaff)
    print(nrow(set))
    SNP.pos <- df.Fst$POS[i]
    tester <- nrow(filter(set, start <= SNP.pos, stop >= SNP.pos))
    if(tester > 0){
      df.Fst$low_complexity[i] <- 1
    }
  }
  return(df.Fst)
}

#Call the function
#OTref_MASTER_raw <- low.complex.index(OTref_MASTER_raw, OTref_repeat_ranges)
#PFref_MASTER_raw <- low.complex.index(PFref_MASTER_raw, PFref_repeat_ranges)
#Write the output
#write.table(OTref_MASTER_raw, file="OTref/OTref_Master.fst", quote=FALSE, row.names=FALSE)
#write.table(PFref_MASTER_raw, file="PFref/PFref_Master.fst", quote=FALSE, row.names=FALSE)
#Read the output back in
OTref_MASTER_raw <- read.table(file="OTref/OTref_Master.fst", header = TRUE, stringsAsFactors = FALSE)
PFref_MASTER_raw <- read.table(file="PFref/PFref_Master.fst", header = TRUE, stringsAsFactors = FALSE)

#2.d Reference genomes have been annotated with gene predictions
# from aligned pooled RNAseq data.  GFF files will be used
# to indicate whether a SNP overlaps with a gene (a binary index column as above)
# with another column indicating the ID of the gene (if any).
# For each reference, this means reading two files - 
# first a GFF file which contains information on where features
# are located (in which features are identified simply
# by a number), and second a csv file produced by GO_FEAT 
# which ties the GFF IDs to real gene names

#Read GFFs
OTref_GFF <- read.table(file = "Annotation/augustus.hints_genes.OT3b.gff3", 
                        header = FALSE, stringsAsFactors = FALSE)
PFref_GFF <- read.table(file = "Annotation/augustus.hints_genes.PF.gff3", 
                        header = FALSE, stringsAsFactors = FALSE)
colnames(OTref_GFF) <- c("SCAFF", "source", "feature", "start", "end", 
                         "score", "strand", "phase", "attributes")
colnames(PFref_GFF) <- c("SCAFF", "source", "feature", "start", "end", 
                         "score", "strand", "phase", "attributes")

#Read GO_FEAT csvs
OTref_GO_FEAT <- read.table(file = "Annotation/GO_FEAT_OT3b_all.csv", 
                            header = FALSE, stringsAsFactors = FALSE, 
                            sep = ";", fill = TRUE, quote = "") %>% select(2:11)
PFref_GO_FEAT <- read.table(file = "Annotation/GO_FEAT_PF_all.csv", 
                            header = FALSE, stringsAsFactors = FALSE, 
                            sep = ";", fill = TRUE, quote = "") %>% select(2:11)
csv_cols <- c("Locus tag", "Length", "Product", "Completeness", 
              "Gene onthology", "SEED", "Protein databases", 
              "Genome annotation databases", "Family and domain databases", 
              "Crossreferences")
colnames(OTref_GO_FEAT) <- csv_cols
colnames(PFref_GO_FEAT) <- csv_cols

#Filter GFFs to just rows with "gene" in the "feature" column
#(Each gene can have several features [exon, intron, etc], which are nested - 
#"gene" is the top-level feature which has start and end positions which
#cover the whole gene)
OTref_GFF <- filter(OTref_GFF, feature == "gene")
PFref_GFF <- filter(PFref_GFF, feature == "gene")

#Reformat "attributes" column of GFFs.  Currently they contain
#data in the form of (e.g.) "ID=g1;" - this needs to be changed to
#"g1" in order to match up with the IDs specified in the csv files
OTref_GFF <- tibble::as_tibble(OTref_GFF) %>% 
  mutate(temp1 = str_replace_all(attributes, 
                                 pattern = "ID=", 
                                 replacement = "")) %>% 
  mutate(geneID = str_replace_all(temp1, 
                                  pattern = ";", 
                                  replacement = "")) %>% 
  select(-all_of(c("attributes", "temp1")))
PFref_GFF <- tibble::as_tibble(PFref_GFF) %>% 
  mutate(temp1 = str_replace_all(attributes, 
                                 pattern = "ID=", 
                                 replacement = "")) %>% 
  mutate(geneID = str_replace_all(temp1, 
                                  pattern = ";", 
                                  replacement = "")) %>% 
  select(-all_of(c("attributes", "temp1")))

#Split the "Locus tag" column in the GFFs into two.
#Currently data there is formated like (e.g.) "g7239.t1" - 
#this will be split into two columns - a geneID column
#and a transcript column (.t1 = transcript 1, .t2 = transcript 2, etc)
OTref_GO_FEAT <- OTref_GO_FEAT %>% separate(`Locus tag`, into = c("geneID", "transcript"))
PFref_GO_FEAT <- PFref_GO_FEAT %>% separate(`Locus tag`, into = c("geneID", "transcript"))
#Drop everything that is not the primary transcript
OTref_GO_FEAT <- filter(OTref_GO_FEAT, transcript == "t1")
PFref_GO_FEAT <- filter(PFref_GO_FEAT, transcript == "t1")

#Sanity check (should return TRUE)
nrow(OTref_GFF) == nrow(OTref_GO_FEAT) & nrow(PFref_GFF) == nrow(PFref_GO_FEAT)

#Join each reference's GFF and GO_FEAT csv files together
OTref_annotations <- inner_join(OTref_GFF, OTref_GO_FEAT)
PFref_annotations <- inner_join(PFref_GFF, PFref_GO_FEAT)
#Drop the empty "SEED" column
OTref_annotations <- select(OTref_annotations, -c(SEED))
PFref_annotations <- select(PFref_annotations, -c(SEED))
#And write to file
write.table(OTref_annotations, file = "Annotation/OTref_annotations.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PFref_annotations, file = "Annotation/PFref_annotations.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)

#A function to index SNPs 1 (intersects with a gene) or 0 (does not intersect with a gene)
gene.index <- function(df.Fst, df.Annotation){
  #Create all 0 column (specific rows will be set to 1 as conditions are satisfied below)
  df.Fst$gene.index <- 0
  #Create an all "NA" column for gene ID
  df.Fst$gene.ID <- "NA"
  #Filter annotations to just those relating to scaffolds in df.Fst
  df.Annotation <- filter(df.Annotation, SCAFF %in% df.Fst$CHROM)
  #Set a range to iterate over
  range <- 1:nrow(df.Fst)
  #Iterate over range
  for(i in range){
    #Identify the scaffold
    scaff <- df.Fst$CHROM[i]
    print(c(i, scaff))
    #Subset annotations to relevant scaffold
    set <- filter(df.Annotation, SCAFF == scaff)
    print(nrow(set))
    #Get the SNP position
    SNP.pos <- df.Fst$POS[i]
    #Define a test - will be 0 if there are no intersections with a gene range
    tester <- nrow(filter(set, start <= SNP.pos, end >= SNP.pos))
    #A condition of tester's value to determine whether the gene.index needs to be set to 1
    if(tester > 0){
      #Set to 1 if condition is violated
      df.Fst$gene.index[i] <- 1
      #Store the rows of set which violate the condition
      set.subset <- filter(set, start <= SNP.pos, end >= SNP.pos)
      #Record all gene IDs
      print(set.subset$geneID)
      #Store as a list in the relevant cell of the relevant column
      df.Fst$gene.ID[i] <- as.list(set.subset$geneID)
    }
  }
  return(df.Fst)
}

#Call the funciton (HORRIBLY SLOW)
#OTref_MASTER_raw_Fst_ann <- gene.index(OTref_MASTER_raw, OTref_annotations)
#PFref_MASTER_raw_Fst_ann <- gene.index(PFref_MASTER_raw, PFref_annotations)
#OTref_MASTER_raw_Fst_ann %<>% mutate_if(is.list, as.character)
#PFref_MASTER_raw_Fst_ann %<>% mutate_if(is.list, as.character)
#Write to file
#write.table(OTref_MASTER_raw_Fst_ann, file = "OTref/OTref_MASTER_raw_Fst_annotations.tsv", 
#            quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(PFref_MASTER_raw_Fst_ann, file = "PFref/PFref_MASTER_raw_Fst_annotations.tsv", 
#            quote = FALSE, sep = "\t", row.names = FALSE)

#Read back in (if necessary)
OTref_MASTER_raw <- read.table(file="OTref/OTref_MASTER_raw_Fst_annotations.tsv", header = TRUE, stringsAsFactors = FALSE)
PFref_MASTER_raw <- read.table(file="PFref/PFref_MASTER_raw_Fst_annotations.tsv", header = TRUE, stringsAsFactors = FALSE)

#2.e A sanity check - have the last two steps worked properly?
# This can be checked by looking at data for one scaffold
# per reference and manually checking whether SNPs are
# correctly identified as being in a masked region, or
# as intersecting a gene

#Get a pair of scaffolds to test
OT_test_scaff <- OTref_MASTER_raw$CHROM[1]
PF_test_scaff <- PFref_MASTER_raw$CHROM[1]
#Get the matching subset of all data for the test scaffolds
OT_Fst_test_scaff <- OTref_MASTER_raw[which(OTref_MASTER_raw$CHROM == OT_test_scaff), ]
PF_Fst_test_scaff <- PFref_MASTER_raw[which(PFref_MASTER_raw$CHROM == PF_test_scaff), ]
OT_anno_test_scaff <- OTref_annotations[which(OTref_annotations$SCAFF == OT_test_scaff), ]
PF_anno_test_scaff <- PFref_annotations[which(PFref_annotations$SCAFF == PF_test_scaff), ]
OT_repeats_test_scaff <- OTref_repeat_ranges[which(OTref_repeat_ranges$CHROM == OT_test_scaff), ]
PF_repeats_test_scaff <- PFref_repeat_ranges[which(PFref_repeat_ranges$CHROM == PF_test_scaff), ]

#Manually inspect the test_scaff subsets.... they look good!

#Perform a final join to include other annotation data where possible
OTref_MASTER_raw <- left_join(OTref_MASTER_raw, OTref_annotations, by = c("gene.ID" = "geneID"))
PFref_MASTER_raw <- left_join(PFref_MASTER_raw, PFref_annotations, by = c("gene.ID" = "geneID"))
#Write to file
write.table(OTref_MASTER_raw, file = "OTref/OTref_MASTER_raw_Fst_annotations.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PFref_MASTER_raw, file = "PFref/PFref_MASTER_raw_Fst_annotations.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)



#######################
#Part 2 - analyse data#
#######################

#1. Read in data
setwd("/Users/u1693640/Documents/PhD/2020/bcftools_Fst")
#Read in combined Fst and annotations
OTref_MASTER_raw <- read.table(file="OTref/OTref_MASTER_raw_Fst_annotations.tsv", 
                               header = TRUE, stringsAsFactors = FALSE, sep = "\t", fill = TRUE, 
                               quote = "")
PFref_MASTER_raw <- read.table(file="PFref/PFref_MASTER_raw_Fst_annotations.tsv", 
                               header = TRUE, stringsAsFactors = FALSE, sep = "\t", fill = TRUE, 
                               quote = "")

#And the scaffold length data
OT_scaffold_lengths <- read.table("OT_included_scaffolds_lengths.txt", 
                                  sep = " ", stringsAsFactors = FALSE)
PF_scaffold_lengths <- read.table("PF_included_scaffolds_lengths.txt", 
                                  sep = " ", stringsAsFactors = FALSE)
colnames(OT_scaffold_lengths) <- c("CHROM", "LENGTH")
colnames(PF_scaffold_lengths) <- c("CHROM", "LENGTH")

#Change any Fst values <0 to 0
OTref_MASTER_raw[,3:12][OTref_MASTER_raw[, 3:12] < 0] <- 0.0
PFref_MASTER_raw[,3:12][PFref_MASTER_raw[, 3:12] < 0] <- 0.0

#Drop all columns containing Fst data pertaining to a comparison against NR
OTref_MASTER_raw <- OTref_MASTER_raw %>% select(!contains("NR"))
PFref_MASTER_raw <- PFref_MASTER_raw %>% select(!contains("NR"))

#2. Establish Fst thresholds to score outliers, using the 
# empirical distribution of Fst in each pairwise comparison

#A function to return a vector of quantiles in 5% steps
get_quantiles <- function(x){
  my_quantiles <- quantile(as.numeric(x), probs = seq(0,1,0.05))
  return(my_quantiles)
}
#OTref
OTref_qs_OT_L <- get_quantiles(OTref_MASTER_raw$OT_L)
OTref_qs_OT_PF <- get_quantiles(OTref_MASTER_raw$OT_PF)
OTref_qs_OT_V <- get_quantiles(OTref_MASTER_raw$OT_V)
OTref_qs_PF_L <- get_quantiles(OTref_MASTER_raw$PF_L)
OTref_qs_PF_V <- get_quantiles(OTref_MASTER_raw$PF_V)
OTref_qs_V_L <- get_quantiles(OTref_MASTER_raw$V_L)
#PFref
PFref_qs_OT_L <- get_quantiles(PFref_MASTER_raw$OT_L)
PFref_qs_OT_PF <- get_quantiles(PFref_MASTER_raw$OT_PF)
PFref_qs_OT_V <- get_quantiles(PFref_MASTER_raw$OT_V)
PFref_qs_PF_L <- get_quantiles(PFref_MASTER_raw$PF_L)
PFref_qs_PF_V <- get_quantiles(PFref_MASTER_raw$PF_V)
PFref_qs_V_L <- get_quantiles(PFref_MASTER_raw$V_L)
#The threshold values are the 13th (60%) and 20th (95%)
#values in the above vectors - these will be stored
upper_q <- 20
lower_q <- 13
#OTref
OTref_qUpper_OT_L <- OTref_qs_OT_L[[upper_q]]
OTref_qLower_OT_L <- OTref_qs_OT_L[[lower_q]]
OTref_qUpper_OT_PF <- OTref_qs_OT_PF[[upper_q]]
OTref_qLower_OT_PF <- OTref_qs_OT_PF[[lower_q]]
OTref_qUpper_OT_V <- OTref_qs_OT_V[[upper_q]]
OTref_qLower_OT_V <- OTref_qs_OT_V[[lower_q]]
OTref_qUpper_PF_L <- OTref_qs_PF_L[[upper_q]]
OTref_qLower_PF_L <- OTref_qs_PF_L[[lower_q]]
OTref_qUpper_PF_V <- OTref_qs_PF_V[[upper_q]]
OTref_qLower_PF_V <- OTref_qs_PF_V[[lower_q]]
OTref_qUpper_V_L <- OTref_qs_V_L[[upper_q]]
OTref_qLower_V_L <- OTref_qs_V_L[[lower_q]]
#PFref
PFref_qUpper_OT_L <- PFref_qs_OT_L[[upper_q]]
PFref_qLower_OT_L <- PFref_qs_OT_L[[lower_q]]
PFref_qUpper_OT_PF <- PFref_qs_OT_PF[[upper_q]]
PFref_qLower_OT_PF <- PFref_qs_OT_PF[[lower_q]]
PFref_qUpper_OT_V <- PFref_qs_OT_V[[upper_q]]
PFref_qLower_OT_V <- PFref_qs_OT_V[[lower_q]]
PFref_qUpper_PF_L <- PFref_qs_PF_L[[upper_q]]
PFref_qLower_PF_L <- PFref_qs_PF_L[[lower_q]]
PFref_qUpper_PF_V <- PFref_qs_PF_V[[upper_q]]
PFref_qLower_PF_V <- PFref_qs_PF_V[[lower_q]]
PFref_qUpper_V_L <- PFref_qs_V_L[[upper_q]]
PFref_qLower_V_L <- PFref_qs_V_L[[lower_q]]

#3. Add a new column to each master dataframe
# to indicate whether a SNP satisfies all
# conditions to be scored as an outlier
#
#OTref
OTref_MASTER_raw <- OTref_MASTER_raw %>% 
  mutate(outlier = ifelse(OT_L > OTref_qUpper_OT_L & 
                            OT_PF > OTref_qUpper_OT_PF & 
                            PF_V > OTref_qUpper_PF_V & 
                            V_L > OTref_qUpper_V_L & 
                            OT_V < OTref_qLower_OT_V & 
                            PF_L < OTref_qLower_PF_L, 
                          "outlier", "background"))
#PFref
PFref_MASTER_raw <- PFref_MASTER_raw %>% 
  mutate(outlier = ifelse(OT_L > PFref_qUpper_OT_L & 
                            OT_PF > PFref_qUpper_OT_PF & 
                            PF_V > PFref_qUpper_PF_V & 
                            V_L > PFref_qUpper_V_L & 
                            OT_V < PFref_qLower_OT_V & 
                            PF_L < PFref_qLower_PF_L, 
                          "outlier", "background"))
#How many outliers in each reference?
sum(OTref_MASTER_raw$outlier == "outlier")
sum(PFref_MASTER_raw$outlier == "outlier")
#Get just the outliers
OTref_outliers <- OTref_MASTER_raw[which(OTref_MASTER_raw$outlier == "outlier"),]
PFref_outliers <- PFref_MASTER_raw[which(PFref_MASTER_raw$outlier == "outlier"),]
#Count outliers in scaffolds
OTref_outlier_counts_per_scaffold <- as.data.frame(table(OTref_outliers$CHROM)) %>% arrange(-Freq)
colnames(OTref_outlier_counts_per_scaffold) <- c("CHROM", "outlier SNPs")
PFref_outlier_counts_per_scaffold <- as.data.frame(table(PFref_outliers$CHROM)) %>% arrange(-Freq)
colnames(PFref_outlier_counts_per_scaffold) <- c("CHROM", "outlier SNPs")
#Add scaffold length
OTref_outlier_counts_per_scaffold <- OTref_outlier_counts_per_scaffold %>% left_join(OT_scaffold_lengths)
PFref_outlier_counts_per_scaffold <- PFref_outlier_counts_per_scaffold %>% left_join(PF_scaffold_lengths)

#4. Repeating above but selecting the 99% quantile as the upper threshold
get_quantiles_2 <- function(x){
  my_quantiles <- quantile(as.numeric(x), probs = seq(0,1,0.01))
  return(my_quantiles)
}
#OTref
OTref_qs2_OT_L <- get_quantiles_2(OTref_MASTER_raw$OT_L)
OTref_qs2_OT_PF <- get_quantiles_2(OTref_MASTER_raw$OT_PF)
OTref_qs2_OT_V <- get_quantiles_2(OTref_MASTER_raw$OT_V)
OTref_qs2_PF_L <- get_quantiles_2(OTref_MASTER_raw$PF_L)
OTref_qs2_PF_V <- get_quantiles_2(OTref_MASTER_raw$PF_V)
OTref_qs2_V_L <- get_quantiles_2(OTref_MASTER_raw$V_L)
#PFref
PFref_qs2_OT_L <- get_quantiles_2(PFref_MASTER_raw$OT_L)
PFref_qs2_OT_PF <- get_quantiles_2(PFref_MASTER_raw$OT_PF)
PFref_qs2_OT_V <- get_quantiles_2(PFref_MASTER_raw$OT_V)
PFref_qs2_PF_L <- get_quantiles_2(PFref_MASTER_raw$PF_L)
PFref_qs2_PF_V <- get_quantiles_2(PFref_MASTER_raw$PF_V)
PFref_qs2_V_L <- get_quantiles_2(PFref_MASTER_raw$V_L)

#The threshold values are the 61st (60%) and 100th (99%)
#values in the above vectors - these will be stored
upper_q <- 100
lower_q <- 61
#OTref
OTref_qUpper2_OT_L <- OTref_qs2_OT_L[[upper_q]]
OTref_qLower2_OT_L <- OTref_qs2_OT_L[[lower_q]]
OTref_qUpper2_OT_PF <- OTref_qs2_OT_PF[[upper_q]]
OTref_qLower2_OT_PF <- OTref_qs2_OT_PF[[lower_q]]
OTref_qUpper2_OT_V <- OTref_qs2_OT_V[[upper_q]]
OTref_qLower2_OT_V <- OTref_qs2_OT_V[[lower_q]]
OTref_qUpper2_PF_L <- OTref_qs2_PF_L[[upper_q]]
OTref_qLower2_PF_L <- OTref_qs2_PF_L[[lower_q]]
OTref_qUpper2_PF_V <- OTref_qs2_PF_V[[upper_q]]
OTref_qLower2_PF_V <- OTref_qs2_PF_V[[lower_q]]
OTref_qUpper2_V_L <- OTref_qs2_V_L[[upper_q]]
OTref_qLower2_V_L <- OTref_qs2_V_L[[lower_q]]
#PFref
PFref_qUpper2_OT_L <- PFref_qs2_OT_L[[upper_q]]
PFref_qLower2_OT_L <- PFref_qs2_OT_L[[lower_q]]
PFref_qUpper2_OT_PF <- PFref_qs2_OT_PF[[upper_q]]
PFref_qLower2_OT_PF <- PFref_qs2_OT_PF[[lower_q]]
PFref_qUpper2_OT_V <- PFref_qs2_OT_V[[upper_q]]
PFref_qLower2_OT_V <- PFref_qs2_OT_V[[lower_q]]
PFref_qUpper2_PF_L <- PFref_qs2_PF_L[[upper_q]]
PFref_qLower2_PF_L <- PFref_qs2_PF_L[[lower_q]]
PFref_qUpper2_PF_V <- PFref_qs2_PF_V[[upper_q]]
PFref_qLower2_PF_V <- PFref_qs2_PF_V[[lower_q]]
PFref_qUpper2_V_L <- PFref_qs2_V_L[[upper_q]]
PFref_qLower2_V_L <- PFref_qs2_V_L[[lower_q]]

#4.b Add a new column to each master dataframe
# to indicate whether a SNP satisfies all
# conditions to be scored as an outlier
#
#OTref
OTref_MASTER_raw <- OTref_MASTER_raw %>% 
  mutate(outlier2 = ifelse(OT_L > OTref_qUpper2_OT_L & 
                             OT_PF > OTref_qUpper2_OT_PF & 
                             PF_V > OTref_qUpper2_PF_V & 
                             V_L > OTref_qUpper2_V_L & 
                             OT_V < OTref_qLower2_OT_V & 
                             PF_L < OTref_qLower2_PF_L, 
                           "outlier", "background"))
#PFref
PFref_MASTER_raw <- PFref_MASTER_raw %>% 
  mutate(outlier2 = ifelse(OT_L > PFref_qUpper2_OT_L & 
                             OT_PF > PFref_qUpper2_OT_PF & 
                             PF_V > PFref_qUpper2_PF_V & 
                             V_L > PFref_qUpper2_V_L & 
                             OT_V < PFref_qLower2_OT_V & 
                             PF_L < PFref_qLower2_PF_L, 
                           "outlier", "background"))
#How many outliers in each reference?
sum(OTref_MASTER_raw$outlier2 == "outlier")
sum(PFref_MASTER_raw$outlier2 == "outlier")
#Get just the outliers
OTref_outliers2 <- OTref_MASTER_raw[which(OTref_MASTER_raw$outlier2 == "outlier"),]
PFref_outliers2 <- PFref_MASTER_raw[which(PFref_MASTER_raw$outlier2 == "outlier"),]
#Count outliers in scaffolds
OTref_outlier2_counts_per_scaffold <- as.data.frame(table(OTref_outliers2$CHROM)) %>% arrange(-Freq)
colnames(OTref_outlier2_counts_per_scaffold) <- c("CHROM", "outlier SNPs")
PFref_outlier2_counts_per_scaffold <- as.data.frame(table(PFref_outliers2$CHROM)) %>% arrange(-Freq)
colnames(PFref_outlier2_counts_per_scaffold) <- c("CHROM", "outlier SNPs")
#Add scaffold length
OTref_outlier2_counts_per_scaffold <- OTref_outlier2_counts_per_scaffold %>% left_join(OT_scaffold_lengths)
PFref_outlier2_counts_per_scaffold <- PFref_outlier2_counts_per_scaffold %>% left_join(PF_scaffold_lengths)

#A function to add total SNP count per scaffold as a column 
snp.counter <- function(perScaff.fst, full.fst){
  #Set all 0
  perScaff.fst$SNP.count <- 0
  #Iterate over scaffolds
  for(i in 1:nrow(perScaff.fst)){
    #Replace 0s with SNP counts
    perScaff.fst$SNP.count[i] <- nrow(full.fst[which(full.fst$CHROM == perScaff.fst$CHROM[i]), ])
  }
  #Return result
  return(perScaff.fst)
}
#Call the function
OTref_outlier_counts_per_scaffold <- snp.counter(OTref_outlier_counts_per_scaffold, OTref_MASTER_raw)
OTref_outlier2_counts_per_scaffold <- snp.counter(OTref_outlier2_counts_per_scaffold, OTref_MASTER_raw)
PFref_outlier_counts_per_scaffold <- snp.counter(PFref_outlier_counts_per_scaffold, PFref_MASTER_raw)
PFref_outlier2_counts_per_scaffold <- snp.counter(PFref_outlier2_counts_per_scaffold, PFref_MASTER_raw)
#Outliers as a % of all SNPs per scaffold
OTref_outlier_counts_per_scaffold <- mutate(OTref_outlier_counts_per_scaffold, perc.outliers = (`outlier SNPs`/SNP.count)*100)
OTref_outlier2_counts_per_scaffold <- mutate(OTref_outlier2_counts_per_scaffold, perc.outliers = (`outlier SNPs`/SNP.count)*100)
PFref_outlier_counts_per_scaffold <- mutate(PFref_outlier_counts_per_scaffold, perc.outliers = (`outlier SNPs`/SNP.count)*100)
PFref_outlier2_counts_per_scaffold <- mutate(PFref_outlier2_counts_per_scaffold, perc.outliers = (`outlier SNPs`/SNP.count)*100)

# 4.c Median within- and between-phenotype per-SNP Fst
#This is a useful summary statistic and is good for 
#plotting (two points plotted per SNP as opposed to 6)
OTref_MASTER_raw <- OTref_MASTER_raw %>% rowwise() %>% mutate(median.within = median(c(OT_V, PF_L)), median.between = median(c(OT_L, OT_PF, PF_V, V_L)))
PFref_MASTER_raw <- PFref_MASTER_raw %>% rowwise() %>% mutate(median.within = median(c(OT_V, PF_L)), median.between = median(c(OT_L, OT_PF, PF_V, V_L)))

#5. #Get outlier SNPs which intersect with genes
#OTref
OTref_95_genes <- OTref_MASTER_raw[which(OTref_MASTER_raw$outlier == "outlier" & 
                                           OTref_MASTER_raw$gene.index == 1), ] %>% 
  select(., c(CHROM, POS, LENGTH, low_complexity, gene.ID, source, feature, start, end, score, strand, 
              transcript, Length, Product, Completeness, Gene.onthology, Protein.databases, Family.and.domain.databases, 
              Crossreferences, OT_L, OT_PF, OT_V, PF_L, PF_V, V_L, 
              median.within, median.between))
OTref_99_genes <- OTref_MASTER_raw[which(OTref_MASTER_raw$outlier2 == "outlier" & 
                                           OTref_MASTER_raw$gene.index == 1), ] %>% 
  select(., c(CHROM, POS, LENGTH, low_complexity, gene.ID, source, feature, start, end, score, strand, 
              transcript, Length, Product, Completeness, Gene.onthology, Protein.databases, Family.and.domain.databases, 
              Crossreferences, OT_L, OT_PF, OT_V, PF_L, PF_V, V_L, 
              median.within, median.between))
#PFref
PFref_95_genes <- PFref_MASTER_raw[which(PFref_MASTER_raw$outlier == "outlier" & 
                                           PFref_MASTER_raw$gene.index == 1), ] %>% 
  select(., c(CHROM, POS, LENGTH, low_complexity, gene.ID, source, feature, start, end, score, strand, 
              transcript, Length, Product, Completeness, Gene.onthology, Protein.databases, Family.and.domain.databases, 
              Crossreferences, OT_L, OT_PF, OT_V, PF_L, PF_V, V_L, 
              median.within, median.between))
PFref_99_genes <- PFref_MASTER_raw[which(PFref_MASTER_raw$outlier2 == "outlier" & 
                                           PFref_MASTER_raw$gene.index == 1), ] %>% 
  select(., c(CHROM, POS, LENGTH, low_complexity, gene.ID, source, feature, start, end, score, strand, 
              transcript, Length, Product, Completeness, Gene.onthology, Protein.databases, Family.and.domain.databases, 
              Crossreferences, OT_L, OT_PF, OT_V, PF_L, PF_V, V_L, 
              median.within, median.between))
#Re-order based on the frequency of each scaffold ("CHROM")
OTref_95_genes <- OTref_95_genes %>% add_count(CHROM) %>% 
  arrange(desc(n)) %>% 
  select(-n)
OTref_99_genes <- OTref_99_genes %>% add_count(CHROM) %>% 
  arrange(desc(n)) %>% 
  select(-n)
PFref_95_genes <- PFref_95_genes %>% add_count(CHROM) %>% 
  arrange(desc(n)) %>% 
  select(-n)
PFref_99_genes <- PFref_99_genes %>% add_count(CHROM) %>% 
  arrange(desc(n)) %>% 
  select(-n)

#6. Write some results files
write.table(OTref_MASTER_raw, file = "Results/OTref_MASTER_Fst_results_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PFref_MASTER_raw, file = "Results/PFref_MASTER_Fst_results_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(OTref_outlier_counts_per_scaffold, file = "Results/OTref_Fst_per_scaff_95_60_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PFref_outlier_counts_per_scaffold, file = "Results/PFref_Fst_per_scaff_95_60_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(OTref_outlier2_counts_per_scaffold, file = "Results/OTref_Fst_per_scaff_99_60_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PFref_outlier2_counts_per_scaffold, file = "Results/PFref_Fst_per_scaff_99_60_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(OTref_95_genes, file = "Results/OTref_95_genes_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(OTref_99_genes, file = "Results/OTref_99_genes_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PFref_95_genes, file = "Results/PFref_95_genes_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(PFref_99_genes, file = "Results/PFref_99_genes_NO_NR.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE)