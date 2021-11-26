############################################################################
################## Sliding Window Diversity Statistics #####################
############################################################################

#Adapted from https://markravinet.github.io/inferring_selection_part1.html
#and https://markravinet.github.io/Chapter8.html
#Using scripts from https://github.com/simonhmartin/genomics_general
#
#Input data: VCFs used in Fst outlier analysis (representative diploids, biallalic SNPs)
#
#VCFs conversion to geno:
#vcf=OTref.biallelic.DP3g1maf05.recode.vcf.gz
#python parseVCF.py \
#-i $vcf | bgzip > OTref.biallelic.DP3g1maf05.recode.geno.gz
#
#Samples and populations files (PFref as example):
#bcftools query -l PFref.biallelic.DP3g1maf05.recode.vcf.gz > PFref.samples
#grep "^NR18" PFref.samples > PFref.samples.NR
#grep "^PF18" PFref.samples > PFref.samples.PF
#grep "^V18" PFref.samples > PFref.samples.V
#grep "^L18" PFref.samples > PFref.samples.L
#grep "^OT18" PFref.samples > PFref.samples.OT
#grep "^OT_18" PFref.samples >> PFref.samples.OT
#awk '{print $1"\tNR"}' PFref.samples.NR > PFref.pop_file
#awk '{print $1"\tPF"}' PFref.samples.PF >> PFref.pop_file
#awk '{print $1"\tV"}' PFref.samples.V >> PFref.pop_file
#awk '{print $1"\tL"}' PFref.samples.L >> PFref.pop_file
#awk '{print $1"\tOT"}' PFref.samples.OT >> PFref.pop_file
#
#Sliding window analysis:
#python popgenWindows.py \
#-g $infile -f phased -w 100000 -m 10 -s 25000 \
#-p NR -p PF -p V -p L -p OT \
#--popsFile $popfile \
#--writeFailedWindow \
#-T 4

############################################################################
############################################################################
############################################################################

#0. Environment setup
rm(list = ls())
#install.packages("gridExtra")
#install.packages("ggridges")
#install.packages("RColorBrewer")
#install.packages("grid")
library(grid)
library(ggridges)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
maindir <- "/Users/u1693640/Documents/PhD/2020/bcftools_Fst/SlidingWindowDiversity"
setwd(maindir)

############################################################################

#1. Read data and prepare for plotting
OT_infile <- "./Data/1a_sliding_window_OT.out"
PF_infile <- "./Data/1b_sliding_window_PF.out"
OTdata <- read_csv(OT_infile, col_types = "cnnnnddddddddddddddddddddddddd")
PFdata <- read_csv(PF_infile, col_types = "cnnnnddddddddddddddddddddddddd")
#Both of these files contain 1 window which contains no data 
#These must be removed
OTdata <- OTdata[-18981, ]
PFdata <- PFdata[-8966, ]

############################################################################

#2. Plotting functions
#
#Required data
btwn_Fst <- c("Fst_NR_PF", "Fst_NR_L", "Fst_PF_V", "Fst_PF_OT", "Fst_V_L", "Fst_L_OT")
wthn_Fst <- c("Fst_NR_V", "Fst_NR_OT", "Fst_PF_L", "Fst_V_OT")
btwn_dxy <- c("dxy_NR_PF", "dxy_NR_L", "dxy_PF_V", "dxy_PF_OT", "dxy_V_L", "dxy_L_OT")
wthn_dxy <- c("dxy_NR_V", "dxy_NR_OT", "dxy_PF_L", "dxy_V_OT")
pi_FM <- c("pi_NR", "pi_V", "pi_OT")
pi_P <- c("pi_PF", "pi_L")
#
#Some more required data (outlier SNP counts per scaffold from Fst outlier analysis)
OT_counts_per_scaffold <- read.table(file="./Data/OTref_Fst_per_scaff_95_60.tsv", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_counts_per_scaffold <- read.table(file="./Data/PFref_Fst_per_scaff_95_60.tsv", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#
#Subset to specified scaffold (base for subsequent subset functions)
scaffold.subset <- function(x, target.scaffold){
  return(x %>% filter(scaffold == target.scaffold) %>% select(-c(scaffold, start, end, sites)))
}
#Subset to Fst
fst.subset <- function(x){
  return(x %>% select(mid, contains("Fst")) %>% 
           gather(-mid, key = "stat", value = "value") %>% 
           mutate(
             type = case_when(
               stat %in% btwn_Fst ~ "between", 
               stat %in% wthn_Fst ~ "within"
             )
           ))
}
#Subset to dXY
dxy.subset <- function(x){
  return(x %>% select(mid, contains("dxy")) %>% 
           gather(-mid, key = "stat", value = "value") %>% 
           mutate(
             type = case_when(
               stat %in% btwn_dxy ~ "between", 
               stat %in% wthn_dxy ~ "within"
             )
           ))
}
#Subset to Pi
pi.subset <- function(x){
  return(x %>% select(mid, contains("pi")) %>% 
           gather(-mid, key = "stat", value = "value") %>% 
           mutate(
             type=case_when(
               stat %in% pi_FM ~ "FM", 
               stat %in% pi_P ~ "P"
             )
           ))
}
#
#Combine and plot
plot.scaffold.diversity <- function(data, target.scaffold, filename, title, reference.name){
  #Subset to scaffold of interest
  data=scaffold.subset(data, target.scaffold)
  #Create further subsets
  fst.data=fst.subset(data)
  dxy.data=dxy.subset(data)
  pi.data=pi.subset(data)
  #Open pdf for plotting
  pdf(file = filename, paper = "a4", title = title, pointsize = 6, width = 14, height = 14)
  #Plot...
  #Fst:
  grid.draw(rbind(
    ggplotGrob(ggplot(fst.data %>% filter(type == "between"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle(paste(reference.name, target.scaffold, sep=" ")) + 
                 scale_color_brewer(palette = "Paired") + 
                 ylab(expression(italic(F)[ST])) + 
                 ylim(c(min(fst.data$value)), max(fst.data$value)) + 
                 theme(legend.title=element_blank())), 
    ggplotGrob(ggplot(fst.data %>% filter(type == "within"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "RdYlBu") + 
                 ylab(expression(italic(F)[ST])) + 
                 ylim(c(min(fst.data$value)), max(fst.data$value)) + 
                 theme(legend.title=element_blank())), 
    #dXY:
    ggplotGrob(ggplot(dxy.data %>% filter(type == "between"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "Paired") + 
                 ylab("dXY") + 
                 ylim(c(min(dxy.data$value)), max(dxy.data$value)) + 
                 theme(legend.title=element_blank())), 
    ggplotGrob(ggplot(dxy.data %>% filter(type == "within"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "RdYlBu") + 
                 ylab("dXY") + 
                 ylim(c(min(dxy.data$value)), max(dxy.data$value)) + 
                 theme(legend.title=element_blank())), 
    #Pi:
    ggplotGrob(ggplot(pi.data %>% filter(type == "P"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "Paired") + 
                 ylab("Pi") + 
                 ylim(c(min(pi.data$value)), max(pi.data$value)) + 
                 theme(legend.title=element_blank())), 
    ggplotGrob(ggplot(pi.data %>% filter(type == "FM"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("Position (Mb)") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "RdYlBu") + 
                 ylab("Pi") + 
                 ylim(c(min(pi.data$value)), max(pi.data$value)) + 
                 theme(legend.title=element_blank())), 
    size="first")
  )
  #Close the file connection
  dev.off()
}
#
#This doesn't show position along the X axis, but still provides
#a good visual comparison between populations and across all stats
stat.ridges <- function(x, target.scaffold, reference.name){
  ggplot(scaffold.subset(x, target.scaffold = target.scaffold) %>% 
           gather(-mid, key="stat", value="value"), 
         aes(x=value, y=stat, fill=stat)) + 
    geom_density_ridges() + 
    ggtitle(paste(reference.name, target.scaffold, sep=" "))
}

############################################################################

#3. Plot
#First just some tests...
#OT 62955
plot.scaffold.diversity(
  data = OTdata, 
  target.scaffold = "62955", 
  filename = "./Plots/OTref_62955_diversity.pdf", 
  title = "OTref_62955", 
  reference.name = "OTref"
)
#PF 346
plot.scaffold.diversity(
  data = PFdata, 
  target.scaffold = "346", 
  filename = "./Plots/PFref_346_diversity.pdf", 
  title = "PFref_346", 
  reference.name = "PFref"
)
# For some reason stat.ridges() fails when wrapped inside a function to write to pdf...
pdf(file = "./Plots/OTref_62955_ridges.pdf", paper = "a4", title = "title", pointsize = 6, width = 14, height = 14)
stat.ridges(
  x=OTdata, 
  target.scaffold="62955", 
  reference.name = "OTref")
dev.off()
pdf(file = "./Plots/PFref_346_ridges.pdf", paper = "a4", title = "title", pointsize = 6, width = 14, height = 14)
stat.ridges(
  x=PFdata, 
  target.scaffold="346", 
  reference.name = "PFref")
dev.off()
#
#These appear to be in order - next all scaffolds with >=2 marker SNPs (from the
#Fst outlier analysis) in each assembly are looped over to plot
#
#Prepare vector of scaffold IDs - keep only IDs with >=2 markers
#also only keep IDs of scaffolds which are >=200kb (i.e. scaffolds with enough data to plot!)
OT_marker_scaffolds <- OT_counts_per_scaffold[which(OT_counts_per_scaffold$outlier.SNPs >= 2 & 
                                                      OT_counts_per_scaffold$LENGTH >= 200000), ]
PF_marker_scaffolds <- PF_counts_per_scaffold[which(PF_counts_per_scaffold$outlier.SNPs >= 2 & 
                                                      PF_counts_per_scaffold$LENGTH >= 200000), ]
#
#The plotting function requires a minor change to produce more meaningful titles
#(specifically only available when looping across multiple scaffolds)
#Unfortunately the simplest way of doing this is to copy and paste the existing function
#and edit the ggtitle() instruction
#(Note that the ggtitle() command given to the first plotted object is "pre_title4" - 
#this is actually defined in the subsequent loop which calls this function - a bit confusing 
#but it does work).
plot.scaffold.diversity2 <- function(data, target.scaffold, filename, title, reference.name){
  #Subset to scaffold of interest
  data=scaffold.subset(data, target.scaffold)
  #Create further subsets
  fst.data=fst.subset(data)
  dxy.data=dxy.subset(data)
  pi.data=pi.subset(data)
  #Open pdf for plotting
  pdf(file = filename, paper = "a4", title = title, pointsize = 6, width = 14, height = 14)
  #Plot...
  #Fst:
  grid.draw(rbind(
    ggplotGrob(ggplot(fst.data %>% filter(type == "between"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle(pre_title4) + 
                 scale_color_brewer(palette = "Paired") + 
                 ylab(expression(italic(F)[ST])) + 
                 ylim(c(min(fst.data$value)), max(fst.data$value)) + 
                 theme(legend.title=element_blank())), 
    ggplotGrob(ggplot(fst.data %>% filter(type == "within"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "RdYlBu") + 
                 ylab(expression(italic(F)[ST])) + 
                 ylim(c(min(fst.data$value)), max(fst.data$value)) + 
                 theme(legend.title=element_blank())), 
    #dXY:
    ggplotGrob(ggplot(dxy.data %>% filter(type == "between"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "Paired") + 
                 ylab("dXY") + 
                 ylim(c(min(dxy.data$value)), max(dxy.data$value)) + 
                 theme(legend.title=element_blank())), 
    ggplotGrob(ggplot(dxy.data %>% filter(type == "within"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "RdYlBu") + 
                 ylab("dXY") + 
                 ylim(c(min(dxy.data$value)), max(dxy.data$value)) + 
                 theme(legend.title=element_blank())), 
    #Pi:
    ggplotGrob(ggplot(pi.data %>% filter(type == "P"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "Paired") + 
                 ylab("Pi") + 
                 ylim(c(min(pi.data$value)), max(pi.data$value)) + 
                 theme(legend.title=element_blank())), 
    ggplotGrob(ggplot(pi.data %>% filter(type == "FM"), aes(mid/10^6, value, color = stat)) + 
                 geom_line() + 
                 xlab("Position (Mb)") + 
                 ggtitle("") + 
                 scale_color_brewer(palette = "RdYlBu") + 
                 ylab("Pi") + 
                 ylim(c(min(pi.data$value)), max(pi.data$value)) + 
                 theme(legend.title=element_blank())), 
    size="first")
  )
  #Close the file connection
  dev.off()
}


#
#Loop (OTref)...
for(i in 1:nrow(OT_marker_scaffolds)){
  target=OT_marker_scaffolds$CHROM[i]
  SNPcount=OT_marker_scaffolds$SNP.count[i]
  markercount=OT_marker_scaffolds$outlier.SNPs[i]
  percentoutliers=signif(OT_marker_scaffolds$perc.outliers[i], 4)
  pre_title = paste("OTref_", target, " (", sep="")
  pre_title2 = paste(pre_title, SNPcount, " SNPs,", sep = "")
  pre_title3 = paste(pre_title2, markercount, "outliers [", sep=" ")
  pre_title4 = paste(pre_title3, percentoutliers, "%])", sep="")
  plot.scaffold.diversity2(
    data = OTdata, 
    target.scaffold = target, 
    filename = paste("./Plots/OT/OTref_", target, "_diversity.pdf", sep=""), 
    title = pre_title4, 
    reference.name = "OTref"
  )
}
#Loop (PFref)...
for(i in 1:nrow(PF_marker_scaffolds)){
  target=PF_marker_scaffolds$CHROM[i]
  SNPcount=PF_marker_scaffolds$SNP.count[i]
  markercount=PF_marker_scaffolds$outlier.SNPs[i]
  percentoutliers=signif(PF_marker_scaffolds$perc.outliers[i], 4)
  pre_title = paste("PFref_", target, " (", sep="")
  pre_title2 = paste(pre_title, SNPcount, " SNPs,", sep = "")
  pre_title3 = paste(pre_title2, markercount, "outliers [", sep=" ")
  pre_title4 = paste(pre_title3, percentoutliers, "%])", sep="")
  plot.scaffold.diversity2(
    data = PFdata, 
    target.scaffold = target, 
    filename = paste("./Plots/PF/PFref_", target, "_diversity.pdf", sep=""), 
    title = pre_title4, 
    reference.name = "PFref"
  )
}
#
#Ridgeline plots....
#
#OT:
pdf(file="./Plots/OT/OTref_all_stat_ridges.pdf", paper="a4", title="OTref_all_stat_ridges", pointsize=6, width=14, height=14)
for(i in 1:nrow(OT_marker_scaffolds)){
  target=OT_marker_scaffolds$CHROM[i]
  print(stat.ridges(
    x=OTdata, 
    target.scaffold=target,
    reference.name="OTref"
  ))
}
dev.off()
#
#PF:
pdf(file="./Plots/PF/PFref_all_stat_ridges.pdf", paper="a4", title="PFref_all_stat_ridges", pointsize=6, width=14, height=14)
for(i in 1:nrow(PF_marker_scaffolds)){
  target=PF_marker_scaffolds$CHROM[i]
  print(stat.ridges(
    x=PFdata, 
    target.scaffold=target,
    reference.name="PFref"
  ))
}
dev.off()


