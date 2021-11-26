#Estimating genome size with findGSE

#Note - 
#This analysis was conducted using data from two individuals sequenced in 2017 for genome assembly.
#Sequencing was performed by Genewiz (NJ, USA) but the resulting data were insufficient
#to assemble a contiguous genome (best N50 achieved was ~50kb).
#Both individuals were males.
#One individual was from the functionally monogynous OT (Orihuela del Tremedal) population
#in spain.
#The other individual was from the polygynous SD (Santon Downham) population
#in the UK.  
#Frequency of heterozygous sites when aligning the SD reads back to the assembly created
#from them indicated the individual was likely a diploid male.

#findGSE uses the output of Jellyfish (count, then histo)
#to estimate genome size from k-mer distribution

#This work was conducted using the settings in 
#findGSE's example on their github (https://github.com/schneebergerlab/findGSE)

#Install findGSE and dependencies
install.packages("pracma")
install.packages("fGarch")
q("no")
install.packages("devtools")
devtools::install_github("schneebergerlab/findGSE")

#To prepare input data, jellyfish was run (externally):
#zcat *.fastq.gz | jellyfish count /dev/fd/0 -C -o test_21mer -m 21 -t 1 -s 5G
#jellyfish histo -h 3000000 -o test_21mer.histo test_21mer


####OT####

#1. Setwd and read histo file
setwd("/Users/u1693640/Documents/PhD/2018/findGSE/OT")
histo21 <- read.table(file = "test_21mer.histo", stringsAsFactors = FALSE)

#2. Run findGSE
library("findGSE")
findGSE(histo="test_21mer.histo", sizek=21, outdir="hom_test_21mer")


####SD####

#1. Setwd
setwd("/Users/u1693640/Documents/PhD/2018/findGSE/SD")

#2. Run findGSE
#(VCFtools indicates 0.3% of heterozygous sites)
findGSE(histo="test_21mer.histo", sizek=21, exp_hom=99.7,  outdir="het_test_21mer")






