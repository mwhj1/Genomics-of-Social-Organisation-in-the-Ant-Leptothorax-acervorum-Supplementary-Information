### All 1Mb+ and 100kb+ scafolds - rehh analysis ###

#Following tutorial at https://speciationgenomics.github.io/haplotypes/ 
#and the rehh vignette at https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html

#Package installation
#library(withr)
#with_makevars(c(PKG_CFLAGS = "-std=c99"), 
#              install.packages("tidyverse"), 
#              assignment = "+=")
#with_makevars(c(PKG_CFLAGS = "-std=c99"), 
#              install.packages("rehh"), 
#              assignment = "+=")
#install.packages("rlist")
#install.packages("tidyverse")
#with_makevars(c(PKG_CFLAGS = "-std=c99"), 
#              devtools::install_github("tidyverse/tidyverse"), 
#              assignment = "+=")
#install.packages("remotes")
#with_makevars(c(PKG_CFLAGS = "-std=c99"), 
#              remotes::install_github("r-lib/tidyselect"), 
#              assignment = "+=")
#with_makevars(c(PKG_CFLAGS = "-std=c99"), 
#              install.packages("dlookr"), 
#              assignment = "+=")
#install.packages("ggcorrplot")
#devtools::install_github("laresbernardo/lares")

################################### L analysis ###################################

#Clean up the environment
rm(list = ls())

#Load libraries
library(rehh)
library(tidyverse)

#1Mb+ scaffolds

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
OTref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"
PFref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"

#For reference-

#All OTref scaffold IDs (1Mb plus):
#62662  62734  62783  62812  62829  62846  62861  62869	62897  62919  62946  62959  62997  63061  63103
#62706  62740  62795  62814  62835  62847  62864  62877	62899  62922  62952  62962  63005  63069  
#62716  62760  62807  62816  62836  62849  62865  62890	62901  62925  62953  62978  63007  63083  
#62724  62770  62809  62818  62842  62858  62867  62894	62917  62940  62955  62988  63038  63101

#All PFref scaffold IDs (1Mb plus):
#120  139  147  232  250  294  353  363	429  77  80  83  87  90  94  
#124  144  160  246  273  296  354  365	467  78  81  85  88  91  95  
#135  146  176  248  291  346  356  370	67   79  82  86  89  92  96

#Convert to vectors of IDs
OTref_scaffold_IDs_1Mbplus <- c("62662", "62734", "62783", "62812", "62829", "62846", "62861", "62869", "62897", "62919", "62946", "62959", "62997", "63061", "63103", 
                                "62706", "62740", "62795", "62814", "62835", "62847", "62864", "62877", "62899", "62922", "62952", "62962", "63005", "63069", 
                                "62716", "62760", "62807", "62816", "62836", "62849", "62865", "62890", "62901", "62925", "62953", "62978", "63007", "63083", 
                                "62724", "62770", "62809", "62818", "62842", "62858", "62867", "62894", "62917", "62940", "62955", "62988", "63038", "63101")
PFref_scaffold_IDs_1Mbplus <- c("120", "139", "147", "232", "250", "294", "353", "363", "429", "77", "80", "83", "87", "90", "94", 
                                "124", "144", "160", "246", "273", "296", "354", "365", "467", "78", "81", "85", "88", "91", "95", 
                                "135", "146", "176", "248", "291", "346", "356", "370", "67", "79", "82", "86", "89", "92", "96")

#Function to create whole-genome dataset from individual scaffold VCFs (performs iHH calculation)
#L version:
get_wgscan <- function(master_directory, scaffold_vector){
  #Setwd to master directory
  setwd(master_directory)
  for(i in 1:length(scaffold_vector)){
    #Get target scaffold name
    target_scaffold = scaffold_vector[i]
    #Move to scaffold directory
    setwd(paste("./", target_scaffold, sep=""))
    #Print progress
    print(getwd())
    #Get VCF files
    vcf_files = dir(pattern='\\.vcf.gz$')
    #Remove first VCF file (not population-specific)
    vcf_files = vcf_files[2:6]
    #Read data
    hh = data2haplohh(hap_file = vcf_files[1], 
                      polarize_vcf = FALSE)
    #Filter on MAF = 0.05
    hh = subset(hh, min_maf = 0.05)
    #Perform scan of single scaffold
    scan = scan_hh(hh, polarized=FALSE)
    #Concatenate scaffold data to whole genome data
    if(i == 1){
      wgscan = scan
    }
    else{
      wgscan = rbind(wgscan, scan)
    }
    #Return to master directory
    setwd(master_directory)
  }
  #Return wgscan
  return(wgscan)
}

#Call the function and read the 1Mb+ scaffold data...
#OTref:
wgs.L.1Mb.OTref <- get_wgscan(OTref_dir_1Mbplus, OTref_scaffold_IDs_1Mbplus)
#PFref:
wgs.L.1Mb.PFref <- get_wgscan(PFref_dir_1Mbplus, PFref_scaffold_IDs_1Mbplus)

#100kb+ scaffolds

#Set scaffold directories
OTref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb"
PFref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb"

#OTref scaffolds:
OTref_scaffold_IDs_100kbplus <- c("39131", "61229", "61581", "61856", "61938", "61995", "62068", "62168", "62259", "62330", "62396", "62442", "62457", "62497", "62569", "62613", "62681", "62848", 
                                  "58028", "61235", "61680", "61882", "61961", "61999", "62075", "62176", "62266", "62336", "62399", "62443", "62458", "62498", "62571", "62621", "62683", "62851", 
                                  "60249", "61334", "61694", "61893", "61975", "62030", "62076", "62184", "62267", "62348", "62400", "62444", "62461", "62510", "62572", "62626", "62684", "63", 
                                  "60458", "61395", "61909", "61983", "62040", "62105", "62196", "62269", "62349", "62401", "62446", "62462", "62527", "62587", "62638", "62698", "63045", 
                                  "60771", "61436", "61704", "61911", "61984", "62042", "62107", "62208", "62270", "62358", "62416", "62449", "62465", "62529", "62590", "62646", "62732", "63049", 
                                  "60779", "61543", "61705", "61913", "61990", "62043", "62113", "62239", "62286", "62366", "62422", "62450", "62478", "62530", "62593", "62647", "62756", "63066", 
                                  "60880", "61551", "61826", "61928", "61992", "62054", "62123", "62243", "62297", "62367", "62425", "62451", "62479", "62550", "62603", "62658", "62785", 
                                  "60971", "61570", "61832", "61930", "61994", "62056", "62133", "62246", "62305", "62372", "62440", "62456", "62492", "62565", "62609", "62672", "62790")

#PFref scaffolds:
PFref_scaffold_IDs_100kbplus <- c("100563", "101478", "102614", "102907", "103389", "103756", "103902", "149", "186", "226", "238", "264", "290", "357", "405", "421", "462", "516", 
                                  "101076", "101853", "102641", "102951", "103416", "103763", "104138", "138999", "152", "189", "228", "239", "268", "302", "362", "407", "422", "465", "84", 
                                  "101196", "101952", "102650", "103187", "103625", "103772", "104162", "139916", "155", "191", "229", "242", "272", "316", "368", "408", "432", "489", "93", 
                                  "101206", "102071", "102697", "103293", "103669", "103807", "105", "139917", "172", "218", "231", "253", "274", "319", "369", "410", "444", "503", 
                                  "101260", "102114", "102735", "103306", "103672", "103828", "107", "142", "174", "219", "234", "254", "280", "341", "375", "412", "456", "504", "97", 
                                  "101397", "102330", "102801", "103355", "103716", "110", "178", "222", "235", "257", "281", "348", "392", "413", "457", "514", "97840", 
                                  "101444", "102442", "103357", "103755", "103879", "112", "143431", "179", "225", "237", "263", "287", "352", "404", "420", "461", "515")

#Call the function and read the 100kb+ scaffold data...
#OTref:
wgs.L.100kb.OTref <- get_wgscan(OTref_dir_100kbplus, OTref_scaffold_IDs_100kbplus)
#PFref:
wgs.L.100kb.PFref <- get_wgscan(PFref_dir_100kbplus, PFref_scaffold_IDs_100kbplus)

#Join 1Mb+ and 100kb+ dataframes together...
#OTref:
wgs.L.ALL.OTref <- rbind(wgs.L.1Mb.OTref, wgs.L.100kb.OTref)
#PFref:
wgs.L.ALL.PFref <- rbind(wgs.L.1Mb.PFref, wgs.L.100kb.PFref)

#Write to file for posterity...
setwd(results_dir)
#OTref:
write.table(wgs.L.ALL.OTref, 
            file="wgs.L.ALL.OTref.scan", quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(wgs.L.ALL.PFref, 
            file="wgs.L.ALL.PFref.scan", quote=FALSE, sep="\t", row.names=FALSE)

##################################################################################

################################### NR analysis ##################################

#Clean up the environment
rm(list = ls())

#Load libraries
library(rehh)
library(tidyverse)

#1Mb+ scaffolds

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
OTref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"
PFref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"

#For reference-

#All OTref scaffold IDs (1Mb plus):
#62662  62734  62783  62812  62829  62846  62861  62869	62897  62919  62946  62959  62997  63061  63103
#62706  62740  62795  62814  62835  62847  62864  62877	62899  62922  62952  62962  63005  63069  
#62716  62760  62807  62816  62836  62849  62865  62890	62901  62925  62953  62978  63007  63083  
#62724  62770  62809  62818  62842  62858  62867  62894	62917  62940  62955  62988  63038  63101

#All PFref scaffold IDs (1Mb plus):
#120  139  147  232  250  294  353  363	429  77  80  83  87  90  94  
#124  144  160  246  273  296  354  365	467  78  81  85  88  91  95  
#135  146  176  248  291  346  356  370	67   79  82  86  89  92  96

#Convert to vectors of IDs
OTref_scaffold_IDs_1Mbplus <- c("62662", "62734", "62783", "62812", "62829", "62846", "62861", "62869", "62897", "62919", "62946", "62959", "62997", "63061", "63103", 
                                "62706", "62740", "62795", "62814", "62835", "62847", "62864", "62877", "62899", "62922", "62952", "62962", "63005", "63069", 
                                "62716", "62760", "62807", "62816", "62836", "62849", "62865", "62890", "62901", "62925", "62953", "62978", "63007", "63083", 
                                "62724", "62770", "62809", "62818", "62842", "62858", "62867", "62894", "62917", "62940", "62955", "62988", "63038", "63101")
PFref_scaffold_IDs_1Mbplus <- c("120", "139", "147", "232", "250", "294", "353", "363", "429", "77", "80", "83", "87", "90", "94", 
                                "124", "144", "160", "246", "273", "296", "354", "365", "467", "78", "81", "85", "88", "91", "95", 
                                "135", "146", "176", "248", "291", "346", "356", "370", "67", "79", "82", "86", "89", "92", "96")

#Function to create whole-genome dataset from individual scaffold VCFs (performs iHH calculation)
#NR version:
get_wgscan <- function(master_directory, scaffold_vector){
  #Setwd to master directory
  setwd(master_directory)
  for(i in 1:length(scaffold_vector)){
    #Get target scaffold name
    target_scaffold = scaffold_vector[i]
    #Move to scaffold directory
    setwd(paste("./", target_scaffold, sep=""))
    #Print progress
    print(getwd())
    #Get VCF files
    vcf_files = dir(pattern='\\.vcf.gz$')
    #Remove first VCF file (not population-specific)
    vcf_files = vcf_files[2:6]
    #Read data
    hh = data2haplohh(hap_file = vcf_files[2], 
                      polarize_vcf = FALSE)
    #Filter on MAF = 0.05
    hh = subset(hh, min_maf = 0.05)
    #Perform scan of single scaffold
    scan = scan_hh(hh, polarized=FALSE)
    #Concatenate scaffold data to whole genome data
    if(i == 1){
      wgscan = scan
    }
    else{
      wgscan = rbind(wgscan, scan)
    }
    #Return to master directory
    setwd(master_directory)
  }
  #Return wgscan
  return(wgscan)
}

#Call the function and read the 1Mb+ scaffold data...
#OTref:
wgs.NR.1Mb.OTref <- get_wgscan(OTref_dir_1Mbplus, OTref_scaffold_IDs_1Mbplus)
#PFref:
wgs.NR.1Mb.PFref <- get_wgscan(PFref_dir_1Mbplus, PFref_scaffold_IDs_1Mbplus)

#100kb+ scaffolds

#Set scaffold directories
OTref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb"
PFref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb"

#OTref scaffolds:
OTref_scaffold_IDs_100kbplus <- c("39131", "61229", "61581", "61856", "61938", "61995", "62068", "62168", "62259", "62330", "62396", "62442", "62457", "62497", "62569", "62613", "62681", "62848", 
                                  "58028", "61235", "61680", "61882", "61961", "61999", "62075", "62176", "62266", "62336", "62399", "62443", "62458", "62498", "62571", "62621", "62683", "62851", 
                                  "60249", "61334", "61694", "61893", "61975", "62030", "62076", "62184", "62267", "62348", "62400", "62444", "62461", "62510", "62572", "62626", "62684", "63", 
                                  "60458", "61395", "61909", "61983", "62040", "62105", "62196", "62269", "62349", "62401", "62446", "62462", "62527", "62587", "62638", "62698", "63045", 
                                  "60771", "61436", "61704", "61911", "61984", "62042", "62107", "62208", "62270", "62358", "62416", "62449", "62465", "62529", "62590", "62646", "62732", "63049", 
                                  "60779", "61543", "61705", "61913", "61990", "62043", "62113", "62239", "62286", "62366", "62422", "62450", "62478", "62530", "62593", "62647", "62756", "63066", 
                                  "60880", "61551", "61826", "61928", "61992", "62054", "62123", "62243", "62297", "62367", "62425", "62451", "62479", "62550", "62603", "62658", "62785", 
                                  "60971", "61570", "61832", "61930", "61994", "62056", "62133", "62246", "62305", "62372", "62440", "62456", "62492", "62565", "62609", "62672", "62790")

#PFref scaffolds:
PFref_scaffold_IDs_100kbplus <- c("100563", "101478", "102614", "102907", "103389", "103756", "103902", "149", "186", "226", "238", "264", "290", "357", "405", "421", "462", "516", 
                                  "101076", "101853", "102641", "102951", "103416", "103763", "104138", "138999", "152", "189", "228", "239", "268", "302", "362", "407", "422", "465", "84", 
                                  "101196", "101952", "102650", "103187", "103625", "103772", "104162", "139916", "155", "191", "229", "242", "272", "316", "368", "408", "432", "489", "93", 
                                  "101206", "102071", "102697", "103293", "103669", "103807", "105", "139917", "172", "218", "231", "253", "274", "319", "369", "410", "444", "503", 
                                  "101260", "102114", "102735", "103306", "103672", "103828", "107", "142", "174", "219", "234", "254", "280", "341", "375", "412", "456", "504", "97", 
                                  "101397", "102330", "102801", "103355", "103716", "110", "178", "222", "235", "257", "281", "348", "392", "413", "457", "514", "97840", 
                                  "101444", "102442", "103357", "103755", "103879", "112", "143431", "179", "225", "237", "263", "287", "352", "404", "420", "461", "515")

#Call the function and read the 100kb+ scaffold data...
#OTref:
wgs.NR.100kb.OTref <- get_wgscan(OTref_dir_100kbplus, OTref_scaffold_IDs_100kbplus)
#PFref:
wgs.NR.100kb.PFref <- get_wgscan(PFref_dir_100kbplus, PFref_scaffold_IDs_100kbplus)

#Join 1Mb+ and 100kb+ dataframes together...
#OTref:
wgs.NR.ALL.OTref <- rbind(wgs.NR.1Mb.OTref, wgs.NR.100kb.OTref)
#PFref:
wgs.NR.ALL.PFref <- rbind(wgs.NR.1Mb.PFref, wgs.NR.100kb.PFref)

#Write to file for posterity...
setwd(results_dir)
#OTref:
write.table(wgs.NR.ALL.OTref, 
            file="wgs.NR.ALL.OTref.scan", quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(wgs.NR.ALL.PFref, 
            file="wgs.NR.ALL.PFref.scan", quote=FALSE, sep="\t", row.names=FALSE)

##################################################################################

################################### OT analysis ##################################

#Clean up the environment
rm(list = ls())

#Load libraries
library(rehh)
library(tidyverse)

#1Mb+ scaffolds

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
OTref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"
PFref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"

#For reference-

#All OTref scaffold IDs (1Mb plus):
#62662  62734  62783  62812  62829  62846  62861  62869	62897  62919  62946  62959  62997  63061  63103
#62706  62740  62795  62814  62835  62847  62864  62877	62899  62922  62952  62962  63005  63069  
#62716  62760  62807  62816  62836  62849  62865  62890	62901  62925  62953  62978  63007  63083  
#62724  62770  62809  62818  62842  62858  62867  62894	62917  62940  62955  62988  63038  63101

#All PFref scaffold IDs (1Mb plus):
#120  139  147  232  250  294  353  363	429  77  80  83  87  90  94  
#124  144  160  246  273  296  354  365	467  78  81  85  88  91  95  
#135  146  176  248  291  346  356  370	67   79  82  86  89  92  96

#Convert to vectors of IDs
OTref_scaffold_IDs_1Mbplus <- c("62662", "62734", "62783", "62812", "62829", "62846", "62861", "62869", "62897", "62919", "62946", "62959", "62997", "63061", "63103", 
                                "62706", "62740", "62795", "62814", "62835", "62847", "62864", "62877", "62899", "62922", "62952", "62962", "63005", "63069", 
                                "62716", "62760", "62807", "62816", "62836", "62849", "62865", "62890", "62901", "62925", "62953", "62978", "63007", "63083", 
                                "62724", "62770", "62809", "62818", "62842", "62858", "62867", "62894", "62917", "62940", "62955", "62988", "63038", "63101")
PFref_scaffold_IDs_1Mbplus <- c("120", "139", "147", "232", "250", "294", "353", "363", "429", "77", "80", "83", "87", "90", "94", 
                                "124", "144", "160", "246", "273", "296", "354", "365", "467", "78", "81", "85", "88", "91", "95", 
                                "135", "146", "176", "248", "291", "346", "356", "370", "67", "79", "82", "86", "89", "92", "96")

#Function to create whole-genome dataset from individual scaffold VCFs (performs iHH calculation)
#OT version:
get_wgscan <- function(master_directory, scaffold_vector){
  #Setwd to master directory
  setwd(master_directory)
  for(i in 1:length(scaffold_vector)){
    #Get target scaffold name
    target_scaffold = scaffold_vector[i]
    #Move to scaffold directory
    setwd(paste("./", target_scaffold, sep=""))
    #Print progress
    print(getwd())
    #Get VCF files
    vcf_files = dir(pattern='\\.vcf.gz$')
    #Remove first VCF file (not population-specific)
    vcf_files = vcf_files[2:6]
    #Read data
    hh = data2haplohh(hap_file = vcf_files[3], 
                      polarize_vcf = FALSE)
    #Filter on MAF = 0.05
    hh = subset(hh, min_maf = 0.05)
    #Perform scan of single scaffold
    scan = scan_hh(hh, polarized=FALSE)
    #Concatenate scaffold data to whole genome data
    if(i == 1){
      wgscan = scan
    }
    else{
      wgscan = rbind(wgscan, scan)
    }
    #Return to master directory
    setwd(master_directory)
  }
  #Return wgscan
  return(wgscan)
}

#Call the function and read the 1Mb+ scaffold data...
#OTref:
wgs.OT.1Mb.OTref <- get_wgscan(OTref_dir_1Mbplus, OTref_scaffold_IDs_1Mbplus)
#PFref:
wgs.OT.1Mb.PFref <- get_wgscan(PFref_dir_1Mbplus, PFref_scaffold_IDs_1Mbplus)

#100kb+ scaffolds

#Set scaffold directories
OTref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb"
PFref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb"

#OTref scaffolds:
OTref_scaffold_IDs_100kbplus <- c("39131", "61229", "61581", "61856", "61938", "61995", "62068", "62168", "62259", "62330", "62396", "62442", "62457", "62497", "62569", "62613", "62681", "62848", 
                                  "58028", "61235", "61680", "61882", "61961", "61999", "62075", "62176", "62266", "62336", "62399", "62443", "62458", "62498", "62571", "62621", "62683", "62851", 
                                  "60249", "61334", "61694", "61893", "61975", "62030", "62076", "62184", "62267", "62348", "62400", "62444", "62461", "62510", "62572", "62626", "62684", "63", 
                                  "60458", "61395", "61909", "61983", "62040", "62105", "62196", "62269", "62349", "62401", "62446", "62462", "62527", "62587", "62638", "62698", "63045", 
                                  "60771", "61436", "61704", "61911", "61984", "62042", "62107", "62208", "62270", "62358", "62416", "62449", "62465", "62529", "62590", "62646", "62732", "63049", 
                                  "60779", "61543", "61705", "61913", "61990", "62043", "62113", "62239", "62286", "62366", "62422", "62450", "62478", "62530", "62593", "62647", "62756", "63066", 
                                  "60880", "61551", "61826", "61928", "61992", "62054", "62123", "62243", "62297", "62367", "62425", "62451", "62479", "62550", "62603", "62658", "62785", 
                                  "60971", "61570", "61832", "61930", "61994", "62056", "62133", "62246", "62305", "62372", "62440", "62456", "62492", "62565", "62609", "62672", "62790")

#PFref scaffolds:
PFref_scaffold_IDs_100kbplus <- c("100563", "101478", "102614", "102907", "103389", "103756", "103902", "149", "186", "226", "238", "264", "290", "357", "405", "421", "462", "516", 
                                  "101076", "101853", "102641", "102951", "103416", "103763", "104138", "138999", "152", "189", "228", "239", "268", "302", "362", "407", "422", "465", "84", 
                                  "101196", "101952", "102650", "103187", "103625", "103772", "104162", "139916", "155", "191", "229", "242", "272", "316", "368", "408", "432", "489", "93", 
                                  "101206", "102071", "102697", "103293", "103669", "103807", "105", "139917", "172", "218", "231", "253", "274", "319", "369", "410", "444", "503", 
                                  "101260", "102114", "102735", "103306", "103672", "103828", "107", "142", "174", "219", "234", "254", "280", "341", "375", "412", "456", "504", "97", 
                                  "101397", "102330", "102801", "103355", "103716", "110", "178", "222", "235", "257", "281", "348", "392", "413", "457", "514", "97840", 
                                  "101444", "102442", "103357", "103755", "103879", "112", "143431", "179", "225", "237", "263", "287", "352", "404", "420", "461", "515")

#Call the function and read the 100kb+ scaffold data...
#OTref:
wgs.OT.100kb.OTref <- get_wgscan(OTref_dir_100kbplus, OTref_scaffold_IDs_100kbplus)
#PFref:
wgs.OT.100kb.PFref <- get_wgscan(PFref_dir_100kbplus, PFref_scaffold_IDs_100kbplus)

#Join 1Mb+ and 100kb+ dataframes together...
#OTref:
wgs.OT.ALL.OTref <- rbind(wgs.OT.1Mb.OTref, wgs.OT.100kb.OTref)
#PFref:
wgs.OT.ALL.PFref <- rbind(wgs.OT.1Mb.PFref, wgs.OT.100kb.PFref)

#Write to file for posterity...
setwd(results_dir)
#OTref:
write.table(wgs.OT.ALL.OTref, 
            file="wgs.OT.ALL.OTref.scan", quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(wgs.OT.ALL.PFref, 
            file="wgs.OT.ALL.PFref.scan", quote=FALSE, sep="\t", row.names=FALSE)

##################################################################################

################################### PF analysis ##################################

#Clean up the environment
rm(list = ls())

#Load libraries
library(rehh)
library(tidyverse)

#1Mb+ scaffolds

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
OTref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"
PFref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"

#For reference-

#All OTref scaffold IDs (1Mb plus):
#62662  62734  62783  62812  62829  62846  62861  62869	62897  62919  62946  62959  62997  63061  63103
#62706  62740  62795  62814  62835  62847  62864  62877	62899  62922  62952  62962  63005  63069  
#62716  62760  62807  62816  62836  62849  62865  62890	62901  62925  62953  62978  63007  63083  
#62724  62770  62809  62818  62842  62858  62867  62894	62917  62940  62955  62988  63038  63101

#All PFref scaffold IDs (1Mb plus):
#120  139  147  232  250  294  353  363	429  77  80  83  87  90  94  
#124  144  160  246  273  296  354  365	467  78  81  85  88  91  95  
#135  146  176  248  291  346  356  370	67   79  82  86  89  92  96

#Convert to vectors of IDs
OTref_scaffold_IDs_1Mbplus <- c("62662", "62734", "62783", "62812", "62829", "62846", "62861", "62869", "62897", "62919", "62946", "62959", "62997", "63061", "63103", 
                                "62706", "62740", "62795", "62814", "62835", "62847", "62864", "62877", "62899", "62922", "62952", "62962", "63005", "63069", 
                                "62716", "62760", "62807", "62816", "62836", "62849", "62865", "62890", "62901", "62925", "62953", "62978", "63007", "63083", 
                                "62724", "62770", "62809", "62818", "62842", "62858", "62867", "62894", "62917", "62940", "62955", "62988", "63038", "63101")
PFref_scaffold_IDs_1Mbplus <- c("120", "139", "147", "232", "250", "294", "353", "363", "429", "77", "80", "83", "87", "90", "94", 
                                "124", "144", "160", "246", "273", "296", "354", "365", "467", "78", "81", "85", "88", "91", "95", 
                                "135", "146", "176", "248", "291", "346", "356", "370", "67", "79", "82", "86", "89", "92", "96")

#Function to create whole-genome dataset from individual scaffold VCFs (performs iHH calculation)
#PF version:
get_wgscan <- function(master_directory, scaffold_vector){
  #Setwd to master directory
  setwd(master_directory)
  for(i in 1:length(scaffold_vector)){
    #Get target scaffold name
    target_scaffold = scaffold_vector[i]
    #Move to scaffold directory
    setwd(paste("./", target_scaffold, sep=""))
    #Print progress
    print(getwd())
    #Get VCF files
    vcf_files = dir(pattern='\\.vcf.gz$')
    #Remove first VCF file (not population-specific)
    vcf_files = vcf_files[2:6]
    #Read data
    hh = data2haplohh(hap_file = vcf_files[4], 
                      polarize_vcf = FALSE)
    #Filter on MAF = 0.05
    hh = subset(hh, min_maf = 0.05)
    #Perform scan of single scaffold
    scan = scan_hh(hh, polarized=FALSE)
    #Concatenate scaffold data to whole genome data
    if(i == 1){
      wgscan = scan
    }
    else{
      wgscan = rbind(wgscan, scan)
    }
    #Return to master directory
    setwd(master_directory)
  }
  #Return wgscan
  return(wgscan)
}

#Call the function and read the 1Mb+ scaffold data...
#OTref:
wgs.PF.1Mb.OTref <- get_wgscan(OTref_dir_1Mbplus, OTref_scaffold_IDs_1Mbplus)
#PFref:
wgs.PF.1Mb.PFref <- get_wgscan(PFref_dir_1Mbplus, PFref_scaffold_IDs_1Mbplus)

#100kb+ scaffolds

#Set scaffold directories
OTref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb"
PFref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb"

#OTref scaffolds:
OTref_scaffold_IDs_100kbplus <- c("39131", "61229", "61581", "61856", "61938", "61995", "62068", "62168", "62259", "62330", "62396", "62442", "62457", "62497", "62569", "62613", "62681", "62848", 
                                  "58028", "61235", "61680", "61882", "61961", "61999", "62075", "62176", "62266", "62336", "62399", "62443", "62458", "62498", "62571", "62621", "62683", "62851", 
                                  "60249", "61334", "61694", "61893", "61975", "62030", "62076", "62184", "62267", "62348", "62400", "62444", "62461", "62510", "62572", "62626", "62684", "63", 
                                  "60458", "61395", "61909", "61983", "62040", "62105", "62196", "62269", "62349", "62401", "62446", "62462", "62527", "62587", "62638", "62698", "63045", 
                                  "60771", "61436", "61704", "61911", "61984", "62042", "62107", "62208", "62270", "62358", "62416", "62449", "62465", "62529", "62590", "62646", "62732", "63049", 
                                  "60779", "61543", "61705", "61913", "61990", "62043", "62113", "62239", "62286", "62366", "62422", "62450", "62478", "62530", "62593", "62647", "62756", "63066", 
                                  "60880", "61551", "61826", "61928", "61992", "62054", "62123", "62243", "62297", "62367", "62425", "62451", "62479", "62550", "62603", "62658", "62785", 
                                  "60971", "61570", "61832", "61930", "61994", "62056", "62133", "62246", "62305", "62372", "62440", "62456", "62492", "62565", "62609", "62672", "62790")

#PFref scaffolds:
PFref_scaffold_IDs_100kbplus <- c("100563", "101478", "102614", "102907", "103389", "103756", "103902", "149", "186", "226", "238", "264", "290", "357", "405", "421", "462", "516", 
                                  "101076", "101853", "102641", "102951", "103416", "103763", "104138", "138999", "152", "189", "228", "239", "268", "302", "362", "407", "422", "465", "84", 
                                  "101196", "101952", "102650", "103187", "103625", "103772", "104162", "139916", "155", "191", "229", "242", "272", "316", "368", "408", "432", "489", "93", 
                                  "101206", "102071", "102697", "103293", "103669", "103807", "105", "139917", "172", "218", "231", "253", "274", "319", "369", "410", "444", "503", 
                                  "101260", "102114", "102735", "103306", "103672", "103828", "107", "142", "174", "219", "234", "254", "280", "341", "375", "412", "456", "504", "97", 
                                  "101397", "102330", "102801", "103355", "103716", "110", "178", "222", "235", "257", "281", "348", "392", "413", "457", "514", "97840", 
                                  "101444", "102442", "103357", "103755", "103879", "112", "143431", "179", "225", "237", "263", "287", "352", "404", "420", "461", "515")

#Call the function and read the 100kb+ scaffold data...
#OTref:
wgs.PF.100kb.OTref <- get_wgscan(OTref_dir_100kbplus, OTref_scaffold_IDs_100kbplus)
#PFref:
wgs.PF.100kb.PFref <- get_wgscan(PFref_dir_100kbplus, PFref_scaffold_IDs_100kbplus)

#Join 1Mb+ and 100kb+ dataframes together...
#OTref:
wgs.PF.ALL.OTref <- rbind(wgs.PF.1Mb.OTref, wgs.PF.100kb.OTref)
#PFref:
wgs.PF.ALL.PFref <- rbind(wgs.PF.1Mb.PFref, wgs.PF.100kb.PFref)

#Write to file for posterity...
setwd(results_dir)
#OTref:
write.table(wgs.PF.ALL.OTref, 
            file="wgs.PF.ALL.OTref.scan", quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(wgs.PF.ALL.PFref, 
            file="wgs.PF.ALL.PFref.scan", quote=FALSE, sep="\t", row.names=FALSE)

##################################################################################

################################### V analysis ##################################

#Clean up the environment
rm(list = ls())

#Load libraries
library(rehh)
library(tidyverse)

#1Mb+ scaffolds

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
OTref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"
PFref_dir_1Mbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb/Scaffolds_over_1Mb"

#For reference-

#All OTref scaffold IDs (1Mb plus):
#62662  62734  62783  62812  62829  62846  62861  62869	62897  62919  62946  62959  62997  63061  63103
#62706  62740  62795  62814  62835  62847  62864  62877	62899  62922  62952  62962  63005  63069  
#62716  62760  62807  62816  62836  62849  62865  62890	62901  62925  62953  62978  63007  63083  
#62724  62770  62809  62818  62842  62858  62867  62894	62917  62940  62955  62988  63038  63101

#All PFref scaffold IDs (1Mb plus):
#120  139  147  232  250  294  353  363	429  77  80  83  87  90  94  
#124  144  160  246  273  296  354  365	467  78  81  85  88  91  95  
#135  146  176  248  291  346  356  370	67   79  82  86  89  92  96

#Convert to vectors of IDs
OTref_scaffold_IDs_1Mbplus <- c("62662", "62734", "62783", "62812", "62829", "62846", "62861", "62869", "62897", "62919", "62946", "62959", "62997", "63061", "63103", 
                                "62706", "62740", "62795", "62814", "62835", "62847", "62864", "62877", "62899", "62922", "62952", "62962", "63005", "63069", 
                                "62716", "62760", "62807", "62816", "62836", "62849", "62865", "62890", "62901", "62925", "62953", "62978", "63007", "63083", 
                                "62724", "62770", "62809", "62818", "62842", "62858", "62867", "62894", "62917", "62940", "62955", "62988", "63038", "63101")
PFref_scaffold_IDs_1Mbplus <- c("120", "139", "147", "232", "250", "294", "353", "363", "429", "77", "80", "83", "87", "90", "94", 
                                "124", "144", "160", "246", "273", "296", "354", "365", "467", "78", "81", "85", "88", "91", "95", 
                                "135", "146", "176", "248", "291", "346", "356", "370", "67", "79", "82", "86", "89", "92", "96")

#Function to create whole-genome dataset from individual scaffold VCFs (performs iHH calculation)
#V version:
get_wgscan <- function(master_directory, scaffold_vector){
  #Setwd to master directory
  setwd(master_directory)
  for(i in 1:length(scaffold_vector)){
    #Get target scaffold name
    target_scaffold = scaffold_vector[i]
    #Move to scaffold directory
    setwd(paste("./", target_scaffold, sep=""))
    #Print progress
    print(getwd())
    #Get VCF files
    vcf_files = dir(pattern='\\.vcf.gz$')
    #Remove first VCF file (not population-specific)
    vcf_files = vcf_files[2:6]
    #Read data
    hh = data2haplohh(hap_file = vcf_files[5], 
                      polarize_vcf = FALSE)
    #Filter on MAF = 0.05
    hh = subset(hh, min_maf = 0.05)
    #Perform scan of single scaffold
    scan = scan_hh(hh, polarized=FALSE)
    #Concatenate scaffold data to whole genome data
    if(i == 1){
      wgscan = scan
    }
    else{
      wgscan = rbind(wgscan, scan)
    }
    #Return to master directory
    setwd(master_directory)
  }
  #Return wgscan
  return(wgscan)
}

#Call the function and read the 1Mb+ scaffold data...
#OTref:
wgs.V.1Mb.OTref <- get_wgscan(OTref_dir_1Mbplus, OTref_scaffold_IDs_1Mbplus)
#PFref:
wgs.V.1Mb.PFref <- get_wgscan(PFref_dir_1Mbplus, PFref_scaffold_IDs_1Mbplus)

#100kb+ scaffolds

#Set scaffold directories
OTref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/OTref/Scaffolds_over_10kb/Scaffolds_over_100kb"
PFref_dir_100kbplus <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/PFref/Scaffolds_over_10kb/Scaffolds_over_100kb"

#OTref scaffolds:
OTref_scaffold_IDs_100kbplus <- c("39131", "61229", "61581", "61856", "61938", "61995", "62068", "62168", "62259", "62330", "62396", "62442", "62457", "62497", "62569", "62613", "62681", "62848", 
                                  "58028", "61235", "61680", "61882", "61961", "61999", "62075", "62176", "62266", "62336", "62399", "62443", "62458", "62498", "62571", "62621", "62683", "62851", 
                                  "60249", "61334", "61694", "61893", "61975", "62030", "62076", "62184", "62267", "62348", "62400", "62444", "62461", "62510", "62572", "62626", "62684", "63", 
                                  "60458", "61395", "61909", "61983", "62040", "62105", "62196", "62269", "62349", "62401", "62446", "62462", "62527", "62587", "62638", "62698", "63045", 
                                  "60771", "61436", "61704", "61911", "61984", "62042", "62107", "62208", "62270", "62358", "62416", "62449", "62465", "62529", "62590", "62646", "62732", "63049", 
                                  "60779", "61543", "61705", "61913", "61990", "62043", "62113", "62239", "62286", "62366", "62422", "62450", "62478", "62530", "62593", "62647", "62756", "63066", 
                                  "60880", "61551", "61826", "61928", "61992", "62054", "62123", "62243", "62297", "62367", "62425", "62451", "62479", "62550", "62603", "62658", "62785", 
                                  "60971", "61570", "61832", "61930", "61994", "62056", "62133", "62246", "62305", "62372", "62440", "62456", "62492", "62565", "62609", "62672", "62790")

#PFref scaffolds:
PFref_scaffold_IDs_100kbplus <- c("100563", "101478", "102614", "102907", "103389", "103756", "103902", "149", "186", "226", "238", "264", "290", "357", "405", "421", "462", "516", 
                                  "101076", "101853", "102641", "102951", "103416", "103763", "104138", "138999", "152", "189", "228", "239", "268", "302", "362", "407", "422", "465", "84", 
                                  "101196", "101952", "102650", "103187", "103625", "103772", "104162", "139916", "155", "191", "229", "242", "272", "316", "368", "408", "432", "489", "93", 
                                  "101206", "102071", "102697", "103293", "103669", "103807", "105", "139917", "172", "218", "231", "253", "274", "319", "369", "410", "444", "503", 
                                  "101260", "102114", "102735", "103306", "103672", "103828", "107", "142", "174", "219", "234", "254", "280", "341", "375", "412", "456", "504", "97", 
                                  "101397", "102330", "102801", "103355", "103716", "110", "178", "222", "235", "257", "281", "348", "392", "413", "457", "514", "97840", 
                                  "101444", "102442", "103357", "103755", "103879", "112", "143431", "179", "225", "237", "263", "287", "352", "404", "420", "461", "515")

#Call the function and read the 100kb+ scaffold data...
#OTref:
wgs.V.100kb.OTref <- get_wgscan(OTref_dir_100kbplus, OTref_scaffold_IDs_100kbplus)
#PFref:
wgs.V.100kb.PFref <- get_wgscan(PFref_dir_100kbplus, PFref_scaffold_IDs_100kbplus)

#Join 1Mb+ and 100kb+ dataframes together...
#OTref:
wgs.V.ALL.OTref <- rbind(wgs.V.1Mb.OTref, wgs.V.100kb.OTref)
#PFref:
wgs.V.ALL.PFref <- rbind(wgs.V.1Mb.PFref, wgs.V.100kb.PFref)

#Write to file for posterity...
setwd(results_dir)
#OTref:
write.table(wgs.V.ALL.OTref, 
            file="wgs.V.ALL.OTref.scan", quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(wgs.V.ALL.PFref, 
            file="wgs.V.ALL.PFref.scan", quote=FALSE, sep="\t", row.names=FALSE)

##################################################################################

############################### Full wgs analysis ################################

#Clean up
rm(list = ls())

#Load libraries
library(rehh)
library(tidyverse)
library(dplyr)
library(rlist)
library(dlookr)
library(ggcorrplot)
library(ggplot2)
library(lares)

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
setwd(results_dir)

#Read back in wgs scans...
#OTref:
wgs.L.ALL.OTref <- data.table::fread(file="wgs.L.ALL.OTref.scan")
wgs.NR.ALL.OTref <- data.table::fread(file="wgs.NR.ALL.OTref.scan")
wgs.OT.ALL.OTref <- data.table::fread(file="wgs.OT.ALL.OTref.scan")
wgs.PF.ALL.OTref <- data.table::fread(file="wgs.PF.ALL.OTref.scan")
wgs.V.ALL.OTref <- data.table::fread(file="wgs.V.ALL.OTref.scan")
#PFref:
wgs.L.ALL.PFref <- data.table::fread(file="wgs.L.ALL.PFref.scan")
wgs.NR.ALL.PFref <- data.table::fread(file="wgs.NR.ALL.PFref.scan")
wgs.OT.ALL.PFref <- data.table::fread(file="wgs.OT.ALL.PFref.scan")
wgs.PF.ALL.PFref <- data.table::fread(file="wgs.PF.ALL.PFref.scan")
wgs.V.ALL.PFref <- data.table::fread(file="wgs.V.ALL.PFref.scan")

### Calculation of selection statistics ###

#Perform iHS calculations with FDR p-value correction...
#OTref:
L.ALL.OTref.iHS <- ihh2ihs(wgs.L.ALL.OTref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)
NR.ALL.OTref.iHS <- ihh2ihs(wgs.NR.ALL.OTref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)
OT.ALL.OTref.iHS <- ihh2ihs(wgs.OT.ALL.OTref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)
PF.ALL.OTref.iHS <- ihh2ihs(wgs.PF.ALL.OTref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)
V.ALL.OTref.iHS <- ihh2ihs(wgs.V.ALL.OTref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)
#PFref:
L.ALL.PFref.iHS <- ihh2ihs(wgs.L.ALL.PFref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)
NR.ALL.PFref.iHS <- ihh2ihs(wgs.NR.ALL.PFref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)
OT.ALL.PFref.iHS <- ihh2ihs(wgs.OT.ALL.PFref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)
PF.ALL.PFref.iHS <- ihh2ihs(wgs.PF.ALL.PFref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)
V.ALL.PFref.iHS <- ihh2ihs(wgs.V.ALL.PFref, freqbin = 1, p.adjust.method="BY", include_freq=TRUE)

#Write to file for posterity (main iHS data)...
#OTref:
write.table(as.data.frame(L.ALL.OTref.iHS[1]), 
            file="L.ALL.OTref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(NR.ALL.OTref.iHS[1]), 
            file="NR.ALL.OTref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(OT.ALL.OTref.iHS[1]), 
            file="OT.ALL.OTref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(PF.ALL.OTref.iHS[1]), 
            file="PF.ALL.OTref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(V.ALL.OTref.iHS[1]), 
            file="V.ALL.OTref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(as.data.frame(L.ALL.PFref.iHS[1]), 
            file="L.ALL.PFref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(NR.ALL.PFref.iHS[1]), 
            file="NR.ALL.PFref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(OT.ALL.PFref.iHS[1]), 
            file="OT.ALL.PFref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(PF.ALL.PFref.iHS[1]), 
            file="PF.ALL.PFref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(V.ALL.PFref.iHS[1]), 
            file="V.ALL.PFref.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Write to file for posterity (summary iHS data)...
#OTref:
write.table(as.data.frame(L.ALL.OTref.iHS[2]), 
            file="L.ALL.OTref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(NR.ALL.OTref.iHS[2]), 
            file="NR.ALL.OTref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(OT.ALL.OTref.iHS[2]), 
            file="OT.ALL.OTref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(PF.ALL.OTref.iHS[2]), 
            file="PF.ALL.OTref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(V.ALL.OTref.iHS[2]), 
            file="V.ALL.OTref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(as.data.frame(L.ALL.PFref.iHS[2]), 
            file="L.ALL.PFref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(NR.ALL.PFref.iHS[2]), 
            file="NR.ALL.PFref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(OT.ALL.PFref.iHS[2]), 
            file="OT.ALL.PFref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(PF.ALL.PFref.iHS[2]), 
            file="PF.ALL.PFref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(V.ALL.PFref.iHS[2]), 
            file="V.ALL.PFref.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Perform iHS calculations **without** FDR p-value correction...
#OTref:
L.ALL.OTref.iHS.noFDR <- ihh2ihs(wgs.L.ALL.OTref, freqbin = 1, include_freq=TRUE)
NR.ALL.OTref.iHS.noFDR <- ihh2ihs(wgs.NR.ALL.OTref, freqbin = 1, include_freq=TRUE)
OT.ALL.OTref.iHS.noFDR <- ihh2ihs(wgs.OT.ALL.OTref, freqbin = 1, include_freq=TRUE)
PF.ALL.OTref.iHS.noFDR <- ihh2ihs(wgs.PF.ALL.OTref, freqbin = 1, include_freq=TRUE)
V.ALL.OTref.iHS.noFDR <- ihh2ihs(wgs.V.ALL.OTref, freqbin = 1, include_freq=TRUE)
#PFref:
L.ALL.PFref.iHS.noFDR <- ihh2ihs(wgs.L.ALL.PFref, freqbin = 1, include_freq=TRUE)
NR.ALL.PFref.iHS.noFDR <- ihh2ihs(wgs.NR.ALL.PFref, freqbin = 1, include_freq=TRUE)
OT.ALL.PFref.iHS.noFDR <- ihh2ihs(wgs.OT.ALL.PFref, freqbin = 1, include_freq=TRUE)
PF.ALL.PFref.iHS.noFDR <- ihh2ihs(wgs.PF.ALL.PFref, freqbin = 1, include_freq=TRUE)
V.ALL.PFref.iHS.noFDR <- ihh2ihs(wgs.V.ALL.PFref, freqbin = 1, include_freq=TRUE)

#Write to file for posterity (main iHS data)...
#OTref:
write.table(as.data.frame(L.ALL.OTref.iHS.noFDR[1]), 
            file="L.ALL.OTref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(NR.ALL.OTref.iHS.noFDR[1]), 
            file="NR.ALL.OTref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(OT.ALL.OTref.iHS.noFDR[1]), 
            file="OT.ALL.OTref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(PF.ALL.OTref.iHS.noFDR[1]), 
            file="PF.ALL.OTref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(V.ALL.OTref.iHS.noFDR[1]), 
            file="V.ALL.OTref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(as.data.frame(L.ALL.PFref.iHS.noFDR[1]), 
            file="L.ALL.PFref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(NR.ALL.PFref.iHS.noFDR[1]), 
            file="NR.ALL.PFref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(OT.ALL.PFref.iHS.noFDR[1]), 
            file="OT.ALL.PFref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(PF.ALL.PFref.iHS.noFDR[1]), 
            file="PF.ALL.PFref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(V.ALL.PFref.iHS.noFDR[1]), 
            file="V.ALL.PFref.noPvalueCorrection.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Write to file for posterity (summary iHS data)...
#OTref:
write.table(as.data.frame(L.ALL.OTref.iHS.noFDR[2]), 
            file="L.ALL.OTref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(NR.ALL.OTref.iHS.noFDR[2]), 
            file="NR.ALL.OTref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(OT.ALL.OTref.iHS.noFDR[2]), 
            file="OT.ALL.OTref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(PF.ALL.OTref.iHS.noFDR[2]), 
            file="PF.ALL.OTref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(V.ALL.OTref.iHS.noFDR[2]), 
            file="V.ALL.OTref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(as.data.frame(L.ALL.PFref.iHS.noFDR[2]), 
            file="L.ALL.PFref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(NR.ALL.PFref.iHS.noFDR[2]), 
            file="NR.ALL.PFref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(OT.ALL.PFref.iHS.noFDR[2]), 
            file="OT.ALL.PFref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(PF.ALL.PFref.iHS.noFDR[2]), 
            file="PF.ALL.PFref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(V.ALL.PFref.iHS.noFDR[2]), 
            file="V.ALL.PFref.noPvalueCorrection.iHS.summary", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Perform xpEHH calculations with FDR p-value correction...
#OTref:
OT_NR_ALL_OTref <- ies2xpehh(wgs.OT.ALL.OTref, wgs.NR.ALL.OTref,
                             popname1 = "OT", popname2 = "NR",
                             include_freq = T, p.adjust.method="BY")
OT_PF_ALL_OTref <- ies2xpehh(wgs.OT.ALL.OTref, wgs.PF.ALL.OTref,
                             popname1 = "OT", popname2 = "PF",
                             include_freq = T, p.adjust.method="BY")
OT_V_ALL_OTref <- ies2xpehh(wgs.OT.ALL.OTref, wgs.V.ALL.OTref,
                            popname1 = "OT", popname2 = "V",
                            include_freq = T, p.adjust.method="BY")
OT_L_ALL_OTref <- ies2xpehh(wgs.OT.ALL.OTref, wgs.L.ALL.OTref,
                            popname1 = "OT", popname2 = "L",
                            include_freq = T, p.adjust.method="BY")
NR_PF_ALL_OTref <- ies2xpehh(wgs.NR.ALL.OTref, wgs.PF.ALL.OTref,
                             popname1 = "NR", popname2 = "PF",
                             include_freq = T, p.adjust.method="BY")
NR_V_ALL_OTref <- ies2xpehh(wgs.NR.ALL.OTref, wgs.V.ALL.OTref,
                            popname1 = "NR", popname2 = "V",
                            include_freq = T, p.adjust.method="BY")
NR_L_ALL_OTref <- ies2xpehh(wgs.NR.ALL.OTref, wgs.L.ALL.OTref,
                            popname1 = "NR", popname2 = "L",
                            include_freq = T, p.adjust.method="BY")
PF_V_ALL_OTref <- ies2xpehh(wgs.PF.ALL.OTref, wgs.V.ALL.OTref,
                            popname1 = "PF", popname2 = "V",
                            include_freq = T, p.adjust.method="BY")
PF_L_ALL_OTref <- ies2xpehh(wgs.PF.ALL.OTref, wgs.L.ALL.OTref,
                            popname1 = "PF", popname2 = "L",
                            include_freq = T, p.adjust.method="BY")
V_L_ALL_OTref <- ies2xpehh(wgs.V.ALL.OTref, wgs.L.ALL.OTref,
                           popname1 = "V", popname2 = "L",
                           include_freq = T, p.adjust.method="BY")
#PFref:
OT_NR_ALL_PFref <- ies2xpehh(wgs.OT.ALL.PFref, wgs.NR.ALL.PFref,
                             popname1 = "OT", popname2 = "NR",
                             include_freq = T, p.adjust.method="BY")
OT_PF_ALL_PFref <- ies2xpehh(wgs.OT.ALL.PFref, wgs.PF.ALL.PFref,
                             popname1 = "OT", popname2 = "PF",
                             include_freq = T, p.adjust.method="BY")
OT_V_ALL_PFref <- ies2xpehh(wgs.OT.ALL.PFref, wgs.V.ALL.PFref,
                            popname1 = "OT", popname2 = "V",
                            include_freq = T, p.adjust.method="BY")
OT_L_ALL_PFref <- ies2xpehh(wgs.OT.ALL.PFref, wgs.L.ALL.PFref,
                            popname1 = "OT", popname2 = "L",
                            include_freq = T, p.adjust.method="BY")
NR_PF_ALL_PFref <- ies2xpehh(wgs.NR.ALL.PFref, wgs.PF.ALL.PFref,
                             popname1 = "NR", popname2 = "PF",
                             include_freq = T, p.adjust.method="BY")
NR_V_ALL_PFref <- ies2xpehh(wgs.NR.ALL.PFref, wgs.V.ALL.PFref,
                            popname1 = "NR", popname2 = "V",
                            include_freq = T, p.adjust.method="BY")
NR_L_ALL_PFref <- ies2xpehh(wgs.NR.ALL.PFref, wgs.L.ALL.PFref,
                            popname1 = "NR", popname2 = "L",
                            include_freq = T, p.adjust.method="BY")
PF_V_ALL_PFref <- ies2xpehh(wgs.PF.ALL.PFref, wgs.V.ALL.PFref,
                            popname1 = "PF", popname2 = "V",
                            include_freq = T, p.adjust.method="BY")
PF_L_ALL_PFref <- ies2xpehh(wgs.PF.ALL.PFref, wgs.L.ALL.PFref,
                            popname1 = "PF", popname2 = "L",
                            include_freq = T, p.adjust.method="BY")
V_L_ALL_PFref <- ies2xpehh(wgs.V.ALL.PFref, wgs.L.ALL.PFref,
                           popname1 = "V", popname2 = "L",
                           include_freq = T, p.adjust.method="BY")

#Write to file for posterity...
#OTref:
write.table(OT_NR_ALL_OTref, 
            file="OT_NR_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_PF_ALL_OTref, 
            file="OT_PF_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_V_ALL_OTref, 
            file="OT_V_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_L_ALL_OTref, 
            file="OT_L_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_PF_ALL_OTref, 
            file="NR_PF_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_V_ALL_OTref, 
            file="NR_V_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_L_ALL_OTref, 
            file="NR_L_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_V_ALL_OTref, 
            file="PF_V_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_L_ALL_OTref, 
            file="PF_L_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(V_L_ALL_OTref, 
            file="V_L_ALL_OTref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(OT_NR_ALL_PFref, 
            file="OT_NR_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_PF_ALL_PFref, 
            file="OT_PF_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_V_ALL_PFref, 
            file="OT_V_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_L_ALL_PFref, 
            file="OT_L_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_PF_ALL_PFref, 
            file="NR_PF_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_V_ALL_PFref, 
            file="NR_V_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_L_ALL_PFref, 
            file="NR_L_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_V_ALL_PFref, 
            file="PF_V_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_L_ALL_PFref, 
            file="PF_L_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(V_L_ALL_PFref, 
            file="V_L_ALL_PFref.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Perform xpEHH calculations **without** FDR p-value correction...
#OTref:
OT_NR_ALL_OTref.noFDR <- ies2xpehh(wgs.OT.ALL.OTref, wgs.NR.ALL.OTref,
                                   popname1 = "OT", popname2 = "NR",
                                   include_freq = T)
OT_PF_ALL_OTref.noFDR <- ies2xpehh(wgs.OT.ALL.OTref, wgs.PF.ALL.OTref,
                                   popname1 = "OT", popname2 = "PF",
                                   include_freq = T)
OT_V_ALL_OTref.noFDR <- ies2xpehh(wgs.OT.ALL.OTref, wgs.V.ALL.OTref,
                                  popname1 = "OT", popname2 = "V",
                                  include_freq = T)
OT_L_ALL_OTref.noFDR <- ies2xpehh(wgs.OT.ALL.OTref, wgs.L.ALL.OTref,
                                  popname1 = "OT", popname2 = "L",
                                  include_freq = T)
NR_PF_ALL_OTref.noFDR <- ies2xpehh(wgs.NR.ALL.OTref, wgs.PF.ALL.OTref,
                                   popname1 = "NR", popname2 = "PF",
                                   include_freq = T)
NR_V_ALL_OTref.noFDR <- ies2xpehh(wgs.NR.ALL.OTref, wgs.V.ALL.OTref,
                                  popname1 = "NR", popname2 = "V",
                                  include_freq = T)
NR_L_ALL_OTref.noFDR <- ies2xpehh(wgs.NR.ALL.OTref, wgs.L.ALL.OTref,
                                  popname1 = "NR", popname2 = "L",
                                  include_freq = T)
PF_V_ALL_OTref.noFDR <- ies2xpehh(wgs.PF.ALL.OTref, wgs.V.ALL.OTref,
                                  popname1 = "PF", popname2 = "V",
                                  include_freq = T)
PF_L_ALL_OTref.noFDR <- ies2xpehh(wgs.PF.ALL.OTref, wgs.L.ALL.OTref,
                                  popname1 = "PF", popname2 = "L",
                                  include_freq = T)
V_L_ALL_OTref.noFDR <- ies2xpehh(wgs.V.ALL.OTref, wgs.L.ALL.OTref,
                                 popname1 = "V", popname2 = "L",
                                 include_freq = T)
#PFref:
OT_NR_ALL_PFref.noFDR <- ies2xpehh(wgs.OT.ALL.PFref, wgs.NR.ALL.PFref,
                                   popname1 = "OT", popname2 = "NR",
                                   include_freq = T)
OT_PF_ALL_PFref.noFDR <- ies2xpehh(wgs.OT.ALL.PFref, wgs.PF.ALL.PFref,
                                   popname1 = "OT", popname2 = "PF",
                                   include_freq = T)
OT_V_ALL_PFref.noFDR <- ies2xpehh(wgs.OT.ALL.PFref, wgs.V.ALL.PFref,
                                  popname1 = "OT", popname2 = "V",
                                  include_freq = T)
OT_L_ALL_PFref.noFDR <- ies2xpehh(wgs.OT.ALL.PFref, wgs.L.ALL.PFref,
                                  popname1 = "OT", popname2 = "L",
                                  include_freq = T)
NR_PF_ALL_PFref.noFDR <- ies2xpehh(wgs.NR.ALL.PFref, wgs.PF.ALL.PFref,
                                   popname1 = "NR", popname2 = "PF",
                                   include_freq = T)
NR_V_ALL_PFref.noFDR <- ies2xpehh(wgs.NR.ALL.PFref, wgs.V.ALL.PFref,
                                  popname1 = "NR", popname2 = "V",
                                  include_freq = T)
NR_L_ALL_PFref.noFDR <- ies2xpehh(wgs.NR.ALL.PFref, wgs.L.ALL.PFref,
                                  popname1 = "NR", popname2 = "L",
                                  include_freq = T)
PF_V_ALL_PFref.noFDR <- ies2xpehh(wgs.PF.ALL.PFref, wgs.V.ALL.PFref,
                                  popname1 = "PF", popname2 = "V",
                                  include_freq = T)
PF_L_ALL_PFref.noFDR <- ies2xpehh(wgs.PF.ALL.PFref, wgs.L.ALL.PFref,
                                  popname1 = "PF", popname2 = "L",
                                  include_freq = T)
V_L_ALL_PFref.noFDR <- ies2xpehh(wgs.V.ALL.PFref, wgs.L.ALL.PFref,
                                 popname1 = "V", popname2 = "L",
                                 include_freq = T)

#Write to file for posterity...
#OTref:
write.table(OT_NR_ALL_OTref.noFDR, 
            file="OT_NR_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_PF_ALL_OTref.noFDR, 
            file="OT_PF_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_V_ALL_OTref.noFDR, 
            file="OT_V_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_L_ALL_OTref.noFDR, 
            file="OT_L_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_PF_ALL_OTref.noFDR, 
            file="NR_PF_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_V_ALL_OTref.noFDR, 
            file="NR_V_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_L_ALL_OTref.noFDR, 
            file="NR_L_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_V_ALL_OTref.noFDR, 
            file="PF_V_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_L_ALL_OTref.noFDR, 
            file="PF_L_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(V_L_ALL_OTref.noFDR, 
            file="V_L_ALL_OTref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(OT_NR_ALL_PFref.noFDR, 
            file="OT_NR_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_PF_ALL_PFref.noFDR, 
            file="OT_PF_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_V_ALL_PFref.noFDR, 
            file="OT_V_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_L_ALL_PFref.noFDR, 
            file="OT_L_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_PF_ALL_PFref.noFDR, 
            file="NR_PF_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_V_ALL_PFref.noFDR, 
            file="NR_V_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_L_ALL_PFref.noFDR, 
            file="NR_L_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_V_ALL_PFref.noFDR, 
            file="PF_V_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_L_ALL_PFref.noFDR, 
            file="PF_L_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(V_L_ALL_PFref.noFDR, 
            file="V_L_ALL_PFref.noPvalueCorrection.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Perform Rsb calculations with FDR p-value correction...
#OTref:
OT_NR_ALL_OTref_Rsb <- ines2rsb(wgs.OT.ALL.OTref, wgs.NR.ALL.OTref,
                                popname1 = "OT", popname2 = "NR",
                                include_freq = T, p.adjust.method="BY")
OT_PF_ALL_OTref_Rsb <- ines2rsb(wgs.OT.ALL.OTref, wgs.PF.ALL.OTref,
                                popname1 = "OT", popname2 = "PF",
                                include_freq = T, p.adjust.method="BY")
OT_V_ALL_OTref_Rsb <- ines2rsb(wgs.OT.ALL.OTref, wgs.V.ALL.OTref,
                               popname1 = "OT", popname2 = "V",
                               include_freq = T, p.adjust.method="BY")
OT_L_ALL_OTref_Rsb <- ines2rsb(wgs.OT.ALL.OTref, wgs.L.ALL.OTref,
                               popname1 = "OT", popname2 = "L",
                               include_freq = T, p.adjust.method="BY")
NR_PF_ALL_OTref_Rsb <- ines2rsb(wgs.NR.ALL.OTref, wgs.PF.ALL.OTref,
                                popname1 = "NR", popname2 = "PF",
                                include_freq = T, p.adjust.method="BY")
NR_V_ALL_OTref_Rsb <- ines2rsb(wgs.NR.ALL.OTref, wgs.V.ALL.OTref,
                               popname1 = "NR", popname2 = "V",
                               include_freq = T, p.adjust.method="BY")
NR_L_ALL_OTref_Rsb <- ines2rsb(wgs.NR.ALL.OTref, wgs.L.ALL.OTref,
                               popname1 = "NR", popname2 = "L",
                               include_freq = T, p.adjust.method="BY")
PF_V_ALL_OTref_Rsb <- ines2rsb(wgs.PF.ALL.OTref, wgs.V.ALL.OTref,
                               popname1 = "PF", popname2 = "V",
                               include_freq = T, p.adjust.method="BY")
PF_L_ALL_OTref_Rsb <- ines2rsb(wgs.PF.ALL.OTref, wgs.L.ALL.OTref,
                               popname1 = "PF", popname2 = "L",
                               include_freq = T, p.adjust.method="BY")
V_L_ALL_OTref_Rsb <- ines2rsb(wgs.V.ALL.OTref, wgs.L.ALL.OTref,
                              popname1 = "V", popname2 = "L",
                              include_freq = T, p.adjust.method="BY")
#PFref:
OT_NR_ALL_PFref_Rsb <- ines2rsb(wgs.OT.ALL.PFref, wgs.NR.ALL.PFref,
                                popname1 = "OT", popname2 = "NR",
                                include_freq = T, p.adjust.method="BY")
OT_PF_ALL_PFref_Rsb <- ines2rsb(wgs.OT.ALL.PFref, wgs.PF.ALL.PFref,
                                popname1 = "OT", popname2 = "PF",
                                include_freq = T, p.adjust.method="BY")
OT_V_ALL_PFref_Rsb <- ines2rsb(wgs.OT.ALL.PFref, wgs.V.ALL.PFref,
                               popname1 = "OT", popname2 = "V",
                               include_freq = T, p.adjust.method="BY")
OT_L_ALL_PFref_Rsb <- ines2rsb(wgs.OT.ALL.PFref, wgs.L.ALL.PFref,
                               popname1 = "OT", popname2 = "L",
                               include_freq = T, p.adjust.method="BY")
NR_PF_ALL_PFref_Rsb <- ines2rsb(wgs.NR.ALL.PFref, wgs.PF.ALL.PFref,
                                popname1 = "NR", popname2 = "PF",
                                include_freq = T, p.adjust.method="BY")
NR_V_ALL_PFref_Rsb <- ines2rsb(wgs.NR.ALL.PFref, wgs.V.ALL.PFref,
                               popname1 = "NR", popname2 = "V",
                               include_freq = T, p.adjust.method="BY")
NR_L_ALL_PFref_Rsb <- ines2rsb(wgs.NR.ALL.PFref, wgs.L.ALL.PFref,
                               popname1 = "NR", popname2 = "L",
                               include_freq = T, p.adjust.method="BY")
PF_V_ALL_PFref_Rsb <- ines2rsb(wgs.PF.ALL.PFref, wgs.V.ALL.PFref,
                               popname1 = "PF", popname2 = "V",
                               include_freq = T, p.adjust.method="BY")
PF_L_ALL_PFref_Rsb <- ines2rsb(wgs.PF.ALL.PFref, wgs.L.ALL.PFref,
                               popname1 = "PF", popname2 = "L",
                               include_freq = T, p.adjust.method="BY")
V_L_ALL_PFref_Rsb <- ines2rsb(wgs.V.ALL.PFref, wgs.L.ALL.PFref,
                              popname1 = "V", popname2 = "L",
                              include_freq = T, p.adjust.method="BY")

#Write to file for posterity...
#OTref:
write.table(OT_NR_ALL_OTref_Rsb, 
            file="OT_NR_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_PF_ALL_OTref_Rsb, 
            file="OT_PF_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_V_ALL_OTref_Rsb, 
            file="OT_V_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_L_ALL_OTref_Rsb, 
            file="OT_L_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_PF_ALL_OTref_Rsb, 
            file="NR_PF_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_V_ALL_OTref_Rsb, 
            file="NR_V_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_L_ALL_OTref_Rsb, 
            file="NR_L_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_V_ALL_OTref_Rsb, 
            file="PF_V_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_L_ALL_OTref_Rsb, 
            file="PF_L_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(V_L_ALL_OTref_Rsb, 
            file="V_L_ALL_OTref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(OT_NR_ALL_PFref_Rsb, 
            file="OT_NR_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_PF_ALL_PFref_Rsb, 
            file="OT_PF_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_V_ALL_PFref_Rsb, 
            file="OT_V_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_L_ALL_PFref_Rsb, 
            file="OT_L_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_PF_ALL_PFref_Rsb, 
            file="NR_PF_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_V_ALL_PFref_Rsb, 
            file="NR_V_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_L_ALL_PFref_Rsb, 
            file="NR_L_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_V_ALL_PFref_Rsb, 
            file="PF_V_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_L_ALL_PFref_Rsb, 
            file="PF_L_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(V_L_ALL_PFref_Rsb, 
            file="V_L_ALL_PFref.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Perform Rsb calculations **without** FDR p-value correction...
#OTref:
OT_NR_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.OT.ALL.OTref, wgs.NR.ALL.OTref,
                                      popname1 = "OT", popname2 = "NR",
                                      include_freq = T)
OT_PF_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.OT.ALL.OTref, wgs.PF.ALL.OTref,
                                      popname1 = "OT", popname2 = "PF",
                                      include_freq = T)
OT_V_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.OT.ALL.OTref, wgs.V.ALL.OTref,
                                     popname1 = "OT", popname2 = "V",
                                     include_freq = T)
OT_L_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.OT.ALL.OTref, wgs.L.ALL.OTref,
                                     popname1 = "OT", popname2 = "L",
                                     include_freq = T)
NR_PF_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.NR.ALL.OTref, wgs.PF.ALL.OTref,
                                      popname1 = "NR", popname2 = "PF",
                                      include_freq = T)
NR_V_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.NR.ALL.OTref, wgs.V.ALL.OTref,
                                     popname1 = "NR", popname2 = "V",
                                     include_freq = T)
NR_L_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.NR.ALL.OTref, wgs.L.ALL.OTref,
                                     popname1 = "NR", popname2 = "L",
                                     include_freq = T)
PF_V_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.PF.ALL.OTref, wgs.V.ALL.OTref,
                                     popname1 = "PF", popname2 = "V",
                                     include_freq = T)
PF_L_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.PF.ALL.OTref, wgs.L.ALL.OTref,
                                     popname1 = "PF", popname2 = "L",
                                     include_freq = T)
V_L_ALL_OTref_Rsb.noFDR <- ines2rsb(wgs.V.ALL.OTref, wgs.L.ALL.OTref,
                                    popname1 = "V", popname2 = "L",
                                    include_freq = T)
#PFref:
OT_NR_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.OT.ALL.PFref, wgs.NR.ALL.PFref,
                                      popname1 = "OT", popname2 = "NR",
                                      include_freq = T)
OT_PF_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.OT.ALL.PFref, wgs.PF.ALL.PFref,
                                      popname1 = "OT", popname2 = "PF",
                                      include_freq = T)
OT_V_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.OT.ALL.PFref, wgs.V.ALL.PFref,
                                     popname1 = "OT", popname2 = "V",
                                     include_freq = T)
OT_L_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.OT.ALL.PFref, wgs.L.ALL.PFref,
                                     popname1 = "OT", popname2 = "L",
                                     include_freq = T)
NR_PF_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.NR.ALL.PFref, wgs.PF.ALL.PFref,
                                      popname1 = "NR", popname2 = "PF",
                                      include_freq = T)
NR_V_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.NR.ALL.PFref, wgs.V.ALL.PFref,
                                     popname1 = "NR", popname2 = "V",
                                     include_freq = T)
NR_L_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.NR.ALL.PFref, wgs.L.ALL.PFref,
                                     popname1 = "NR", popname2 = "L",
                                     include_freq = T)
PF_V_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.PF.ALL.PFref, wgs.V.ALL.PFref,
                                     popname1 = "PF", popname2 = "V",
                                     include_freq = T)
PF_L_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.PF.ALL.PFref, wgs.L.ALL.PFref,
                                     popname1 = "PF", popname2 = "L",
                                     include_freq = T)
V_L_ALL_PFref_Rsb.noFDR <- ines2rsb(wgs.V.ALL.PFref, wgs.L.ALL.PFref,
                                    popname1 = "V", popname2 = "L",
                                    include_freq = T)

#Write to file for posterity...
#OTref:
write.table(OT_NR_ALL_OTref_Rsb.noFDR, 
            file="OT_NR_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_PF_ALL_OTref_Rsb.noFDR, 
            file="OT_PF_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_V_ALL_OTref_Rsb.noFDR, 
            file="OT_V_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_L_ALL_OTref_Rsb.noFDR, 
            file="OT_L_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_PF_ALL_OTref_Rsb.noFDR, 
            file="NR_PF_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_V_ALL_OTref_Rsb.noFDR, 
            file="NR_V_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_L_ALL_OTref_Rsb.noFDR, 
            file="NR_L_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_V_ALL_OTref_Rsb.noFDR, 
            file="PF_V_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_L_ALL_OTref_Rsb.noFDR, 
            file="PF_L_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(V_L_ALL_OTref_Rsb.noFDR, 
            file="V_L_ALL_OTref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(OT_NR_ALL_PFref_Rsb.noFDR, 
            file="OT_NR_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_PF_ALL_PFref_Rsb.noFDR, 
            file="OT_PF_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_V_ALL_PFref_Rsb.noFDR, 
            file="OT_V_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OT_L_ALL_PFref_Rsb.noFDR, 
            file="OT_L_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_PF_ALL_PFref_Rsb.noFDR, 
            file="NR_PF_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_V_ALL_PFref_Rsb.noFDR, 
            file="NR_V_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(NR_L_ALL_PFref_Rsb.noFDR, 
            file="NR_L_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_V_ALL_PFref_Rsb.noFDR, 
            file="PF_V_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_L_ALL_PFref_Rsb.noFDR, 
            file="PF_L_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(V_L_ALL_PFref_Rsb.noFDR, 
            file="V_L_ALL_PFref.noPvalueCorrection.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)

### Calculation of candidate regions ###

#iHS - without FDR correction
#Calculate candidate regions in 50kb windows, 0 overalp, threshold = 6, no joining of neighbouring regions (from iHS P values)...
#OTref:
cr.L.OTref <- calc_candidate_regions(L.ALL.OTref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR.OTref <- calc_candidate_regions(NR.ALL.OTref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT.OTref <- calc_candidate_regions(OT.ALL.OTref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF.OTref <- calc_candidate_regions(PF.ALL.OTref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V.OTref <- calc_candidate_regions(V.ALL.OTref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
#PFref:
cr.L.PFref <- calc_candidate_regions(L.ALL.PFref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR.PFref <- calc_candidate_regions(NR.ALL.PFref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT.PFref <- calc_candidate_regions(OT.ALL.PFref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF.PFref <- calc_candidate_regions(PF.ALL.PFref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V.PFref <- calc_candidate_regions(V.ALL.PFref.iHS.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)

#xpEHH - without FDR correction
#Calculate candidate regions in 50kb windows, 0 overalp, threshold = 6, (from xpEHH P values)...
#OTref:
cr.OT_NR.OTref <- calc_candidate_regions(OT_NR_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_PF.OTref <- calc_candidate_regions(OT_PF_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_V.OTref <- calc_candidate_regions(OT_V_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_L.OTref <- calc_candidate_regions(OT_L_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_PF.OTref <- calc_candidate_regions(NR_PF_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_V.OTref <- calc_candidate_regions(NR_V_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_L.OTref <- calc_candidate_regions(NR_L_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_V.OTref <- calc_candidate_regions(PF_V_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_L.OTref <- calc_candidate_regions(PF_L_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V_L.OTref <- calc_candidate_regions(V_L_ALL_OTref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
#PFref:
cr.OT_NR.PFref <- calc_candidate_regions(OT_NR_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_PF.PFref <- calc_candidate_regions(OT_PF_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_V.PFref <- calc_candidate_regions(OT_V_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_L.PFref <- calc_candidate_regions(OT_L_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_PF.PFref <- calc_candidate_regions(NR_PF_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_V.PFref <- calc_candidate_regions(NR_V_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_L.PFref <- calc_candidate_regions(NR_L_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_V.PFref <- calc_candidate_regions(PF_V_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_L.PFref <- calc_candidate_regions(PF_L_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V_L.PFref <- calc_candidate_regions(V_L_ALL_PFref.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)

#Rsb - without FDR correction
#Calculate candidate regions in 50kb windows, 0 overalp, threshold = 6, (from Rsb P values)...
#OTref:
cr.OT_NR.Rsb.OTref <- calc_candidate_regions(OT_NR_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_PF.Rsb.OTref <- calc_candidate_regions(OT_PF_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_V.Rsb.OTref <- calc_candidate_regions(OT_V_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_L.Rsb.OTref <- calc_candidate_regions(OT_L_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_PF.Rsb.OTref <- calc_candidate_regions(NR_PF_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_V.Rsb.OTref <- calc_candidate_regions(NR_V_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_L.Rsb.OTref <- calc_candidate_regions(NR_L_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_V.Rsb.OTref <- calc_candidate_regions(PF_V_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_L.Rsb.OTref <- calc_candidate_regions(PF_L_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V_L.Rsb.OTref <- calc_candidate_regions(V_L_ALL_OTref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
#PFref:
cr.OT_NR.Rsb.PFref <- calc_candidate_regions(OT_NR_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_PF.Rsb.PFref <- calc_candidate_regions(OT_PF_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_V.Rsb.PFref <- calc_candidate_regions(OT_V_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_L.Rsb.PFref <- calc_candidate_regions(OT_L_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_PF.Rsb.PFref <- calc_candidate_regions(NR_PF_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_V.Rsb.PFref <- calc_candidate_regions(NR_V_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_L.Rsb.PFref <- calc_candidate_regions(NR_L_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_V.Rsb.PFref <- calc_candidate_regions(PF_V_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_L.Rsb.PFref <- calc_candidate_regions(PF_L_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V_L.Rsb.PFref <- calc_candidate_regions(V_L_ALL_PFref_Rsb.noFDR, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)

#iHS - with FDR correction
#Calculate candidate regions in 50kb windows, 0 overalp, threshold = 6, (from iHS P values)...
#OTref:
cr.L.OTref.withFDR <- calc_candidate_regions(L.ALL.OTref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR.OTref.withFDR <- calc_candidate_regions(NR.ALL.OTref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT.OTref.withFDR <- calc_candidate_regions(OT.ALL.OTref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF.OTref.withFDR <- calc_candidate_regions(PF.ALL.OTref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V.OTref.withFDR <- calc_candidate_regions(V.ALL.OTref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
#PFref:
cr.L.PFref.withFDR <- calc_candidate_regions(L.ALL.PFref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR.PFref.withFDR <- calc_candidate_regions(NR.ALL.PFref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT.PFref.withFDR <- calc_candidate_regions(OT.ALL.PFref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF.PFref.withFDR <- calc_candidate_regions(PF.ALL.PFref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V.PFref.withFDR <- calc_candidate_regions(V.ALL.PFref.iHS, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)

#xpEHH - with FDR correction
#Calculate candidate regions in 50kb windows, 0 overalp, threshold = 6, (from xpEHH P values)...
#OTref:
cr.OT_NR.OTref.withFDR <- calc_candidate_regions(OT_NR_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_PF.OTref.withFDR <- calc_candidate_regions(OT_PF_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_V.OTref.withFDR <- calc_candidate_regions(OT_V_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_L.OTref.withFDR <- calc_candidate_regions(OT_L_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_PF.OTref.withFDR <- calc_candidate_regions(NR_PF_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_V.OTref.withFDR <- calc_candidate_regions(NR_V_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_L.OTref.withFDR <- calc_candidate_regions(NR_L_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_V.OTref.withFDR <- calc_candidate_regions(PF_V_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_L.OTref.withFDR <- calc_candidate_regions(PF_L_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V_L.OTref.withFDR <- calc_candidate_regions(V_L_ALL_OTref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
#PFref:
cr.OT_NR.PFref.withFDR <- calc_candidate_regions(OT_NR_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_PF.PFref.withFDR <- calc_candidate_regions(OT_PF_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_V.PFref.withFDR <- calc_candidate_regions(OT_V_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_L.PFref.withFDR <- calc_candidate_regions(OT_L_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_PF.PFref.withFDR <- calc_candidate_regions(NR_PF_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_V.PFref.withFDR <- calc_candidate_regions(NR_V_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_L.PFref.withFDR <- calc_candidate_regions(NR_L_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_V.PFref.withFDR <- calc_candidate_regions(PF_V_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_L.PFref.withFDR <- calc_candidate_regions(PF_L_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V_L.PFref.withFDR <- calc_candidate_regions(V_L_ALL_PFref, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)

#Rsb - with FDR correction
#Calculate candidate regions in 50kb windows, 0 overalp, threshold = 6, (from Rsb P values)...
#OTref:
cr.OT_NR.Rsb.OTref.withFDR <- calc_candidate_regions(OT_NR_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_PF.Rsb.OTref.withFDR <- calc_candidate_regions(OT_PF_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_V.Rsb.OTref.withFDR <- calc_candidate_regions(OT_V_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_L.Rsb.OTref.withFDR <- calc_candidate_regions(OT_L_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_PF.Rsb.OTref.withFDR <- calc_candidate_regions(NR_PF_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_V.Rsb.OTref.withFDR <- calc_candidate_regions(NR_V_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_L.Rsb.OTref.withFDR <- calc_candidate_regions(NR_L_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_V.Rsb.OTref.withFDR <- calc_candidate_regions(PF_V_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_L.Rsb.OTref.withFDR <- calc_candidate_regions(PF_L_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V_L.Rsb.OTref.withFDR <- calc_candidate_regions(V_L_ALL_OTref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
#PFref:
cr.OT_NR.Rsb.PFref.withFDR <- calc_candidate_regions(OT_NR_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_PF.Rsb.PFref.withFDR <- calc_candidate_regions(OT_PF_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_V.Rsb.PFref.withFDR <- calc_candidate_regions(OT_V_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.OT_L.Rsb.PFref.withFDR <- calc_candidate_regions(OT_L_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_PF.Rsb.PFref.withFDR <- calc_candidate_regions(NR_PF_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_V.Rsb.PFref.withFDR <- calc_candidate_regions(NR_V_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.NR_L.Rsb.PFref.withFDR <- calc_candidate_regions(NR_L_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_V.Rsb.PFref.withFDR <- calc_candidate_regions(PF_V_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.PF_L.Rsb.PFref.withFDR <- calc_candidate_regions(PF_L_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)
cr.V_L.Rsb.PFref.withFDR <- calc_candidate_regions(V_L_ALL_PFref_Rsb, window_size=5E4, overlap=0, threshold=6, min_n_extr_mrk=1, pval=TRUE, join_neighbors=FALSE)


#Write to file for posterity

#iHS (no FDR)...
#OTref:
write.table(cr.L.OTref, file="cr.L.OTref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR.OTref, file="cr.NR.OTref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT.OTref, file="cr.OT.OTref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF.OTref, file="cr.PF.OTref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V.OTref, file="cr.V.OTref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(cr.L.PFref, file="cr.L.PFref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR.PFref, file="cr.NR.PFref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT.PFref, file="cr.OT.PFref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF.PFref, file="cr.PF.PFref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V.PFref, file="cr.V.PFref.noPvalueCorrection.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)

#xpEHH (no FDR)...
#OTref:
write.table(cr.OT_NR.OTref, file="cr.OT_NR.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.OTref, file="cr.OT_PF.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.OTref, file="cr.OT_V.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.OTref, file="cr.OT_L.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.OTref, file="cr.NR_PF.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.OTref, file="cr.NR_V.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.OTref, file="cr.NR_L.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.OTref, file="cr.PF_V.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.OTref, file="cr.PF_L.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.OTref, file="cr.V_L.OTref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(cr.OT_NR.PFref, file="cr.OT_NR.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.PFref, file="cr.OT_PF.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.PFref, file="cr.OT_V.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.PFref, file="cr.OT_L.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.PFref, file="cr.NR_PF.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.PFref, file="cr.NR_V.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.PFref, file="cr.NR_L.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.PFref, file="cr.PF_V.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.PFref, file="cr.PF_L.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.PFref, file="cr.V_L.PFref.noPvalueCorrection.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Rsb (no FDR)...
#OTref:
write.table(cr.OT_NR.Rsb.OTref, file="cr.OT_NR.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.Rsb.OTref, file="cr.OT_PF.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.Rsb.OTref, file="cr.OT_V.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.Rsb.OTref, file="cr.OT_L.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.Rsb.OTref, file="cr.NR_PF.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.Rsb.OTref, file="cr.NR_V.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.Rsb.OTref, file="cr.NR_L.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.Rsb.OTref, file="cr.PF_V.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.Rsb.OTref, file="cr.PF_L.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.Rsb.OTref, file="cr.V_L.OTref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(cr.OT_NR.Rsb.PFref, file="cr.OT_NR.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.Rsb.PFref, file="cr.OT_PF.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.Rsb.PFref, file="cr.OT_V.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.Rsb.PFref, file="cr.OT_L.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.Rsb.PFref, file="cr.NR_PF.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.Rsb.PFref, file="cr.NR_V.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.Rsb.PFref, file="cr.NR_L.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.Rsb.PFref, file="cr.PF_V.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.Rsb.PFref, file="cr.PF_L.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.Rsb.PFref, file="cr.V_L.PFref.noPvalueCorrection.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)

#iHS (with FDR)...
#OTref:
write.table(cr.L.OTref.withFDR, file="cr.L.OTref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR.OTref.withFDR, file="cr.NR.OTref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT.OTref.withFDR, file="cr.OT.OTref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF.OTref.withFDR, file="cr.PF.OTref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V.OTref.withFDR, file="cr.V.OTref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(cr.L.PFref.withFDR, file="cr.L.PFref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR.PFref.withFDR, file="cr.NR.PFref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT.PFref.withFDR, file="cr.OT.PFref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF.PFref.withFDR, file="cr.PF.PFref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V.PFref.withFDR, file="cr.V.PFref.noOverlap.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)

#xpEHH (with FDR)...
#OTref:
write.table(cr.OT_NR.OTref.withFDR, file="cr.OT_NR.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.OTref.withFDR, file="cr.OT_PF.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.OTref.withFDR, file="cr.OT_V.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.OTref.withFDR, file="cr.OT_L.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.OTref.withFDR, file="cr.NR_PF.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.OTref.withFDR, file="cr.NR_V.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.OTref.withFDR, file="cr.NR_L.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.OTref.withFDR, file="cr.PF_V.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.OTref.withFDR, file="cr.PF_L.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.OTref.withFDR, file="cr.V_L.OTref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(cr.OT_NR.PFref.withFDR, file="cr.OT_NR.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.PFref.withFDR, file="cr.OT_PF.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.PFref.withFDR, file="cr.OT_V.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.PFref.withFDR, file="cr.OT_L.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.PFref.withFDR, file="cr.NR_PF.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.PFref.withFDR, file="cr.NR_V.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.PFref.withFDR, file="cr.NR_L.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.PFref.withFDR, file="cr.PF_V.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.PFref.withFDR, file="cr.PF_L.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.PFref.withFDR, file="cr.V_L.PFref.noOverlap.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Rsb (with FDR)...
#OTref:
write.table(cr.OT_NR.Rsb.OTref.withFDR, file="cr.OT_NR.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.Rsb.OTref.withFDR, file="cr.OT_PF.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.Rsb.OTref.withFDR, file="cr.OT_V.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.Rsb.OTref.withFDR, file="cr.OT_L.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.Rsb.OTref.withFDR, file="cr.NR_PF.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.Rsb.OTref.withFDR, file="cr.NR_V.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.Rsb.OTref.withFDR, file="cr.NR_L.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.Rsb.OTref.withFDR, file="cr.PF_V.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.Rsb.OTref.withFDR, file="cr.PF_L.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.Rsb.OTref.withFDR, file="cr.V_L.OTref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(cr.OT_NR.Rsb.PFref.withFDR, file="cr.OT_NR.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.Rsb.PFref.withFDR, file="cr.OT_PF.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.Rsb.PFref.withFDR, file="cr.OT_V.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.Rsb.PFref.withFDR, file="cr.OT_L.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.Rsb.PFref.withFDR, file="cr.NR_PF.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.Rsb.PFref.withFDR, file="cr.NR_V.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.Rsb.PFref.withFDR, file="cr.NR_L.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.Rsb.PFref.withFDR, file="cr.PF_V.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.Rsb.PFref.withFDR, file="cr.PF_L.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.Rsb.PFref.withFDR, file="cr.V_L.PFref.noOverlap.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Amalgamate candidate regions...

#First read scaffold length data
OTref_scaffold_lengths <- read.table(file="/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/Additional_files/OT3b_lengths.txt", 
                                     quote="", stringsAsFactors=FALSE)
colnames(OTref_scaffold_lengths) <- c("CHR", "LENGTH")
PFref_scaffold_lengths <- read.table(file="/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/Additional_files/PF_lengths.txt", 
                                     quote="", stringsAsFactors=FALSE)
colnames(PFref_scaffold_lengths) <- c("CHR", "LENGTH")

#Next create a dummy data frame for each reference, splitting each scaffold into 50kb ranges
split.scaffolds <- function(length.df){
  #Create empty result data frame
  main.result.df = data.frame(CHR=numeric(), START=numeric(), END=numeric())
  #Remove scaffolds < 50kb
  length.df = length.df[which(length.df$LENGTH >= 50000), ]
  #Iterate over scaffolds
  for(i in 1:nrow(length.df)){
    #Store scaffold
    scaffold = length.df$CHR[i]
    #Store length
    scafflen = length.df$LENGTH[i]
    #Make vector of 50kb positions, starting from 0, along the ith scaffold's length
    posvec = seq(from=0, to=scafflen, by=50000)
    #Create CHR, START and END vectors for this scaffold
    startvec=posvec[1:length(posvec)-1]
    endvec=posvec[2:length(posvec)]
    chrvec=rep(scaffold, length(posvec)-1)
    #Bind together in a temporary result data frame
    tmp.df = data.frame(CHR=chrvec, START=startvec, END=endvec)
    #rbind to main result data frame
    main.result.df = rbind(main.result.df, tmp.df)
  }
  #Return main result data frame
  return(main.result.df)
}

#Run the function on each reference
OTref_in_50kb <- split.scaffolds(OTref_scaffold_lengths)
PFref_in_50kb <- split.scaffolds(PFref_scaffold_lengths)

#A function to create a combined position column (used in the following function)
poscolumn <- function(df){
  #Do the pasting
  df$POS <- paste(df$CHR, format(df$START, scientific=F), sep="-")
  df$POS <- paste(df$POS, format(df$END, scientific=F), sep=":")
  #Remove white space which appears
  df <- df %>% mutate_if(is.character, str_squish)
  df <- df %>% mutate(POS = str_replace_all(POS, " ", ""))
}

#All-in-one amalgamation of candidate regions function
#--lots of inputs, must be supplied in correct order
amalgamate.cr <- function(scaffolds_in_50kb, 
                          iHS.cr.PF, iHS.cr.L, iHS.cr.OT, iHS.cr.V, iHS.cr.NR, 
                          iHS.cr.PF.FDR, iHS.cr.L.FDR, iHS.cr.OT.FDR, iHS.cr.V.FDR, iHS.cr.NR.FDR, 
                          cr.xpEHH.OT_NR, cr.xpEHH.OT_PF, cr.xpEHH.OT_V, cr.xpEHH.OT_L, cr.xpEHH.NR_PF, 
                          cr.xpEHH.NR_V, cr.xpEHH.NR_L, cr.xpEHH.PF_V, cr.xpEHH.PF_L, cr.xpEHH.V_L, 
                          cr.xpEHH.OT_NR.FDR, cr.xpEHH.OT_PF.FDR, cr.xpEHH.OT_V.FDR, cr.xpEHH.OT_L.FDR, cr.xpEHH.NR_PF.FDR, 
                          cr.xpEHH.NR_V.FDR, cr.xpEHH.NR_L.FDR, cr.xpEHH.PF_V.FDR, cr.xpEHH.PF_L.FDR, cr.xpEHH.V_L.FDR, 
                          cr.Rsb.OT_NR, cr.Rsb.OT_PF, cr.Rsb.OT_V, cr.Rsb.OT_L, cr.Rsb.NR_PF, 
                          cr.Rsb.NR_V, cr.Rsb.NR_L, cr.Rsb.PF_V, cr.Rsb.PF_L, cr.Rsb.V_L, 
                          cr.Rsb.OT_NR.FDR, cr.Rsb.OT_PF.FDR, cr.Rsb.OT_V.FDR, cr.Rsb.OT_L.FDR, cr.Rsb.NR_PF.FDR, 
                          cr.Rsb.NR_V.FDR, cr.Rsb.NR_L.FDR, cr.Rsb.PF_V.FDR, cr.Rsb.PF_L.FDR, cr.Rsb.V_L.FDR){
  #Add position column to each cr data frame
  scaffolds_in_50kb = poscolumn(scaffolds_in_50kb)
  iHS.cr.PF = poscolumn(iHS.cr.PF)
  iHS.cr.L = poscolumn(iHS.cr.L)
  iHS.cr.OT = poscolumn(iHS.cr.OT)
  iHS.cr.V = poscolumn(iHS.cr.V)
  iHS.cr.NR = poscolumn(iHS.cr.NR)
  cr.xpEHH.OT_NR = poscolumn(cr.xpEHH.OT_NR)
  cr.xpEHH.OT_PF = poscolumn(cr.xpEHH.OT_PF)
  cr.xpEHH.OT_V = poscolumn(cr.xpEHH.OT_V)
  cr.xpEHH.OT_L = poscolumn(cr.xpEHH.OT_L)
  cr.xpEHH.NR_PF = poscolumn(cr.xpEHH.NR_PF)
  cr.xpEHH.NR_V = poscolumn(cr.xpEHH.NR_V)
  cr.xpEHH.NR_L = poscolumn(cr.xpEHH.NR_L)
  cr.xpEHH.PF_V = poscolumn(cr.xpEHH.PF_V)
  cr.xpEHH.PF_L = poscolumn(cr.xpEHH.PF_L)
  cr.xpEHH.V_L = poscolumn(cr.xpEHH.V_L)
  cr.Rsb.OT_NR = poscolumn(cr.Rsb.OT_NR)
  cr.Rsb.OT_PF = poscolumn(cr.Rsb.OT_PF)
  cr.Rsb.OT_V = poscolumn(cr.Rsb.OT_V)
  cr.Rsb.OT_L = poscolumn(cr.Rsb.OT_L)
  cr.Rsb.NR_PF = poscolumn(cr.Rsb.NR_PF)
  cr.Rsb.NR_V = poscolumn(cr.Rsb.NR_V)
  cr.Rsb.NR_L = poscolumn(cr.Rsb.NR_L)
  cr.Rsb.PF_V = poscolumn(cr.Rsb.PF_V)
  cr.Rsb.PF_L = poscolumn(cr.Rsb.PF_L)
  cr.Rsb.V_L = poscolumn(cr.Rsb.V_L)
  iHS.cr.PF.FDR = poscolumn(iHS.cr.PF.FDR)
  iHS.cr.L.FDR = poscolumn(iHS.cr.L.FDR)
  iHS.cr.OT.FDR = poscolumn(iHS.cr.OT.FDR)
  iHS.cr.V.FDR = poscolumn(iHS.cr.V.FDR)
  iHS.cr.NR.FDR = poscolumn(iHS.cr.NR.FDR)
  cr.xpEHH.OT_NR.FDR = poscolumn(cr.xpEHH.OT_NR.FDR)
  cr.xpEHH.OT_PF.FDR = poscolumn(cr.xpEHH.OT_PF.FDR)
  cr.xpEHH.OT_V.FDR = poscolumn(cr.xpEHH.OT_V.FDR)
  cr.xpEHH.OT_L.FDR = poscolumn(cr.xpEHH.OT_L.FDR)
  cr.xpEHH.NR_PF.FDR = poscolumn(cr.xpEHH.NR_PF.FDR)
  cr.xpEHH.NR_V.FDR = poscolumn(cr.xpEHH.NR_V.FDR)
  cr.xpEHH.NR_L.FDR = poscolumn(cr.xpEHH.NR_L.FDR)
  cr.xpEHH.PF_V.FDR = poscolumn(cr.xpEHH.PF_V.FDR)
  cr.xpEHH.PF_L.FDR = poscolumn(cr.xpEHH.PF_L.FDR)
  cr.xpEHH.V_L.FDR = poscolumn(cr.xpEHH.V_L.FDR)
  cr.Rsb.OT_NR.FDR = poscolumn(cr.Rsb.OT_NR.FDR)
  cr.Rsb.OT_PF.FDR = poscolumn(cr.Rsb.OT_PF.FDR)
  cr.Rsb.OT_V.FDR = poscolumn(cr.Rsb.OT_V.FDR)
  cr.Rsb.OT_L.FDR = poscolumn(cr.Rsb.OT_L.FDR)
  cr.Rsb.NR_PF.FDR = poscolumn(cr.Rsb.NR_PF.FDR)
  cr.Rsb.NR_V.FDR = poscolumn(cr.Rsb.NR_V.FDR)
  cr.Rsb.NR_L.FDR = poscolumn(cr.Rsb.NR_L.FDR)
  cr.Rsb.PF_V.FDR = poscolumn(cr.Rsb.PF_V.FDR)
  cr.Rsb.PF_L.FDR = poscolumn(cr.Rsb.PF_L.FDR)
  cr.Rsb.V_L.FDR = poscolumn(cr.Rsb.V_L.FDR)
  #Subset scaffolds in 50kb data frame to scaffolds contained in candidate region data frames
  #and store in main result data frame
  scaffold_list = unique(c(iHS.cr.PF$CHR, iHS.cr.L$CHR, iHS.cr.OT$CHR, iHS.cr.V$CHR, iHS.cr.NR$CHR, 
                           iHS.cr.PF.FDR$CHR, iHS.cr.L.FDR$CHR, iHS.cr.OT.FDR$CHR, iHS.cr.V.FDR$CHR, iHS.cr.NR.FDR$CHR, 
                           cr.xpEHH.OT_NR$CHR, cr.xpEHH.OT_PF$CHR, cr.xpEHH.OT_V$CHR, cr.xpEHH.OT_L$CHR, cr.xpEHH.NR_PF$CHR, 
                           cr.xpEHH.NR_V$CHR, cr.xpEHH.NR_L$CHR, cr.xpEHH.PF_V$CHR, cr.xpEHH.PF_L$CHR, cr.xpEHH.V_L$CHR, 
                           cr.xpEHH.OT_NR.FDR$CHR, cr.xpEHH.OT_PF.FDR$CHR, cr.xpEHH.OT_V.FDR$CHR, cr.xpEHH.OT_L.FDR$CHR, cr.xpEHH.NR_PF.FDR$CHR, 
                           cr.xpEHH.NR_V.FDR$CHR, cr.xpEHH.NR_L.FDR$CHR, cr.xpEHH.PF_V.FDR$CHR, cr.xpEHH.PF_L.FDR$CHR, cr.xpEHH.V_L.FDR$CHR, 
                           cr.Rsb.OT_NR$CHR, cr.Rsb.OT_PF$CHR, cr.Rsb.OT_V$CHR, cr.Rsb.OT_L$CHR, cr.Rsb.NR_PF$CHR, 
                           cr.Rsb.NR_V$CHR, cr.Rsb.NR_L$CHR, cr.Rsb.PF_V$CHR, cr.Rsb.PF_L$CHR, cr.Rsb.V_L$CHR, 
                           cr.Rsb.OT_NR.FDR$CHR, cr.Rsb.OT_PF.FDR$CHR, cr.Rsb.OT_V.FDR$CHR, cr.Rsb.OT_L.FDR$CHR, cr.Rsb.NR_PF.FDR$CHR, 
                           cr.Rsb.NR_V.FDR$CHR, cr.Rsb.NR_L.FDR$CHR, cr.Rsb.PF_V.FDR$CHR, cr.Rsb.PF_L.FDR$CHR, cr.Rsb.V_L.FDR$CHR))
  main.result = scaffolds_in_50kb[which(scaffolds_in_50kb$CHR %in% scaffold_list), ]
  #Create empty columns for main result data frame
  main.result$iHS.PF <- 0
  main.result$iHS.PF.FDR <- 0
  main.result$iHS.L <- 0
  main.result$iHS.L.FDR <- 0
  main.result$iHS.OT <- 0
  main.result$iHS.OT.FDR <- 0
  main.result$iHS.V <- 0
  main.result$iHS.V.FDR <- 0
  main.result$iHS.NR <- 0
  main.result$iHS.NR.FDR <- 0
  main.result$xpEHH.OT_NR <- 0
  main.result$xpEHH.OT_NR.FDR <- 0
  main.result$xpEHH.OT_PF <- 0
  main.result$xpEHH.OT_PF.FDR <- 0
  main.result$xpEHH.OT_V <- 0
  main.result$xpEHH.OT_V.FDR <- 0
  main.result$xpEHH.OT_L <- 0
  main.result$xpEHH.OT_L.FDR <- 0
  main.result$xpEHH.NR_PF <- 0
  main.result$xpEHH.NR_PF.FDR <- 0
  main.result$xpEHH.NR_V <- 0
  main.result$xpEHH.NR_V.FDR <- 0
  main.result$xpEHH.NR_L <- 0
  main.result$xpEHH.NR_L.FDR <- 0
  main.result$xpEHH.PF_V <- 0
  main.result$xpEHH.PF_V.FDR <- 0
  main.result$xpEHH.PF_L <- 0
  main.result$xpEHH.PF_L.FDR <- 0
  main.result$xpEHH.V_L <- 0
  main.result$xpEHH.V_L.FDR <- 0
  main.result$Rsb.OT_NR <- 0
  main.result$Rsb.OT_NR.FDR <- 0
  main.result$Rsb.OT_PF <- 0
  main.result$Rsb.OT_PF.FDR <- 0
  main.result$Rsb.OT_V <- 0
  main.result$Rsb.OT_V.FDR <- 0
  main.result$Rsb.OT_L <- 0
  main.result$Rsb.OT_L.FDR <- 0
  main.result$Rsb.NR_PF <- 0
  main.result$Rsb.NR_PF.FDR <- 0
  main.result$Rsb.NR_V <- 0
  main.result$Rsb.NR_V.FDR <- 0
  main.result$Rsb.NR_L <- 0
  main.result$Rsb.NR_L.FDR <- 0
  main.result$Rsb.PF_V <- 0
  main.result$Rsb.PF_V.FDR <- 0
  main.result$Rsb.PF_L <- 0
  main.result$Rsb.PF_L.FDR <- 0
  main.result$Rsb.V_L <- 0
  main.result$Rsb.V_L.FDR <- 0
  #Iterate over scaffold positions
  for(i in 1:nrow(main.result)){
    #Get position
    testpos = main.result$POS[i]
    #Populate empty columns...
    #iHS (no FDR)
    if(nrow(iHS.cr.PF[which(iHS.cr.PF$POS == testpos), ]) >= 1){
      main.result$iHS.PF[i] <- iHS.cr.PF[which(iHS.cr.PF$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(iHS.cr.L[which(iHS.cr.L$POS == testpos), ]) >= 1){
      main.result$iHS.L[i] <- iHS.cr.L[which(iHS.cr.L$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(iHS.cr.OT[which(iHS.cr.OT$POS == testpos), ]) >= 1){
      main.result$iHS.OT[i] <- iHS.cr.OT[which(iHS.cr.OT$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(iHS.cr.V[which(iHS.cr.V$POS == testpos), ]) >= 1){
      main.result$iHS.V[i] <- iHS.cr.V[which(iHS.cr.V$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(iHS.cr.NR[which(iHS.cr.NR$POS == testpos), ]) >= 1){
      main.result$iHS.NR[i] <- iHS.cr.NR[which(iHS.cr.NR$POS == testpos), ]$N_EXTR_MRK
    }
    #iHS (with FDR)
    if(nrow(iHS.cr.PF.FDR[which(iHS.cr.PF.FDR$POS == testpos), ]) >= 1){
      main.result$iHS.PF.FDR[i] <- iHS.cr.PF.FDR[which(iHS.cr.PF.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(iHS.cr.L.FDR[which(iHS.cr.L.FDR$POS == testpos), ]) >= 1){
      main.result$iHS.L.FDR[i] <- iHS.cr.L.FDR[which(iHS.cr.L.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(iHS.cr.OT.FDR[which(iHS.cr.OT.FDR$POS == testpos), ]) >= 1){
      main.result$iHS.OT.FDR[i] <- iHS.cr.OT.FDR[which(iHS.cr.OT.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(iHS.cr.V.FDR[which(iHS.cr.V.FDR$POS == testpos), ]) >= 1){
      main.result$iHS.V.FDR[i] <- iHS.cr.V.FDR[which(iHS.cr.V.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(iHS.cr.NR.FDR[which(iHS.cr.NR.FDR$POS == testpos), ]) >= 1){
      main.result$iHS.NR.FDR[i] <- iHS.cr.NR.FDR[which(iHS.cr.NR.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    #xpEHH (no FDR)
    if(nrow(cr.xpEHH.OT_NR[which(cr.xpEHH.OT_NR$POS == testpos), ] >= 1)){
      main.result$xpEHH.OT_NR[i] <- cr.xpEHH.OT_NR[which(cr.xpEHH.OT_NR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.OT_PF[which(cr.xpEHH.OT_PF$POS == testpos), ] >= 1)){
      main.result$xpEHH.OT_PF[i] <- cr.xpEHH.OT_PF[which(cr.xpEHH.OT_PF$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.OT_V[which(cr.xpEHH.OT_V$POS == testpos), ] >= 1)){
      main.result$xpEHH.OT_V[i] <- cr.xpEHH.OT_V[which(cr.xpEHH.OT_V$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.OT_L[which(cr.xpEHH.OT_L$POS == testpos), ] >= 1)){
      main.result$xpEHH.OT_L[i] <- cr.xpEHH.OT_L[which(cr.xpEHH.OT_L$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.NR_PF[which(cr.xpEHH.NR_PF$POS == testpos), ] >= 1)){
      main.result$xpEHH.NR_PF[i] <- cr.xpEHH.NR_PF[which(cr.xpEHH.NR_PF$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.NR_V[which(cr.xpEHH.NR_V$POS == testpos), ] >= 1)){
      main.result$xpEHH.NR_V[i] <- cr.xpEHH.NR_V[which(cr.xpEHH.NR_V$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.NR_L[which(cr.xpEHH.NR_L$POS == testpos), ] >= 1)){
      main.result$xpEHH.NR_L[i] <- cr.xpEHH.NR_L[which(cr.xpEHH.NR_L$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.PF_V[which(cr.xpEHH.PF_V$POS == testpos), ] >= 1)){
      main.result$xpEHH.PF_V[i] <- cr.xpEHH.PF_V[which(cr.xpEHH.PF_V$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.PF_L[which(cr.xpEHH.PF_L$POS == testpos), ] >= 1)){
      main.result$xpEHH.PF_L[i] <- cr.xpEHH.PF_L[which(cr.xpEHH.PF_L$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.V_L[which(cr.xpEHH.V_L$POS == testpos), ] >= 1)){
      main.result$xpEHH.V_L[i] <- cr.xpEHH.V_L[which(cr.xpEHH.V_L$POS == testpos), ]$N_EXTR_MRK
    }
    #xpEHH (with FDR)
    if(nrow(cr.xpEHH.OT_NR.FDR[which(cr.xpEHH.OT_NR.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.OT_NR.FDR[i] <- cr.xpEHH.OT_NR.FDR[which(cr.xpEHH.OT_NR.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.OT_PF.FDR[which(cr.xpEHH.OT_PF.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.OT_PF.FDR[i] <- cr.xpEHH.OT_PF.FDR[which(cr.xpEHH.OT_PF.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.OT_V.FDR[which(cr.xpEHH.OT_V.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.OT_V.FDR[i] <- cr.xpEHH.OT_V.FDR[which(cr.xpEHH.OT_V.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.OT_L.FDR[which(cr.xpEHH.OT_L.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.OT_L.FDR[i] <- cr.xpEHH.OT_L.FDR[which(cr.xpEHH.OT_L.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.NR_PF.FDR[which(cr.xpEHH.NR_PF.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.NR_PF.FDR[i] <- cr.xpEHH.NR_PF.FDR[which(cr.xpEHH.NR_PF.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.NR_V.FDR[which(cr.xpEHH.NR_V.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.NR_V.FDR[i] <- cr.xpEHH.NR_V.FDR[which(cr.xpEHH.NR_V.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.NR_L.FDR[which(cr.xpEHH.NR_L.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.NR_L.FDR[i] <- cr.xpEHH.NR_L.FDR[which(cr.xpEHH.NR_L.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.PF_V.FDR[which(cr.xpEHH.PF_V.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.PF_V.FDR[i] <- cr.xpEHH.PF_V.FDR[which(cr.xpEHH.PF_V.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.PF_L.FDR[which(cr.xpEHH.PF_L.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.PF_L.FDR[i] <- cr.xpEHH.PF_L.FDR[which(cr.xpEHH.PF_L.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.xpEHH.V_L.FDR[which(cr.xpEHH.V_L.FDR$POS == testpos), ] >= 1)){
      main.result$xpEHH.V_L.FDR[i] <- cr.xpEHH.V_L.FDR[which(cr.xpEHH.V_L.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    #Rsb (no FDR)
    if(nrow(cr.Rsb.OT_NR[which(cr.Rsb.OT_NR$POS == testpos), ] >= 1)){
      main.result$Rsb.OT_NR[i] <- cr.Rsb.OT_NR[which(cr.Rsb.OT_NR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.OT_PF[which(cr.Rsb.OT_PF$POS == testpos), ] >= 1)){
      main.result$Rsb.OT_PF[i] <- cr.Rsb.OT_PF[which(cr.Rsb.OT_PF$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.OT_V[which(cr.Rsb.OT_V$POS == testpos), ] >= 1)){
      main.result$Rsb.OT_V[i] <- cr.Rsb.OT_V[which(cr.Rsb.OT_V$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.OT_L[which(cr.Rsb.OT_L$POS == testpos), ] >= 1)){
      main.result$Rsb.OT_L[i] <- cr.Rsb.OT_L[which(cr.Rsb.OT_L$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.NR_PF[which(cr.Rsb.NR_PF$POS == testpos), ] >= 1)){
      main.result$Rsb.NR_PF[i] <- cr.Rsb.NR_PF[which(cr.Rsb.NR_PF$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.NR_V[which(cr.Rsb.NR_V$POS == testpos), ] >= 1)){
      main.result$Rsb.NR_V[i] <- cr.Rsb.NR_V[which(cr.Rsb.NR_V$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.NR_L[which(cr.Rsb.NR_L$POS == testpos), ] >= 1)){
      main.result$Rsb.NR_L[i] <- cr.Rsb.NR_L[which(cr.Rsb.NR_L$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.PF_V[which(cr.Rsb.PF_V$POS == testpos), ] >= 1)){
      main.result$Rsb.PF_V[i] <- cr.Rsb.PF_V[which(cr.Rsb.PF_V$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.PF_L[which(cr.Rsb.PF_L$POS == testpos), ] >= 1)){
      main.result$Rsb.PF_L[i] <- cr.Rsb.PF_L[which(cr.Rsb.PF_L$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.V_L[which(cr.Rsb.V_L$POS == testpos), ] >= 1)){
      main.result$Rsb.V_L[i] <- cr.Rsb.V_L[which(cr.Rsb.V_L$POS == testpos), ]$N_EXTR_MRK
    }
    #Rsb (with FDR)
    if(nrow(cr.Rsb.OT_NR.FDR[which(cr.Rsb.OT_NR.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.OT_NR.FDR[i] <- cr.Rsb.OT_NR.FDR[which(cr.Rsb.OT_NR.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.OT_PF.FDR[which(cr.Rsb.OT_PF.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.OT_PF.FDR[i] <- cr.Rsb.OT_PF.FDR[which(cr.Rsb.OT_PF.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.OT_V.FDR[which(cr.Rsb.OT_V.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.OT_V.FDR[i] <- cr.Rsb.OT_V.FDR[which(cr.Rsb.OT_V.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.OT_L.FDR[which(cr.Rsb.OT_L.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.OT_L.FDR[i] <- cr.Rsb.OT_L.FDR[which(cr.Rsb.OT_L.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.NR_PF.FDR[which(cr.Rsb.NR_PF.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.NR_PF.FDR[i] <- cr.Rsb.NR_PF.FDR[which(cr.Rsb.NR_PF.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.NR_V.FDR[which(cr.Rsb.NR_V.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.NR_V.FDR[i] <- cr.Rsb.NR_V.FDR[which(cr.Rsb.NR_V.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.NR_L.FDR[which(cr.Rsb.NR_L.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.NR_L.FDR[i] <- cr.Rsb.NR_L.FDR[which(cr.Rsb.NR_L.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.PF_V.FDR[which(cr.Rsb.PF_V.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.PF_V.FDR[i] <- cr.Rsb.PF_V.FDR[which(cr.Rsb.PF_V.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.PF_L.FDR[which(cr.Rsb.PF_L.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.PF_L.FDR[i] <- cr.Rsb.PF_L.FDR[which(cr.Rsb.PF_L.FDR$POS == testpos), ]$N_EXTR_MRK
    }
    if(nrow(cr.Rsb.V_L.FDR[which(cr.Rsb.V_L.FDR$POS == testpos), ] >= 1)){
      main.result$Rsb.V_L.FDR[i] <- cr.Rsb.V_L.FDR[which(cr.Rsb.V_L.FDR$POS == testpos), ]$N_EXTR_MRK
    }
  }
  return(main.result)
}

#Run the function
OTref_all_candidate_regions <- amalgamate.cr(OTref_in_50kb, 
                                             cr.PF.OTref, cr.L.OTref, cr.OT.OTref, cr.V.OTref, cr.NR.OTref, 
                                             cr.PF.OTref.withFDR, cr.L.OTref.withFDR, cr.OT.OTref.withFDR, cr.V.OTref.withFDR, cr.NR_L.OTref.withFDR, 
                                             cr.OT_NR.OTref, 
                                             cr.OT_PF.OTref, 
                                             cr.OT_V.OTref, 
                                             cr.OT_L.OTref, 
                                             cr.NR_PF.OTref, 
                                             cr.NR_V.OTref, 
                                             cr.NR_L.OTref, 
                                             cr.PF_V.OTref, 
                                             cr.PF_L.OTref, 
                                             cr.V_L.OTref, 
                                             cr.OT_NR.OTref.withFDR, 
                                             cr.OT_PF.OTref.withFDR, 
                                             cr.OT_V.OTref.withFDR, 
                                             cr.OT_L.OTref.withFDR, 
                                             cr.NR_PF.OTref.withFDR, 
                                             cr.NR_V.OTref.withFDR, 
                                             cr.NR_L.OTref.withFDR, 
                                             cr.PF_V.OTref.withFDR, 
                                             cr.PF_L.OTref.withFDR, 
                                             cr.V_L.OTref.withFDR, 
                                             cr.OT_NR.Rsb.OTref, 
                                             cr.OT_PF.Rsb.OTref, 
                                             cr.OT_V.Rsb.OTref, 
                                             cr.OT_L.Rsb.OTref, 
                                             cr.NR_PF.Rsb.OTref, 
                                             cr.NR_V.Rsb.OTref, 
                                             cr.NR_L.Rsb.OTref, 
                                             cr.PF_V.Rsb.OTref, 
                                             cr.PF_L.Rsb.OTref, 
                                             cr.V_L.Rsb.OTref, 
                                             cr.OT_NR.Rsb.OTref.withFDR, 
                                             cr.OT_PF.Rsb.OTref.withFDR, 
                                             cr.OT_V.Rsb.OTref.withFDR, 
                                             cr.OT_L.Rsb.OTref.withFDR, 
                                             cr.NR_PF.Rsb.OTref.withFDR, 
                                             cr.NR_V.Rsb.OTref.withFDR, 
                                             cr.NR_L.Rsb.OTref.withFDR, 
                                             cr.PF_V.Rsb.OTref.withFDR, 
                                             cr.PF_L.Rsb.OTref.withFDR, 
                                             cr.V_L.Rsb.OTref.withFDR)
PFref_all_candidate_regions <- amalgamate.cr(PFref_in_50kb, 
                                             cr.PF.PFref, cr.L.PFref, cr.OT.PFref, cr.V.PFref, cr.NR.PFref, 
                                             cr.PF.PFref.withFDR, cr.L.PFref.withFDR, cr.OT.PFref.withFDR, cr.V.PFref.withFDR, cr.NR_L.PFref.withFDR, 
                                             cr.OT_NR.PFref, 
                                             cr.OT_PF.PFref, 
                                             cr.OT_V.PFref, 
                                             cr.OT_L.PFref, 
                                             cr.NR_PF.PFref, 
                                             cr.NR_V.PFref, 
                                             cr.NR_L.PFref, 
                                             cr.PF_V.PFref, 
                                             cr.PF_L.PFref, 
                                             cr.V_L.PFref, 
                                             cr.OT_NR.PFref.withFDR, 
                                             cr.OT_PF.PFref.withFDR, 
                                             cr.OT_V.PFref.withFDR, 
                                             cr.OT_L.PFref.withFDR, 
                                             cr.NR_PF.PFref.withFDR, 
                                             cr.NR_V.PFref.withFDR, 
                                             cr.NR_L.PFref.withFDR, 
                                             cr.PF_V.PFref.withFDR, 
                                             cr.PF_L.PFref.withFDR, 
                                             cr.V_L.PFref.withFDR, 
                                             cr.OT_NR.Rsb.PFref, 
                                             cr.OT_PF.Rsb.PFref, 
                                             cr.OT_V.Rsb.PFref, 
                                             cr.OT_L.Rsb.PFref, 
                                             cr.NR_PF.Rsb.PFref, 
                                             cr.NR_V.Rsb.PFref, 
                                             cr.NR_L.Rsb.PFref, 
                                             cr.PF_V.Rsb.PFref, 
                                             cr.PF_L.Rsb.PFref, 
                                             cr.V_L.Rsb.PFref, 
                                             cr.OT_NR.Rsb.PFref.withFDR, 
                                             cr.OT_PF.Rsb.PFref.withFDR, 
                                             cr.OT_V.Rsb.PFref.withFDR, 
                                             cr.OT_L.Rsb.PFref.withFDR, 
                                             cr.NR_PF.Rsb.PFref.withFDR, 
                                             cr.NR_V.Rsb.PFref.withFDR, 
                                             cr.NR_L.Rsb.PFref.withFDR, 
                                             cr.PF_V.Rsb.PFref.withFDR, 
                                             cr.PF_L.Rsb.PFref.withFDR, 
                                             cr.V_L.Rsb.PFref.withFDR)

#Add scaffold length column and re-order
OTref_all_candidate_regions <- OTref_all_candidate_regions %>% left_join(OTref_scaffold_lengths)
OTref_all_candidate_regions <- OTref_all_candidate_regions[c(1,2,3,55,seq(4,54))]
PFref_all_candidate_regions <- PFref_all_candidate_regions %>% left_join(PFref_scaffold_lengths)
PFref_all_candidate_regions <- PFref_all_candidate_regions[c(1,2,3,55,seq(4,54))]

#Write to file
write.table(OTref_all_candidate_regions, 
            file="OTref_all_candidate_regions.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_all_candidate_regions, 
            file="PFref_all_candidate_regions.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Automated EDA reports with dlookr
OTref_all_candidate_regions %>% 
  eda_report(output_format="html", 
             output_dir=getwd(), 
             output_file="OTref_all_candidate_regions_EDA_report.html")
PFref_all_candidate_regions %>% 
  eda_report(output_format="html", 
             output_dir=getwd(), 
             output_file="PFref_all_candidate_regions_EDA_report.html")

#Correlation matrices of candidate regions across different tests

#First just select the FDR-corrected data
OTref_corr_data <- OTref_all_candidate_regions %>% select(contains("FDR"))
PFref_corr_data <- PFref_all_candidate_regions %>% select(contains("FDR"))
#Create correlation plots
OTref_FDR_cr_corrplot <- ggcorrplot(cor(OTref_corr_data), tl.cex = 8, tl.srt = 90)
PFref_FDR_cr_corrplot <- ggcorrplot(cor(PFref_corr_data), tl.cex = 8, tl.srt = 90)
#Save
setwd("/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Correlation_plots")
ggsave("OTref_FDR_cr_corrplot.pdf", OTref_FDR_cr_corrplot, width=11, height=8.5)
ggsave("PFref_FDR_cr_corrplot.pdf", PFref_FDR_cr_corrplot, width=11, height=8.5)

#Plots showing ranked correlations between pairs of variables...
#All tests:
OTref_all_corr_cross <- corr_cross(OTref_corr_data)
PFref_all_corr_cross <- corr_cross(PFref_corr_data)
#xpEHH:
OTref_xpEHH_corr_cross <- corr_cross(OTref_corr_data %>% select(contains("xpEHH")))
PFref_xpEHH_corr_cross <- corr_cross(PFref_corr_data %>% select(contains("xpEHH")))
#Rsb:
OTref_Rsb_corr_cross <- corr_cross(OTref_corr_data %>% select(contains("Rsb")))
PFref_Rsb_corr_cross <- corr_cross(PFref_corr_data %>% select(contains("Rsb")))
#Save
ggsave("OTref_ranked_correlation_candidate_regions_all_tests.pdf", 
       OTref_all_corr_cross, width=11, height=8.5)
ggsave("PFref_ranked_correlation_candidate_regions_all_tests.pdf", 
       PFref_all_corr_cross, width=11, height=8.5)
ggsave("OTref_ranked_correlation_candidate_regions_xpEHH.pdf", 
       OTref_xpEHH_corr_cross, width=11, height=8.5)
ggsave("PFref_ranked_correlation_candidate_regions_xpEHH.pdf", 
       PFref_xpEHH_corr_cross, width=11, height=8.5)
ggsave("OTref_ranked_correlation_candidate_regions_Rsb.pdf", 
       OTref_Rsb_corr_cross, width=11, height=8.5)
ggsave("PFref_ranked_correlation_candidate_regions_Rsb.pdf", 
       PFref_Rsb_corr_cross, width=11, height=8.5)

#A function to group a master candidate regions dataframe 
#by scaffold and count the number of outlier markers
#for each test
scaffold_group_and_sum <- function(cr.df){
  #Select only FDR columns
  data = cr.df %>%  select(c(CHR, contains("FDR")))
  #Convert CHR column to factor, for grouping
  data$CHR <- as.factor(as.character(data$CHR))
  #Group by CHR and sum the number of outlier markers in each CHR
  data = data %>% group_by(CHR) %>% 
    summarise(sum.iHS.PF.FDR = sum(iHS.PF.FDR), 
              sum.iHS.L.FDR = sum(iHS.L.FDR), 
              sum.iHS.OT.FDR = sum(iHS.OT.FDR), 
              sum.iHS.V.FDR = sum(iHS.V.FDR), 
              sum.iHS.NR.FDR = sum(iHS.NR.FDR), 
              sum.xpEHH.OT_NR.FDR = sum(xpEHH.OT_NR.FDR), 
              sum.xpEHH.OT_PF.FDR = sum(xpEHH.OT_PF.FDR), 
              sum.xpEHH.OT_V.FDR = sum(xpEHH.OT_V.FDR), 
              sum.xpEHH.OT_L.FDR = sum(xpEHH.OT_L.FDR), 
              sum.xpEHH.NR_PF.FDR = sum(xpEHH.NR_PF.FDR), 
              sum.xpEHH.NR_V.FDR = sum(xpEHH.NR_V.FDR), 
              sum.xpEHH.NR_L.FDR = sum(xpEHH.NR_L.FDR), 
              sum.xpEHH.PF_V.FDR = sum(xpEHH.PF_V.FDR), 
              sum.xpEHH.PF_L.FDR = sum(xpEHH.PF_L.FDR), 
              sum.xpEHH.V_L.FDR = sum(xpEHH.V_L.FDR), 
              sum.Rsb.OT_NR.FDR = sum(Rsb.OT_NR.FDR), 
              sum.Rsb.OT_PF.FDR = sum(Rsb.OT_PF.FDR), 
              sum.Rsb.OT_V.FDR = sum(Rsb.OT_V.FDR), 
              sum.Rsb.OT_L.FDR = sum(Rsb.OT_L.FDR), 
              sum.Rsb.NR_PF.FDR = sum(Rsb.NR_PF.FDR), 
              sum.Rsb.NR_V.FDR = sum(Rsb.NR_V.FDR), 
              sum.Rsb.NR_L.FDR = sum(Rsb.NR_L.FDR), 
              sum.Rsb.PF_V.FDR = sum(Rsb.PF_V.FDR), 
              sum.Rsb.PF_L.FDR = sum(Rsb.PF_L.FDR), 
              sum.Rsb.V_L.FDR = sum(Rsb.V_L.FDR))
  #Create a temporary rowwise sum column for ordering
  data$sum <- apply(data[,c(2:26)], 1, sum)
  #Order by new sum column, and remove sum column
  data = data %>% arrange(desc(sum)) %>% select(-(sum))
  #Convert CHR column back to integer
  data$CHR <- as.integer(as.character(data$CHR))
  #Return result
  return(data)
}
#Call the function
OTref_all_candidate_regions_group_sum <- scaffold_group_and_sum(OTref_all_candidate_regions)
PFref_all_candidate_regions_group_sum <- scaffold_group_and_sum(PFref_all_candidate_regions)
#Write to file
write.table(OTref_all_candidate_regions_group_sum, file="OTref_candidate_scaffolds_extreme_marker_counts.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_all_candidate_regions_group_sum, file="PFref_candidate_scaffolds_extreme_marker_counts.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)


##################################################################################

##################### Visualisation - Manhattan plots ############################

#Clean up
rm(list = ls())

#Load libraries
library(rehh)
library(tidyverse)
library(rlist)
library(gridExtra)
library(ggpubr)

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
setwd(results_dir)

#Read data back in (only FDR-corrected data)

#iHS (OTref):
PF_iHS_OTref.FDR <- read.table(file="PF.ALL.OTref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
L_iHS_OTref.FDR <- read.table(file="L.ALL.OTref.iHS", 
                              header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_iHS_OTref.FDR <- read.table(file="OT.ALL.OTref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_iHS_OTref.FDR <- read.table(file="V.ALL.OTref.iHS", 
                              header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_iHS_OTref.FDR <- read.table(file="NR.ALL.OTref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#iHS (PFref):
PF_iHS_PFref.FDR <- read.table(file="PF.ALL.PFref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
L_iHS_PFref.FDR <- read.table(file="L.ALL.PFref.iHS", 
                              header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_iHS_PFref.FDR <- read.table(file="OT.ALL.PFref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_iHS_PFref.FDR <- read.table(file="V.ALL.PFref.iHS", 
                              header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_iHS_PFref.FDR <- read.table(file="NR.ALL.PFref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#xpEHH (OTref):
OT_NR_ALL_OTref.xpEHH.FDR <- read.table(file="OT_NR_ALL_OTref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_OTref.xpEHH.FDR <- read.table(file="OT_PF_ALL_OTref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_OTref.xpEHH.FDR <- read.table(file="OT_V_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_OTref.xpEHH.FDR <- read.table(file="OT_L_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_OTref.xpEHH.FDR <- read.table(file="NR_PF_ALL_OTref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_OTref.xpEHH.FDR <- read.table(file="NR_V_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_OTref.xpEHH.FDR <- read.table(file="NR_L_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_OTref.xpEHH.FDR <- read.table(file="PF_V_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_OTref.xpEHH.FDR <- read.table(file="PF_L_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_OTref.xpEHH.FDR <- read.table(file="V_L_ALL_OTref.xpEHH", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#xpEHH (PFref):
OT_NR_ALL_PFref.xpEHH.FDR <- read.table(file="OT_NR_ALL_PFref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_PFref.xpEHH.FDR <- read.table(file="OT_PF_ALL_PFref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_PFref.xpEHH.FDR <- read.table(file="OT_V_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_PFref.xpEHH.FDR <- read.table(file="OT_L_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_PFref.xpEHH.FDR <- read.table(file="NR_PF_ALL_PFref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_PFref.xpEHH.FDR <- read.table(file="NR_V_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_PFref.xpEHH.FDR <- read.table(file="NR_L_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_PFref.xpEHH.FDR <- read.table(file="PF_V_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_PFref.xpEHH.FDR <- read.table(file="PF_L_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_PFref.xpEHH.FDR <- read.table(file="V_L_ALL_PFref.xpEHH", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Rsb (OTref):
OT_NR_ALL_OTref.Rsb.FDR <- read.table(file="OT_NR_ALL_OTref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_OTref.Rsb.FDR <- read.table(file="OT_PF_ALL_OTref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_OTref.Rsb.FDR <- read.table(file="OT_V_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_OTref.Rsb.FDR <- read.table(file="OT_L_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_OTref.Rsb.FDR <- read.table(file="NR_PF_ALL_OTref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_OTref.Rsb.FDR <- read.table(file="NR_V_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_OTref.Rsb.FDR <- read.table(file="NR_L_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_OTref.Rsb.FDR <- read.table(file="PF_V_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_OTref.Rsb.FDR <- read.table(file="PF_L_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_OTref.Rsb.FDR <- read.table(file="V_L_ALL_OTref.Rsb", 
                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Rsb (PFref):
OT_NR_ALL_PFref.Rsb.FDR <- read.table(file="OT_NR_ALL_PFref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_PFref.Rsb.FDR <- read.table(file="OT_PF_ALL_PFref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_PFref.Rsb.FDR <- read.table(file="OT_V_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_PFref.Rsb.FDR <- read.table(file="OT_L_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_PFref.Rsb.FDR <- read.table(file="NR_PF_ALL_PFref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_PFref.Rsb.FDR <- read.table(file="NR_V_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_PFref.Rsb.FDR <- read.table(file="NR_L_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_PFref.Rsb.FDR <- read.table(file="PF_V_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_PFref.Rsb.FDR <- read.table(file="PF_L_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_PFref.Rsb.FDR <- read.table(file="V_L_ALL_PFref.Rsb", 
                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#Candidate regions master files
OTref_all_candidate_regions <- read.table(file="OTref_all_candidate_regions.tsv", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_all_candidate_regions <- read.table(file="PFref_all_candidate_regions.tsv", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#Scaffold length files
OTref_scaffold_lengths <- read.table(file="/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/Additional_files/OT3b_lengths.txt", 
                                     quote="", stringsAsFactors=FALSE)
colnames(OTref_scaffold_lengths) <- c("CHR", "LENGTH")
PFref_scaffold_lengths <- read.table(file="/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/Additional_files/PF_lengths.txt", 
                                     quote="", stringsAsFactors=FALSE)
colnames(PFref_scaffold_lengths) <- c("CHR", "LENGTH")




#Following https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

#Before starting, scaffold length data need to be incorporated into
#each data frame for plotting - these then need to be sorted by length
#with the longest scaffolds at the top, and then filtered to remove
#scaffolds <1Mb

#iHS function:
reduce.by.length.iHS <- function(df.iHS, df.length){
  #Change ihs.CHR to character
  df.iHS$ihs.CHR <- as.character(df.iHS$ihs.CHR)
  #Join, sort, and filter
  result = left_join(df.iHS, df.length, by=c("ihs.CHR" = "CHR")) %>% 
    arrange(desc(LENGTH)) %>% 
    filter(LENGTH >= 1000000)
  #Return the result
  return(result)
}
#cross-population version
reduce.by.length.xp <- function(df.xp, df.length){
  #Change CHR to character
  df.xp$CHR <- as.character(df.xp$CHR)
  #Join, sort and filter
  result = left_join(df.xp, df.length) %>% 
    arrange(desc(LENGTH)) %>% 
    filter(LENGTH >= 1000000)
  #Return the result
  return(result)
}

#Apply the functions...

#iHS (OTref):
PF_iHS_OTref.FDR <- reduce.by.length.iHS(PF_iHS_OTref.FDR, OTref_scaffold_lengths)
L_iHS_OTref.FDR <- reduce.by.length.iHS(L_iHS_OTref.FDR, OTref_scaffold_lengths)
OT_iHS_OTref.FDR <- reduce.by.length.iHS(OT_iHS_OTref.FDR, OTref_scaffold_lengths)
V_iHS_OTref.FDR <- reduce.by.length.iHS(V_iHS_OTref.FDR, OTref_scaffold_lengths)
NR_iHS_OTref.FDR <- reduce.by.length.iHS(NR_iHS_OTref.FDR, OTref_scaffold_lengths)
#iHS (PFref):
PF_iHS_PFref.FDR <- reduce.by.length.iHS(PF_iHS_PFref.FDR, PFref_scaffold_lengths)
L_iHS_PFref.FDR <- reduce.by.length.iHS(L_iHS_PFref.FDR, PFref_scaffold_lengths)
OT_iHS_PFref.FDR <- reduce.by.length.iHS(OT_iHS_PFref.FDR, PFref_scaffold_lengths)
V_iHS_PFref.FDR <- reduce.by.length.iHS(V_iHS_PFref.FDR, PFref_scaffold_lengths)
NR_iHS_PFref.FDR <- reduce.by.length.iHS(NR_iHS_PFref.FDR, PFref_scaffold_lengths)
#xpEHH (OTref):
OT_NR_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(OT_NR_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
OT_PF_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(OT_PF_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
OT_V_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(OT_V_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
OT_L_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(OT_L_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
NR_PF_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(NR_PF_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
NR_V_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(NR_V_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
NR_L_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(NR_L_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
PF_V_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(PF_V_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
PF_L_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(PF_L_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
V_L_ALL_OTref.xpEHH.FDR <- reduce.by.length.xp(V_L_ALL_OTref.xpEHH.FDR, OTref_scaffold_lengths)
#xpEHH (PFref):
OT_NR_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(OT_NR_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
OT_PF_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(OT_PF_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
OT_V_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(OT_V_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
OT_L_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(OT_L_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
NR_PF_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(NR_PF_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
NR_V_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(NR_V_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
NR_L_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(NR_L_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
PF_V_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(PF_V_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
PF_L_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(PF_L_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
V_L_ALL_PFref.xpEHH.FDR <- reduce.by.length.xp(V_L_ALL_PFref.xpEHH.FDR, PFref_scaffold_lengths)
#Rsb (OTref):
OT_NR_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(OT_NR_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
OT_PF_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(OT_PF_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
OT_V_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(OT_V_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
OT_L_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(OT_L_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
NR_PF_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(NR_PF_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
NR_V_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(NR_V_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
NR_L_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(NR_L_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
PF_V_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(PF_V_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
PF_L_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(PF_L_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
V_L_ALL_OTref.Rsb.FDR <- reduce.by.length.xp(V_L_ALL_OTref.Rsb.FDR, OTref_scaffold_lengths)
#Rsb (PFref):
OT_NR_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(OT_NR_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)
OT_PF_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(OT_PF_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)
OT_V_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(OT_V_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)
OT_L_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(OT_L_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)
NR_PF_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(NR_PF_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)
NR_V_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(NR_V_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)
NR_L_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(NR_L_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)
PF_V_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(PF_V_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)
PF_L_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(PF_L_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)
V_L_ALL_PFref.Rsb.FDR <- reduce.by.length.xp(V_L_ALL_PFref.Rsb.FDR, PFref_scaffold_lengths)

#To make the data easier to plot, the first step is to
#reduce the number of non-significant points to  plot.
#A function to accomplish this, using the threshold
#employed above (logP >= 6):
#cross-population version:
reduce.noise <- function(df){
  #Get the significant data
  sig.df = df %>% subset(LOGPVALUE >= 6)
  #Get the non-significant data and reduce
  notsig.df = df %>% subset(LOGPVALUE < 6) %>% 
    slice(sample(nrow(.), nrow(.) / 5))
  #Bind back together
  result = rbind(sig.df, notsig.df)
  #Sort by descending length and return
  result = result %>% arrange(desc(LENGTH), POSITION)
  return(result)
}
#iHS version:
reduce.noise.ihs <- function(df){
  #Get the significant data
  sig.df = df %>% subset(ihs.LOGPVALUE >= 6)
  #Get the non-significant data and reduce
  notsig.df = df %>% subset(ihs.LOGPVALUE < 6) %>% 
    slice(sample(nrow(.), nrow(.) / 5))
  #Bind back together
  result = rbind(sig.df, notsig.df)
  #Sort by descending length and return
  result = result %>% arrange(desc(LENGTH), ihs.POSITION)
  return(result)
}

#A further reduction in the amount of data plotted 
#is achieved by only plotting scaffolds which contain
#candidate selected regions.  
#A function to accomplish this
#(df.cr is a master candidate region data frame)...
#cross-population version:
get.cr.scaffolds <- function(df.logP, df.cr){
  result = df.logP[which(df.logP$CHR %in% unique(df.cr$CHR)), ]
  result = result %>% arrange(desc(LENGTH), POSITION)
  return(result)
}
#iHS version
get.cr.scaffolds.iHS <- function(df.logP, df.cr){
  result = df.logP[which(df.logP$ihs.CHR %in% unique(df.cr$CHR)), ]
  result = result %>% arrange(desc(LENGTH), ihs.POSITION)
  return(result)
}

#Next, a function add a column of cumulative position in bp...
#cross-population version:
get.bpcum <- function(gwas.dat){
  nCHR <- length(unique(gwas.dat$CHR))
  gwas.dat$BPcum <- 0
  s <- 0
  nbp <- c()
  gwas.dat <- gwas.dat[with(gwas.dat, order(CHR, POSITION)), ]
  for(i in 1:length(unique(gwas.dat$CHR))){
    target = unique(gwas.dat$CHR)[i]
    nbp[i] = max(gwas.dat[which(gwas.dat$CHR == target), ]$POSITION, na.rm=TRUE)
    gwas.dat[gwas.dat$CHR == target,"BPcum"] <- gwas.dat[gwas.dat$CHR == target,"POSITION"] + s
    s <- s + nbp[i]
  }
  return(gwas.dat)
}
#iHS version:
get.bpcum.iHS <- function(gwas.dat){
  nCHR <- length(unique(gwas.dat$ihs.CHR))
  gwas.dat$BPcum <- NA
  s <- 0
  nbp <- c()
  gwas.dat <- gwas.dat[with(gwas.dat, order(ihs.CHR, ihs.POSITION)), ]
  for(i in 1:length(unique(gwas.dat$ihs.CHR))){
    target = unique(gwas.dat$ihs.CHR)[i]
    nbp[i] = max(gwas.dat[which(gwas.dat$ihs.CHR == target), ]$ihs.POSITION, na.rm=TRUE)
    gwas.dat[gwas.dat$ihs.CHR == target,"BPcum"] <- gwas.dat[gwas.dat$ihs.CHR == target,"ihs.POSITION"] + s
    s <- s + nbp[i]
  }
  return(gwas.dat)
}

#Run these functions on all the data to be plotted.

#iHS...
#OTref:
PF.iHS.OTref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(PF_iHS_OTref.FDR, 
                                                                                 OTref_all_candidate_regions)))
L.iHS.OTref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(L_iHS_OTref.FDR, 
                                                                                OTref_all_candidate_regions)))
OT.iHS.OTref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(OT_iHS_OTref.FDR, 
                                                                                 OTref_all_candidate_regions)))
V.iHS.OTref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(V_iHS_OTref.FDR, 
                                                                                OTref_all_candidate_regions)))
NR.iHS.OTref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(NR_iHS_OTref.FDR, 
                                                                                 OTref_all_candidate_regions)))
#PFref:
PF.iHS.PFref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(PF_iHS_PFref.FDR, 
                                                                                 PFref_all_candidate_regions)))
L.iHS.PFref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(L_iHS_PFref.FDR, 
                                                                                PFref_all_candidate_regions)))
OT.iHS.PFref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(OT_iHS_PFref.FDR, 
                                                                                 PFref_all_candidate_regions)))
V.iHS.PFref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(V_iHS_PFref.FDR, 
                                                                                PFref_all_candidate_regions)))
NR.iHS.PFref_for_plotting <- get.bpcum.iHS(reduce.noise.ihs(get.cr.scaffolds.iHS(NR_iHS_PFref.FDR, 
                                                                                 PFref_all_candidate_regions)))
#xpEHH...
#OTref:
OT_NR.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_NR_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
OT_PF.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_PF_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
OT_V.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_V_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
OT_L.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_L_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
NR_PF.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_PF_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
NR_V.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_V_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
NR_L.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_L_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
PF_V.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(PF_V_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
PF_L.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(PF_L_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
V_L.xpEHH.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(V_L_ALL_OTref.xpEHH.FDR, OTref_all_candidate_regions)))
#PFref:
OT_NR.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_NR_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))
OT_PF.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_PF_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))
OT_V.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_V_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))
OT_L.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_L_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))
NR_PF.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_PF_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))
NR_V.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_V_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))
NR_L.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_L_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))
PF_V.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(PF_V_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))
PF_L.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(PF_L_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))
V_L.xpEHH.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(V_L_ALL_PFref.xpEHH.FDR, PFref_all_candidate_regions)))

#Rsb...
#OTref:
OT_NR.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_NR_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
OT_PF.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_PF_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
OT_V.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_V_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
OT_L.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_L_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
NR_PF.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_PF_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
NR_V.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_V_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
NR_L.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_L_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
PF_V.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(PF_V_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
PF_L.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(PF_L_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
V_L.Rsb.OTref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(V_L_ALL_OTref.Rsb.FDR, OTref_all_candidate_regions)))
#PFref:
OT_NR.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_NR_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))
OT_PF.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_PF_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))
OT_V.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_V_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))
OT_L.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(OT_L_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))
NR_PF.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_PF_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))
NR_V.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_V_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))
NR_L.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(NR_L_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))
PF_V.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(PF_V_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))
PF_L.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(PF_L_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))
V_L.Rsb.PFref_for_plotting <- get.bpcum(reduce.noise(get.cr.scaffolds(V_L_ALL_PFref.Rsb.FDR, PFref_all_candidate_regions)))

#Next, a central position for each scaffold in each dataframe needs to be calculated and stored.
#Fuctions for doing this...
#cross-population version:
get.center <- function(df){
  result <- df %>% group_by(CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  return(result)
}
#iHS version:
get.center.ihs <- function(df){
  result <- df %>% group_by(ihs.CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  return(result)
}

#Run this function...
#iHs - OTref:
L.iHS.OTref_centers <- get.center.ihs(L.iHS.OTref_for_plotting)
NR.iHS.OTref_centers <- get.center.ihs(NR.iHS.OTref_for_plotting)
OT.iHS.OTref_centers <- get.center.ihs(OT.iHS.OTref_for_plotting)
PF.iHS.OTref_centers <- get.center.ihs(PF.iHS.OTref_for_plotting)
V.iHS.OTref_centers <- get.center.ihs(V.iHS.OTref_for_plotting)
#iHs - PFref:
L.iHS.PFref_centers <- get.center.ihs(L.iHS.PFref_for_plotting)
NR.iHS.PFref_centers <- get.center.ihs(NR.iHS.PFref_for_plotting)
OT.iHS.PFref_centers <- get.center.ihs(OT.iHS.PFref_for_plotting)
PF.iHS.PFref_centers <- get.center.ihs(PF.iHS.PFref_for_plotting)
V.iHS.PFref_centers <- get.center.ihs(V.iHS.PFref_for_plotting)
#xpEHH - OTref:
OT_NR.xpEHH.OTref_centers <- get.center(OT_NR.xpEHH.OTref_for_plotting)
OT_PF.xpEHH.OTref_centers <- get.center(OT_PF.xpEHH.OTref_for_plotting)
OT_V.xpEHH.OTref_centers <- get.center(OT_V.xpEHH.OTref_for_plotting)
OT_L.xpEHH.OTref_centers <- get.center(OT_L.xpEHH.OTref_for_plotting)
NR_PF.xpEHH.OTref_centers <- get.center(NR_PF.xpEHH.OTref_for_plotting)
NR_V.xpEHH.OTref_centers <- get.center(NR_V.xpEHH.OTref_for_plotting)
NR_L.xpEHH.OTref_centers <- get.center(NR_L.xpEHH.OTref_for_plotting)
PF_V.xpEHH.OTref_centers <- get.center(PF_V.xpEHH.OTref_for_plotting)
PF_L.xpEHH.OTref_centers <- get.center(PF_L.xpEHH.OTref_for_plotting)
V_L.xpEHH.OTref_centers <- get.center(V_L.xpEHH.OTref_for_plotting)
#xpEHH - PFref:
OT_NR.xpEHH.PFref_centers <- get.center(OT_NR.xpEHH.PFref_for_plotting)
OT_PF.xpEHH.PFref_centers <- get.center(OT_PF.xpEHH.PFref_for_plotting)
OT_V.xpEHH.PFref_centers <- get.center(OT_V.xpEHH.PFref_for_plotting)
OT_L.xpEHH.PFref_centers <- get.center(OT_L.xpEHH.PFref_for_plotting)
NR_PF.xpEHH.PFref_centers <- get.center(NR_PF.xpEHH.PFref_for_plotting)
NR_V.xpEHH.PFref_centers <- get.center(NR_V.xpEHH.PFref_for_plotting)
NR_L.xpEHH.PFref_centers <- get.center(NR_L.xpEHH.PFref_for_plotting)
PF_V.xpEHH.PFref_centers <- get.center(PF_V.xpEHH.PFref_for_plotting)
PF_L.xpEHH.PFref_centers <- get.center(PF_L.xpEHH.PFref_for_plotting)
V_L.xpEHH.PFref_centers <- get.center(V_L.xpEHH.PFref_for_plotting)
#Rsb - OTref:
OT_NR.Rsb.OTref_centers <- get.center(OT_NR.Rsb.OTref_for_plotting)
OT_PF.Rsb.OTref_centers <- get.center(OT_PF.Rsb.OTref_for_plotting)
OT_V.Rsb.OTref_centers <- get.center(OT_V.Rsb.OTref_for_plotting)
OT_L.Rsb.OTref_centers <- get.center(OT_L.Rsb.OTref_for_plotting)
NR_PF.Rsb.OTref_centers <- get.center(NR_PF.Rsb.OTref_for_plotting)
NR_V.Rsb.OTref_centers <- get.center(NR_V.Rsb.OTref_for_plotting)
NR_L.Rsb.OTref_centers <- get.center(NR_L.Rsb.OTref_for_plotting)
PF_V.Rsb.OTref_centers <- get.center(PF_V.Rsb.OTref_for_plotting)
PF_L.Rsb.OTref_centers <- get.center(PF_L.Rsb.OTref_for_plotting)
V_L.Rsb.OTref_centers <- get.center(V_L.Rsb.OTref_for_plotting)
#Rsb - PFref:
OT_NR.Rsb.PFref_centers <- get.center(OT_NR.Rsb.PFref_for_plotting)
OT_PF.Rsb.PFref_centers <- get.center(OT_PF.Rsb.PFref_for_plotting)
OT_V.Rsb.PFref_centers <- get.center(OT_V.Rsb.PFref_for_plotting)
OT_L.Rsb.PFref_centers <- get.center(OT_L.Rsb.PFref_for_plotting)
NR_PF.Rsb.PFref_centers <- get.center(NR_PF.Rsb.PFref_for_plotting)
NR_V.Rsb.PFref_centers <- get.center(NR_V.Rsb.PFref_for_plotting)
NR_L.Rsb.PFref_centers <- get.center(NR_L.Rsb.PFref_for_plotting)
PF_V.Rsb.PFref_centers <- get.center(PF_V.Rsb.PFref_for_plotting)
PF_L.Rsb.PFref_centers <- get.center(PF_L.Rsb.PFref_for_plotting)
V_L.Rsb.PFref_centers <- get.center(V_L.Rsb.PFref_for_plotting)

#Final preparation for plotting - set a y-axis limit for 
#each statistic / reference assembly pair...
#iHS - OTref:
ylim.iHS.OTref <- max(ceiling(c(L.iHS.OTref_for_plotting$ihs.LOGPVALUE, 
                                NR.iHS.OTref_for_plotting$ihs.LOGPVALUE, 
                                OT.iHS.OTref_for_plotting$ihs.LOGPVALUE, 
                                PF.iHS.OTref_for_plotting$ihs.LOGPVALUE, 
                                V.iHS.OTref_for_plotting$ihs.LOGPVALUE))) + 2
#iHS - PFref:
ylim.iHS.PFref <- max(ceiling(c(L.iHS.PFref_for_plotting$ihs.LOGPVALUE, 
                                NR.iHS.PFref_for_plotting$ihs.LOGPVALUE, 
                                OT.iHS.PFref_for_plotting$ihs.LOGPVALUE, 
                                PF.iHS.PFref_for_plotting$ihs.LOGPVALUE, 
                                V.iHS.PFref_for_plotting$ihs.LOGPVALUE))) + 2
#xpEHH - OTref:
ylim.xpEHH.OTref <- max(ceiling(c(OT_NR.xpEHH.OTref_for_plotting$LOGPVALUE, 
                                  OT_PF.xpEHH.OTref_for_plotting$LOGPVALUE, 
                                  OT_V.xpEHH.OTref_for_plotting$LOGPVALUE, 
                                  OT_L.xpEHH.OTref_for_plotting$LOGPVALUE, 
                                  NR_PF.xpEHH.OTref_for_plotting$LOGPVALUE, 
                                  NR_V.xpEHH.OTref_for_plotting$LOGPVALUE, 
                                  NR_L.xpEHH.OTref_for_plotting$LOGPVALUE, 
                                  PF_V.xpEHH.OTref_for_plotting$LOGPVALUE, 
                                  PF_L.xpEHH.OTref_for_plotting$LOGPVALUE, 
                                  V_L.xpEHH.OTref_for_plotting$LOGPVALUE))) + 2
#xpEHH - PFref:
ylim.xpEHH.PFref <- max(ceiling(c(OT_NR.xpEHH.PFref_for_plotting$LOGPVALUE, 
                                  OT_PF.xpEHH.PFref_for_plotting$LOGPVALUE, 
                                  OT_V.xpEHH.PFref_for_plotting$LOGPVALUE, 
                                  OT_L.xpEHH.PFref_for_plotting$LOGPVALUE, 
                                  NR_PF.xpEHH.PFref_for_plotting$LOGPVALUE, 
                                  NR_V.xpEHH.PFref_for_plotting$LOGPVALUE, 
                                  NR_L.xpEHH.PFref_for_plotting$LOGPVALUE, 
                                  PF_V.xpEHH.PFref_for_plotting$LOGPVALUE, 
                                  PF_L.xpEHH.PFref_for_plotting$LOGPVALUE, 
                                  V_L.xpEHH.PFref_for_plotting$LOGPVALUE))) + 2
#Rsb - OTref:
ylim.Rsb.OTref <- max(ceiling(c(OT_NR.Rsb.OTref_for_plotting$LOGPVALUE, 
                                OT_PF.Rsb.OTref_for_plotting$LOGPVALUE, 
                                OT_V.Rsb.OTref_for_plotting$LOGPVALUE, 
                                OT_L.Rsb.OTref_for_plotting$LOGPVALUE, 
                                NR_PF.Rsb.OTref_for_plotting$LOGPVALUE, 
                                NR_V.Rsb.OTref_for_plotting$LOGPVALUE, 
                                NR_L.Rsb.OTref_for_plotting$LOGPVALUE, 
                                PF_V.Rsb.OTref_for_plotting$LOGPVALUE, 
                                PF_L.Rsb.OTref_for_plotting$LOGPVALUE, 
                                V_L.Rsb.OTref_for_plotting$LOGPVALUE))) + 2
#Rsb - PFref:
ylim.Rsb.PFref <- max(ceiling(c(OT_NR.Rsb.PFref_for_plotting$LOGPVALUE, 
                                OT_PF.Rsb.PFref_for_plotting$LOGPVALUE, 
                                OT_V.Rsb.PFref_for_plotting$LOGPVALUE, 
                                OT_L.Rsb.PFref_for_plotting$LOGPVALUE, 
                                NR_PF.Rsb.PFref_for_plotting$LOGPVALUE, 
                                NR_V.Rsb.PFref_for_plotting$LOGPVALUE, 
                                NR_L.Rsb.PFref_for_plotting$LOGPVALUE, 
                                PF_V.Rsb.PFref_for_plotting$LOGPVALUE, 
                                PF_L.Rsb.PFref_for_plotting$LOGPVALUE, 
                                V_L.Rsb.PFref_for_plotting$LOGPVALUE))) + 2

#And store the logpvalue significance threshold
sig <- 6

#Manhattan plotting function (cross-population version):
manhplot <- function(gwas.dat, axis.set, ylimit, title){
  p <- ggplot(gwas.dat, aes(x = BPcum, y = LOGPVALUE, 
                            color = as.factor(CHR), size = LOGPVALUE)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = sig, color = "grey40", linetype = "dashed") + 
    scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylimit)) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"), length(unique(gwas.dat$CHR)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log10(p)") + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5)
    ) + 
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

#Manhattan plotting function (iHS version):
manhplot.iHS <- function(gwas.dat, axis.set, ylimit, title){
  p <- ggplot(gwas.dat, aes(x = BPcum, y = ihs.LOGPVALUE, 
                            color = as.factor(ihs.CHR), size = ihs.LOGPVALUE)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = sig, color = "grey40", linetype = "dashed") + 
    scale_x_continuous(label = axis.set$ihs.CHR, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylimit)) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"), length(unique(gwas.dat$ihs.CHR)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log10(p)") + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5)
    ) + 
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

#Create individual plots...

#iHS (OTref):
p.OT.iHS.OTref <- manhplot.iHS(OT.iHS.OTref_for_plotting, OT.iHS.OTref_centers, ylim.iHS.OTref, "OT")
p.V.iHS.OTref <- manhplot.iHS(V.iHS.OTref_for_plotting, V.iHS.OTref_centers, ylim.iHS.OTref, "V")
p.NR.iHS.OTref <- manhplot.iHS(NR.iHS.OTref_for_plotting, NR.iHS.OTref_centers, ylim.iHS.OTref, "NR")
p.PF.iHS.OTref <- manhplot.iHS(PF.iHS.OTref_for_plotting, PF.iHS.OTref_centers, ylim.iHS.OTref, "PF")
p.L.iHS.OTref <- manhplot.iHS(L.iHS.OTref_for_plotting, L.iHS.OTref_centers, ylim.iHS.OTref, "L")
#iHS (PFref):
p.OT.iHS.PFref <- manhplot.iHS(OT.iHS.PFref_for_plotting, OT.iHS.PFref_centers, ylim.iHS.PFref, "OT")
p.V.iHS.PFref <- manhplot.iHS(V.iHS.PFref_for_plotting, V.iHS.PFref_centers, ylim.iHS.PFref, "V")
p.NR.iHS.PFref <- manhplot.iHS(NR.iHS.PFref_for_plotting, NR.iHS.PFref_centers, ylim.iHS.PFref, "NR")
p.PF.iHS.PFref <- manhplot.iHS(PF.iHS.PFref_for_plotting, PF.iHS.PFref_centers, ylim.iHS.PFref, "PF")
p.L.iHS.PFref <- manhplot.iHS(L.iHS.PFref_for_plotting, L.iHS.PFref_centers, ylim.iHS.PFref, "L")

#xpEHH (OTref):
p.OT_PF.xpEHH.OTref <- manhplot(OT_PF.xpEHH.OTref_for_plotting, OT_PF.xpEHH.OTref_centers, ylim.xpEHH.OTref, "OT_PF")
p.OT_L.xpEHH.OTref <- manhplot(OT_L.xpEHH.OTref_for_plotting, OT_L.xpEHH.OTref_centers, ylim.xpEHH.OTref, "OT_L")
p.PF_V.xpEHH.OTref <- manhplot(PF_V.xpEHH.OTref_for_plotting, PF_V.xpEHH.OTref_centers, ylim.xpEHH.OTref, "PF_V")
p.V_L.xpEHH.OTref <- manhplot(V_L.xpEHH.OTref_for_plotting, V_L.xpEHH.OTref_centers, ylim.xpEHH.OTref, "V_L")
p.OT_V.xpEHH.OTref <- manhplot(OT_V.xpEHH.OTref_for_plotting, OT_V.xpEHH.OTref_centers, ylim.xpEHH.OTref, "OT_V")
p.PF_L.xpEHH.OTref <- manhplot(PF_L.xpEHH.OTref_for_plotting, PF_L.xpEHH.OTref_centers, ylim.xpEHH.OTref, "PF_L")
p.OT_NR.xpEHH.OTref <- manhplot(OT_NR.xpEHH.OTref_for_plotting, OT_NR.xpEHH.OTref_centers, ylim.xpEHH.OTref, "OT_NR")
p.NR_V.xpEHH.OTref <- manhplot(NR_V.xpEHH.OTref_for_plotting, NR_V.xpEHH.OTref_centers, ylim.xpEHH.OTref, "NR_V")
p.NR_PF.xpEHH.OTref <- manhplot(NR_PF.xpEHH.OTref_for_plotting, NR_PF.xpEHH.OTref_centers, ylim.xpEHH.OTref, "NR_PF")
p.NR_L.xpEHH.OTref <- manhplot(NR_L.xpEHH.OTref_for_plotting, NR_L.xpEHH.OTref_centers, ylim.xpEHH.OTref, "NR_L")
#xpEHH (PFref):
p.OT_PF.xpEHH.PFref <- manhplot(OT_PF.xpEHH.PFref_for_plotting, OT_PF.xpEHH.PFref_centers, ylim.xpEHH.PFref, "OT_PF")
p.OT_L.xpEHH.PFref <- manhplot(OT_L.xpEHH.PFref_for_plotting, OT_L.xpEHH.PFref_centers, ylim.xpEHH.PFref, "OT_L")
p.PF_V.xpEHH.PFref <- manhplot(PF_V.xpEHH.PFref_for_plotting, PF_V.xpEHH.PFref_centers, ylim.xpEHH.PFref, "PF_V")
p.V_L.xpEHH.PFref <- manhplot(V_L.xpEHH.PFref_for_plotting, V_L.xpEHH.PFref_centers, ylim.xpEHH.PFref, "V_L")
p.OT_V.xpEHH.PFref <- manhplot(OT_V.xpEHH.PFref_for_plotting, OT_V.xpEHH.PFref_centers, ylim.xpEHH.PFref, "OT_V")
p.PF_L.xpEHH.PFref <- manhplot(PF_L.xpEHH.PFref_for_plotting, PF_L.xpEHH.PFref_centers, ylim.xpEHH.PFref, "PF_L")
p.OT_NR.xpEHH.PFref <- manhplot(OT_NR.xpEHH.PFref_for_plotting, OT_NR.xpEHH.PFref_centers, ylim.xpEHH.PFref, "OT_NR")
p.NR_V.xpEHH.PFref <- manhplot(NR_V.xpEHH.PFref_for_plotting, NR_V.xpEHH.PFref_centers, ylim.xpEHH.PFref, "NR_V")
p.NR_PF.xpEHH.PFref <- manhplot(NR_PF.xpEHH.PFref_for_plotting, NR_PF.xpEHH.PFref_centers, ylim.xpEHH.PFref, "NR_PF")
p.NR_L.xpEHH.PFref <- manhplot(NR_L.xpEHH.PFref_for_plotting, NR_L.xpEHH.PFref_centers, ylim.xpEHH.PFref, "NR_L")

#Rsb (OTref):
p.OT_PF.Rsb.OTref <- manhplot(OT_PF.Rsb.OTref_for_plotting, OT_PF.Rsb.OTref_centers, ylim.Rsb.OTref, "OT_PF")
p.OT_L.Rsb.OTref <- manhplot(OT_L.Rsb.OTref_for_plotting, OT_L.Rsb.OTref_centers, ylim.Rsb.OTref, "OT_L")
p.PF_V.Rsb.OTref <- manhplot(PF_V.Rsb.OTref_for_plotting, PF_V.Rsb.OTref_centers, ylim.Rsb.OTref, "PF_V")
p.V_L.Rsb.OTref <- manhplot(V_L.Rsb.OTref_for_plotting, V_L.Rsb.OTref_centers, ylim.Rsb.OTref, "V_L")
p.OT_V.Rsb.OTref <- manhplot(OT_V.Rsb.OTref_for_plotting, OT_V.Rsb.OTref_centers, ylim.Rsb.OTref, "OT_V")
p.PF_L.Rsb.OTref <- manhplot(PF_L.Rsb.OTref_for_plotting, PF_L.Rsb.OTref_centers, ylim.Rsb.OTref, "PF_L")
p.OT_NR.Rsb.OTref <- manhplot(OT_NR.Rsb.OTref_for_plotting, OT_NR.Rsb.OTref_centers, ylim.Rsb.OTref, "OT_NR")
p.NR_V.Rsb.OTref <- manhplot(NR_V.Rsb.OTref_for_plotting, NR_V.Rsb.OTref_centers, ylim.Rsb.OTref, "NR_V")
p.NR_PF.Rsb.OTref <- manhplot(NR_PF.Rsb.OTref_for_plotting, NR_PF.Rsb.OTref_centers, ylim.Rsb.OTref, "NR_PF")
p.NR_L.Rsb.OTref <- manhplot(NR_L.Rsb.OTref_for_plotting, NR_L.Rsb.OTref_centers, ylim.Rsb.OTref, "NR_L")
#Rsb (PFref):
p.OT_PF.Rsb.PFref <- manhplot(OT_PF.Rsb.PFref_for_plotting, OT_PF.Rsb.PFref_centers, ylim.Rsb.PFref, "OT_PF")
p.OT_L.Rsb.PFref <- manhplot(OT_L.Rsb.PFref_for_plotting, OT_L.Rsb.PFref_centers, ylim.Rsb.PFref, "OT_L")
p.PF_V.Rsb.PFref <- manhplot(PF_V.Rsb.PFref_for_plotting, PF_V.Rsb.PFref_centers, ylim.Rsb.PFref, "PF_V")
p.V_L.Rsb.PFref <- manhplot(V_L.Rsb.PFref_for_plotting, V_L.Rsb.PFref_centers, ylim.Rsb.PFref, "V_L")
p.OT_V.Rsb.PFref <- manhplot(OT_V.Rsb.PFref_for_plotting, OT_V.Rsb.PFref_centers, ylim.Rsb.PFref, "OT_V")
p.PF_L.Rsb.PFref <- manhplot(PF_L.Rsb.PFref_for_plotting, PF_L.Rsb.PFref_centers, ylim.Rsb.PFref, "PF_L")
p.OT_NR.Rsb.PFref <- manhplot(OT_NR.Rsb.PFref_for_plotting, OT_NR.Rsb.PFref_centers, ylim.Rsb.PFref, "OT_NR")
p.NR_V.Rsb.PFref <- manhplot(NR_V.Rsb.PFref_for_plotting, NR_V.Rsb.PFref_centers, ylim.Rsb.PFref, "NR_V")
p.NR_PF.Rsb.PFref <- manhplot(NR_PF.Rsb.PFref_for_plotting, NR_PF.Rsb.PFref_centers, ylim.Rsb.PFref, "NR_PF")
p.NR_L.Rsb.PFref <- manhplot(NR_L.Rsb.PFref_for_plotting, NR_L.Rsb.PFref_centers, ylim.Rsb.PFref, "NR_L")

#Arrange each set of plots in a list of pages...

#iHS (OTref):
pagelist.iHS.OTref <- vector("list", 1)
pagelist.iHS.OTref[[1]] <- ggarrange(p.OT.iHS.OTref, p.V.iHS.OTref, 
                                     ncol = 2, nrow = 1)
pagelist.iHS.OTref[[2]] <- ggarrange(p.PF.iHS.OTref, p.L.iHS.OTref, 
                                     ncol = 2, nrow = 1)
pagelist.iHS.OTref[[3]] <- ggarrange(p.NR.iHS.OTref, 
                                     ncol = 2, nrow = 1)
#iHS (PFref):
pagelist.iHS.PFref <- vector("list", 1)
pagelist.iHS.PFref[[1]] <- ggarrange(p.OT.iHS.PFref, p.V.iHS.PFref, 
                                     ncol = 2, nrow = 1)
pagelist.iHS.PFref[[2]] <- ggarrange(p.PF.iHS.PFref, p.L.iHS.PFref, 
                                     ncol = 2, nrow = 1)
pagelist.iHS.PFref[[3]] <- ggarrange(p.NR.iHS.PFref, 
                                     ncol = 2, nrow = 1)

#xpEHH (OTref):
pagelist.xpEHH.OTref <- vector("list", 5)
pagelist.xpEHH.OTref[[1]] <- ggarrange(p.OT_PF.xpEHH.OTref, p.OT_L.xpEHH.OTref, 
                                       ncol = 2, nrow = 1)
pagelist.xpEHH.OTref[[2]] <- ggarrange(p.PF_V.xpEHH.OTref, p.V_L.xpEHH.OTref, 
                                       ncol = 2, nrow = 1)
pagelist.xpEHH.OTref[[3]] <- ggarrange(p.OT_V.xpEHH.OTref, p.PF_L.xpEHH.OTref, 
                                       ncol = 2, nrow = 1)
pagelist.xpEHH.OTref[[4]] <- ggarrange(p.OT_NR.xpEHH.OTref, p.NR_V.xpEHH.OTref, 
                                       ncol = 2, nrow = 1)
pagelist.xpEHH.OTref[[5]] <- ggarrange(p.NR_PF.xpEHH.OTref, p.NR_L.xpEHH.OTref, 
                                       ncol = 2, nrow = 1)
#xpEHH (PFref):
pagelist.xpEHH.PFref <- vector("list", 5)
pagelist.xpEHH.PFref[[1]] <- ggarrange(p.OT_PF.xpEHH.PFref, p.OT_L.xpEHH.PFref, 
                                       ncol = 2, nrow = 1)
pagelist.xpEHH.PFref[[2]] <- ggarrange(p.PF_V.xpEHH.PFref, p.V_L.xpEHH.PFref, 
                                       ncol = 2, nrow = 1)
pagelist.xpEHH.PFref[[3]] <- ggarrange(p.OT_V.xpEHH.PFref, p.PF_L.xpEHH.PFref, 
                                       ncol = 2, nrow = 1)
pagelist.xpEHH.PFref[[4]] <- ggarrange(p.OT_NR.xpEHH.PFref, p.NR_V.xpEHH.PFref, 
                                       ncol = 2, nrow = 1)
pagelist.xpEHH.PFref[[5]] <- ggarrange(p.NR_PF.xpEHH.PFref, p.NR_L.xpEHH.PFref, 
                                       ncol = 2, nrow = 1)

#Rsb (OTref):
pagelist.Rsb.OTref <- vector("list", 5)
pagelist.Rsb.OTref[[1]] <- ggarrange(p.OT_PF.Rsb.OTref, p.OT_L.Rsb.OTref, 
                                     ncol = 2, nrow = 1)
pagelist.Rsb.OTref[[2]] <- ggarrange(p.PF_V.Rsb.OTref, p.V_L.Rsb.OTref, 
                                     ncol = 2, nrow = 1)
pagelist.Rsb.OTref[[3]] <- ggarrange(p.OT_V.Rsb.OTref, p.PF_L.Rsb.OTref, 
                                     ncol = 2, nrow = 1)
pagelist.Rsb.OTref[[4]] <- ggarrange(p.OT_NR.Rsb.OTref, p.NR_V.Rsb.OTref, 
                                     ncol = 2, nrow = 1)
pagelist.Rsb.OTref[[5]] <- ggarrange(p.NR_PF.Rsb.OTref, p.NR_L.Rsb.OTref, 
                                     ncol = 2, nrow = 1)
#Rsb (PFref):
pagelist.Rsb.PFref <- vector("list", 5)
pagelist.Rsb.PFref[[1]] <- ggarrange(p.OT_PF.Rsb.PFref, p.OT_L.Rsb.PFref, 
                                     ncol = 2, nrow = 1)
pagelist.Rsb.PFref[[2]] <- ggarrange(p.PF_V.Rsb.PFref, p.V_L.Rsb.PFref, 
                                     ncol = 2, nrow = 1)
pagelist.Rsb.PFref[[3]] <- ggarrange(p.OT_V.Rsb.PFref, p.PF_L.Rsb.PFref, 
                                     ncol = 2, nrow = 1)
pagelist.Rsb.PFref[[4]] <- ggarrange(p.OT_NR.Rsb.PFref, p.NR_V.Rsb.PFref, 
                                     ncol = 2, nrow = 1)
pagelist.Rsb.PFref[[5]] <- ggarrange(p.NR_PF.Rsb.PFref, p.NR_L.Rsb.PFref, 
                                     ncol = 2, nrow = 1)

#Save manhattan plots
manhplotdir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/ManhattanPlots"
setwd(manhplotdir)

#iHS (OTref):
ggsave("OTref_iHS_candidate_scaffolds_Manhattan.pdf", marrangeGrob(grobs=pagelist.iHS.OTref, nrow = 1, ncol = 1), width=11, height=8.5)
#iHS (PFref):
ggsave("PFref_iHS_candidate_scaffolds_Manhattan.pdf", marrangeGrob(grobs=pagelist.iHS.PFref, nrow = 1, ncol = 1), width=11, height=8.5)

#xpEHH (OTref):
ggsave("OTref_xpEHH_candidate_scaffolds_Manhattan.pdf", marrangeGrob(grobs=pagelist.xpEHH.OTref, nrow = 1, ncol = 1), width=11, height=8.5)
#xpEHH (PFref):
ggsave("PFref_xpEHH_candidate_scaffolds_Manhattan.pdf", marrangeGrob(grobs=pagelist.xpEHH.PFref, nrow = 1, ncol = 1), width=11, height=8.5)

#Rsb (OTref):
ggsave("OTref_Rsb_candidate_scaffolds_Manhattan.pdf", marrangeGrob(grobs=pagelist.Rsb.OTref, nrow = 1, ncol = 1), width=11, height=8.5)
#Rsb (PFref):
ggsave("PFref_Rsb_candidate_scaffolds_Manhattan.pdf", marrangeGrob(grobs=pagelist.Rsb.PFref, nrow = 1, ncol = 1), width=11, height=8.5)

##################################################################################

##################### Candidate Genes from Annotations ###########################

#Adapted from https://speciationgenomics.github.io/candidate_genes/

#Clean up
rm(list = ls())

#Load libraries
library(tidyverse)
library(ape)

#Set directories
gofeat_dir <- "/data/genomicsocorg/mwhj1/Data/Annotation/Supernova/GO_FEAT/"
OTgff_dir <- "/data/genomicsocorg/mwhj1/Data/Annotation/Supernova/OT3b_Results.dir/"
PFgff_dir <- "/data/genomicsocorg/mwhj1/Data/Annotation/Supernova/PF_Results.dir/"
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"

#Read annotation data
OT_gff <- read.gff(file=paste(OTgff_dir, "augustus.hints_genes.OT3b.gff3", sep=""))
PF_gff <- read.gff(file=paste(PFgff_dir, "augustus.hints_genes.PF.gff3", sep=""))
OT_gofeat <- read.delim(file=paste(gofeat_dir, "GO_FEAT_OT3b_all.csv", sep=""), 
                        sep=";", fill=TRUE, comment.char="")
PF_gofeat <- read.delim(file=paste(gofeat_dir, "GO_FEAT_PF_all.csv", sep=""), 
                        sep=";", fill=TRUE, comment.char="")

#Prepare gff data...
#Select genes only
OT_gff_genes <- OT_gff %>% filter(type == "gene")
PF_gff_genes <- PF_gff %>% filter(type == "gene")
#Reformat attribute column to give gene ID in gff_genes
OT_gff_genes$geneID <- sapply(strsplit(OT_gff_genes$attributes, split='=', fixed=TRUE), function(x) (x[2]))
OT_gff_genes$geneID <- sapply(strsplit(OT_gff_genes$geneID, split=';', fixed=TRUE), function(x) (x[1]))
PF_gff_genes$geneID <- sapply(strsplit(PF_gff_genes$attributes, split='=', fixed=TRUE), function(x) (x[2]))
PF_gff_genes$geneID <- sapply(strsplit(PF_gff_genes$geneID, split=';', fixed=TRUE), function(x) (x[1]))
#Reformat Locus.tag column in gofeat results to give gene ID
OT_gofeat$geneID <- sapply(strsplit(OT_gofeat$Locus.tag, split=".", fixed=TRUE), function(x) x[1])
PF_gofeat$geneID <- sapply(strsplit(PF_gofeat$Locus.tag, split=".", fixed=TRUE), function(x) x[1])
#Join gffs and gofeat results to get master annotation dataframes
OT_annotations <- left_join(OT_gff_genes, OT_gofeat, by="geneID") %>% select(-(contains("X")))
PF_annotations <- left_join(PF_gff_genes, PF_gofeat, by="geneID") %>% select(-(contains("X")))
#Create a mid-point column
OT_annotations <- OT_annotations %>% mutate(mid = start + (end-start)/2)
PF_annotations <- PF_annotations %>% mutate(mid = start + (end-start)/2)
#Remove weird length 0 annotations
OT_annotations <- OT_annotations %>% filter(Length > 0)
PF_annotations <- PF_annotations %>% filter(Length > 0)
#Convert seqids (scaffolds) to integers
OT_annotations$seqid <- as.integer(as.character(OT_annotations$seqid))
PF_annotations$seqid <- as.integer(as.character(PF_annotations$seqid))
#Write to file for posterity
write.table(OT_annotations, file="OT3b_gene_annotations_master.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PF_annotations, file="PF_gene_annotations_master.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Read data back in
setwd(results_dir)
#Candidate regions master files
OTref_all_candidate_regions <- read.table(file="OTref_all_candidate_regions.tsv", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_all_candidate_regions <- read.table(file="PFref_all_candidate_regions.tsv", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#Individual candidate regions results (only FDR-corrected data)..
options(scipen=999)
#iHS (with FDR)...
#OTref:
cr.L.OTref.withFDR <- read.table(file="cr.L.OTref.noOverlap.iHS", 
                                 quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR.OTref.withFDR <- read.table(file="cr.NR.OTref.noOverlap.iHS", 
                                  quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT.OTref.withFDR <- read.table(file="cr.OT.OTref.noOverlap.iHS", 
                                  quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF.OTref.withFDR <- read.table(file="cr.PF.OTref.noOverlap.iHS", 
                                  quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.V.OTref.withFDR <- read.table(file="cr.V.OTref.noOverlap.iHS", 
                                 quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
#PFref:
cr.L.PFref.withFDR <- read.table(file="cr.L.PFref.noOverlap.iHS", 
                                 quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR.PFref.withFDR <- read.table(file="cr.NR.PFref.noOverlap.iHS", 
                                  quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT.PFref.withFDR <- read.table(file="cr.OT.PFref.noOverlap.iHS", 
                                  quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF.PFref.withFDR <- read.table(file="cr.PF.PFref.noOverlap.iHS", 
                                  quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.V.PFref.withFDR <- read.table(file="cr.V.PFref.noOverlap.iHS", 
                                 quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)

#xpEHH (with FDR)...
#OTref:
cr.OT_NR.OTref.xpEHH.withFDR <- read.table(file="cr.OT_NR.OTref.noOverlap.xpEHH", 
                                           quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_PF.OTref.xpEHH.withFDR <- read.table(file="cr.OT_PF.OTref.noOverlap.xpEHH", 
                                           quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_V.OTref.xpEHH.withFDR <- read.table(file="cr.OT_V.OTref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_L.OTref.xpEHH.withFDR <- read.table(file="cr.OT_L.OTref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_PF.OTref.xpEHH.withFDR <- read.table(file="cr.NR_PF.OTref.noOverlap.xpEHH", 
                                           quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_V.OTref.xpEHH.withFDR <- read.table(file="cr.NR_V.OTref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_L.OTref.xpEHH.withFDR <- read.table(file="cr.NR_L.OTref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF_V.OTref.xpEHH.withFDR <- read.table(file="cr.PF_V.OTref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF_L.OTref.xpEHH.withFDR <- read.table(file="cr.PF_L.OTref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.V_L.OTref.xpEHH.withFDR <- read.table(file="cr.V_L.OTref.noOverlap.xpEHH", 
                                         quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
#PFref:
cr.OT_NR.PFref.xpEHH.withFDR <- read.table(file="cr.OT_NR.PFref.noOverlap.xpEHH", 
                                           quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_PF.PFref.xpEHH.withFDR <- read.table(file="cr.OT_PF.PFref.noOverlap.xpEHH", 
                                           quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_V.PFref.xpEHH.withFDR <- read.table(file="cr.OT_V.PFref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_L.PFref.xpEHH.withFDR <- read.table(file="cr.OT_L.PFref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_PF.PFref.xpEHH.withFDR <- read.table(file="cr.NR_PF.PFref.noOverlap.xpEHH", 
                                           quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_V.PFref.xpEHH.withFDR <- read.table(file="cr.NR_V.PFref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_L.PFref.xpEHH.withFDR <- read.table(file="cr.NR_L.PFref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF_V.PFref.xpEHH.withFDR <- read.table(file="cr.PF_V.PFref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF_L.PFref.xpEHH.withFDR <- read.table(file="cr.PF_L.PFref.noOverlap.xpEHH", 
                                          quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.V_L.PFref.xpEHH.withFDR <- read.table(file="cr.V_L.PFref.noOverlap.xpEHH", 
                                         quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)

#Rsb (with FDR)...
#Rsb (with FDR)...
#OTref:
cr.OT_NR.OTref.Rsb.withFDR <- read.table(file="cr.OT_NR.OTref.noOverlap.Rsb", 
                                         quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_PF.OTref.Rsb.withFDR <- read.table(file="cr.OT_PF.OTref.noOverlap.Rsb", 
                                         quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_V.OTref.Rsb.withFDR <- read.table(file="cr.OT_V.OTref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_L.OTref.Rsb.withFDR <- read.table(file="cr.OT_L.OTref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_PF.OTref.Rsb.withFDR <- read.table(file="cr.NR_PF.OTref.noOverlap.Rsb", 
                                         quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_V.OTref.Rsb.withFDR <- read.table(file="cr.NR_V.OTref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_L.OTref.Rsb.withFDR <- read.table(file="cr.NR_L.OTref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF_V.OTref.Rsb.withFDR <- read.table(file="cr.PF_V.OTref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF_L.OTref.Rsb.withFDR <- read.table(file="cr.PF_L.OTref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.V_L.OTref.Rsb.withFDR <- read.table(file="cr.V_L.OTref.noOverlap.Rsb", 
                                       quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
#PFref:
cr.OT_NR.PFref.Rsb.withFDR <- read.table(file="cr.OT_NR.PFref.noOverlap.Rsb", 
                                         quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_PF.PFref.Rsb.withFDR <- read.table(file="cr.OT_PF.PFref.noOverlap.Rsb", 
                                         quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_V.PFref.Rsb.withFDR <- read.table(file="cr.OT_V.PFref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.OT_L.PFref.Rsb.withFDR <- read.table(file="cr.OT_L.PFref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_PF.PFref.Rsb.withFDR <- read.table(file="cr.NR_PF.PFref.noOverlap.Rsb", 
                                         quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_V.PFref.Rsb.withFDR <- read.table(file="cr.NR_V.PFref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.NR_L.PFref.Rsb.withFDR <- read.table(file="cr.NR_L.PFref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF_V.PFref.Rsb.withFDR <- read.table(file="cr.PF_V.PFref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.PF_L.PFref.Rsb.withFDR <- read.table(file="cr.PF_L.PFref.noOverlap.Rsb", 
                                        quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
cr.V_L.PFref.Rsb.withFDR <- read.table(file="cr.V_L.PFref.noOverlap.Rsb", 
                                       quote="", sep="\t",stringsAsFactors=FALSE, fill=TRUE, header=TRUE)

#A function to find annotated genes in master candidate regions dataframes:
find.candidate.genes <- function(annotation.df, cr.df){
  #Create empty column in cr.df to hold gene IDs of any intersecting genes
  cr.df$geneIDs <- as.character("NA")
  #And columns for gene products and GO terms
  cr.df$gene.products <- as.character("NA")
  cr.df$gene.GO.terms <- as.character("NA")
  #And a mid-point column
  cr.df = cr.df %>% mutate(mid = START + (END-START)/2)
  #Iterate over cr.df rows
  for(i in 1:nrow(cr.df)){
    #Get mid-point of cr and scaffold
    midpoint = cr.df$mid[i]
    scaffold = cr.df$CHR[i]
    #Subset annotation to genes with mid-point within 50kb of selected mid-point
    annotation.subset.chr = annotation.df[which(annotation.df$seqid == scaffold), ]
    annotation.subset = annotation.subset.chr[which(abs(annotation.subset.chr$mid - midpoint) <= 50000), ]
    #If there are annotations found...
    if(nrow(annotation.subset) > 0){
      #Populate new columns
      cr.df$geneIDs[i] <- str_c(annotation.subset$geneID, collapse = ":")
      cr.df$gene.products[i] <- str_c(annotation.subset$Product, collapse=":")
      cr.df$gene.GO.terms[i] <- str_c(annotation.subset$Gene.onthology, collapse=":")
    }
  }
  #Return cr.df with correctly ordered columns
  result = cr.df %>% select(-mid)
  result = result[c(1:5, 56:58, 6:55)]
  return(result)
}

#Call the function...
OTref_candidate_regions_with_annotations <- find.candidate.genes(OT_annotations, OTref_all_candidate_regions)
PFref_candidate_regions_with_annotations <- find.candidate.genes(PF_annotations, PFref_all_candidate_regions)

#A function to find annotated genes in cr dataframes for individual tests:
find.candidate.genes.test <- function(annotation.df, cr.df){
  #Create empty column in cr.df to hold gene IDs of any intersecting genes
  cr.df$geneIDs <- as.character("NA")
  #And columns for gene products and GO terms
  cr.df$gene.products <- as.character("NA")
  cr.df$gene.GO.terms <- as.character("NA")
  #And a mid-point column
  cr.df = cr.df %>% mutate(mid = START + (END-START)/2)
  #Iterate over cr.df rows
  for(i in 1:nrow(cr.df)){
    #Get mid-point of cr and scaffold
    midpoint = cr.df$mid[i]
    scaffold = cr.df$CHR[i]
    #Subset annotation to genes with mid-point within 50kb of selected mid-point
    annotation.subset.chr = annotation.df[which(annotation.df$seqid == scaffold), ]
    annotation.subset = annotation.subset.chr[which(abs(annotation.subset.chr$mid - midpoint) <= 50000), ]
    #If there are annotations found...
    if(nrow(annotation.subset) > 0){
      #Populate new columns
      cr.df$geneIDs[i] <- str_c(annotation.subset$geneID, collapse = ":")
      cr.df$gene.products[i] <- str_c(annotation.subset$Product, collapse=":")
      cr.df$gene.GO.terms[i] <- str_c(annotation.subset$Gene.onthology, collapse=":")
    }
  }
  #Return cr.df without mid column
  result = cr.df %>% select(-mid)
  return(result)
}

#And on each FDR-corrected test-specific cr dataframe...
#iHS (OTref):
cr.L.OTref.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.L.OTref.withFDR)
cr.NR.OTref.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.NR.OTref.withFDR)
cr.OT.OTref.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.OT.OTref.withFDR)
cr.PF.OTref.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.PF.OTref.withFDR)
cr.V.OTref.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.V.OTref.withFDR)
#iHS (PFref):
cr.L.PFref.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.L.PFref.withFDR) #no FDR-corrected crs
cr.NR.PFref.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.NR.PFref.withFDR)
cr.OT.PFref.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.OT.PFref.withFDR) #no FDR-corrected crs
cr.PF.PFref.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.PF.PFref.withFDR) #no FDR-corrected crs
cr.V.PFref.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.V.PFref.withFDR)
#xpEHH (OTref):
cr.OT_NR.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.OT_NR.OTref.xpEHH.withFDR)
cr.OT_PF.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.OT_PF.OTref.xpEHH.withFDR)
cr.OT_V.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.OT_V.OTref.xpEHH.withFDR)
cr.OT_L.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.OT_L.OTref.xpEHH.withFDR)
cr.NR_PF.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.NR_PF.OTref.xpEHH.withFDR)
cr.NR_V.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.NR_V.OTref.xpEHH.withFDR)
cr.NR_L.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.NR_L.OTref.xpEHH.withFDR)
cr.PF_V.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.PF_V.OTref.xpEHH.withFDR)
cr.PF_L.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.PF_L.OTref.xpEHH.withFDR)
cr.V_L.OTref.xpEHH.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.V_L.OTref.xpEHH.withFDR)
#xpEHH (PFref):
cr.OT_NR.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.OT_NR.PFref.xpEHH.withFDR)
cr.OT_PF.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.OT_PF.PFref.xpEHH.withFDR) #no FDR-corrected crs
cr.OT_V.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.OT_V.PFref.xpEHH.withFDR)
cr.OT_L.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.OT_L.PFref.xpEHH.withFDR)
cr.NR_PF.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.NR_PF.PFref.xpEHH.withFDR)
cr.NR_V.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.NR_V.PFref.xpEHH.withFDR)
cr.NR_L.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.NR_L.PFref.xpEHH.withFDR)
cr.PF_V.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.PF_V.PFref.xpEHH.withFDR)
cr.PF_L.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.PF_L.PFref.xpEHH.withFDR)
cr.V_L.PFref.xpEHH.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.V_L.PFref.xpEHH.withFDR)
#Rsb (OTref):
cr.OT_NR.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.OT_NR.OTref.Rsb.withFDR)
cr.OT_PF.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.OT_PF.OTref.Rsb.withFDR)
cr.OT_V.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.OT_V.OTref.Rsb.withFDR)
cr.OT_L.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.OT_L.OTref.Rsb.withFDR)
cr.NR_PF.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.NR_PF.OTref.Rsb.withFDR)
cr.NR_V.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.NR_V.OTref.Rsb.withFDR)
cr.NR_L.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.NR_L.OTref.Rsb.withFDR)
cr.PF_V.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.PF_V.OTref.Rsb.withFDR)
cr.PF_L.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.PF_L.OTref.Rsb.withFDR)
cr.V_L.OTref.Rsb.withFDR.annotations <- find.candidate.genes.test(OT_annotations, cr.V_L.OTref.Rsb.withFDR)
#Rsb (PFref):
cr.OT_NR.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.OT_NR.PFref.Rsb.withFDR)
cr.OT_PF.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.OT_PF.PFref.Rsb.withFDR) #no FDR-corrected crs
cr.OT_V.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.OT_V.PFref.Rsb.withFDR)
cr.OT_L.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.OT_L.PFref.Rsb.withFDR)
cr.NR_PF.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.NR_PF.PFref.Rsb.withFDR) #no FDR-corrected crs
cr.NR_V.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.NR_V.PFref.Rsb.withFDR)
cr.NR_L.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.NR_L.PFref.Rsb.withFDR)
cr.PF_V.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.PF_V.PFref.Rsb.withFDR)
cr.PF_L.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.PF_L.PFref.Rsb.withFDR)
cr.V_L.PFref.Rsb.withFDR.annotations <- find.candidate.genes.test(PF_annotations, cr.V_L.PFref.Rsb.withFDR)

#Write to file...
#Master cr dataframes:
write.table(OTref_candidate_regions_with_annotations, file="OTref_all_candidate_regions_annotations.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_candidate_regions_with_annotations, file="PFref_all_candidate_regions_annotations.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
#iHS (OTref):
write.table(cr.L.OTref.withFDR.annotations, file="cr.L.OTref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR.OTref.withFDR.annotations, file="cr.NR.OTref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT.OTref.withFDR.annotations, file="cr.OT.OTref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF.OTref.withFDR.annotations, file="cr.PF.OTref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V.OTref.withFDR.annotations, file="cr.V.OTref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
#iHS (PFref):
write.table(cr.L.PFref.withFDR.annotations, file="cr.L.PFref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR.PFref.withFDR.annotations, file="cr.NR.PFref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT.PFref.withFDR.annotations, file="cr.OT.PFref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF.PFref.withFDR.annotations, file="cr.PF.PFref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V.PFref.withFDR.annotations, file="cr.V.PFref.withFDR.noOverlap.annotations.iHS", 
            quote=FALSE, sep="\t", row.names=FALSE)
#xpEHH (OTref):
write.table(cr.OT_NR.OTref.xpEHH.withFDR.annotations, file="cr.OT_NR.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.OTref.xpEHH.withFDR.annotations, file="cr.OT_PF.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.OTref.xpEHH.withFDR.annotations, file="cr.OT_V.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.OTref.xpEHH.withFDR.annotations, file="cr.OT_L.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.OTref.xpEHH.withFDR.annotations, file="cr.NR_PF.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.OTref.xpEHH.withFDR.annotations, file="cr.NR_V.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.OTref.xpEHH.withFDR.annotations, file="cr.NR_L.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.OTref.xpEHH.withFDR.annotations, file="cr.PF_V.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.OTref.xpEHH.withFDR.annotations, file="cr.PF_L.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.OTref.xpEHH.withFDR.annotations, file="cr.V_L.OTref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
#xpEHH (PFref):
write.table(cr.OT_NR.PFref.xpEHH.withFDR.annotations, file="cr.OT_NR.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.PFref.xpEHH.withFDR.annotations, file="cr.OT_PF.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.PFref.xpEHH.withFDR.annotations, file="cr.OT_V.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.PFref.xpEHH.withFDR.annotations, file="cr.OT_L.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.PFref.xpEHH.withFDR.annotations, file="cr.NR_PF.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.PFref.xpEHH.withFDR.annotations, file="cr.NR_V.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.PFref.xpEHH.withFDR.annotations, file="cr.NR_L.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.PFref.xpEHH.withFDR.annotations, file="cr.PF_V.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.PFref.xpEHH.withFDR.annotations, file="cr.PF_L.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.PFref.xpEHH.withFDR.annotations, file="cr.V_L.PFref.noOverlap.annotations.xpEHH", 
            quote=FALSE, sep="\t", row.names=FALSE)
#Rsb (OTref):
write.table(cr.OT_NR.OTref.Rsb.withFDR.annotations, file="cr.OT_NR.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.OTref.Rsb.withFDR.annotations, file="cr.OT_PF.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.OTref.Rsb.withFDR.annotations, file="cr.OT_V.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.OTref.Rsb.withFDR.annotations, file="cr.OT_L.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.OTref.Rsb.withFDR.annotations, file="cr.NR_PF.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.OTref.Rsb.withFDR.annotations, file="cr.NR_V.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.OTref.Rsb.withFDR.annotations, file="cr.NR_L.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.OTref.Rsb.withFDR.annotations, file="cr.PF_V.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.OTref.Rsb.withFDR.annotations, file="cr.PF_L.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.OTref.Rsb.withFDR.annotations, file="cr.V_L.OTref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
#Rsb (PFref):
write.table(cr.OT_NR.PFref.Rsb.withFDR.annotations, file="cr.OT_NR.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_PF.PFref.Rsb.withFDR.annotations, file="cr.OT_PF.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_V.PFref.Rsb.withFDR.annotations, file="cr.OT_V.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.OT_L.PFref.Rsb.withFDR.annotations, file="cr.OT_L.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_PF.PFref.Rsb.withFDR.annotations, file="cr.NR_PF.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_V.PFref.Rsb.withFDR.annotations, file="cr.NR_V.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.NR_L.PFref.Rsb.withFDR.annotations, file="cr.NR_L.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_V.PFref.Rsb.withFDR.annotations, file="cr.PF_V.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.PF_L.PFref.Rsb.withFDR.annotations, file="cr.PF_L.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(cr.V_L.PFref.Rsb.withFDR.annotations, file="cr.V_L.PFref.noOverlap.annotations.Rsb", 
            quote=FALSE, sep="\t", row.names=FALSE)

################################################################################

######################## Per-scaffold Line graphs ##############################
############### n extreme markers in 50kb sliding windows ######################

#Clean up
rm(list = ls())

#Load libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
setwd(results_dir)

#Read in data...
#Candidate regions master files:
OTref_all_candidate_regions <- read.table(file="OTref_all_candidate_regions.tsv", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_all_candidate_regions <- read.table(file="PFref_all_candidate_regions.tsv", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#Candidate regions per-scaffold counts
OTref_all_candidate_regions_group_sum <- read.table(file="OTref_candidate_scaffolds_extreme_marker_counts.tsv", 
                                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_all_candidate_regions_group_sum <- read.table(file="PFref_candidate_scaffolds_extreme_marker_counts.tsv", 
                                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Required data for plotting function:
FM_cols <- c("iHS.OT", "iHS.OT.FDR", "iHS.V", "iHS.V.FDR")
P_cols <- c("iHS.PF", "iHS.PF.FDR", "iHS.L", "iHS.L.FDR")
NR_cols <- c("iHS.NR", "iHS.NR.FDR")
btwn_cols <- c("xpEHH.OT_PF", "xpEHH.OT_PF.FDR","xpEHH.OT_L", "xpEHH.OT_L.FDR", "xpEHH.PF_V", "xpEHH.PF_V.FDR", "xpEHH.V_L", "xpEHH.V_L.FDR", 
               "Rsb.OT_PF", "Rsb.OT_PF.FDR","Rsb.OT_L", "Rsb.OT_L.FDR", "Rsb.PF_V", "Rsb.PF_V.FDR", "Rsb.V_L", "Rsb.V_L.FDR")
wthn_cols <- c("xpEHH.OT_V", "xpEHH.OT_V.FDR", "xpEHH.PF_L", "xpEHH.PF_L.FDR", 
               "Rsb.OT_V", "Rsb.OT_V.FDR", "Rsb.PF_L", "Rsb.PF_L.FDR")
unkwn_cols <- c("xpEHH.OT_NR", "xpEHH.OT_NR.FDR", "xpEHH.NR_PF", "xpEHH.NR_PF.FDR", "xpEHH.NR_V", "xpEHH.NR_V.FDR", "xpEHH.NR_L", "xpEHH.NR_L.FDR", 
                "Rsb.OT_NR", "Rsb.OT_NR.FDR", "Rsb.NR_PF", "Rsb.NR_PF.FDR", "Rsb.NR_V", "Rsb.NR_V.FDR", "Rsb.NR_L", "Rsb.NR_L.FDR")

#Subsetting functions:
scaffold.subset <- function(cr.df, target.scaffold){
  return(cr.df[which(cr.df$CHR == target.scaffold), ])
}
xpEHH.subset <- function(x){
  return(x %>% select(mid, contains("xpEHH")) %>% 
           gather(-mid, key = "stat", value = "value") %>% 
           mutate(
             type = case_when(
               stat %in% btwn_cols ~ "between", 
               stat %in% wthn_cols ~ "within", 
               stat %in% unkwn_cols ~ "unknown"
             )
           ))
}
Rsb.subset <- function(x){
  return(x %>% select(mid, contains("Rsb")) %>% 
           gather(-mid, key = "stat", value = "value") %>% 
           mutate(
             type = case_when(
               stat %in% btwn_cols ~ "between", 
               stat %in% wthn_cols ~ "within", 
               stat %in% unkwn_cols ~ "unknown"
             )
           ))
}
iHS.subset <- function(x){
  return(x %>% select(mid, contains("iHS")) %>% 
           gather(-mid, key = "stat", value = "value") %>% 
           mutate(
             type = case_when(
               stat %in% FM_cols ~ "FM", 
               stat %in% P_cols ~ "P", 
               stat %in% NR_cols ~ "NR"
             )
           ))
}



#Plotting function
#cr.df=all_candidate_regions; sum.df=candidate_regions_group_sum; 
#number=number of scaffolds to plot; test=test statistic ("iHS", "xpEHH", "Rsb)
#FDR=TRUE or FALSE (plot FDR-corrected or un-corrected data)
per_scaffold_extr_marker_plots <- function(cr.df, sum.df, number, FDR){
  #Create page list
  pagelist <- vector("list", number)
  #Create mid-point column 
  cr.df = cr.df %>% mutate(mid = START + (END-START)/2)
  #Subset according to whether FDR = TRUE or FALSE
  if(FDR == TRUE){
    cr.df = cr.df %>% select(CHR, 6:ncol(cr.df)) %>% 
      select(CHR, mid, contains("FDR"))
  }
  if(FDR == FALSE){
    cr.df = cr.df %>% select(CHR, 6:ncol(cr.df)) %>% 
      select(CHR, mid, !contains("FDR"))
  }
  #Iterate from 1 to number
  for(i in 1:number){
    #Get target scaffold
    scaffold = sum.df$CHR[i]
    #Subset cr.df to target scaffold
    cr.df.scaffold = scaffold.subset(cr.df, scaffold)
    #Get highest value to set a common y limit
    ymax = max(cr.df.scaffold[3:ncol(cr.df.scaffold)])+1
    #Create subsets
    iHS.data = iHS.subset(cr.df.scaffold)
    xpEHH.data = xpEHH.subset(cr.df.scaffold)
    Rsb.data = Rsb.subset(cr.df.scaffold)
    #Create individual plots...
    #xpEHH (between)
    p1 = ggplot(xpEHH.data %>% filter(type=="between"), aes(mid/10^6, value, color = stat)) + 
      geom_line() + 
      xlab("") + 
      ggtitle(as.character(scaffold)) + 
      scale_color_brewer(palette = "Paired") + 
      ylab("") + 
      ylim(c(0, ymax))  + 
      theme(legend.title=element_blank())
    #xpEHH (within)
    p2 = ggplot(xpEHH.data %>% filter(type=="within"), aes(mid/10^6, value, color = stat)) + 
      geom_line() + 
      xlab("") + 
      ggtitle("") + 
      scale_color_brewer(palette = "RdYlBu") + 
      ylab("") + 
      ylim(c(0, ymax))  + 
      theme(legend.title=element_blank())
    #xpEHH (unknown - NR)
    p3 = ggplot(xpEHH.data %>% filter(type=="unknown"), aes(mid/10^6, value, color = stat)) + 
      geom_line() + 
      xlab("") + 
      ggtitle("") + 
      ylab("") + 
      ylim(c(0, ymax))  + 
      theme(legend.title=element_blank())
    #Rsb (between)
    p4 = ggplot(Rsb.data %>% filter(type=="between"), aes(mid/10^6, value, color = stat)) + 
      geom_line() + 
      xlab("") + 
      ggtitle("") + 
      scale_color_brewer(palette = "Paired") + 
      ylab("") + 
      ylim(c(0, ymax))  + 
      theme(legend.title=element_blank())
    #Rsb (within)
    p5 = ggplot(Rsb.data %>% filter(type=="within"), aes(mid/10^6, value, color = stat)) + 
      geom_line() + 
      xlab("") + 
      ggtitle("") + 
      scale_color_brewer(palette = "RdYlBu") + 
      ylab("Outlier SNP count") + 
      ylim(c(0, ymax))  + 
      theme(legend.title=element_blank())
    #Rsb (unknown - NR)
    p6 = ggplot(Rsb.data %>% filter(type=="unknown"), aes(mid/10^6, value, color = stat)) + 
      geom_line() + 
      xlab("") + 
      ggtitle("") + 
      ylab("") + 
      ylim(c(0, ymax))  + 
      theme(legend.title=element_blank())
    #iHS (FM)
    p7 = ggplot(iHS.data %>% filter(type=="FM"), aes(mid/10^6, value, color = stat)) + 
      geom_line() + 
      xlab("") + 
      ggtitle("") + 
      scale_color_brewer(palette = "Paired") + 
      ylab("") + 
      ylim(c(0, ymax))  + 
      theme(legend.title=element_blank())
    #iHS (P)
    p8 = ggplot(iHS.data %>% filter(type=="P"), aes(mid/10^6, value, color = stat)) + 
      geom_line() + 
      xlab("") + 
      ggtitle("") + 
      scale_color_brewer(palette = "RdYlBu") + 
      ylab("") + 
      ylim(c(0, ymax))  + 
      theme(legend.title=element_blank())
    #iHS (NR)
    p9 = ggplot(iHS.data %>% filter(type=="NR"), aes(mid/10^6, value, color = stat)) + 
      geom_line() + 
      xlab("position (Mb)") + 
      ggtitle("") + 
      ylab("") + 
      ylim(c(0, ymax))  + 
      theme(legend.title=element_blank())
    #Arrange plots on page
    pagelist[[i]] <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, 
                               ncol = 1, nrow = 9, align = "v")
  }
  return(pagelist)
}

#Call function and save plots
plotdir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Outlier_scaffold_plots"
setwd(plotdir)
#OTref(FDR):
OTref_FDR_extr_mrkr_plots <- per_scaffold_extr_marker_plots(OTref_all_candidate_regions, OTref_all_candidate_regions_group_sum, nrow(OTref_all_candidate_regions_group_sum), FDR=TRUE)
#PFref(FDR):
PFref_FDR_extr_mrkr_plots <- per_scaffold_extr_marker_plots(PFref_all_candidate_regions, PFref_all_candidate_regions_group_sum, nrow(PFref_all_candidate_regions_group_sum), FDR=TRUE)
#OTref(no FDR):
OTref_No_FDR_extr_mrkr_plots <- per_scaffold_extr_marker_plots(OTref_all_candidate_regions, OTref_all_candidate_regions_group_sum, nrow(OTref_all_candidate_regions_group_sum), FDR=FALSE)
#PFref(no FDR):
PFref_No_FDR_extr_mrkr_plots <- per_scaffold_extr_marker_plots(PFref_all_candidate_regions, PFref_all_candidate_regions_group_sum, nrow(PFref_all_candidate_regions_group_sum), FDR=FALSE)

ggsave("OTref_FDR_extreme_marker_counts_plots.pdf", marrangeGrob(grobs=OTref_FDR_extr_mrkr_plots, nrow=1, ncol=1), width=8.5, height=11)
ggsave("PFref_FDR_extreme_marker_counts_plots.pdf", marrangeGrob(grobs=PFref_FDR_extr_mrkr_plots, nrow=1, ncol=1), width=8.5, height=11)
ggsave("OTref_No_FDR_extreme_marker_counts_plots.pdf", marrangeGrob(grobs=OTref_No_FDR_extr_mrkr_plots, nrow=1, ncol=1), width=8.5, height=11)
ggsave("PFref_No_FDR_extreme_marker_counts_plots.pdf", marrangeGrob(grobs=PFref_No_FDR_extr_mrkr_plots, nrow=1, ncol=1), width=8.5, height=11)


####################################################################################

######################## Identifying candidate regions #############################
################## which fit the expectations of divergence ########################
######################### between social phenotypes ################################

#Clean up
rm(list = ls())

#Load libraries
library(tidyverse)

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
setwd(results_dir)

#Read in data...
#Candidate regions master files (with annotations):
OTref_candidate_regions_with_annotations <- read.table(file="OTref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_candidate_regions_with_annotations <- read.table(file="PFref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Gene annotations
OT_annotations <- read.table(file="OT3b_gene_annotations_master.tsv", 
                             header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_annotations <- read.table(file="PF_gene_annotations_master.tsv", 
                             header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#Required data:
FM_cols <- c("iHS.OT", "iHS.OT.FDR", "iHS.V", "iHS.V.FDR")
P_cols <- c("iHS.PF", "iHS.PF.FDR", "iHS.L", "iHS.L.FDR")
NR_cols <- c("iHS.NR", "iHS.NR.FDR")
btwn_cols <- c("xpEHH.OT_PF", "xpEHH.OT_PF.FDR","xpEHH.OT_L", "xpEHH.OT_L.FDR", "xpEHH.PF_V", "xpEHH.PF_V.FDR", "xpEHH.V_L", "xpEHH.V_L.FDR", 
               "Rsb.OT_PF", "Rsb.OT_PF.FDR","Rsb.OT_L", "Rsb.OT_L.FDR", "Rsb.PF_V", "Rsb.PF_V.FDR", "Rsb.V_L", "Rsb.V_L.FDR")
wthn_cols <- c("xpEHH.OT_V", "xpEHH.OT_V.FDR", "xpEHH.PF_L", "xpEHH.PF_L.FDR", 
               "Rsb.OT_V", "Rsb.OT_V.FDR", "Rsb.PF_L", "Rsb.PF_L.FDR")
unkwn_cols <- c("xpEHH.OT_NR", "xpEHH.OT_NR.FDR", "xpEHH.NR_PF", "xpEHH.NR_PF.FDR", "xpEHH.NR_V", "xpEHH.NR_V.FDR", "xpEHH.NR_L", "xpEHH.NR_L.FDR", 
                "Rsb.OT_NR", "Rsb.OT_NR.FDR", "Rsb.NR_PF", "Rsb.NR_PF.FDR", "Rsb.NR_V", "Rsb.NR_V.FDR", "Rsb.NR_L", "Rsb.NR_L.FDR")

#Function called in main function:
check.conditions <- function(df, vec, index){
  if(ncol(df[which(df[,] >= 1)]) == ncol(df)){
    vec[index] <- 1
  }
  return(vec)
}

#Main function
soc.phen.cr.test <- function(cr.df){
  #Create vectors to hold some results
  is.cr.btwn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.Rsb <- rep(0, nrow(cr.df))
  #Iterate over rows in cr.df
  for(i in 1:nrow(cr.df)){
    #Create subsets
    subset.btwn.FDR <- cr.df[i,] %>% select(all_of(btwn_cols)) %>% select(contains("FDR"))
    subset.wthn.FDR <- cr.df[i,] %>% select(all_of(wthn_cols)) %>% select(contains("FDR"))
    subset.btwn.noFDR <- cr.df[i,] %>% select(all_of(btwn_cols)) %>% select(!contains("FDR"))
    subset.wthn.noFDR <- cr.df[i,] %>% select(all_of(wthn_cols)) %>% select(!contains("FDR"))
    subset.btwn.FDR.xpEHH <- subset.btwn.FDR %>% select(contains("xpEHH"))
    subset.wthn.FDR.xpEHH <- subset.wthn.FDR %>% select(contains("xpEHH"))
    subset.btwn.noFDR.xpEHH <- subset.btwn.noFDR %>% select(contains("xpEHH"))
    subset.wthn.noFDR.xpEHH <- subset.wthn.noFDR %>% select(contains("xpEHH"))
    subset.btwn.FDR.Rsb <- subset.btwn.FDR %>% select(contains("Rsb"))
    subset.wthn.FDR.Rsb <- subset.wthn.FDR %>% select(contains("Rsb"))
    subset.btwn.noFDR.Rsb <- subset.btwn.noFDR %>% select(contains("Rsb"))
    subset.wthn.noFDR.Rsb <- subset.wthn.noFDR %>% select(contains("Rsb"))
    #Change ith index of result vectors to 1 if conditions are satisfied
    is.cr.btwn.FDR.xpEHH <- check.conditions(subset.btwn.FDR.xpEHH, is.cr.btwn.FDR.xpEHH, i)
    is.cr.wthn.FDR.xpEHH <- check.conditions(subset.wthn.FDR.xpEHH, is.cr.wthn.FDR.xpEHH, i)
    is.cr.btwn.noFDR.xpEHH <- check.conditions(subset.btwn.noFDR.xpEHH, is.cr.btwn.noFDR.xpEHH, i)
    is.cr.wthn.noFDR.xpEHH <- check.conditions(subset.wthn.noFDR.xpEHH, is.cr.wthn.noFDR.xpEHH, i)
    is.cr.btwn.FDR.Rsb <- check.conditions(subset.btwn.FDR.Rsb, is.cr.btwn.FDR.Rsb, i)
    is.cr.wthn.FDR.Rsb <- check.conditions(subset.wthn.FDR.Rsb, is.cr.wthn.FDR.Rsb, i)
    is.cr.btwn.noFDR.Rsb <- check.conditions(subset.btwn.noFDR.Rsb, is.cr.btwn.noFDR.Rsb, i)
    is.cr.wthn.noFDR.Rsb <- check.conditions(subset.wthn.noFDR.Rsb, is.cr.wthn.noFDR.Rsb, i)
  }
  #Add result vectors as columns to cr.df
  cr.df$is.cr.btwn.FDR.xpEHH <- is.cr.btwn.FDR.xpEHH
  cr.df$is.cr.wthn.FDR.xpEHH <- is.cr.wthn.FDR.xpEHH
  cr.df$is.cr.btwn.noFDR.xpEHH <- is.cr.btwn.noFDR.xpEHH
  cr.df$is.cr.wthn.noFDR.xpEHH <- is.cr.wthn.noFDR.xpEHH
  cr.df$is.cr.btwn.FDR.Rsb <- is.cr.btwn.FDR.Rsb
  cr.df$is.cr.wthn.FDR.Rsb <- is.cr.wthn.FDR.Rsb
  cr.df$is.cr.btwn.noFDR.Rsb <- is.cr.btwn.noFDR.Rsb
  cr.df$is.cr.wthn.noFDR.Rsb <- is.cr.wthn.noFDR.Rsb
  #Create summary columns
  cr.df$is.soc.cr.FDR.xpEHH <- 0
  cr.df$is.soc.cr.noFDR.xpEHH <- 0
  cr.df$is.soc.cr.FDR.Rsb <- 0
  cr.df$is.soc.cr.noFDR.Rsb <- 0
  #Assign 1 to summary columns when no cr detected for within-phentoype comparisons
  #and a cr is detected for between-phenotype comparisons
  cr.df$is.soc.cr.FDR.xpEHH[which(cr.df$is.cr.wthn.FDR.xpEHH == 0 & cr.df$is.cr.btwn.FDR.xpEHH == 1)] <- 1
  cr.df$is.soc.cr.noFDR.xpEHH[which(cr.df$is.cr.wthn.noFDR.xpEHH == 0 & cr.df$is.cr.btwn.noFDR.xpEHH == 1)] <- 1
  cr.df$is.soc.cr.FDR.Rsb[which(cr.df$is.cr.wthn.FDR.Rsb == 0 & cr.df$is.cr.btwn.FDR.Rsb == 1)] <- 1
  cr.df$is.soc.cr.noFDR.Rsb[which(cr.df$is.cr.wthn.noFDR.Rsb == 0 & cr.df$is.cr.btwn.noFDR.Rsb == 1)] <- 1
  #Sort by length
  cr.df = cr.df %>% arrange(desc(LENGTH), START)
  #Return result
  return(cr.df)
}

#Call the function
OTref_social_candidate_regions <- soc.phen.cr.test(OTref_candidate_regions_with_annotations)
PFref_social_candidate_regions <- soc.phen.cr.test(PFref_candidate_regions_with_annotations)

#Check results
nrow(OTref_social_candidate_regions[which(OTref_social_candidate_regions$is.soc.cr.FDR.xpEHH == 1), ])
nrow(OTref_social_candidate_regions[which(OTref_social_candidate_regions$is.soc.cr.noFDR.xpEHH == 1), ])
nrow(OTref_social_candidate_regions[which(OTref_social_candidate_regions$is.soc.cr.FDR.Rsb == 1), ])
nrow(OTref_social_candidate_regions[which(OTref_social_candidate_regions$is.soc.cr.noFDR.Rsb == 1), ])
nrow(OTref_social_candidate_regions[which(OTref_social_candidate_regions$is.soc.cr.noFDR.xpEHH == 1 & OTref_social_candidate_regions$is.soc.cr.noFDR.Rsb == 1), ])
nrow(PFref_social_candidate_regions[which(PFref_social_candidate_regions$is.soc.cr.FDR.xpEHH == 1), ])
nrow(PFref_social_candidate_regions[which(PFref_social_candidate_regions$is.soc.cr.noFDR.xpEHH == 1), ])
nrow(PFref_social_candidate_regions[which(PFref_social_candidate_regions$is.soc.cr.FDR.Rsb == 1), ])
nrow(PFref_social_candidate_regions[which(PFref_social_candidate_regions$is.soc.cr.noFDR.Rsb == 1), ])
nrow(PFref_social_candidate_regions[which(PFref_social_candidate_regions$is.soc.cr.noFDR.xpEHH == 1 & PFref_social_candidate_regions$is.soc.cr.noFDR.Rsb == 1), ])

#Create subsets...
#OTref:
OTref_social_candidate_regions.FDR.xpEHH <- OTref_social_candidate_regions[which(OTref_social_candidate_regions$is.soc.cr.FDR.xpEHH == 1), ]
OTref_social_candidate_regions.noFDR.xpEHH <- OTref_social_candidate_regions[which(OTref_social_candidate_regions$is.soc.cr.noFDR.xpEHH == 1), ]
OTref_social_candidate_regions.FDR.Rsb <- OTref_social_candidate_regions[which(OTref_social_candidate_regions$is.soc.cr.FDR.Rsb == 1), ]
OTref_social_candidate_regions.noFDR.Rsb <- OTref_social_candidate_regions[which(OTref_social_candidate_regions$is.soc.cr.noFDR.Rsb == 1), ]
#PFref:
PFref_social_candidate_regions.FDR.xpEHH <- PFref_social_candidate_regions[which(PFref_social_candidate_regions$is.soc.cr.FDR.xpEHH == 1), ]
PFref_social_candidate_regions.noFDR.xpEHH <- PFref_social_candidate_regions[which(PFref_social_candidate_regions$is.soc.cr.noFDR.xpEHH == 1), ]
PFref_social_candidate_regions.FDR.Rsb <- PFref_social_candidate_regions[which(PFref_social_candidate_regions$is.soc.cr.FDR.Rsb == 1), ]
PFref_social_candidate_regions.noFDR.Rsb <- PFref_social_candidate_regions[which(PFref_social_candidate_regions$is.soc.cr.noFDR.Rsb == 1), ]

#Bind together...
#OTref:
OTref_social_candidate_regions_master <- data.table::rbindlist(list(OTref_social_candidate_regions.FDR.xpEHH, 
                                                                    OTref_social_candidate_regions.noFDR.xpEHH, 
                                                                    OTref_social_candidate_regions.FDR.Rsb, 
                                                                    OTref_social_candidate_regions.noFDR.Rsb)) %>% 
  distinct()
#PFref:
PFref_social_candidate_regions_master <- data.table::rbindlist(list(PFref_social_candidate_regions.FDR.xpEHH, 
                                                                    PFref_social_candidate_regions.noFDR.xpEHH, 
                                                                    PFref_social_candidate_regions.FDR.Rsb, 
                                                                    PFref_social_candidate_regions.noFDR.Rsb)) %>% 
  distinct()

#Write to file...
#Master sets:
write.table(OTref_social_candidate_regions_master, 
            file="OTref_social_candidate_regions_master.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_social_candidate_regions_master, 
            file="PFref_social_candidate_regions_master.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
#Subsets...
#OTref:
write.table(OTref_social_candidate_regions.FDR.xpEHH, 
            file="OTref_social_candidate_regions.FDR.xpEHH.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OTref_social_candidate_regions.noFDR.xpEHH, 
            file="OTref_social_candidate_regions.noFDR.xpEHH.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OTref_social_candidate_regions.FDR.Rsb, 
            file="OTref_social_candidate_regions.FDR.Rsb.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OTref_social_candidate_regions.noFDR.Rsb, 
            file="OTref_social_candidate_regions.noFDR.Rsb.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
#PFref:
write.table(PFref_social_candidate_regions.FDR.xpEHH, 
            file="PFref_social_candidate_regions.FDR.xpEHH.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_social_candidate_regions.noFDR.xpEHH, 
            file="PFref_social_candidate_regions.noFDR.xpEHH.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_social_candidate_regions.FDR.Rsb, 
            file="PFref_social_candidate_regions.FDR.Rsb.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_social_candidate_regions.noFDR.Rsb, 
            file="PFref_social_candidate_regions.noFDR.Rsb.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)

#A function to subset gene annotation files to genes in social candidate regions data frames
get.candidate.genes <- function(cr.df, ann.df){
  #First step is to get the gene IDs from cr.df and format them as a unique list...
  #Get gene ID column
  geneID.column = cr.df$geneIDs
  #Create emptyp vector
  gene.ID.vector = c()
  #Iterate over gene ID column entries
  for(i in 1:length(geneID.column)){
    #Split gene IDs (separated by ":")
    x = strsplit(geneID.column[i], split=":", fixed=TRUE)
    #Add to gene.ID.vector
    gene.ID.vector = c(gene.ID.vector, x)
  }
  #Reformat as vector, remove NAs, and make unique
  gene.ID.vector = unlist(gene.ID.vector)
  gene.ID.vector = gene.ID.vector[!is.na(gene.ID.vector)]
  gene.ID.vector = unique(gene.ID.vector)
  #Create an empty list of the same length as gene.ID.vector
  #(this will be populated by rows from the gene annotation data frame - to be bound later)
  annotation.list <- vector("list", length(gene.ID.vector))
  #Iterate over unique gene IDs
  for(j in 1:length(gene.ID.vector)){
    #Get jth gene ID
    ID = gene.ID.vector[j]
    #Make jth index of annotation.list the corresponding subset of the annotation data frame
    annotation.list[[j]] = ann.df %>% filter(geneID == ID)
  }
  #Create result
  result = do.call(rbind, annotation.list) %>% distinct()
  #Return result
  return(result)
}

#Call the function...
#Master versions (all social candidate regions detected by any test, FDR-corrected and -uncorrected)
OTref_social_candidate_genes <- get.candidate.genes(OTref_social_candidate_regions_master, 
                                                    OT_annotations)
PFref_social_candidate_genes <- get.candidate.genes(PFref_social_candidate_regions_master, 
                                                    PF_annotations)
#Test-specific versions...
#xpEHH (FDR-corrected)
OTref_social_candidate_genes.FDR.xpEHH <- get.candidate.genes(OTref_social_candidate_regions.FDR.xpEHH, 
                                                              OT_annotations)
#(nothing detected for PF FDR-corrected xpEHH)
#xpEHH (not FDR-corrected)
OTref_social_candidate_genes.noFDR.xpEHH <- get.candidate.genes(OTref_social_candidate_regions.noFDR.xpEHH, 
                                                                OT_annotations)
PFref_social_candidate_genes.noFDR.xpEHH <- get.candidate.genes(PFref_social_candidate_regions.noFDR.xpEHH, 
                                                                PF_annotations)
#Rsb (FDR-corrected)
OTref_social_candidate_genes.FDR.Rsb <- get.candidate.genes(OTref_social_candidate_regions.FDR.Rsb, 
                                                            OT_annotations)
#(nothing detected for PF FDR-corrected Rsb)
#Rsb (not FDR-corrected)
OTref_social_candidate_genes.noFDR.Rsb <- get.candidate.genes(OTref_social_candidate_regions.noFDR.Rsb, 
                                                              OT_annotations)
PFref_social_candidate_genes.noFDR.Rsb <- get.candidate.genes(PFref_social_candidate_regions.noFDR.Rsb, 
                                                              PF_annotations)

#Check out the number of shared results in the master sets (any test, corrected and un-corrected)
#(Subtract 1 from these because of "Uncharacterized protein")
shared_social_candidate_genes <- intersect(OTref_social_candidate_genes$Product, 
                                           PFref_social_candidate_genes$Product)
length(shared_social_candidate_genes) - 1

#Write to file...
#Master candidate genes (any test, corrected and un-corrected):
write.table(OTref_social_candidate_genes, 
            file="OT3b_social_candidate_genes_with_midpoint_50kb_or_less_from_midpoint_of_50kb_cr_window.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_social_candidate_genes, 
            file="PF_social_candidate_genes_with_midpoint_50kb_or_less_from_midpoint_of_50kb_cr_window.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
#xpEHH:
write.table(OTref_social_candidate_genes.FDR.xpEHH, 
            file="OT3b_social_candidate_genes_with_midpoint_50kb_or_less_from_midpoint_of_50kb_cr_window.xpEHH.FDR.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OTref_social_candidate_genes.noFDR.xpEHH, 
            file="OT3b_social_candidate_genes_with_midpoint_50kb_or_less_from_midpoint_of_50kb_cr_window.xpEHH.noFDR.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_social_candidate_genes.noFDR.xpEHH, 
            file="PF_social_candidate_genes_with_midpoint_50kb_or_less_from_midpoint_of_50kb_cr_window.xpEHH.noFDR.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
#Rsb:
write.table(OTref_social_candidate_genes.FDR.Rsb, 
            file="OT3b_social_candidate_genes_with_midpoint_50kb_or_less_from_midpoint_of_50kb_cr_window.Rsb.FDR.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(OTref_social_candidate_genes.noFDR.Rsb, 
            file="OT3b_social_candidate_genes_with_midpoint_50kb_or_less_from_midpoint_of_50kb_cr_window.Rsb.noFDR.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_social_candidate_genes.noFDR.Rsb, 
            file="PF_social_candidate_genes_with_midpoint_50kb_or_less_from_midpoint_of_50kb_cr_window.Rsb.noFDR.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)



####################################################################################

################################ Summary tables ####################################

#Clean up
rm(list = ls())

#Load libraries
library(tidyverse)

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
setwd(results_dir)

#Read data back in...
#Social candidate regions (any test, corrected and un-corrected):
OTref_social_candidate_regions_master <- read.table(file="OTref_social_candidate_regions_master.tsv", 
                                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_social_candidate_regions_master <- read.table(file="PFref_social_candidate_regions_master.tsv", 
                                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Candidate regions master files (with annotations):
OTref_candidate_regions_with_annotations <- read.table(file="OTref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_candidate_regions_with_annotations <- read.table(file="PFref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#iHS summary table...
#Create iHS cr subsets:
OTref_iHS_crs <- OTref_candidate_regions_with_annotations %>% select(CHR, START, END, LENGTH, contains("iHS"))
PFref_iHS_crs <- PFref_candidate_regions_with_annotations %>% select(CHR, START, END, LENGTH, contains("iHS"))
#Create summary table vectors:
iHS.summary.pop <- rep(c("PF", "L", "OT", "V", "NR"), each=2)
iHS.summary.ref <- rep(c("OT", "PF"), 5)
iHS.n.outliers.noFDR <- as.character(c(sum(OTref_iHS_crs$iHS.PF), sum(PFref_iHS_crs$iHS.PF), 
                                       sum(OTref_iHS_crs$iHS.L), sum(PFref_iHS_crs$iHS.L), 
                                       sum(OTref_iHS_crs$iHS.OT), sum(PFref_iHS_crs$iHS.OT), 
                                       sum(OTref_iHS_crs$iHS.V), sum(PFref_iHS_crs$iHS.V), 
                                       sum(OTref_iHS_crs$iHS.NR), sum(PFref_iHS_crs$iHS.NR)))
iHS.n.outliers.FDR <- as.character(c(sum(OTref_iHS_crs$iHS.PF.FDR), sum(PFref_iHS_crs$iHS.PF.FDR), 
                                     sum(OTref_iHS_crs$iHS.L.FDR), sum(PFref_iHS_crs$iHS.L.FDR), 
                                     sum(OTref_iHS_crs$iHS.OT.FDR), sum(PFref_iHS_crs$iHS.OT.FDR), 
                                     sum(OTref_iHS_crs$iHS.V.FDR), sum(PFref_iHS_crs$iHS.V.FDR), 
                                     sum(OTref_iHS_crs$iHS.NR.FDR), sum(PFref_iHS_crs$iHS.NR.FDR)))
iHS.n.cr.noFDR <- as.character(c(nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.PF != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.PF != 0), ]), 
                                 nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.L != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.L != 0), ]), 
                                 nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.OT != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.OT != 0), ]), 
                                 nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.V != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.V != 0), ]), 
                                 nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.NR != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.NR != 0), ])))
iHS.n.cr.FDR <- as.character(c(nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.PF.FDR != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.PF.FDR != 0), ]), 
                               nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.L.FDR != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.L.FDR != 0), ]), 
                               nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.OT.FDR != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.OT.FDR != 0), ]), 
                               nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.V.FDR != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.V.FDR != 0), ]), 
                               nrow(OTref_iHS_crs[which(OTref_iHS_crs$iHS.NR.FDR != 0), ]), nrow(PFref_iHS_crs[which(PFref_iHS_crs$iHS.NR.FDR != 0), ])))
iHS.n.scaffolds.noFDR <- as.character(c(length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.PF != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.PF != 0)])), 
                                        length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.L != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.L != 0)])), 
                                        length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.OT != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.OT != 0)])), 
                                        length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.V != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.V != 0)])), 
                                        length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.NR != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.NR != 0)]))))
iHS.n.scaffolds.FDR <- as.character(c(length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.PF.FDR != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.PF.FDR != 0)])), 
                                      length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.L.FDR != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.L.FDR != 0)])), 
                                      length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.OT.FDR != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.OT.FDR != 0)])), 
                                      length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.V.FDR != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.V.FDR != 0)])), 
                                      length(unique(OTref_iHS_crs$CHR[which(OTref_iHS_crs$iHS.NR.FDR != 0)])), length(unique(PFref_iHS_crs$CHR[which(PFref_iHS_crs$iHS.NR.FDR != 0)]))))
iHS.mean.outliers.per.cr.noFDR <- as.character(c(round(mean(OTref_iHS_crs$iHS.PF), 3), round(mean(PFref_iHS_crs$iHS.PF), 3), 
                                                 round(mean(OTref_iHS_crs$iHS.L), 3), round(mean(PFref_iHS_crs$iHS.L), 3), 
                                                 round(mean(OTref_iHS_crs$iHS.OT), 3), round(mean(PFref_iHS_crs$iHS.OT), 3), 
                                                 round(mean(OTref_iHS_crs$iHS.V), 3), round(mean(PFref_iHS_crs$iHS.V), 3), 
                                                 round(mean(OTref_iHS_crs$iHS.NR), 3), round(mean(PFref_iHS_crs$iHS.NR), 3)))
iHS.mean.outliers.per.cr.FDR <- as.character(c(round(mean(OTref_iHS_crs$iHS.PF.FDR), 3), round(mean(PFref_iHS_crs$iHS.PF.FDR), 3), 
                                               round(mean(OTref_iHS_crs$iHS.L.FDR), 3), round(mean(PFref_iHS_crs$iHS.L.FDR), 3), 
                                               round(mean(OTref_iHS_crs$iHS.OT.FDR), 3), round(mean(PFref_iHS_crs$iHS.OT.FDR), 3), 
                                               round(mean(OTref_iHS_crs$iHS.V.FDR), 3), round(mean(PFref_iHS_crs$iHS.V.FDR), 3), 
                                               round(mean(OTref_iHS_crs$iHS.NR.FDR), 3), round(mean(PFref_iHS_crs$iHS.NR.FDR), 3)))
iHS.max.N.outliers.in.cr.noFDR <- as.character(c(max(OTref_iHS_crs$iHS.PF), max(PFref_iHS_crs$iHS.PF), 
                                                 max(OTref_iHS_crs$iHS.L), max(PFref_iHS_crs$iHS.L), 
                                                 max(OTref_iHS_crs$iHS.OT), max(PFref_iHS_crs$iHS.OT), 
                                                 max(OTref_iHS_crs$iHS.V), max(PFref_iHS_crs$iHS.V), 
                                                 max(OTref_iHS_crs$iHS.NR), max(PFref_iHS_crs$iHS.NR)))
iHS.max.N.outliers.in.cr.FDR <- as.character(c(max(OTref_iHS_crs$iHS.PF.FDR), max(PFref_iHS_crs$iHS.PF.FDR), 
                                               max(OTref_iHS_crs$iHS.L.FDR), max(PFref_iHS_crs$iHS.L.FDR), 
                                               max(OTref_iHS_crs$iHS.OT.FDR), max(PFref_iHS_crs$iHS.OT.FDR), 
                                               max(OTref_iHS_crs$iHS.V.FDR), max(PFref_iHS_crs$iHS.V.FDR), 
                                               max(OTref_iHS_crs$iHS.NR.FDR), max(PFref_iHS_crs$iHS.NR.FDR)))
#Create summary data frame:
iHS.summary.data.frame <- data.frame(Population = iHS.summary.pop, 
                                     Reference = iHS.summary.ref, 
                                     N.outliers = paste(paste(iHS.n.outliers.noFDR, iHS.n.outliers.FDR, sep=" ("), ")", sep=""), 
                                     N.candidate.regions = paste(paste(iHS.n.cr.noFDR, iHS.n.cr.FDR, sep=" ("), ")", sep=""), 
                                     max.N.outliers.in.cr = paste(paste(iHS.max.N.outliers.in.cr.noFDR, iHS.max.N.outliers.in.cr.FDR, sep=" ("), ")", sep=""), 
                                     mean.outliers.per.cr = paste(paste(iHS.mean.outliers.per.cr.noFDR, iHS.mean.outliers.per.cr.FDR, sep=" ("), ")", sep=""), 
                                     N.scaffolds.with.crs = paste(paste(iHS.n.scaffolds.noFDR, iHS.n.scaffolds.FDR, sep=" ("), ")", sep=""))
#Write to file:
write.table(iHS.summary.data.frame, file="iHS_master_summary_table.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Cross-population summary tables...
#Create xpEHH cr subsets:
OTref_xpEHH_crs <- OTref_candidate_regions_with_annotations %>% select(CHR, START, END, LENGTH, contains("xpEHH"))
PFref_xpEHH_crs <- PFref_candidate_regions_with_annotations %>% select(CHR, START, END, LENGTH, contains("xpEHH"))
#Create Rsb cr subsets:
OTref_Rsb_crs <- OTref_candidate_regions_with_annotations %>% select(CHR, START, END, LENGTH, contains("Rsb"))
PFref_Rsb_crs <- PFref_candidate_regions_with_annotations %>% select(CHR, START, END, LENGTH, contains("Rsb"))
#Create summary table vectors:
xp.pop.pairs <- rep(c("OT_PF", "OT_L", "PF_V", "V_L", 
                      "OT_V", "PF_L", 
                      "OT_NR", "NR_V", "NR_PF", "NR_L"), 
                    each=2)
xp.ref <- rep(c("OT", "PF"), each=10)
xp.comparison <- rep(c(rep("between", 4), rep("within", 2), rep("unknown", 4)), each=2)
#xpEHH-specific:
xp.xpEHH.n.mrkrs.noFDR <- as.character(c(sum(OTref_xpEHH_crs$xpEHH.OT_PF), sum(PFref_xpEHH_crs$xpEHH.OT_PF), 
                                         sum(OTref_xpEHH_crs$xpEHH.OT_L), sum(PFref_xpEHH_crs$xpEHH.OT_L), 
                                         sum(OTref_xpEHH_crs$xpEHH.PF_V), sum(PFref_xpEHH_crs$xpEHH.PF_V), 
                                         sum(OTref_xpEHH_crs$xpEHH.V_L), sum(PFref_xpEHH_crs$xpEHH.V_L), 
                                         sum(OTref_xpEHH_crs$xpEHH.OT_V), sum(PFref_xpEHH_crs$xpEHH.OT_V), 
                                         sum(OTref_xpEHH_crs$xpEHH.PF_L), sum(PFref_xpEHH_crs$xpEHH.PF_L), 
                                         sum(OTref_xpEHH_crs$xpEHH.OT_NR), sum(PFref_xpEHH_crs$xpEHH.OT_NR), 
                                         sum(OTref_xpEHH_crs$xpEHH.NR_V), sum(PFref_xpEHH_crs$xpEHH.NR_V), 
                                         sum(OTref_xpEHH_crs$xpEHH.NR_PF), sum(PFref_xpEHH_crs$xpEHH.NR_PF), 
                                         sum(OTref_xpEHH_crs$xpEHH.NR_L), sum(PFref_xpEHH_crs$xpEHH.NR_L)))
xp.xpEHH.n.mrkrs.FDR <- as.character(c(sum(OTref_xpEHH_crs$xpEHH.OT_PF.FDR), sum(PFref_xpEHH_crs$xpEHH.OT_PF.FDR), 
                                       sum(OTref_xpEHH_crs$xpEHH.OT_L.FDR), sum(PFref_xpEHH_crs$xpEHH.OT_L.FDR), 
                                       sum(OTref_xpEHH_crs$xpEHH.PF_V.FDR), sum(PFref_xpEHH_crs$xpEHH.PF_V.FDR), 
                                       sum(OTref_xpEHH_crs$xpEHH.V_L.FDR), sum(PFref_xpEHH_crs$xpEHH.V_L.FDR), 
                                       sum(OTref_xpEHH_crs$xpEHH.OT_V.FDR), sum(PFref_xpEHH_crs$xpEHH.OT_V.FDR), 
                                       sum(OTref_xpEHH_crs$xpEHH.PF_L.FDR), sum(PFref_xpEHH_crs$xpEHH.PF_L.FDR), 
                                       sum(OTref_xpEHH_crs$xpEHH.OT_NR.FDR), sum(PFref_xpEHH_crs$xpEHH.OT_NR.FDR), 
                                       sum(OTref_xpEHH_crs$xpEHH.NR_V.FDR), sum(PFref_xpEHH_crs$xpEHH.NR_V.FDR), 
                                       sum(OTref_xpEHH_crs$xpEHH.NR_PF.FDR), sum(PFref_xpEHH_crs$xpEHH.NR_PF.FDR), 
                                       sum(OTref_xpEHH_crs$xpEHH.NR_L.FDR), sum(PFref_xpEHH_crs$xpEHH.NR_L.FDR)))
xp.xpEHH.n.cr.noFDR <- as.character(c(nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.OT_PF != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.OT_PF != 0), ]), 
                                      nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.OT_L != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.OT_L != 0), ]), 
                                      nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.PF_V != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.PF_V != 0), ]), 
                                      nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.V_L != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.V_L != 0), ]), 
                                      nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.OT_V != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.OT_V != 0), ]), 
                                      nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.PF_L != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.PF_L != 0), ]), 
                                      nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.OT_NR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.OT_NR != 0), ]), 
                                      nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.NR_V != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.NR_V != 0), ]), 
                                      nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.NR_PF != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.NR_PF != 0), ]), 
                                      nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.NR_L != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.NR_L != 0), ])))
xp.xpEHH.n.cr.FDR <- as.character(c(nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.OT_PF.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.OT_PF.FDR != 0), ]), 
                                    nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.OT_L.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.OT_L.FDR != 0), ]), 
                                    nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.PF_V.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.PF_V.FDR != 0), ]), 
                                    nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.V_L.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.V_L.FDR != 0), ]), 
                                    nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.OT_V.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.OT_V.FDR != 0), ]), 
                                    nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.PF_L.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.PF_L.FDR != 0), ]), 
                                    nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.OT_NR.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.OT_NR.FDR != 0), ]), 
                                    nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.NR_V.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.NR_V.FDR != 0), ]), 
                                    nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.NR_PF.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.NR_PF.FDR != 0), ]), 
                                    nrow(OTref_xpEHH_crs[which(OTref_xpEHH_crs$xpEHH.NR_L.FDR != 0), ]), nrow(PFref_xpEHH_crs[which(PFref_xpEHH_crs$xpEHH.NR_L.FDR != 0), ])))
xp.xpEHH.max.N.outliers.in.cr.noFDR <- as.character(c(max(OTref_xpEHH_crs$xpEHH.OT_PF), max(PFref_xpEHH_crs$xpEHH.OT_PF), 
                                                      max(OTref_xpEHH_crs$xpEHH.OT_L), max(PFref_xpEHH_crs$xpEHH.OT_L), 
                                                      max(OTref_xpEHH_crs$xpEHH.PF_V), max(PFref_xpEHH_crs$xpEHH.PF_V), 
                                                      max(OTref_xpEHH_crs$xpEHH.V_L), max(PFref_xpEHH_crs$xpEHH.V_L), 
                                                      max(OTref_xpEHH_crs$xpEHH.OT_V), max(PFref_xpEHH_crs$xpEHH.OT_V), 
                                                      max(OTref_xpEHH_crs$xpEHH.PF_L), max(PFref_xpEHH_crs$xpEHH.PF_L), 
                                                      max(OTref_xpEHH_crs$xpEHH.OT_NR), max(PFref_xpEHH_crs$xpEHH.OT_NR), 
                                                      max(OTref_xpEHH_crs$xpEHH.NR_V), max(PFref_xpEHH_crs$xpEHH.NR_V), 
                                                      max(OTref_xpEHH_crs$xpEHH.NR_PF), max(PFref_xpEHH_crs$xpEHH.NR_PF), 
                                                      max(OTref_xpEHH_crs$xpEHH.NR_L), max(PFref_xpEHH_crs$xpEHH.NR_L)))
xp.xpEHH.max.N.outliers.in.cr.FDR <- as.character(c(max(OTref_xpEHH_crs$xpEHH.OT_PF.FDR), max(PFref_xpEHH_crs$xpEHH.OT_PF.FDR), 
                                                    max(OTref_xpEHH_crs$xpEHH.OT_L.FDR), max(PFref_xpEHH_crs$xpEHH.OT_L.FDR), 
                                                    max(OTref_xpEHH_crs$xpEHH.PF_V.FDR), max(PFref_xpEHH_crs$xpEHH.PF_V.FDR), 
                                                    max(OTref_xpEHH_crs$xpEHH.V_L.FDR), max(PFref_xpEHH_crs$xpEHH.V_L.FDR), 
                                                    max(OTref_xpEHH_crs$xpEHH.OT_V.FDR), max(PFref_xpEHH_crs$xpEHH.OT_V.FDR), 
                                                    max(OTref_xpEHH_crs$xpEHH.PF_L.FDR), max(PFref_xpEHH_crs$xpEHH.PF_L.FDR), 
                                                    max(OTref_xpEHH_crs$xpEHH.OT_NR.FDR), max(PFref_xpEHH_crs$xpEHH.OT_NR.FDR), 
                                                    max(OTref_xpEHH_crs$xpEHH.NR_V.FDR), max(PFref_xpEHH_crs$xpEHH.NR_V.FDR), 
                                                    max(OTref_xpEHH_crs$xpEHH.NR_PF.FDR), max(PFref_xpEHH_crs$xpEHH.NR_PF.FDR), 
                                                    max(OTref_xpEHH_crs$xpEHH.NR_L.FDR), max(PFref_xpEHH_crs$xpEHH.NR_L.FDR)))
xp.xpEHH.mean.outliers.per.cr.noFDR <- as.character(c(round(mean(OTref_xpEHH_crs$xpEHH.OT_PF), 3), round(mean(PFref_xpEHH_crs$xpEHH.OT_PF), 3), 
                                                      round(mean(OTref_xpEHH_crs$xpEHH.OT_L), 3), round(mean(PFref_xpEHH_crs$xpEHH.OT_L), 3), 
                                                      round(mean(OTref_xpEHH_crs$xpEHH.PF_V), 3), round(mean(PFref_xpEHH_crs$xpEHH.PF_V), 3), 
                                                      round(mean(OTref_xpEHH_crs$xpEHH.V_L), 3), round(mean(PFref_xpEHH_crs$xpEHH.V_L), 3), 
                                                      round(mean(OTref_xpEHH_crs$xpEHH.OT_V), 3), round(mean(PFref_xpEHH_crs$xpEHH.OT_V), 3), 
                                                      round(mean(OTref_xpEHH_crs$xpEHH.PF_L), 3), round(mean(PFref_xpEHH_crs$xpEHH.PF_L), 3), 
                                                      round(mean(OTref_xpEHH_crs$xpEHH.OT_NR), 3), round(mean(PFref_xpEHH_crs$xpEHH.OT_NR), 3), 
                                                      round(mean(OTref_xpEHH_crs$xpEHH.NR_V), 3), round(mean(PFref_xpEHH_crs$xpEHH.NR_V), 3), 
                                                      round(mean(OTref_xpEHH_crs$xpEHH.NR_PF), 3), round(mean(PFref_xpEHH_crs$xpEHH.NR_PF), 3), 
                                                      round(mean(OTref_xpEHH_crs$xpEHH.NR_L), 3), round(mean(PFref_xpEHH_crs$xpEHH.NR_L), 3)))
xp.xpEHH.mean.outliers.per.cr.FDR <- as.character(c(round(mean(OTref_xpEHH_crs$xpEHH.OT_PF.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.OT_PF.FDR), 3), 
                                                    round(mean(OTref_xpEHH_crs$xpEHH.OT_L.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.OT_L.FDR), 3), 
                                                    round(mean(OTref_xpEHH_crs$xpEHH.PF_V.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.PF_V.FDR), 3), 
                                                    round(mean(OTref_xpEHH_crs$xpEHH.V_L.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.V_L.FDR), 3), 
                                                    round(mean(OTref_xpEHH_crs$xpEHH.OT_V.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.OT_V.FDR), 3), 
                                                    round(mean(OTref_xpEHH_crs$xpEHH.PF_L.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.PF_L.FDR), 3), 
                                                    round(mean(OTref_xpEHH_crs$xpEHH.OT_NR.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.OT_NR.FDR), 3), 
                                                    round(mean(OTref_xpEHH_crs$xpEHH.NR_V.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.NR_V.FDR), 3), 
                                                    round(mean(OTref_xpEHH_crs$xpEHH.NR_PF.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.NR_PF.FDR), 3), 
                                                    round(mean(OTref_xpEHH_crs$xpEHH.NR_L.FDR), 3), round(mean(PFref_xpEHH_crs$xpEHH.NR_L.FDR), 3)))
xp.xpEHH.n.scaffolds.noFDR <- as.character(c(length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.OT_PF != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.OT_PF != 0)])), 
                                             length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.OT_L != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.OT_L != 0)])), 
                                             length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.PF_V != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.PF_V != 0)])), 
                                             length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.V_L != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.V_L != 0)])), 
                                             length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.OT_V != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.OT_V != 0)])), 
                                             length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.PF_L != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.PF_L != 0)])), 
                                             length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.OT_NR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.OT_NR != 0)])), 
                                             length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.NR_V != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.NR_V != 0)])), 
                                             length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.NR_PF != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.NR_PF != 0)])), 
                                             length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.NR_L != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.NR_L != 0)]))))
xp.xpEHH.n.scaffolds.FDR <- as.character(c(length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.OT_PF.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.OT_PF.FDR != 0)])), 
                                           length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.OT_L.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.OT_L.FDR != 0)])), 
                                           length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.PF_V.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.PF_V.FDR != 0)])), 
                                           length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.V_L.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.V_L.FDR != 0)])), 
                                           length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.OT_V.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.OT_V.FDR != 0)])), 
                                           length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.PF_L.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.PF_L.FDR != 0)])), 
                                           length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.OT_NR.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.OT_NR.FDR != 0)])), 
                                           length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.NR_V.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.NR_V.FDR != 0)])), 
                                           length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.NR_PF.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.NR_PF.FDR != 0)])), 
                                           length(unique(OTref_xpEHH_crs$CHR[which(OTref_xpEHH_crs$xpEHH.NR_L.FDR != 0)])), length(unique(PFref_xpEHH_crs$CHR[which(PFref_xpEHH_crs$xpEHH.NR_L.FDR != 0)]))))
#Rsb-specific:
xp.Rsb.n.mrkrs.noFDR <- as.character(c(sum(OTref_Rsb_crs$Rsb.OT_PF), sum(PFref_Rsb_crs$Rsb.OT_PF), 
                                       sum(OTref_Rsb_crs$Rsb.OT_L), sum(PFref_Rsb_crs$Rsb.OT_L), 
                                       sum(OTref_Rsb_crs$Rsb.PF_V), sum(PFref_Rsb_crs$Rsb.PF_V), 
                                       sum(OTref_Rsb_crs$Rsb.V_L), sum(PFref_Rsb_crs$Rsb.V_L), 
                                       sum(OTref_Rsb_crs$Rsb.OT_V), sum(PFref_Rsb_crs$Rsb.OT_V), 
                                       sum(OTref_Rsb_crs$Rsb.PF_L), sum(PFref_Rsb_crs$Rsb.PF_L), 
                                       sum(OTref_Rsb_crs$Rsb.OT_NR), sum(PFref_Rsb_crs$Rsb.OT_NR), 
                                       sum(OTref_Rsb_crs$Rsb.NR_V), sum(PFref_Rsb_crs$Rsb.NR_V), 
                                       sum(OTref_Rsb_crs$Rsb.NR_PF), sum(PFref_Rsb_crs$Rsb.NR_PF), 
                                       sum(OTref_Rsb_crs$Rsb.NR_L), sum(PFref_Rsb_crs$Rsb.NR_L)))
xp.Rsb.n.mrkrs.FDR <- as.character(c(sum(OTref_Rsb_crs$Rsb.OT_PF.FDR), sum(PFref_Rsb_crs$Rsb.OT_PF.FDR), 
                                     sum(OTref_Rsb_crs$Rsb.OT_L.FDR), sum(PFref_Rsb_crs$Rsb.OT_L.FDR), 
                                     sum(OTref_Rsb_crs$Rsb.PF_V.FDR), sum(PFref_Rsb_crs$Rsb.PF_V.FDR), 
                                     sum(OTref_Rsb_crs$Rsb.V_L.FDR), sum(PFref_Rsb_crs$Rsb.V_L.FDR), 
                                     sum(OTref_Rsb_crs$Rsb.OT_V.FDR), sum(PFref_Rsb_crs$Rsb.OT_V.FDR), 
                                     sum(OTref_Rsb_crs$Rsb.PF_L.FDR), sum(PFref_Rsb_crs$Rsb.PF_L.FDR), 
                                     sum(OTref_Rsb_crs$Rsb.OT_NR.FDR), sum(PFref_Rsb_crs$Rsb.OT_NR.FDR), 
                                     sum(OTref_Rsb_crs$Rsb.NR_V.FDR), sum(PFref_Rsb_crs$Rsb.NR_V.FDR), 
                                     sum(OTref_Rsb_crs$Rsb.NR_PF.FDR), sum(PFref_Rsb_crs$Rsb.NR_PF.FDR), 
                                     sum(OTref_Rsb_crs$Rsb.NR_L.FDR), sum(PFref_Rsb_crs$Rsb.NR_L.FDR)))
xp.Rsb.n.cr.noFDR <- as.character(c(nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.OT_PF != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.OT_PF != 0), ]), 
                                    nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.OT_L != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.OT_L != 0), ]), 
                                    nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.PF_V != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.PF_V != 0), ]), 
                                    nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.V_L != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.V_L != 0), ]), 
                                    nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.OT_V != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.OT_V != 0), ]), 
                                    nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.PF_L != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.PF_L != 0), ]), 
                                    nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.OT_NR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.OT_NR != 0), ]), 
                                    nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.NR_V != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.NR_V != 0), ]), 
                                    nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.NR_PF != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.NR_PF != 0), ]), 
                                    nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.NR_L != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.NR_L != 0), ])))
xp.Rsb.n.cr.FDR <- as.character(c(nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.OT_PF.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.OT_PF.FDR != 0), ]), 
                                  nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.OT_L.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.OT_L.FDR != 0), ]), 
                                  nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.PF_V.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.PF_V.FDR != 0), ]), 
                                  nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.V_L.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.V_L.FDR != 0), ]), 
                                  nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.OT_V.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.OT_V.FDR != 0), ]), 
                                  nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.PF_L.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.PF_L.FDR != 0), ]), 
                                  nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.OT_NR.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.OT_NR.FDR != 0), ]), 
                                  nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.NR_V.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.NR_V.FDR != 0), ]), 
                                  nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.NR_PF.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.NR_PF.FDR != 0), ]), 
                                  nrow(OTref_Rsb_crs[which(OTref_Rsb_crs$Rsb.NR_L.FDR != 0), ]), nrow(PFref_Rsb_crs[which(PFref_Rsb_crs$Rsb.NR_L.FDR != 0), ])))
xp.Rsb.max.N.outliers.in.cr.noFDR <- as.character(c(max(OTref_Rsb_crs$Rsb.OT_PF), max(PFref_Rsb_crs$Rsb.OT_PF), 
                                                    max(OTref_Rsb_crs$Rsb.OT_L), max(PFref_Rsb_crs$Rsb.OT_L), 
                                                    max(OTref_Rsb_crs$Rsb.PF_V), max(PFref_Rsb_crs$Rsb.PF_V), 
                                                    max(OTref_Rsb_crs$Rsb.V_L), max(PFref_Rsb_crs$Rsb.V_L), 
                                                    max(OTref_Rsb_crs$Rsb.OT_V), max(PFref_Rsb_crs$Rsb.OT_V), 
                                                    max(OTref_Rsb_crs$Rsb.PF_L), max(PFref_Rsb_crs$Rsb.PF_L), 
                                                    max(OTref_Rsb_crs$Rsb.OT_NR), max(PFref_Rsb_crs$Rsb.OT_NR), 
                                                    max(OTref_Rsb_crs$Rsb.NR_V), max(PFref_Rsb_crs$Rsb.NR_V), 
                                                    max(OTref_Rsb_crs$Rsb.NR_PF), max(PFref_Rsb_crs$Rsb.NR_PF), 
                                                    max(OTref_Rsb_crs$Rsb.NR_L), max(PFref_Rsb_crs$Rsb.NR_L)))
xp.Rsb.max.N.outliers.in.cr.FDR <- as.character(c(max(OTref_Rsb_crs$Rsb.OT_PF.FDR), max(PFref_Rsb_crs$Rsb.OT_PF.FDR), 
                                                  max(OTref_Rsb_crs$Rsb.OT_L.FDR), max(PFref_Rsb_crs$Rsb.OT_L.FDR), 
                                                  max(OTref_Rsb_crs$Rsb.PF_V.FDR), max(PFref_Rsb_crs$Rsb.PF_V.FDR), 
                                                  max(OTref_Rsb_crs$Rsb.V_L.FDR), max(PFref_Rsb_crs$Rsb.V_L.FDR), 
                                                  max(OTref_Rsb_crs$Rsb.OT_V.FDR), max(PFref_Rsb_crs$Rsb.OT_V.FDR), 
                                                  max(OTref_Rsb_crs$Rsb.PF_L.FDR), max(PFref_Rsb_crs$Rsb.PF_L.FDR), 
                                                  max(OTref_Rsb_crs$Rsb.OT_NR.FDR), max(PFref_Rsb_crs$Rsb.OT_NR.FDR), 
                                                  max(OTref_Rsb_crs$Rsb.NR_V.FDR), max(PFref_Rsb_crs$Rsb.NR_V.FDR), 
                                                  max(OTref_Rsb_crs$Rsb.NR_PF.FDR), max(PFref_Rsb_crs$Rsb.NR_PF.FDR), 
                                                  max(OTref_Rsb_crs$Rsb.NR_L.FDR), max(PFref_Rsb_crs$Rsb.NR_L.FDR)))
xp.Rsb.mean.outliers.per.cr.noFDR <- as.character(c(round(mean(OTref_Rsb_crs$Rsb.OT_PF), 3), round(mean(PFref_Rsb_crs$Rsb.OT_PF), 3), 
                                                    round(mean(OTref_Rsb_crs$Rsb.OT_L), 3), round(mean(PFref_Rsb_crs$Rsb.OT_L), 3), 
                                                    round(mean(OTref_Rsb_crs$Rsb.PF_V), 3), round(mean(PFref_Rsb_crs$Rsb.PF_V), 3), 
                                                    round(mean(OTref_Rsb_crs$Rsb.V_L), 3), round(mean(PFref_Rsb_crs$Rsb.V_L), 3), 
                                                    round(mean(OTref_Rsb_crs$Rsb.OT_V), 3), round(mean(PFref_Rsb_crs$Rsb.OT_V), 3), 
                                                    round(mean(OTref_Rsb_crs$Rsb.PF_L), 3), round(mean(PFref_Rsb_crs$Rsb.PF_L), 3), 
                                                    round(mean(OTref_Rsb_crs$Rsb.OT_NR), 3), round(mean(PFref_Rsb_crs$Rsb.OT_NR), 3), 
                                                    round(mean(OTref_Rsb_crs$Rsb.NR_V), 3), round(mean(PFref_Rsb_crs$Rsb.NR_V), 3), 
                                                    round(mean(OTref_Rsb_crs$Rsb.NR_PF), 3), round(mean(PFref_Rsb_crs$Rsb.NR_PF), 3), 
                                                    round(mean(OTref_Rsb_crs$Rsb.NR_L), 3), round(mean(PFref_Rsb_crs$Rsb.NR_L), 3)))
xp.Rsb.mean.outliers.per.cr.FDR <- as.character(c(round(mean(OTref_Rsb_crs$Rsb.OT_PF.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.OT_PF.FDR), 3), 
                                                  round(mean(OTref_Rsb_crs$Rsb.OT_L.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.OT_L.FDR), 3), 
                                                  round(mean(OTref_Rsb_crs$Rsb.PF_V.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.PF_V.FDR), 3), 
                                                  round(mean(OTref_Rsb_crs$Rsb.V_L.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.V_L.FDR), 3), 
                                                  round(mean(OTref_Rsb_crs$Rsb.OT_V.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.OT_V.FDR), 3), 
                                                  round(mean(OTref_Rsb_crs$Rsb.PF_L.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.PF_L.FDR), 3), 
                                                  round(mean(OTref_Rsb_crs$Rsb.OT_NR.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.OT_NR.FDR), 3), 
                                                  round(mean(OTref_Rsb_crs$Rsb.NR_V.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.NR_V.FDR), 3), 
                                                  round(mean(OTref_Rsb_crs$Rsb.NR_PF.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.NR_PF.FDR), 3), 
                                                  round(mean(OTref_Rsb_crs$Rsb.NR_L.FDR), 3), round(mean(PFref_Rsb_crs$Rsb.NR_L.FDR), 3)))
xp.Rsb.n.scaffolds.noFDR <- as.character(c(length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.OT_PF != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.OT_PF != 0)])), 
                                           length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.OT_L != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.OT_L != 0)])), 
                                           length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.PF_V != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.PF_V != 0)])), 
                                           length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.V_L != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.V_L != 0)])), 
                                           length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.OT_V != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.OT_V != 0)])), 
                                           length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.PF_L != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.PF_L != 0)])), 
                                           length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.OT_NR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.OT_NR != 0)])), 
                                           length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.NR_V != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.NR_V != 0)])), 
                                           length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.NR_PF != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.NR_PF != 0)])), 
                                           length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.NR_L != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.NR_L != 0)]))))
xp.Rsb.n.scaffolds.FDR <- as.character(c(length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.OT_PF.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.OT_PF.FDR != 0)])), 
                                         length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.OT_L.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.OT_L.FDR != 0)])), 
                                         length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.PF_V.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.PF_V.FDR != 0)])), 
                                         length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.V_L.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.V_L.FDR != 0)])), 
                                         length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.OT_V.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.OT_V.FDR != 0)])), 
                                         length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.PF_L.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.PF_L.FDR != 0)])), 
                                         length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.OT_NR.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.OT_NR.FDR != 0)])), 
                                         length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.NR_V.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.NR_V.FDR != 0)])), 
                                         length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.NR_PF.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.NR_PF.FDR != 0)])), 
                                         length(unique(OTref_Rsb_crs$CHR[which(OTref_Rsb_crs$Rsb.NR_L.FDR != 0)])), length(unique(PFref_Rsb_crs$CHR[which(PFref_Rsb_crs$Rsb.NR_L.FDR != 0)]))))
#Create xpEHH summary table:
xpEHH.summary.data.frame <- data.frame(Population.pair = xp.pop.pairs, 
                                       Reference = xp.ref, 
                                       Comparison.type = xp.comparison, 
                                       N.outliers = paste(paste(xp.xpEHH.n.mrkrs.noFDR, xp.xpEHH.n.mrkrs.FDR, sep=" ("), ")", sep=""), 
                                       N.candidate.regions = paste(paste(xp.xpEHH.n.cr.noFDR, xp.xpEHH.n.cr.FDR, sep=" ("), ")", sep=""), 
                                       max.N.outliers.in.cr = paste(paste(xp.xpEHH.max.N.outliers.in.cr.noFDR, xp.xpEHH.max.N.outliers.in.cr.FDR, sep=" ("), ")", sep=""), 
                                       mean.outliers.per.cr = paste(paste(xp.xpEHH.mean.outliers.per.cr.noFDR, xp.xpEHH.mean.outliers.per.cr.FDR, sep=" ("), ")", sep=""), 
                                       N.scaffolds.with.crs = paste(paste(xp.xpEHH.n.scaffolds.noFDR, xp.xpEHH.n.scaffolds.FDR, sep=" ("), ")", sep=""))
#Create Rsb summary table:
Rsb.summary.data.frame <- data.frame(Population.pair = xp.pop.pairs, 
                                     Reference = xp.ref, 
                                     Comparison.type = xp.comparison, 
                                     N.outliers = paste(paste(xp.Rsb.n.mrkrs.noFDR, xp.Rsb.n.mrkrs.FDR, sep=" ("), ")", sep=""), 
                                     N.candidate.regions = paste(paste(xp.Rsb.n.cr.noFDR, xp.Rsb.n.cr.FDR, sep=" ("), ")", sep=""), 
                                     max.N.outliers.in.cr = paste(paste(xp.Rsb.max.N.outliers.in.cr.noFDR, xp.Rsb.max.N.outliers.in.cr.FDR, sep=" ("), ")", sep=""), 
                                     mean.outliers.per.cr = paste(paste(xp.Rsb.mean.outliers.per.cr.noFDR, xp.Rsb.mean.outliers.per.cr.FDR, sep=" ("), ")", sep=""), 
                                     N.scaffolds.with.crs = paste(paste(xp.Rsb.n.scaffolds.noFDR, xp.Rsb.n.scaffolds.FDR, sep=" ("), ")", sep=""))
#Write to file:
write.table(xpEHH.summary.data.frame, file="xpEHH_master_summary_table.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(Rsb.summary.data.frame, file="Rsb_master_summary_table.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)

#Social candidate regions summary table...

#First create a subset of each social candidate regions master data frame
#containing only FDR-corrected results:
OTref_social_candidate_regions_master.FDR <- OTref_social_candidate_regions_master %>% 
  select(CHR, contains("is.soc")) %>% 
  select(CHR, contains(".FDR.")) %>%
  filter(is.soc.cr.FDR.xpEHH == 1 | is.soc.cr.FDR.Rsb == 1)
PFref_social_candidate_regions_master.FDR <- PFref_social_candidate_regions_master %>% 
  select(CHR, contains("is.soc")) %>% 
  select(CHR, contains(".FDR.")) %>%
  filter(is.soc.cr.FDR.xpEHH == 1 | is.soc.cr.FDR.Rsb == 1)
#Some more subsets:
OTref_social_candidate_regions_xpEHH <- OTref_social_candidate_regions_master %>% 
  select(CHR, contains("is.soc")) %>% 
  select(CHR, contains("xpEHH")) %>% 
  filter(is.soc.cr.noFDR.xpEHH == 1)
PFref_social_candidate_regions_xpEHH <- PFref_social_candidate_regions_master %>% 
  select(CHR, contains("is.soc")) %>% 
  select(CHR, contains("xpEHH")) %>% 
  filter(is.soc.cr.noFDR.xpEHH == 1)
OTref_social_candidate_regions_Rsb <- OTref_social_candidate_regions_master %>% 
  select(CHR, contains("is.soc")) %>% 
  select(CHR, contains("Rsb")) %>% 
  filter(is.soc.cr.noFDR.Rsb == 1)
PFref_social_candidate_regions_Rsb <- PFref_social_candidate_regions_master %>% 
  select(CHR, contains("is.soc")) %>% 
  select(CHR, contains("Rsb")) %>% 
  filter(is.soc.cr.noFDR.Rsb == 1)
OTref_social_candidate_regions_xpEHH.FDR <- OTref_social_candidate_regions_xpEHH %>% 
  filter(is.soc.cr.FDR.xpEHH == 1)
PFref_social_candidate_regions_xpEHH.FDR <- PFref_social_candidate_regions_xpEHH %>% 
  filter(is.soc.cr.FDR.xpEHH == 1)
OTref_social_candidate_regions_Rsb.FDR <- OTref_social_candidate_regions_Rsb %>% 
  filter(is.soc.cr.FDR.Rsb == 1)
PFref_social_candidate_regions_Rsb.FDR <- PFref_social_candidate_regions_Rsb %>% 
  filter(is.soc.cr.FDR.Rsb == 1)
#Create columns:
soc.cr.ref <- c("OT", "PF")
soc.cr.totals <- c(paste(paste(as.character(nrow(OTref_social_candidate_regions_master)), as.character(nrow(OTref_social_candidate_regions_master.FDR)), sep=" ("), ")", sep=""), 
                   paste(paste(as.character(nrow(PFref_social_candidate_regions_master)), as.character(nrow(PFref_social_candidate_regions_master.FDR)), sep=" ("), ")", sep=""))
soc.cr.total.scaffolds <- c(paste(paste(as.character(length(unique(OTref_social_candidate_regions_master$CHR))), as.character(length(unique(OTref_social_candidate_regions_master.FDR$CHR))), sep=" ("), ")", sep=""), 
                            paste(paste(as.character(length(unique(PFref_social_candidate_regions_master$CHR))), as.character(length(unique(PFref_social_candidate_regions_master.FDR$CHR))), sep=" ("), ")", sep=""))
soc.cr.n.xpEHH <- c(paste(paste(as.character(sum(OTref_social_candidate_regions_master$is.soc.cr.noFDR.xpEHH)), as.character(sum(OTref_social_candidate_regions_master$is.soc.cr.FDR.xpEHH)), sep=" ("), ")", sep=""), 
                    paste(paste(as.character(sum(PFref_social_candidate_regions_master$is.soc.cr.noFDR.xpEHH)), as.character(sum(PFref_social_candidate_regions_master$is.soc.cr.FDR.xpEHH)), sep=" ("), ")", sep=""))
soc.cr.nScaff.xpEHH <- c(paste(paste(as.character(length(unique(OTref_social_candidate_regions_master$CHR[which(OTref_social_candidate_regions_master$is.soc.cr.noFDR.xpEHH >= 1)]))), as.character(length(unique(OTref_social_candidate_regions_master.FDR$CHR[which(OTref_social_candidate_regions_master.FDR$is.soc.cr.FDR.xpEHH >= 1)]))), sep=" ("), ")", sep=""), 
                         paste(paste(as.character(length(unique(PFref_social_candidate_regions_master$CHR[which(PFref_social_candidate_regions_master$is.soc.cr.noFDR.xpEHH >= 1)]))), as.character(length(unique(PFref_social_candidate_regions_master.FDR$CHR[which(PFref_social_candidate_regions_master.FDR$is.soc.cr.FDR.xpEHH >= 1)]))), sep=" ("), ")", sep=""))
soc.cr.n.Rsb <- c(paste(paste(as.character(sum(OTref_social_candidate_regions_master$is.soc.cr.noFDR.Rsb)), as.character(sum(OTref_social_candidate_regions_master$is.soc.cr.FDR.Rsb)), sep=" ("), ")", sep=""), 
                  paste(paste(as.character(sum(PFref_social_candidate_regions_master$is.soc.cr.noFDR.Rsb)), as.character(sum(PFref_social_candidate_regions_master$is.soc.cr.FDR.Rsb)), sep=" ("), ")", sep=""))
soc.cr.nScaff.Rsb <- c(paste(paste(as.character(length(unique(OTref_social_candidate_regions_master$CHR[which(OTref_social_candidate_regions_master$is.soc.cr.noFDR.Rsb >= 1)]))), as.character(length(unique(OTref_social_candidate_regions_master.FDR$CHR[which(OTref_social_candidate_regions_master.FDR$is.soc.cr.FDR.Rsb >= 1)]))), sep=" ("), ")", sep=""), 
                       paste(paste(as.character(length(unique(PFref_social_candidate_regions_master$CHR[which(PFref_social_candidate_regions_master$is.soc.cr.noFDR.Rsb >= 1)]))), as.character(length(unique(PFref_social_candidate_regions_master.FDR$CHR[which(PFref_social_candidate_regions_master.FDR$is.soc.cr.FDR.Rsb >= 1)]))), sep=" ("), ")", sep=""))
soc.cr.xpEHH.nScaff.inRsb <- c(paste(paste(as.character(length(OTref_social_candidate_regions_xpEHH$CHR %in% OTref_social_candidate_regions_Rsb$CHR)), as.character(length(OTref_social_candidate_regions_xpEHH.FDR$CHR %in% OTref_social_candidate_regions_Rsb.FDR$CHR)), sep=" ("), ")", sep=""), 
                               paste(paste(as.character(length(PFref_social_candidate_regions_xpEHH$CHR %in% PFref_social_candidate_regions_Rsb$CHR)), as.character(length(PFref_social_candidate_regions_xpEHH.FDR$CHR %in% PFref_social_candidate_regions_Rsb.FDR$CHR)), sep=" ("), ")", sep=""))
soc.cr.Rsb.nScaff.inxpEHH <- c(paste(paste(as.character(length(OTref_social_candidate_regions_Rsb$CHR %in% OTref_social_candidate_regions_xpEHH$CHR)), as.character(length(OTref_social_candidate_regions_Rsb.FDR$CHR %in% OTref_social_candidate_regions_xpEHH.FDR$CHR)), sep=" ("), ")", sep=""), 
                               paste(paste(as.character(length(PFref_social_candidate_regions_Rsb$CHR %in% PFref_social_candidate_regions_xpEHH$CHR)), as.character(length(PFref_social_candidate_regions_Rsb.FDR$CHR %in% PFref_social_candidate_regions_xpEHH.FDR$CHR)), sep=" ("), ")", sep=""))
#Cbind to create summary table:
soc.cr.summary.table <- as.data.frame(cbind(soc.cr.ref, 
                                            soc.cr.totals, 
                                            soc.cr.total.scaffolds, 
                                            soc.cr.n.xpEHH, 
                                            soc.cr.nScaff.xpEHH, 
                                            soc.cr.n.Rsb, 
                                            soc.cr.nScaff.Rsb, 
                                            soc.cr.xpEHH.nScaff.inRsb, 
                                            soc.cr.Rsb.nScaff.inxpEHH))
#Write to file
write.table(soc.cr.summary.table, file="social_candidate_regions_summary_table.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)


##### Rough work
OTref.FDR.crs <- OTref_candidate_regions_with_annotations %>% select(contains("FDR"))
PFref.FDR.crs <- PFref_candidate_regions_with_annotations %>% select(contains("FDR"))
OTref.FDR.crs$row.sum <- rowSums(OTref.FDR.crs[,])
PFref.FDR.crs$row.sum <- rowSums(PFref.FDR.crs[,])
nrow(OTref.FDR.crs[which(OTref.FDR.crs$row.sum != 0), ])
nrow(PFref.FDR.crs[which(PFref.FDR.crs$row.sum != 0), ])

OTref.soc.FDR.crs <- OTref_social_candidate_regions_master %>% select(contains("is.soc")) %>% select(!contains("noFDR"))
PFref.soc.FDR.crs <- PFref_social_candidate_regions_master %>% select(contains("is.soc")) %>% select(!contains("noFDR"))
OTref.soc.FDR.crs$row.sum <- rowSums(OTref.soc.FDR.crs[,])
PFref.soc.FDR.crs$row.sum <- rowSums(PFref.soc.FDR.crs[,])
nrow(OTref.soc.FDR.crs[which(OTref.soc.FDR.crs$row.sum != 0), ])
nrow(PFref.soc.FDR.crs[which(PFref.soc.FDR.crs$row.sum != 0), ])


#####################################################################################

###################### Validating social candidate regions ##########################

#The test for social candidate regions relies on the assumption that social phenotype,
#and not chance, explains the diversity of regions identified.  This assumption can
#be tested by randomising social phenotype assignments and checking how many (if any)
#regions are identified using the same technique as was employed above.

#The rough approach is as follows...

#Original test:
#FM:
#  OT V
#P:
#  PF L
#
#soc.cr=cr(OT_PF, OT_L, V_PF, V_L), !cr(OT_V, PF_L)
#
#Rearrangements:
#
#Rearrange 1...
#
#FM:
#  OT PF
#P:
#  V L
#unsoc.cr=cr(OT_V, OT_L, PF_V, PF_L), !cr(OT_PF, V_L)
#
#
#Rearrange 2...
#
#FM:
#  V PF
#P:
#  OT L
#unsoc.cr=cr(V_OT, V_L, PF_OT, PF_L), !cr(V_PF, OT_L)
#
#
#Rearrange 3...
#
#FM:
#  V L
#P:
#  OT PF
#unsoc.cr=cr(V_OT, V_PF, L_OT, L_PF), !cr(V_L, OT_PF)
#
#Rearrange 4...
#
#FM:
#  OT L
#P:
#  V PF
#unsoc.cr=cr(OT_V, OT_PF, L_V, L_PF), !cr(OT_L, V_PF)
#
#... a 5th re-assignment would result in the original combinations
#... of within- and between-phenotype comparisons

#If none or only a very few regions are identified in all of these 
#re-arrangements this would strongly discount chance as the explanatory
#factor for the regions identified with the correct assignments.

#Clean up
rm(list = ls())

#Load libraries
library(tidyverse)

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
setwd(results_dir)

#Read data back in...
#Candidate regions master files (with annotations):
OTref_candidate_regions_with_annotations <- read.table(file="OTref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_candidate_regions_with_annotations <- read.table(file="PFref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#Reassigned columns...
#Original (correct):
btwn_cols <- c("xpEHH.OT_PF", "xpEHH.OT_PF.FDR","xpEHH.OT_L", "xpEHH.OT_L.FDR", "xpEHH.PF_V", "xpEHH.PF_V.FDR", "xpEHH.V_L", "xpEHH.V_L.FDR", 
               "Rsb.OT_PF", "Rsb.OT_PF.FDR","Rsb.OT_L", "Rsb.OT_L.FDR", "Rsb.PF_V", "Rsb.PF_V.FDR", "Rsb.V_L", "Rsb.V_L.FDR")
wthn_cols <- c("xpEHH.OT_V", "xpEHH.OT_V.FDR", "xpEHH.PF_L", "xpEHH.PF_L.FDR", 
               "Rsb.OT_V", "Rsb.OT_V.FDR", "Rsb.PF_L", "Rsb.PF_L.FDR")
#Reassignment 1:
btwn_cols.r1 <- c("xpEHH.OT_V", "xpEHH.OT_V.FDR","xpEHH.OT_L", "xpEHH.OT_L.FDR", "xpEHH.PF_V", "xpEHH.PF_V.FDR", "xpEHH.PF_L", "xpEHH.PF_L.FDR", 
                  "Rsb.OT_V", "Rsb.OT_V.FDR","Rsb.OT_L", "Rsb.OT_L.FDR", "Rsb.PF_V", "Rsb.PF_V.FDR", "Rsb.PF_L", "Rsb.PF_L.FDR")
wthn_cols.r1 <- c("xpEHH.OT_PF", "xpEHH.OT_PF.FDR", "xpEHH.V_L", "xpEHH.V_L.FDR", 
                  "Rsb.OT_PF", "Rsb.OT_PF.FDR", "Rsb.V_L", "Rsb.V_L.FDR")
#Reassignment 2:
btwn_cols.r2 <- c("xpEHH.OT_V", "xpEHH.OT_V.FDR","xpEHH.V_L", "xpEHH.V_L.FDR", "xpEHH.OT_PF", "xpEHH.OT_PF.FDR", "xpEHH.PF_L", "xpEHH.PF_L.FDR", 
                  "Rsb.OT_V", "Rsb.OT_V.FDR","Rsb.V_L", "Rsb.V_L.FDR", "Rsb.OT_PF", "Rsb.OT_PF.FDR", "Rsb.PF_L", "Rsb.PF_L.FDR")
wthn_cols.r2 <- c("xpEHH.PF_V", "xpEHH.PF_V.FDR", "xpEHH.OT_L", "xpEHH.OT_L.FDR", 
                  "Rsb.PF_V", "Rsb.PF_V.FDR", "Rsb.OT_L", "Rsb.OT_L.FDR")
#Reassignment 3:
btwn_cols.r3 <- c("xpEHH.OT_V", "xpEHH.OT_V.FDR","xpEHH.OT_L", "xpEHH.OT_L.FDR", "xpEHH.PF_V", "xpEHH.PF_V.FDR", "xpEHH.PF_L", "xpEHH.PF_L.FDR", 
                  "Rsb.OT_V", "Rsb.OT_V.FDR","Rsb.OT_L", "Rsb.OT_L.FDR", "Rsb.PF_V", "Rsb.PF_V.FDR", "Rsb.PF_L", "Rsb.PF_L.FDR")
wthn_cols.r3 <- c("xpEHH.OT_PF", "xpEHH.OT_PF.FDR", "xpEHH.V_L", "xpEHH.V_L.FDR", 
                  "Rsb.OT_PF", "Rsb.OT_PF.FDR", "Rsb.V_L", "Rsb.V_L.FDR")
#Reassignment 4:
btwn_cols.r4 <- c("xpEHH.OT_PF", "xpEHH.OT_PF.FDR","xpEHH.OT_V", "xpEHH.OT_V.FDR", "xpEHH.PF_L", "xpEHH.PF_L.FDR", "xpEHH.V_L", "xpEHH.V_L.FDR", 
                  "Rsb.OT_PF", "Rsb.OT_PF.FDR","Rsb.OT_V", "Rsb.OT_V.FDR", "Rsb.PF_L", "Rsb.PF_L.FDR", "Rsb.V_L", "Rsb.V_L.FDR")
wthn_cols.r4 <- c("xpEHH.OT_L", "xpEHH.OT_L.FDR", "xpEHH.PF_V", "xpEHH.PF_V.FDR", 
                  "Rsb.OT_L", "Rsb.OT_L.FDR", "Rsb.PF_V", "Rsb.PF_V.FDR")



#Function called in main function:
check.conditions <- function(df, vec, index){
  if(ncol(df[which(df[,] >= 1)]) == ncol(df)){
    vec[index] <- 1
  }
  return(vec)
}

#Main functions...
#This is messy, but the simplest way to adapt the existing function
#from above is to create a separate version for each reassignment

#Reassignment 1:
unsoc.phen.cr.test.1 <- function(cr.df){
  #Create vectors to hold some results
  is.cr.btwn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.Rsb <- rep(0, nrow(cr.df))
  #Iterate over rows in cr.df
  for(i in 1:nrow(cr.df)){
    #Create subsets
    subset.btwn.FDR <- cr.df[i,] %>% select(all_of(btwn_cols.r1)) %>% select(contains("FDR"))
    subset.wthn.FDR <- cr.df[i,] %>% select(all_of(wthn_cols.r1)) %>% select(contains("FDR"))
    subset.btwn.noFDR <- cr.df[i,] %>% select(all_of(btwn_cols.r1)) %>% select(!contains("FDR"))
    subset.wthn.noFDR <- cr.df[i,] %>% select(all_of(wthn_cols.r1)) %>% select(!contains("FDR"))
    subset.btwn.FDR.xpEHH <- subset.btwn.FDR %>% select(contains("xpEHH"))
    subset.wthn.FDR.xpEHH <- subset.wthn.FDR %>% select(contains("xpEHH"))
    subset.btwn.noFDR.xpEHH <- subset.btwn.noFDR %>% select(contains("xpEHH"))
    subset.wthn.noFDR.xpEHH <- subset.wthn.noFDR %>% select(contains("xpEHH"))
    subset.btwn.FDR.Rsb <- subset.btwn.FDR %>% select(contains("Rsb"))
    subset.wthn.FDR.Rsb <- subset.wthn.FDR %>% select(contains("Rsb"))
    subset.btwn.noFDR.Rsb <- subset.btwn.noFDR %>% select(contains("Rsb"))
    subset.wthn.noFDR.Rsb <- subset.wthn.noFDR %>% select(contains("Rsb"))
    #Change ith index of result vectors to 1 if conditions are satisfied
    is.cr.btwn.FDR.xpEHH <- check.conditions(subset.btwn.FDR.xpEHH, is.cr.btwn.FDR.xpEHH, i)
    is.cr.wthn.FDR.xpEHH <- check.conditions(subset.wthn.FDR.xpEHH, is.cr.wthn.FDR.xpEHH, i)
    is.cr.btwn.noFDR.xpEHH <- check.conditions(subset.btwn.noFDR.xpEHH, is.cr.btwn.noFDR.xpEHH, i)
    is.cr.wthn.noFDR.xpEHH <- check.conditions(subset.wthn.noFDR.xpEHH, is.cr.wthn.noFDR.xpEHH, i)
    is.cr.btwn.FDR.Rsb <- check.conditions(subset.btwn.FDR.Rsb, is.cr.btwn.FDR.Rsb, i)
    is.cr.wthn.FDR.Rsb <- check.conditions(subset.wthn.FDR.Rsb, is.cr.wthn.FDR.Rsb, i)
    is.cr.btwn.noFDR.Rsb <- check.conditions(subset.btwn.noFDR.Rsb, is.cr.btwn.noFDR.Rsb, i)
    is.cr.wthn.noFDR.Rsb <- check.conditions(subset.wthn.noFDR.Rsb, is.cr.wthn.noFDR.Rsb, i)
  }
  #Add result vectors as columns to cr.df
  cr.df$is.cr.r1.btwn.FDR.xpEHH <- is.cr.btwn.FDR.xpEHH
  cr.df$is.cr.r1.wthn.FDR.xpEHH <- is.cr.wthn.FDR.xpEHH
  cr.df$is.cr.r1.btwn.noFDR.xpEHH <- is.cr.btwn.noFDR.xpEHH
  cr.df$is.cr.r1.wthn.noFDR.xpEHH <- is.cr.wthn.noFDR.xpEHH
  cr.df$is.cr.r1.btwn.FDR.Rsb <- is.cr.btwn.FDR.Rsb
  cr.df$is.cr.r1.wthn.FDR.Rsb <- is.cr.wthn.FDR.Rsb
  cr.df$is.cr.r1.btwn.noFDR.Rsb <- is.cr.btwn.noFDR.Rsb
  cr.df$is.cr.r1.wthn.noFDR.Rsb <- is.cr.wthn.noFDR.Rsb
  #Create summary columns
  cr.df$is.unsoc.cr.r1.FDR.xpEHH <- 0
  cr.df$is.unsoc.cr.r1.noFDR.xpEHH <- 0
  cr.df$is.unsoc.cr.r1.FDR.Rsb <- 0
  cr.df$is.unsoc.cr.r1.noFDR.Rsb <- 0
  #Assign 1 to summary columns when no cr detected for within-phentoype comparisons
  #and a cr is detected for between-phenotype comparisons
  cr.df$is.unsoc.cr.r1.FDR.xpEHH[which(cr.df$is.cr.r1.wthn.FDR.xpEHH == 0 & cr.df$is.cr.r1.btwn.FDR.xpEHH == 1)] <- 1
  cr.df$is.unsoc.cr.r1.noFDR.xpEHH[which(cr.df$is.cr.r1.wthn.noFDR.xpEHH == 0 & cr.df$is.cr.r1.btwn.noFDR.xpEHH == 1)] <- 1
  cr.df$is.unsoc.cr.r1.FDR.Rsb[which(cr.df$is.cr.r1.wthn.FDR.Rsb == 0 & cr.df$is.cr.r1.btwn.FDR.Rsb == 1)] <- 1
  cr.df$is.unsoc.cr.r1.noFDR.Rsb[which(cr.df$is.cr.r1.wthn.noFDR.Rsb == 0 & cr.df$is.cr.r1.btwn.noFDR.Rsb == 1)] <- 1
  #Sort by length
  cr.df = cr.df %>% arrange(desc(LENGTH), START)
  #Return result
  return(cr.df)
}
#Reassignment 2:
unsoc.phen.cr.test.2 <- function(cr.df){
  #Create vectors to hold some results
  is.cr.btwn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.Rsb <- rep(0, nrow(cr.df))
  #Iterate over rows in cr.df
  for(i in 1:nrow(cr.df)){
    #Create subsets
    subset.btwn.FDR <- cr.df[i,] %>% select(all_of(btwn_cols.r2)) %>% select(contains("FDR"))
    subset.wthn.FDR <- cr.df[i,] %>% select(all_of(wthn_cols.r2)) %>% select(contains("FDR"))
    subset.btwn.noFDR <- cr.df[i,] %>% select(all_of(btwn_cols.r2)) %>% select(!contains("FDR"))
    subset.wthn.noFDR <- cr.df[i,] %>% select(all_of(wthn_cols.r2)) %>% select(!contains("FDR"))
    subset.btwn.FDR.xpEHH <- subset.btwn.FDR %>% select(contains("xpEHH"))
    subset.wthn.FDR.xpEHH <- subset.wthn.FDR %>% select(contains("xpEHH"))
    subset.btwn.noFDR.xpEHH <- subset.btwn.noFDR %>% select(contains("xpEHH"))
    subset.wthn.noFDR.xpEHH <- subset.wthn.noFDR %>% select(contains("xpEHH"))
    subset.btwn.FDR.Rsb <- subset.btwn.FDR %>% select(contains("Rsb"))
    subset.wthn.FDR.Rsb <- subset.wthn.FDR %>% select(contains("Rsb"))
    subset.btwn.noFDR.Rsb <- subset.btwn.noFDR %>% select(contains("Rsb"))
    subset.wthn.noFDR.Rsb <- subset.wthn.noFDR %>% select(contains("Rsb"))
    #Change ith index of result vectors to 1 if conditions are satisfied
    is.cr.btwn.FDR.xpEHH <- check.conditions(subset.btwn.FDR.xpEHH, is.cr.btwn.FDR.xpEHH, i)
    is.cr.wthn.FDR.xpEHH <- check.conditions(subset.wthn.FDR.xpEHH, is.cr.wthn.FDR.xpEHH, i)
    is.cr.btwn.noFDR.xpEHH <- check.conditions(subset.btwn.noFDR.xpEHH, is.cr.btwn.noFDR.xpEHH, i)
    is.cr.wthn.noFDR.xpEHH <- check.conditions(subset.wthn.noFDR.xpEHH, is.cr.wthn.noFDR.xpEHH, i)
    is.cr.btwn.FDR.Rsb <- check.conditions(subset.btwn.FDR.Rsb, is.cr.btwn.FDR.Rsb, i)
    is.cr.wthn.FDR.Rsb <- check.conditions(subset.wthn.FDR.Rsb, is.cr.wthn.FDR.Rsb, i)
    is.cr.btwn.noFDR.Rsb <- check.conditions(subset.btwn.noFDR.Rsb, is.cr.btwn.noFDR.Rsb, i)
    is.cr.wthn.noFDR.Rsb <- check.conditions(subset.wthn.noFDR.Rsb, is.cr.wthn.noFDR.Rsb, i)
  }
  #Add result vectors as columns to cr.df
  cr.df$is.cr.r2.btwn.FDR.xpEHH <- is.cr.btwn.FDR.xpEHH
  cr.df$is.cr.r2.wthn.FDR.xpEHH <- is.cr.wthn.FDR.xpEHH
  cr.df$is.cr.r2.btwn.noFDR.xpEHH <- is.cr.btwn.noFDR.xpEHH
  cr.df$is.cr.r2.wthn.noFDR.xpEHH <- is.cr.wthn.noFDR.xpEHH
  cr.df$is.cr.r2.btwn.FDR.Rsb <- is.cr.btwn.FDR.Rsb
  cr.df$is.cr.r2.wthn.FDR.Rsb <- is.cr.wthn.FDR.Rsb
  cr.df$is.cr.r2.btwn.noFDR.Rsb <- is.cr.btwn.noFDR.Rsb
  cr.df$is.cr.r2.wthn.noFDR.Rsb <- is.cr.wthn.noFDR.Rsb
  #Create summary columns
  cr.df$is.unsoc.cr.r2.FDR.xpEHH <- 0
  cr.df$is.unsoc.cr.r2.noFDR.xpEHH <- 0
  cr.df$is.unsoc.cr.r2.FDR.Rsb <- 0
  cr.df$is.unsoc.cr.r2.noFDR.Rsb <- 0
  #Assign 1 to summary columns when no cr detected for within-phentoype comparisons
  #and a cr is detected for between-phenotype comparisons
  cr.df$is.unsoc.cr.r2.FDR.xpEHH[which(cr.df$is.cr.r2.wthn.FDR.xpEHH == 0 & cr.df$is.cr.r2.btwn.FDR.xpEHH == 1)] <- 1
  cr.df$is.unsoc.cr.r2.noFDR.xpEHH[which(cr.df$is.cr.r2.wthn.noFDR.xpEHH == 0 & cr.df$is.cr.r2.btwn.noFDR.xpEHH == 1)] <- 1
  cr.df$is.unsoc.cr.r2.FDR.Rsb[which(cr.df$is.cr.r2.wthn.FDR.Rsb == 0 & cr.df$is.cr.r2.btwn.FDR.Rsb == 1)] <- 1
  cr.df$is.unsoc.cr.r2.noFDR.Rsb[which(cr.df$is.cr.r2.wthn.noFDR.Rsb == 0 & cr.df$is.cr.r2.btwn.noFDR.Rsb == 1)] <- 1
  #Sort by length
  cr.df = cr.df %>% arrange(desc(LENGTH), START)
  #Return result
  return(cr.df)
}
#Reassignment 3:
unsoc.phen.cr.test.3 <- function(cr.df){
  #Create vectors to hold some results
  is.cr.btwn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.Rsb <- rep(0, nrow(cr.df))
  #Iterate over rows in cr.df
  for(i in 1:nrow(cr.df)){
    #Create subsets
    subset.btwn.FDR <- cr.df[i,] %>% select(all_of(btwn_cols.r3)) %>% select(contains("FDR"))
    subset.wthn.FDR <- cr.df[i,] %>% select(all_of(wthn_cols.r3)) %>% select(contains("FDR"))
    subset.btwn.noFDR <- cr.df[i,] %>% select(all_of(btwn_cols.r3)) %>% select(!contains("FDR"))
    subset.wthn.noFDR <- cr.df[i,] %>% select(all_of(wthn_cols.r3)) %>% select(!contains("FDR"))
    subset.btwn.FDR.xpEHH <- subset.btwn.FDR %>% select(contains("xpEHH"))
    subset.wthn.FDR.xpEHH <- subset.wthn.FDR %>% select(contains("xpEHH"))
    subset.btwn.noFDR.xpEHH <- subset.btwn.noFDR %>% select(contains("xpEHH"))
    subset.wthn.noFDR.xpEHH <- subset.wthn.noFDR %>% select(contains("xpEHH"))
    subset.btwn.FDR.Rsb <- subset.btwn.FDR %>% select(contains("Rsb"))
    subset.wthn.FDR.Rsb <- subset.wthn.FDR %>% select(contains("Rsb"))
    subset.btwn.noFDR.Rsb <- subset.btwn.noFDR %>% select(contains("Rsb"))
    subset.wthn.noFDR.Rsb <- subset.wthn.noFDR %>% select(contains("Rsb"))
    #Change ith index of result vectors to 1 if conditions are satisfied
    is.cr.btwn.FDR.xpEHH <- check.conditions(subset.btwn.FDR.xpEHH, is.cr.btwn.FDR.xpEHH, i)
    is.cr.wthn.FDR.xpEHH <- check.conditions(subset.wthn.FDR.xpEHH, is.cr.wthn.FDR.xpEHH, i)
    is.cr.btwn.noFDR.xpEHH <- check.conditions(subset.btwn.noFDR.xpEHH, is.cr.btwn.noFDR.xpEHH, i)
    is.cr.wthn.noFDR.xpEHH <- check.conditions(subset.wthn.noFDR.xpEHH, is.cr.wthn.noFDR.xpEHH, i)
    is.cr.btwn.FDR.Rsb <- check.conditions(subset.btwn.FDR.Rsb, is.cr.btwn.FDR.Rsb, i)
    is.cr.wthn.FDR.Rsb <- check.conditions(subset.wthn.FDR.Rsb, is.cr.wthn.FDR.Rsb, i)
    is.cr.btwn.noFDR.Rsb <- check.conditions(subset.btwn.noFDR.Rsb, is.cr.btwn.noFDR.Rsb, i)
    is.cr.wthn.noFDR.Rsb <- check.conditions(subset.wthn.noFDR.Rsb, is.cr.wthn.noFDR.Rsb, i)
  }
  #Add result vectors as columns to cr.df
  cr.df$is.cr.r3.btwn.FDR.xpEHH <- is.cr.btwn.FDR.xpEHH
  cr.df$is.cr.r3.wthn.FDR.xpEHH <- is.cr.wthn.FDR.xpEHH
  cr.df$is.cr.r3.btwn.noFDR.xpEHH <- is.cr.btwn.noFDR.xpEHH
  cr.df$is.cr.r3.wthn.noFDR.xpEHH <- is.cr.wthn.noFDR.xpEHH
  cr.df$is.cr.r3.btwn.FDR.Rsb <- is.cr.btwn.FDR.Rsb
  cr.df$is.cr.r3.wthn.FDR.Rsb <- is.cr.wthn.FDR.Rsb
  cr.df$is.cr.r3.btwn.noFDR.Rsb <- is.cr.btwn.noFDR.Rsb
  cr.df$is.cr.r3.wthn.noFDR.Rsb <- is.cr.wthn.noFDR.Rsb
  #Create summary columns
  cr.df$is.unsoc.cr.r3.FDR.xpEHH <- 0
  cr.df$is.unsoc.cr.r3.noFDR.xpEHH <- 0
  cr.df$is.unsoc.cr.r3.FDR.Rsb <- 0
  cr.df$is.unsoc.cr.r3.noFDR.Rsb <- 0
  #Assign 1 to summary columns when no cr detected for within-phentoype comparisons
  #and a cr is detected for between-phenotype comparisons
  cr.df$is.unsoc.cr.r3.FDR.xpEHH[which(cr.df$is.cr.r3.wthn.FDR.xpEHH == 0 & cr.df$is.cr.r3.btwn.FDR.xpEHH == 1)] <- 1
  cr.df$is.unsoc.cr.r3.noFDR.xpEHH[which(cr.df$is.cr.r3.wthn.noFDR.xpEHH == 0 & cr.df$is.cr.r3.btwn.noFDR.xpEHH == 1)] <- 1
  cr.df$is.unsoc.cr.r3.FDR.Rsb[which(cr.df$is.cr.r3.wthn.FDR.Rsb == 0 & cr.df$is.cr.r3.btwn.FDR.Rsb == 1)] <- 1
  cr.df$is.unsoc.cr.r3.noFDR.Rsb[which(cr.df$is.cr.r3.wthn.noFDR.Rsb == 0 & cr.df$is.cr.r3.btwn.noFDR.Rsb == 1)] <- 1
  #Sort by length
  cr.df = cr.df %>% arrange(desc(LENGTH), START)
  #Return result
  return(cr.df)
}
#Reassignment 4:
unsoc.phen.cr.test.4 <- function(cr.df){
  #Create vectors to hold some results
  is.cr.btwn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.xpEHH <- rep(0, nrow(cr.df))
  is.cr.btwn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.FDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.btwn.noFDR.Rsb <- rep(0, nrow(cr.df))
  is.cr.wthn.noFDR.Rsb <- rep(0, nrow(cr.df))
  #Iterate over rows in cr.df
  for(i in 1:nrow(cr.df)){
    #Create subsets
    subset.btwn.FDR <- cr.df[i,] %>% select(all_of(btwn_cols.r4)) %>% select(contains("FDR"))
    subset.wthn.FDR <- cr.df[i,] %>% select(all_of(wthn_cols.r4)) %>% select(contains("FDR"))
    subset.btwn.noFDR <- cr.df[i,] %>% select(all_of(btwn_cols.r4)) %>% select(!contains("FDR"))
    subset.wthn.noFDR <- cr.df[i,] %>% select(all_of(wthn_cols.r4)) %>% select(!contains("FDR"))
    subset.btwn.FDR.xpEHH <- subset.btwn.FDR %>% select(contains("xpEHH"))
    subset.wthn.FDR.xpEHH <- subset.wthn.FDR %>% select(contains("xpEHH"))
    subset.btwn.noFDR.xpEHH <- subset.btwn.noFDR %>% select(contains("xpEHH"))
    subset.wthn.noFDR.xpEHH <- subset.wthn.noFDR %>% select(contains("xpEHH"))
    subset.btwn.FDR.Rsb <- subset.btwn.FDR %>% select(contains("Rsb"))
    subset.wthn.FDR.Rsb <- subset.wthn.FDR %>% select(contains("Rsb"))
    subset.btwn.noFDR.Rsb <- subset.btwn.noFDR %>% select(contains("Rsb"))
    subset.wthn.noFDR.Rsb <- subset.wthn.noFDR %>% select(contains("Rsb"))
    #Change ith index of result vectors to 1 if conditions are satisfied
    is.cr.btwn.FDR.xpEHH <- check.conditions(subset.btwn.FDR.xpEHH, is.cr.btwn.FDR.xpEHH, i)
    is.cr.wthn.FDR.xpEHH <- check.conditions(subset.wthn.FDR.xpEHH, is.cr.wthn.FDR.xpEHH, i)
    is.cr.btwn.noFDR.xpEHH <- check.conditions(subset.btwn.noFDR.xpEHH, is.cr.btwn.noFDR.xpEHH, i)
    is.cr.wthn.noFDR.xpEHH <- check.conditions(subset.wthn.noFDR.xpEHH, is.cr.wthn.noFDR.xpEHH, i)
    is.cr.btwn.FDR.Rsb <- check.conditions(subset.btwn.FDR.Rsb, is.cr.btwn.FDR.Rsb, i)
    is.cr.wthn.FDR.Rsb <- check.conditions(subset.wthn.FDR.Rsb, is.cr.wthn.FDR.Rsb, i)
    is.cr.btwn.noFDR.Rsb <- check.conditions(subset.btwn.noFDR.Rsb, is.cr.btwn.noFDR.Rsb, i)
    is.cr.wthn.noFDR.Rsb <- check.conditions(subset.wthn.noFDR.Rsb, is.cr.wthn.noFDR.Rsb, i)
  }
  #Add result vectors as columns to cr.df
  cr.df$is.cr.r4.btwn.FDR.xpEHH <- is.cr.btwn.FDR.xpEHH
  cr.df$is.cr.r4.wthn.FDR.xpEHH <- is.cr.wthn.FDR.xpEHH
  cr.df$is.cr.r4.btwn.noFDR.xpEHH <- is.cr.btwn.noFDR.xpEHH
  cr.df$is.cr.r4.wthn.noFDR.xpEHH <- is.cr.wthn.noFDR.xpEHH
  cr.df$is.cr.r4.btwn.FDR.Rsb <- is.cr.btwn.FDR.Rsb
  cr.df$is.cr.r4.wthn.FDR.Rsb <- is.cr.wthn.FDR.Rsb
  cr.df$is.cr.r4.btwn.noFDR.Rsb <- is.cr.btwn.noFDR.Rsb
  cr.df$is.cr.r4.wthn.noFDR.Rsb <- is.cr.wthn.noFDR.Rsb
  #Create summary columns
  cr.df$is.unsoc.cr.r4.FDR.xpEHH <- 0
  cr.df$is.unsoc.cr.r4.noFDR.xpEHH <- 0
  cr.df$is.unsoc.cr.r4.FDR.Rsb <- 0
  cr.df$is.unsoc.cr.r4.noFDR.Rsb <- 0
  #Assign 1 to summary columns when no cr detected for within-phentoype comparisons
  #and a cr is detected for between-phenotype comparisons
  cr.df$is.unsoc.cr.r4.FDR.xpEHH[which(cr.df$is.cr.r4.wthn.FDR.xpEHH == 0 & cr.df$is.cr.r4.btwn.FDR.xpEHH == 1)] <- 1
  cr.df$is.unsoc.cr.r4.noFDR.xpEHH[which(cr.df$is.cr.r4.wthn.noFDR.xpEHH == 0 & cr.df$is.cr.r4.btwn.noFDR.xpEHH == 1)] <- 1
  cr.df$is.unsoc.cr.r4.FDR.Rsb[which(cr.df$is.cr.r4.wthn.FDR.Rsb == 0 & cr.df$is.cr.r4.btwn.FDR.Rsb == 1)] <- 1
  cr.df$is.unsoc.cr.r4.noFDR.Rsb[which(cr.df$is.cr.r4.wthn.noFDR.Rsb == 0 & cr.df$is.cr.r4.btwn.noFDR.Rsb == 1)] <- 1
  #Sort by length
  cr.df = cr.df %>% arrange(desc(LENGTH), START)
  #Return result
  return(cr.df)
}

#Call the functions...
OTref_unsoc_cr <- unsoc.phen.cr.test.4(unsoc.phen.cr.test.3(unsoc.phen.cr.test.2(unsoc.phen.cr.test.1(OTref_candidate_regions_with_annotations))))
PFref_unsoc_cr <- unsoc.phen.cr.test.4(unsoc.phen.cr.test.3(unsoc.phen.cr.test.2(unsoc.phen.cr.test.1(PFref_candidate_regions_with_annotations))))
#Get counts
colSums(OTref_unsoc_cr %>% select(contains("is.unsoc.cr")))
colSums(PFref_unsoc_cr %>% select(contains("is.unsoc.cr")))
#Write to table
write.table(OTref_unsoc_cr, file="OTref_all_candidate_regions_annotations_social_phenotyp_reassignments.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(PFref_unsoc_cr, file="PFref_all_candidate_regions_annotations_social_phenotyp_reassignments.tsv", 
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(as.data.frame(colSums(OTref_unsoc_cr %>% select(contains("is.unsoc.cr")))), 
            file="OTref_candidate_regions_phenotype_reassignments_summary.tsv", 
            quote=FALSE, sep="\t", row.names=TRUE)
write.table(as.data.frame(colSums(PFref_unsoc_cr %>% select(contains("is.unsoc.cr")))), 
            file="PFref_candidate_regions_phenotype_reassignments_summary.tsv", 
            quote=FALSE, sep="\t", row.names=TRUE)

#Rough work
OTref_unsoc_cr <- read.table(file="OTref_all_candidate_regions_annotations_social_phenotyp_reassignments.tsv", 
                             header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_unsoc_cr <- read.table(file="PFref_all_candidate_regions_annotations_social_phenotyp_reassignments.tsv", 
                             header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#####################################################################################

####################### Comparison to Fst social regions ############################

#Clean up
rm(list = ls())

#Load libraries
library(tidyverse)

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
setwd(results_dir)

#Read data back in...
#Candidate regions master files (with annotations):
OTref_candidate_regions_with_annotations <- read.table(file="OTref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_candidate_regions_with_annotations <- read.table(file="PFref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Social candidate regions master files (with annotations):
OTref_social_candidate_regions_master <- read.table(file="OTref_social_candidate_regions_master.tsv", 
                                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_social_candidate_regions_master <- read.table(file="PFref_social_candidate_regions_master.tsv", 
                                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Fst data (per SNP):
OTref_Fst_outliers <- read.table(file="/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/Additional_files/Fst_data/OTref_perSNP_Fst_Just_Outliers_noNR.tsv", 
                                 header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_Fst_outliers <- read.table(file="/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/Additional_files/Fst_data/PFref_perSNP_Fst_Just_Outliers_noNR.tsv", 
                                 header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Fst data (per scaffold):
OTref_Fst_scaffolds <- read.table(file="/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/Additional_files/Fst_data/OTref_per_scaffold_comparison_final.tsv", 
                                  header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_Fst_scaffolds <- read.table(file="/scratch/genomicsocorg/mwhj1/VCF_phasing/Data/Additional_files/Fst_data/PFref_per_scaffold_comparison_final.tsv", 
                                  header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Reformat CHR column in rehh candidate regions data frames from integer to character
#(to match the Fst outliers CHROM column)
OTref_candidate_regions_with_annotations$CHR <- as.character(OTref_candidate_regions_with_annotations$CHR)
PFref_candidate_regions_with_annotations$CHR <- as.character(PFref_candidate_regions_with_annotations$CHR)
OTref_social_candidate_regions_master$CHR <- as.character(OTref_social_candidate_regions_master$CHR)
PFref_social_candidate_regions_master$CHR <- as.character(PFref_social_candidate_regions_master$CHR)

#Rough work
length(unique(OTref_social_candidate_regions_master$CHR))
length(unique(OTref_social_candidate_regions_master$CHR[which(OTref_social_candidate_regions_master$CHR %in% OTref_Fst_outliers$CHROM)]))
length(unique(PFref_social_candidate_regions_master$CHR))
length(unique(PFref_social_candidate_regions_master$CHR[which(PFref_social_candidate_regions_master$CHR %in% PFref_Fst_outliers$CHROM)]))

length(unique(OTref_Fst_outliers$CHROM))
length(unique(OTref_Fst_outliers$CHROM[which(OTref_Fst_outliers$CHROM %in% OTref_social_candidate_regions_master$CHR)]))
length(unique(PFref_Fst_outliers$CHROM))
length(unique(PFref_Fst_outliers$CHROM[which(PFref_Fst_outliers$CHROM %in% PFref_social_candidate_regions_master$CHR)]))

OTref_Fst_scaffolds.filtered <- OTref_Fst_scaffolds %>% 
  mutate(present.in.rehh = ifelse(CHROM %in% OTref_social_candidate_regions_master$CHR, 1, 0)) %>% 
  filter(LENGTH >= 100000) %>% 
  select(CHROM, LENGTH, n.markers.without.NR, present.in.rehh)

PFref_Fst_scaffolds.filtered <- PFref_Fst_scaffolds %>% 
  mutate(present.in.rehh = ifelse(CHROM %in% PFref_social_candidate_regions_master$CHR, 1, 0)) %>% 
  filter(LENGTH >= 100000) %>% 
  select(CHROM, LENGTH, n.markers.without.NR, present.in.rehh)

OTref_Fst_scaffolds_top10 <- OTref_Fst_scaffolds.filtered %>% top_n(10,n.markers.without.NR)
PFref_Fst_scaffolds_top10 <- PFref_Fst_scaffolds.filtered %>% top_n(10,n.markers.without.NR)

OTref_FDR_social_candidates <- OTref_social_candidate_regions_master %>% 
  select(CHR, START, END, LENGTH, contains("is.soc.cr.FDR")) %>% 
  filter(is.soc.cr.FDR.xpEHH == 1 | is.soc.cr.FDR.Rsb == 1)




########################### Detailed Plots of Candidate Regions ###########################

#Clean up
rm(list = ls())

#Load libraries
library(rehh)
library(tidyverse)
library(rlist)
library(gridExtra)
library(ggpubr)

#Set directories
results_dir <- "/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis"
setwd(results_dir)

#Read data back in (only FDR-corrected data)

#iHS (OTref):
PF_iHS_OTref.FDR <- read.table(file="PF.ALL.OTref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
L_iHS_OTref.FDR <- read.table(file="L.ALL.OTref.iHS", 
                              header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_iHS_OTref.FDR <- read.table(file="OT.ALL.OTref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_iHS_OTref.FDR <- read.table(file="V.ALL.OTref.iHS", 
                              header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_iHS_OTref.FDR <- read.table(file="NR.ALL.OTref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#iHS (PFref):
PF_iHS_PFref.FDR <- read.table(file="PF.ALL.PFref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
L_iHS_PFref.FDR <- read.table(file="L.ALL.PFref.iHS", 
                              header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_iHS_PFref.FDR <- read.table(file="OT.ALL.PFref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_iHS_PFref.FDR <- read.table(file="V.ALL.PFref.iHS", 
                              header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_iHS_PFref.FDR <- read.table(file="NR.ALL.PFref.iHS", 
                               header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#xpEHH (OTref):
OT_NR_ALL_OTref.xpEHH.FDR <- read.table(file="OT_NR_ALL_OTref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_OTref.xpEHH.FDR <- read.table(file="OT_PF_ALL_OTref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_OTref.xpEHH.FDR <- read.table(file="OT_V_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_OTref.xpEHH.FDR <- read.table(file="OT_L_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_OTref.xpEHH.FDR <- read.table(file="NR_PF_ALL_OTref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_OTref.xpEHH.FDR <- read.table(file="NR_V_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_OTref.xpEHH.FDR <- read.table(file="NR_L_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_OTref.xpEHH.FDR <- read.table(file="PF_V_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_OTref.xpEHH.FDR <- read.table(file="PF_L_ALL_OTref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_OTref.xpEHH.FDR <- read.table(file="V_L_ALL_OTref.xpEHH", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#xpEHH (PFref):
OT_NR_ALL_PFref.xpEHH.FDR <- read.table(file="OT_NR_ALL_PFref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_PFref.xpEHH.FDR <- read.table(file="OT_PF_ALL_PFref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_PFref.xpEHH.FDR <- read.table(file="OT_V_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_PFref.xpEHH.FDR <- read.table(file="OT_L_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_PFref.xpEHH.FDR <- read.table(file="NR_PF_ALL_PFref.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_PFref.xpEHH.FDR <- read.table(file="NR_V_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_PFref.xpEHH.FDR <- read.table(file="NR_L_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_PFref.xpEHH.FDR <- read.table(file="PF_V_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_PFref.xpEHH.FDR <- read.table(file="PF_L_ALL_PFref.xpEHH", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_PFref.xpEHH.FDR <- read.table(file="V_L_ALL_PFref.xpEHH", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Rsb (OTref):
OT_NR_ALL_OTref.Rsb.FDR <- read.table(file="OT_NR_ALL_OTref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_OTref.Rsb.FDR <- read.table(file="OT_PF_ALL_OTref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_OTref.Rsb.FDR <- read.table(file="OT_V_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_OTref.Rsb.FDR <- read.table(file="OT_L_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_OTref.Rsb.FDR <- read.table(file="NR_PF_ALL_OTref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_OTref.Rsb.FDR <- read.table(file="NR_V_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_OTref.Rsb.FDR <- read.table(file="NR_L_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_OTref.Rsb.FDR <- read.table(file="PF_V_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_OTref.Rsb.FDR <- read.table(file="PF_L_ALL_OTref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_OTref.Rsb.FDR <- read.table(file="V_L_ALL_OTref.Rsb", 
                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Rsb (PFref):
OT_NR_ALL_PFref.Rsb.FDR <- read.table(file="OT_NR_ALL_PFref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_PFref.Rsb.FDR <- read.table(file="OT_PF_ALL_PFref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_PFref.Rsb.FDR <- read.table(file="OT_V_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_PFref.Rsb.FDR <- read.table(file="OT_L_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_PFref.Rsb.FDR <- read.table(file="NR_PF_ALL_PFref.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_PFref.Rsb.FDR <- read.table(file="NR_V_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_PFref.Rsb.FDR <- read.table(file="NR_L_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_PFref.Rsb.FDR <- read.table(file="PF_V_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_PFref.Rsb.FDR <- read.table(file="PF_L_ALL_PFref.Rsb", 
                                     header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_PFref.Rsb.FDR <- read.table(file="V_L_ALL_PFref.Rsb", 
                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#FDR-un-corrected data...

#iHS (OTref):
PF_iHS_OTref.noFDR <- read.table(file="PF.ALL.OTref.noPvalueCorrection.iHS", 
                                 header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
L_iHS_OTref.noFDR <- read.table(file="L.ALL.OTref.noPvalueCorrection.iHS", 
                                header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_iHS_OTref.noFDR <- read.table(file="OT.ALL.OTref.noPvalueCorrection.iHS", 
                                 header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_iHS_OTref.noFDR <- read.table(file="V.ALL.OTref.noPvalueCorrection.iHS", 
                                header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_iHS_OTref.noFDR <- read.table(file="NR.ALL.OTref.noPvalueCorrection.iHS", 
                                 header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#iHS (PFref):
PF_iHS_PFref.noFDR <- read.table(file="PF.ALL.PFref.noPvalueCorrection.iHS", 
                                 header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
L_iHS_PFref.noFDR <- read.table(file="L.ALL.PFref.noPvalueCorrection.iHS", 
                                header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_iHS_PFref.noFDR <- read.table(file="OT.ALL.PFref.noPvalueCorrection.iHS", 
                                 header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_iHS_PFref.noFDR <- read.table(file="V.ALL.PFref.noPvalueCorrection.iHS", 
                                header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_iHS_PFref.noFDR <- read.table(file="NR.ALL.PFref.noPvalueCorrection.iHS", 
                                 header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#xpEHH (OTref):
OT_NR_ALL_OTref.xpEHH.noFDR <- read.table(file="OT_NR_ALL_OTref.noPvalueCorrection.xpEHH", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_OTref.xpEHH.noFDR <- read.table(file="OT_PF_ALL_OTref.noPvalueCorrection.xpEHH", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_OTref.xpEHH.noFDR <- read.table(file="OT_V_ALL_OTref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_OTref.xpEHH.noFDR <- read.table(file="OT_L_ALL_OTref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_OTref.xpEHH.noFDR <- read.table(file="NR_PF_ALL_OTref.noPvalueCorrection.xpEHH", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_OTref.xpEHH.noFDR <- read.table(file="NR_V_ALL_OTref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_OTref.xpEHH.noFDR <- read.table(file="NR_L_ALL_OTref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_OTref.xpEHH.noFDR <- read.table(file="PF_V_ALL_OTref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_OTref.xpEHH.noFDR <- read.table(file="PF_L_ALL_OTref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_OTref.xpEHH.noFDR <- read.table(file="V_L_ALL_OTref.noPvalueCorrection.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#xpEHH (PFref):
OT_NR_ALL_PFref.xpEHH.noFDR <- read.table(file="OT_NR_ALL_PFref.noPvalueCorrection.xpEHH", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_PFref.xpEHH.noFDR <- read.table(file="OT_PF_ALL_PFref.noPvalueCorrection.xpEHH", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_PFref.xpEHH.noFDR <- read.table(file="OT_V_ALL_PFref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_PFref.xpEHH.noFDR <- read.table(file="OT_L_ALL_PFref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_PFref.xpEHH.noFDR <- read.table(file="NR_PF_ALL_PFref.noPvalueCorrection.xpEHH", 
                                          header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_PFref.xpEHH.noFDR <- read.table(file="NR_V_ALL_PFref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_PFref.xpEHH.noFDR <- read.table(file="NR_L_ALL_PFref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_PFref.xpEHH.noFDR <- read.table(file="PF_V_ALL_PFref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_PFref.xpEHH.noFDR <- read.table(file="PF_L_ALL_PFref.noPvalueCorrection.xpEHH", 
                                         header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_PFref.xpEHH.noFDR <- read.table(file="V_L_ALL_PFref.noPvalueCorrection.xpEHH", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Rsb (OTref):
OT_NR_ALL_OTref.Rsb.noFDR <- read.table(file="OT_NR_ALL_OTref.noPvalueCorrection.Rsb", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_OTref.Rsb.noFDR <- read.table(file="OT_PF_ALL_OTref.noPvalueCorrection.Rsb", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_OTref.Rsb.noFDR <- read.table(file="OT_V_ALL_OTref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_OTref.Rsb.noFDR <- read.table(file="OT_L_ALL_OTref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_OTref.Rsb.noFDR <- read.table(file="NR_PF_ALL_OTref.noPvalueCorrection.Rsb", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_OTref.Rsb.noFDR <- read.table(file="NR_V_ALL_OTref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_OTref.Rsb.noFDR <- read.table(file="NR_L_ALL_OTref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_OTref.Rsb.noFDR <- read.table(file="PF_V_ALL_OTref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_OTref.Rsb.noFDR <- read.table(file="PF_L_ALL_OTref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_OTref.Rsb.noFDR <- read.table(file="V_L_ALL_OTref.noPvalueCorrection.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Rsb (PFref):
OT_NR_ALL_PFref.Rsb.noFDR <- read.table(file="OT_NR_ALL_PFref.noPvalueCorrection.Rsb", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_PF_ALL_PFref.Rsb.noFDR <- read.table(file="OT_PF_ALL_PFref.noPvalueCorrection.Rsb", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_V_ALL_PFref.Rsb.noFDR <- read.table(file="OT_V_ALL_PFref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
OT_L_ALL_PFref.Rsb.noFDR <- read.table(file="OT_L_ALL_PFref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_PF_ALL_PFref.Rsb.noFDR <- read.table(file="NR_PF_ALL_PFref.noPvalueCorrection.Rsb", 
                                        header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_V_ALL_PFref.Rsb.noFDR <- read.table(file="NR_V_ALL_PFref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
NR_L_ALL_PFref.Rsb.noFDR <- read.table(file="NR_L_ALL_PFref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_V_ALL_PFref.Rsb.noFDR <- read.table(file="PF_V_ALL_PFref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PF_L_ALL_PFref.Rsb.noFDR <- read.table(file="PF_L_ALL_PFref.noPvalueCorrection.Rsb", 
                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
V_L_ALL_PFref.Rsb.noFDR <- read.table(file="V_L_ALL_PFref.noPvalueCorrection.Rsb", 
                                      header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#Social candidate regions (any test, corrected and un-corrected):
OTref_social_candidate_regions_master <- read.table(file="OTref_social_candidate_regions_master.tsv", 
                                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_social_candidate_regions_master <- read.table(file="PFref_social_candidate_regions_master.tsv", 
                                                    header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
#Candidate regions master files (with annotations):
OTref_candidate_regions_with_annotations <- read.table(file="OTref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
PFref_candidate_regions_with_annotations <- read.table(file="PFref_all_candidate_regions_annotations.tsv", 
                                                       header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)

#One scaffold shows up as containing the majority of characterised genes in social candidate regions:
#OT 62955
#PF 144
#These need to be plotted - both in their entirety and across the subset of their range which contains 
#social candidate regions

#Full-scaffold plotting function (xpEHH and Rsb combined)
#(requires input in the exact order specified unless used with 'argument=data' format)
#*Important - specify plot_dir without a trailing "/":
plot.xpop.scaffold.logP <- function(scaffold.ID, 
                                    filename, 
                                    plot_dir, 
                                    OT_V_xpEHH, 
                                    PF_L_xpEHH, 
                                    OT_V_Rsb, 
                                    PF_L_Rsb, 
                                    OT_PF_xpEHH, 
                                    OT_L_xpEHH, 
                                    OT_PF_Rsb, 
                                    OT_L_Rsb, 
                                    PF_V_xpEHH, 
                                    V_L_xpEHH, 
                                    PF_V_Rsb, 
                                    V_L_Rsb, 
                                    OT_NR_xpEHH, 
                                    NR_V_xpEHH, 
                                    OT_NR_Rsb, 
                                    NR_V_Rsb, 
                                    NR_PF_xpEHH, 
                                    NR_L_xpEHH, 
                                    NR_PF_Rsb, 
                                    NR_L_Rsb){
  #Subset data to matching scaffold
  OT_V_xpEHH <- OT_V_xpEHH[which(OT_V_xpEHH$CHR == scaffold.ID), ]
  PF_L_xpEHH <- PF_L_xpEHH[which(PF_L_xpEHH$CHR == scaffold.ID), ]
  OT_V_Rsb <- OT_V_Rsb[which(OT_V_Rsb$CHR == scaffold.ID), ]
  PF_L_Rsb <- PF_L_Rsb[which(PF_L_Rsb$CHR == scaffold.ID), ]
  OT_PF_xpEHH <- OT_PF_xpEHH[which(OT_PF_xpEHH$CHR == scaffold.ID), ]
  OT_L_xpEHH <- OT_L_xpEHH[which(OT_L_xpEHH$CHR == scaffold.ID), ]
  OT_PF_Rsb <- OT_PF_Rsb[which(OT_PF_Rsb$CHR == scaffold.ID), ]
  OT_L_Rsb <- OT_L_Rsb[which(OT_L_Rsb$CHR == scaffold.ID), ]
  PF_V_xpEHH <- PF_V_xpEHH[which(PF_V_xpEHH$CHR == scaffold.ID), ]
  V_L_xpEHH <- V_L_xpEHH[which(V_L_xpEHH$CHR == scaffold.ID), ]
  PF_V_Rsb <- PF_V_Rsb[which(PF_V_Rsb$CHR == scaffold.ID), ]
  V_L_Rsb <- V_L_Rsb[which(V_L_Rsb$CHR == scaffold.ID), ]
  OT_NR_xpEHH <- OT_NR_xpEHH[which(OT_NR_xpEHH$CHR == scaffold.ID), ]
  NR_V_xpEHH <- NR_V_xpEHH[which(NR_V_xpEHH$CHR == scaffold.ID), ]
  OT_NR_Rsb <- OT_NR_Rsb[which(OT_NR_Rsb$CHR == scaffold.ID), ]
  NR_V_Rsb <- NR_V_Rsb[which(NR_V_Rsb$CHR == scaffold.ID), ]
  NR_PF_xpEHH <- NR_PF_xpEHH[which(NR_PF_xpEHH$CHR == scaffold.ID), ]
  NR_L_xpEHH <- NR_L_xpEHH[which(NR_L_xpEHH$CHR == scaffold.ID), ]
  NR_PF_Rsb <- NR_PF_Rsb[which(NR_PF_Rsb$CHR == scaffold.ID), ]
  NR_L_Rsb <- NR_L_Rsb[which(NR_L_Rsb$CHR == scaffold.ID), ]
  #Get max logP value for xpEHH and Rsb
  logPmax_xpEHH <- ceiling(max(c(OT_V_xpEHH$LOGPVALUE, 
                                 PF_L_xpEHH$LOGPVALUE, 
                                 OT_PF_xpEHH$LOGPVALUE, 
                                 OT_L_xpEHH$LOGPVALUE, 
                                 PF_V_xpEHH$LOGPVALUE, 
                                 V_L_xpEHH$LOGPVALUE, 
                                 OT_NR_xpEHH$LOGPVALUE, 
                                 NR_V_xpEHH$LOGPVALUE, 
                                 NR_PF_xpEHH$LOGPVALUE, 
                                 NR_L_xpEHH$LOGPVALUE), 
                               na.rm=TRUE)) + 1
  logPmax_Rsb <- ceiling(max(c(OT_V_Rsb$LOGPVALUE, 
                               PF_L_Rsb$LOGPVALUE, 
                               OT_PF_Rsb$LOGPVALUE, 
                               OT_L_Rsb$LOGPVALUE, 
                               PF_V_Rsb$LOGPVALUE, 
                               V_L_Rsb$LOGPVALUE, 
                               OT_NR_Rsb$LOGPVALUE, 
                               NR_V_Rsb$LOGPVALUE, 
                               NR_PF_Rsb$LOGPVALUE, 
                               NR_L_Rsb$LOGPVALUE), 
                             na.rm=TRUE)) + 1
  if(logPmax_xpEHH > logPmax_Rsb){
    logPmax_Rsb <- logPmax_xpEHH 
  }
  if(logPmax_Rsb > logPmax_xpEHH){
    logPmax_xpEHH <- logPmax_Rsb
  }
  #Create plots...
  #Page 1:
  p1.plot1 <- ggplot(OT_V_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_V xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p1.plot2 <- ggplot(PF_L_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("PF_L xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p1.plot3 <- ggplot(OT_V_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_V Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p1.plot4 <- ggplot(PF_L_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("PF_L Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Page 2:
  p2.plot1 <- ggplot(OT_PF_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_PF xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p2.plot2 <- ggplot(OT_L_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_L xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p2.plot3 <- ggplot(OT_PF_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_PF Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p2.plot4 <- ggplot(OT_L_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_L Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Page 3:
  p3.plot1 <- ggplot(PF_V_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("PF_V xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p3.plot2 <- ggplot(V_L_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("V_L xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p3.plot3 <- ggplot(PF_V_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("PF_V Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p3.plot4 <- ggplot(V_L_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("V_L Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Page 4:
  p4.plot1 <- ggplot(OT_NR_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_NR xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p4.plot2 <- ggplot(NR_V_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_V xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p4.plot3 <- ggplot(OT_NR_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_NR Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p4.plot4 <- ggplot(NR_V_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_V Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Page 5:
  p5.plot1 <- ggplot(NR_PF_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_PF xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p5.plot2 <- ggplot(NR_L_xpEHH, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_L xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p5.plot3 <- ggplot(NR_PF_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_PF Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p5.plot4 <- ggplot(NR_L_Rsb, aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_L Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Arrange on pages:
  pagelist <- vector("list", 5)
  pagelist[[1]] <- ggarrange(p1.plot1, p1.plot2, p1.plot3, p1.plot4, 
                             ncol = 2, nrow = 2)
  pagelist[[2]] <- ggarrange(p2.plot1, p2.plot2, p2.plot3, p2.plot4, 
                             ncol = 2, nrow = 2)
  pagelist[[3]] <- ggarrange(p3.plot1, p3.plot2, p3.plot3, p3.plot4, 
                             ncol = 2, nrow = 2)
  pagelist[[4]] <- ggarrange(p4.plot1, p4.plot2, p4.plot3, p4.plot4, 
                             ncol = 2, nrow = 2)
  pagelist[[5]] <- ggarrange(p5.plot1, p5.plot2, p5.plot3, p5.plot4, 
                             ncol = 2, nrow = 2)
  #Save to file:
  ggsave(paste(plot_dir, filename, sep="/"), marrangeGrob(grobs=pagelist, nrow = 1, ncol = 1, list(top=NULL)), width=11, height=8.5)
}

#Plot and save OT 62955...
#No FDR:
plot.xpop.scaffold.logP(scaffold.ID=62955, 
                        filename="OT_62955_xpEHH_and_Rsb_logP_noFDR.pdf", 
                        plot_dir="/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Single_scaffold_plots", 
                        OT_V_xpEHH=OT_V_ALL_OTref.xpEHH.noFDR, 
                        PF_L_xpEHH=PF_L_ALL_OTref.xpEHH.noFDR, 
                        OT_V_Rsb=OT_V_ALL_OTref.Rsb.noFDR, 
                        PF_L_Rsb=PF_L_ALL_OTref.Rsb.noFDR, 
                        OT_PF_xpEHH=OT_PF_ALL_OTref.xpEHH.noFDR, 
                        OT_L_xpEHH=OT_L_ALL_OTref.xpEHH.noFDR, 
                        OT_PF_Rsb=OT_PF_ALL_OTref.Rsb.noFDR, 
                        OT_L_Rsb=OT_L_ALL_OTref.Rsb.noFDR, 
                        PF_V_xpEHH=PF_V_ALL_OTref.xpEHH.noFDR, 
                        V_L_xpEHH=V_L_ALL_OTref.xpEHH.noFDR, 
                        PF_V_Rsb=PF_V_ALL_OTref.Rsb.noFDR, 
                        V_L_Rsb=V_L_ALL_OTref.Rsb.noFDR, 
                        OT_NR_xpEHH=OT_NR_ALL_OTref.xpEHH.noFDR, 
                        NR_V_xpEHH=NR_V_ALL_OTref.xpEHH.noFDR, 
                        OT_NR_Rsb=OT_NR_ALL_OTref.Rsb.noFDR, 
                        NR_V_Rsb=NR_V_ALL_OTref.Rsb.noFDR, 
                        NR_PF_xpEHH=NR_PF_ALL_OTref.xpEHH.noFDR, 
                        NR_L_xpEHH=NR_L_ALL_OTref.xpEHH.noFDR, 
                        NR_PF_Rsb=NR_PF_ALL_OTref.Rsb.noFDR, 
                        NR_L_Rsb=NR_L_ALL_OTref.Rsb.noFDR)

#With FDR:
plot.xpop.scaffold.logP(scaffold.ID=62955, 
                        filename="OT_62955_xpEHH_and_Rsb_logP_withFDR.pdf", 
                        plot_dir="/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Single_scaffold_plots", 
                        OT_V_xpEHH=OT_V_ALL_OTref.xpEHH.FDR, 
                        PF_L_xpEHH=PF_L_ALL_OTref.xpEHH.FDR, 
                        OT_V_Rsb=OT_V_ALL_OTref.Rsb.FDR, 
                        PF_L_Rsb=PF_L_ALL_OTref.Rsb.FDR, 
                        OT_PF_xpEHH=OT_PF_ALL_OTref.xpEHH.FDR, 
                        OT_L_xpEHH=OT_L_ALL_OTref.xpEHH.FDR, 
                        OT_PF_Rsb=OT_PF_ALL_OTref.Rsb.FDR, 
                        OT_L_Rsb=OT_L_ALL_OTref.Rsb.FDR, 
                        PF_V_xpEHH=PF_V_ALL_OTref.xpEHH.FDR, 
                        V_L_xpEHH=V_L_ALL_OTref.xpEHH.FDR, 
                        PF_V_Rsb=PF_V_ALL_OTref.Rsb.FDR, 
                        V_L_Rsb=V_L_ALL_OTref.Rsb.FDR, 
                        OT_NR_xpEHH=OT_NR_ALL_OTref.xpEHH.FDR, 
                        NR_V_xpEHH=NR_V_ALL_OTref.xpEHH.FDR, 
                        OT_NR_Rsb=OT_NR_ALL_OTref.Rsb.FDR, 
                        NR_V_Rsb=NR_V_ALL_OTref.Rsb.FDR, 
                        NR_PF_xpEHH=NR_PF_ALL_OTref.xpEHH.FDR, 
                        NR_L_xpEHH=NR_L_ALL_OTref.xpEHH.FDR, 
                        NR_PF_Rsb=NR_PF_ALL_OTref.Rsb.FDR, 
                        NR_L_Rsb=NR_L_ALL_OTref.Rsb.FDR)

#Plot and save PF 144...
#No FDR:
plot.xpop.scaffold.logP(scaffold.ID=144, 
                        filename="PF_144_xpEHH_and_Rsb_logP_noFDR.pdf", 
                        plot_dir="/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Single_scaffold_plots", 
                        OT_V_xpEHH=OT_V_ALL_PFref.xpEHH.noFDR, 
                        PF_L_xpEHH=PF_L_ALL_PFref.xpEHH.noFDR, 
                        OT_V_Rsb=OT_V_ALL_PFref.Rsb.noFDR, 
                        PF_L_Rsb=PF_L_ALL_PFref.Rsb.noFDR, 
                        OT_PF_xpEHH=OT_PF_ALL_PFref.xpEHH.noFDR, 
                        OT_L_xpEHH=OT_L_ALL_PFref.xpEHH.noFDR, 
                        OT_PF_Rsb=OT_PF_ALL_PFref.Rsb.noFDR, 
                        OT_L_Rsb=OT_L_ALL_PFref.Rsb.noFDR, 
                        PF_V_xpEHH=PF_V_ALL_PFref.xpEHH.noFDR, 
                        V_L_xpEHH=V_L_ALL_PFref.xpEHH.noFDR, 
                        PF_V_Rsb=PF_V_ALL_PFref.Rsb.noFDR, 
                        V_L_Rsb=V_L_ALL_PFref.Rsb.noFDR, 
                        OT_NR_xpEHH=OT_NR_ALL_PFref.xpEHH.noFDR, 
                        NR_V_xpEHH=NR_V_ALL_PFref.xpEHH.noFDR, 
                        OT_NR_Rsb=OT_NR_ALL_PFref.Rsb.noFDR, 
                        NR_V_Rsb=NR_V_ALL_PFref.Rsb.noFDR, 
                        NR_PF_xpEHH=NR_PF_ALL_PFref.xpEHH.noFDR, 
                        NR_L_xpEHH=NR_L_ALL_PFref.xpEHH.noFDR, 
                        NR_PF_Rsb=NR_PF_ALL_PFref.Rsb.noFDR, 
                        NR_L_Rsb=NR_L_ALL_PFref.Rsb.noFDR)

#With FDR:
plot.xpop.scaffold.logP(scaffold.ID=144, 
                        filename="PF_144_xpEHH_and_Rsb_logP_withFDR.pdf", 
                        plot_dir="/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Single_scaffold_plots", 
                        OT_V_xpEHH=OT_V_ALL_PFref.xpEHH.FDR, 
                        PF_L_xpEHH=PF_L_ALL_PFref.xpEHH.FDR, 
                        OT_V_Rsb=OT_V_ALL_PFref.Rsb.FDR, 
                        PF_L_Rsb=PF_L_ALL_PFref.Rsb.FDR, 
                        OT_PF_xpEHH=OT_PF_ALL_PFref.xpEHH.FDR, 
                        OT_L_xpEHH=OT_L_ALL_PFref.xpEHH.FDR, 
                        OT_PF_Rsb=OT_PF_ALL_PFref.Rsb.FDR, 
                        OT_L_Rsb=OT_L_ALL_PFref.Rsb.FDR, 
                        PF_V_xpEHH=PF_V_ALL_PFref.xpEHH.FDR, 
                        V_L_xpEHH=V_L_ALL_PFref.xpEHH.FDR, 
                        PF_V_Rsb=PF_V_ALL_PFref.Rsb.FDR, 
                        V_L_Rsb=V_L_ALL_PFref.Rsb.FDR, 
                        OT_NR_xpEHH=OT_NR_ALL_PFref.xpEHH.FDR, 
                        NR_V_xpEHH=NR_V_ALL_PFref.xpEHH.FDR, 
                        OT_NR_Rsb=OT_NR_ALL_PFref.Rsb.FDR, 
                        NR_V_Rsb=NR_V_ALL_PFref.Rsb.FDR, 
                        NR_PF_xpEHH=NR_PF_ALL_PFref.xpEHH.FDR, 
                        NR_L_xpEHH=NR_L_ALL_PFref.xpEHH.FDR, 
                        NR_PF_Rsb=NR_PF_ALL_PFref.Rsb.FDR, 
                        NR_L_Rsb=NR_L_ALL_PFref.Rsb.FDR)





#Scaffold-specific social candidate region plotting function (xpEHH and Rsb combined)
#(requires input in the exact order specified unless used with 'argument=data' format)
#*Important - specify plot_dir without a trailing "/":
plot.xpop.social.region.logP <- function(scaffold.ID, 
                                         filename, 
                                         plot_dir, 
                                         social_candidate_regions, 
                                         OT_V_xpEHH, 
                                         PF_L_xpEHH, 
                                         OT_V_Rsb, 
                                         PF_L_Rsb, 
                                         OT_PF_xpEHH, 
                                         OT_L_xpEHH, 
                                         OT_PF_Rsb, 
                                         OT_L_Rsb, 
                                         PF_V_xpEHH, 
                                         V_L_xpEHH, 
                                         PF_V_Rsb, 
                                         V_L_Rsb, 
                                         OT_NR_xpEHH, 
                                         NR_V_xpEHH, 
                                         OT_NR_Rsb, 
                                         NR_V_Rsb, 
                                         NR_PF_xpEHH, 
                                         NR_L_xpEHH, 
                                         NR_PF_Rsb, 
                                         NR_L_Rsb){
  #Subset social candidate regions to matching scaffold
  social_candidate_regions <- social_candidate_regions[which(social_candidate_regions$CHR == scaffold.ID), ]
  #Get start position of first social candidate window (subtract 50kb as a buffer)
  social_start <- min(social_candidate_regions$START) - 50000
  #Get end position of last social candidate window (add 50kb as a buffer)
  social_end <- max(social_candidate_regions$END) + 50000
  #Subset data to matching scaffold
  OT_V_xpEHH <- OT_V_xpEHH[which(OT_V_xpEHH$CHR == scaffold.ID), ]
  PF_L_xpEHH <- PF_L_xpEHH[which(PF_L_xpEHH$CHR == scaffold.ID), ]
  OT_V_Rsb <- OT_V_Rsb[which(OT_V_Rsb$CHR == scaffold.ID), ]
  PF_L_Rsb <- PF_L_Rsb[which(PF_L_Rsb$CHR == scaffold.ID), ]
  OT_PF_xpEHH <- OT_PF_xpEHH[which(OT_PF_xpEHH$CHR == scaffold.ID), ]
  OT_L_xpEHH <- OT_L_xpEHH[which(OT_L_xpEHH$CHR == scaffold.ID), ]
  OT_PF_Rsb <- OT_PF_Rsb[which(OT_PF_Rsb$CHR == scaffold.ID), ]
  OT_L_Rsb <- OT_L_Rsb[which(OT_L_Rsb$CHR == scaffold.ID), ]
  PF_V_xpEHH <- PF_V_xpEHH[which(PF_V_xpEHH$CHR == scaffold.ID), ]
  V_L_xpEHH <- V_L_xpEHH[which(V_L_xpEHH$CHR == scaffold.ID), ]
  PF_V_Rsb <- PF_V_Rsb[which(PF_V_Rsb$CHR == scaffold.ID), ]
  V_L_Rsb <- V_L_Rsb[which(V_L_Rsb$CHR == scaffold.ID), ]
  OT_NR_xpEHH <- OT_NR_xpEHH[which(OT_NR_xpEHH$CHR == scaffold.ID), ]
  NR_V_xpEHH <- NR_V_xpEHH[which(NR_V_xpEHH$CHR == scaffold.ID), ]
  OT_NR_Rsb <- OT_NR_Rsb[which(OT_NR_Rsb$CHR == scaffold.ID), ]
  NR_V_Rsb <- NR_V_Rsb[which(NR_V_Rsb$CHR == scaffold.ID), ]
  NR_PF_xpEHH <- NR_PF_xpEHH[which(NR_PF_xpEHH$CHR == scaffold.ID), ]
  NR_L_xpEHH <- NR_L_xpEHH[which(NR_L_xpEHH$CHR == scaffold.ID), ]
  NR_PF_Rsb <- NR_PF_Rsb[which(NR_PF_Rsb$CHR == scaffold.ID), ]
  NR_L_Rsb <- NR_L_Rsb[which(NR_L_Rsb$CHR == scaffold.ID), ]
  #Get max logP value for xpEHH and Rsb
  logPmax_xpEHH <- ceiling(max(c(OT_V_xpEHH$LOGPVALUE, 
                                 PF_L_xpEHH$LOGPVALUE, 
                                 OT_PF_xpEHH$LOGPVALUE, 
                                 OT_L_xpEHH$LOGPVALUE, 
                                 PF_V_xpEHH$LOGPVALUE, 
                                 V_L_xpEHH$LOGPVALUE, 
                                 OT_NR_xpEHH$LOGPVALUE, 
                                 NR_V_xpEHH$LOGPVALUE, 
                                 NR_PF_xpEHH$LOGPVALUE, 
                                 NR_L_xpEHH$LOGPVALUE), 
                               na.rm=TRUE)) + 1
  logPmax_Rsb <- ceiling(max(c(OT_V_Rsb$LOGPVALUE, 
                               PF_L_Rsb$LOGPVALUE, 
                               OT_PF_Rsb$LOGPVALUE, 
                               OT_L_Rsb$LOGPVALUE, 
                               PF_V_Rsb$LOGPVALUE, 
                               V_L_Rsb$LOGPVALUE, 
                               OT_NR_Rsb$LOGPVALUE, 
                               NR_V_Rsb$LOGPVALUE, 
                               NR_PF_Rsb$LOGPVALUE, 
                               NR_L_Rsb$LOGPVALUE), 
                             na.rm=TRUE)) + 1
  if(logPmax_xpEHH > logPmax_Rsb){
    logPmax_Rsb <- logPmax_xpEHH 
  }
  if(logPmax_Rsb > logPmax_xpEHH){
    logPmax_xpEHH <- logPmax_Rsb
  }
  #Create plots...
  #Page 1:
  p1.plot1 <- ggplot(subset(OT_V_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_V xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p1.plot2 <- ggplot(subset(PF_L_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("PF_L xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p1.plot3 <- ggplot(subset(OT_V_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_V Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p1.plot4 <- ggplot(subset(PF_L_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("PF_L Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Page 2:
  p2.plot1 <- ggplot(subset(OT_PF_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_PF xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p2.plot2 <- ggplot(subset(OT_L_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_L xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p2.plot3 <- ggplot(subset(OT_PF_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_PF Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p2.plot4 <- ggplot(subset(OT_L_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_L Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Page 3:
  p3.plot1 <- ggplot(subset(PF_V_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("PF_V xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p3.plot2 <- ggplot(subset(V_L_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("V_L xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p3.plot3 <- ggplot(subset(PF_V_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("PF_V Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p3.plot4 <- ggplot(subset(V_L_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("V_L Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Page 4:
  p4.plot1 <- ggplot(subset(OT_NR_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_NR xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p4.plot2 <- ggplot(subset(NR_V_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_V xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p4.plot3 <- ggplot(subset(OT_NR_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("OT_NR Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p4.plot4 <- ggplot(subset(NR_V_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_V Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Page 5:
  p5.plot1 <- ggplot(subset(NR_PF_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_PF xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  p5.plot2 <- ggplot(subset(NR_L_xpEHH, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_L xpEHH") + 
    ylim(NA, logPmax_xpEHH) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p5.plot3 <- ggplot(subset(NR_PF_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_PF Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1))
  p5.plot4 <- ggplot(subset(NR_L_Rsb, POSITION >=social_start & POSITION <= social_end), aes(POSITION, LOGPVALUE)) + 
    geom_point(size=1) + 
    ggtitle("NR_L Rsb") + 
    ylim(NA, logPmax_Rsb) + 
    scale_x_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 8), axis.text.x = element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #Arrange on pages:
  pagelist <- vector("list", 5)
  pagelist[[1]] <- ggarrange(p1.plot1, p1.plot2, p1.plot3, p1.plot4, 
                             ncol = 2, nrow = 2)
  pagelist[[2]] <- ggarrange(p2.plot1, p2.plot2, p2.plot3, p2.plot4, 
                             ncol = 2, nrow = 2)
  pagelist[[3]] <- ggarrange(p3.plot1, p3.plot2, p3.plot3, p3.plot4, 
                             ncol = 2, nrow = 2)
  pagelist[[4]] <- ggarrange(p4.plot1, p4.plot2, p4.plot3, p4.plot4, 
                             ncol = 2, nrow = 2)
  pagelist[[5]] <- ggarrange(p5.plot1, p5.plot2, p5.plot3, p5.plot4, 
                             ncol = 2, nrow = 2)
  #Save to file:
  ggsave(paste(plot_dir, filename, sep="/"), marrangeGrob(grobs=pagelist, nrow = 1, ncol = 1, list(top=NULL)), width=11, height=8.5)
}






#Plot and save OT 62955...
#No FDR:
plot.xpop.social.region.logP(scaffold.ID=62955, 
                             filename="OT_62955_xpEHH_and_Rsb_logP_noFDR_social_region.pdf", 
                             plot_dir="/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Single_scaffold_plots/Social_regions_only", 
                             social_candidate_regions=OTref_social_candidate_regions_master, 
                             OT_V_xpEHH=OT_V_ALL_OTref.xpEHH.noFDR, 
                             PF_L_xpEHH=PF_L_ALL_OTref.xpEHH.noFDR, 
                             OT_V_Rsb=OT_V_ALL_OTref.Rsb.noFDR, 
                             PF_L_Rsb=PF_L_ALL_OTref.Rsb.noFDR, 
                             OT_PF_xpEHH=OT_PF_ALL_OTref.xpEHH.noFDR, 
                             OT_L_xpEHH=OT_L_ALL_OTref.xpEHH.noFDR, 
                             OT_PF_Rsb=OT_PF_ALL_OTref.Rsb.noFDR, 
                             OT_L_Rsb=OT_L_ALL_OTref.Rsb.noFDR, 
                             PF_V_xpEHH=PF_V_ALL_OTref.xpEHH.noFDR, 
                             V_L_xpEHH=V_L_ALL_OTref.xpEHH.noFDR, 
                             PF_V_Rsb=PF_V_ALL_OTref.Rsb.noFDR, 
                             V_L_Rsb=V_L_ALL_OTref.Rsb.noFDR, 
                             OT_NR_xpEHH=OT_NR_ALL_OTref.xpEHH.noFDR, 
                             NR_V_xpEHH=NR_V_ALL_OTref.xpEHH.noFDR, 
                             OT_NR_Rsb=OT_NR_ALL_OTref.Rsb.noFDR, 
                             NR_V_Rsb=NR_V_ALL_OTref.Rsb.noFDR, 
                             NR_PF_xpEHH=NR_PF_ALL_OTref.xpEHH.noFDR, 
                             NR_L_xpEHH=NR_L_ALL_OTref.xpEHH.noFDR, 
                             NR_PF_Rsb=NR_PF_ALL_OTref.Rsb.noFDR, 
                             NR_L_Rsb=NR_L_ALL_OTref.Rsb.noFDR)

#With FDR:
plot.xpop.social.region.logP(scaffold.ID=62955, 
                             filename="OT_62955_xpEHH_and_Rsb_logP_withFDR_social_region.pdf", 
                             plot_dir="/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Single_scaffold_plots/Social_regions_only", 
                             social_candidate_regions=OTref_social_candidate_regions_master, 
                             OT_V_xpEHH=OT_V_ALL_OTref.xpEHH.FDR, 
                             PF_L_xpEHH=PF_L_ALL_OTref.xpEHH.FDR, 
                             OT_V_Rsb=OT_V_ALL_OTref.Rsb.FDR, 
                             PF_L_Rsb=PF_L_ALL_OTref.Rsb.FDR, 
                             OT_PF_xpEHH=OT_PF_ALL_OTref.xpEHH.FDR, 
                             OT_L_xpEHH=OT_L_ALL_OTref.xpEHH.FDR, 
                             OT_PF_Rsb=OT_PF_ALL_OTref.Rsb.FDR, 
                             OT_L_Rsb=OT_L_ALL_OTref.Rsb.FDR, 
                             PF_V_xpEHH=PF_V_ALL_OTref.xpEHH.FDR, 
                             V_L_xpEHH=V_L_ALL_OTref.xpEHH.FDR, 
                             PF_V_Rsb=PF_V_ALL_OTref.Rsb.FDR, 
                             V_L_Rsb=V_L_ALL_OTref.Rsb.FDR, 
                             OT_NR_xpEHH=OT_NR_ALL_OTref.xpEHH.FDR, 
                             NR_V_xpEHH=NR_V_ALL_OTref.xpEHH.FDR, 
                             OT_NR_Rsb=OT_NR_ALL_OTref.Rsb.FDR, 
                             NR_V_Rsb=NR_V_ALL_OTref.Rsb.FDR, 
                             NR_PF_xpEHH=NR_PF_ALL_OTref.xpEHH.FDR, 
                             NR_L_xpEHH=NR_L_ALL_OTref.xpEHH.FDR, 
                             NR_PF_Rsb=NR_PF_ALL_OTref.Rsb.FDR, 
                             NR_L_Rsb=NR_L_ALL_OTref.Rsb.FDR)

#Plot and save PF 144...
#No FDR:
plot.xpop.social.region.logP(scaffold.ID=144, 
                             filename="PF_144_xpEHH_and_Rsb_logP_noFDR_social_region.pdf", 
                             plot_dir="/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Single_scaffold_plots/Social_regions_only", 
                             social_candidate_regions=PFref_social_candidate_regions_master, 
                             OT_V_xpEHH=OT_V_ALL_PFref.xpEHH.noFDR, 
                             PF_L_xpEHH=PF_L_ALL_PFref.xpEHH.noFDR, 
                             OT_V_Rsb=OT_V_ALL_PFref.Rsb.noFDR, 
                             PF_L_Rsb=PF_L_ALL_PFref.Rsb.noFDR, 
                             OT_PF_xpEHH=OT_PF_ALL_PFref.xpEHH.noFDR, 
                             OT_L_xpEHH=OT_L_ALL_PFref.xpEHH.noFDR, 
                             OT_PF_Rsb=OT_PF_ALL_PFref.Rsb.noFDR, 
                             OT_L_Rsb=OT_L_ALL_PFref.Rsb.noFDR, 
                             PF_V_xpEHH=PF_V_ALL_PFref.xpEHH.noFDR, 
                             V_L_xpEHH=V_L_ALL_PFref.xpEHH.noFDR, 
                             PF_V_Rsb=PF_V_ALL_PFref.Rsb.noFDR, 
                             V_L_Rsb=V_L_ALL_PFref.Rsb.noFDR, 
                             OT_NR_xpEHH=OT_NR_ALL_PFref.xpEHH.noFDR, 
                             NR_V_xpEHH=NR_V_ALL_PFref.xpEHH.noFDR, 
                             OT_NR_Rsb=OT_NR_ALL_PFref.Rsb.noFDR, 
                             NR_V_Rsb=NR_V_ALL_PFref.Rsb.noFDR, 
                             NR_PF_xpEHH=NR_PF_ALL_PFref.xpEHH.noFDR, 
                             NR_L_xpEHH=NR_L_ALL_PFref.xpEHH.noFDR, 
                             NR_PF_Rsb=NR_PF_ALL_PFref.Rsb.noFDR, 
                             NR_L_Rsb=NR_L_ALL_PFref.Rsb.noFDR)

#With FDR:
plot.xpop.social.region.logP(scaffold.ID=144, 
                             filename="PF_144_xpEHH_and_Rsb_logP_withFDR_social_region.pdf", 
                             plot_dir="/scratch/genomicsocorg/mwhj1/VCF_phasing/Analysis/Results/FullAnalysis/Plots/Single_scaffold_plots/Social_regions_only", 
                             social_candidate_regions=PFref_social_candidate_regions_master, 
                             OT_V_xpEHH=OT_V_ALL_PFref.xpEHH.FDR, 
                             PF_L_xpEHH=PF_L_ALL_PFref.xpEHH.FDR, 
                             OT_V_Rsb=OT_V_ALL_PFref.Rsb.FDR, 
                             PF_L_Rsb=PF_L_ALL_PFref.Rsb.FDR, 
                             OT_PF_xpEHH=OT_PF_ALL_PFref.xpEHH.FDR, 
                             OT_L_xpEHH=OT_L_ALL_PFref.xpEHH.FDR, 
                             OT_PF_Rsb=OT_PF_ALL_PFref.Rsb.FDR, 
                             OT_L_Rsb=OT_L_ALL_PFref.Rsb.FDR, 
                             PF_V_xpEHH=PF_V_ALL_PFref.xpEHH.FDR, 
                             V_L_xpEHH=V_L_ALL_PFref.xpEHH.FDR, 
                             PF_V_Rsb=PF_V_ALL_PFref.Rsb.FDR, 
                             V_L_Rsb=V_L_ALL_PFref.Rsb.FDR, 
                             OT_NR_xpEHH=OT_NR_ALL_PFref.xpEHH.FDR, 
                             NR_V_xpEHH=NR_V_ALL_PFref.xpEHH.FDR, 
                             OT_NR_Rsb=OT_NR_ALL_PFref.Rsb.FDR, 
                             NR_V_Rsb=NR_V_ALL_PFref.Rsb.FDR, 
                             NR_PF_xpEHH=NR_PF_ALL_PFref.xpEHH.FDR, 
                             NR_L_xpEHH=NR_L_ALL_PFref.xpEHH.FDR, 
                             NR_PF_Rsb=NR_PF_ALL_PFref.Rsb.FDR, 
                             NR_L_Rsb=NR_L_ALL_PFref.Rsb.FDR)

