############ Supernova assemblies - Scaffold Length Distribution ############

# I ran a biopython script to generate sequence lengths from both assemblies:

#import sys
#from Bio import SeqIO
#
#FastaFile = open(sys.argv[1], 'rU')
#
#for rec in SeqIO.parse(FastaFile, 'fasta'):
#  name = rec.id
#seq = rec.seq
#seqLen = len(rec)
#print (name, seqLen)
#
#FastaFile.close()

# I will now plot frequency distributions of these lengths.


#1. Setwd() and read in data

setwd("/Users/u1693640/Documents/PhD/2019/Supernova/Sequence_lengths")
OT_lengths_in <- read.table(file="OT3b_lengths.txt", col.names = c("ID", "Length"))
PF_lengths_in <- read.table(file="PF_lengths.txt", col.names = c("ID", "Length"))


#2. Define data to plot

OT_lengths <- OT_lengths_in$Length
PF_lengths <- PF_lengths_in$Length

#Calculate limits
OT_range <- range(OT_lengths)
PF_range <- range(PF_lengths)
OT_range[1]
PF_range[1]
OT_range[2]
PF_range[2]
#Range is 500 bp to 10028720 bp

#Break the range (+-1 on either side) into non-overlapping bins
breaks <- seq(499, 10028721, by=1000)

#Classify the lengths according to the bins with cut
OT_length.cut <- cut(OT_lengths, breaks, right=FALSE)
PF_length.cut <- cut(PF_lengths, breaks, right=FALSE)

#Compute the frequency of lengths in each bin
OT_length.freq <- table(OT_length.cut)
PF_length.freq <- table(PF_length.cut)

#Format
OT_to_plot <- cbind(OT_length.freq)
PF_to_plot <- cbind(PF_length.freq)


#3. Plot

#Check the plots look right with log scaled x and y axes
plot(PF_to_plot, log="xy", type='h', lwd=1, lend=2, ylab="Frequency", xlab="Scaffold length (kb)", ylim=c(1, 25000), las=1, main="P")
plot(OT_to_plot, log="xy", type='h', lwd=1, lend=2, ylab="Frequency", xlab="Scaffold length (kb)", ylim=c(1, 25000), las=1, main="FM")

#Write to file
pdf(file="Supernova_scaffold_lengths.pdf", paper="a4r", pointsize=10)
par(mfrow=c(2,1))
plot(PF_to_plot, log="xy", type='h', lwd=2, lend=2, ylab="Frequency", xlab="Scaffold length (kb)", ylim=c(1, 25000), las=1, main="P")
plot(OT_to_plot, log="xy", type='h', lwd=2, lend=2, ylab="Frequency", xlab="Scaffold length (kb)", ylim=c(1, 25000), las=1, main="FM")
dev.off()
