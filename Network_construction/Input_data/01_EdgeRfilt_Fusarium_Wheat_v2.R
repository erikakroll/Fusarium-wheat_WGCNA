#########################################################################################################################
#
# SCRIPT CONTRIBUTORS     : E.Kroll <erika.kroll@rothamsted.ac.uk>; #
#
# DESCRIPTION : Script to prepare data input for dual F. graminearum PH1 
#                     and Wheat (Bobwhite) WGCNA Network 
#
# KEY STAGES  :
#        01. Reading in counts, quality and filtering checks
#        02. Normalise/ quality control
#        03. Separating Wheat Data for WGCNA
#        04. Separating PH1 Data for WGCNA
#       
#
# NOTES       :
#
#
#########################################################################################################################
#######################################01. Reading in counts, quality and filtering checks ########################################


################### load libraries

library(edgeR)
library(ggplot2)
library(genefilter)
library(DESeq2)
library(xlsx)
library(WGCNA)
library(tidyverse)
library(gridExtra)


setwd("Y://Fusarium-Wheat_Network//counts")



PH1_7d_SP_A.f = 'PH1_7d_SP_AReadsPerGene.out.tab'
PH1_7d_SP_B.f = 'PH1_7d_SP_BReadsPerGene.out.tab'
PH1_7d_SP_C.f = 'PH1_7d_SP_CReadsPerGene.out.tab'
PH1_3d_SP_A.f = 'PH1_3d_SP_AReadsPerGene.out.tab'
PH1_3d_SP_B.f = 'PH1_3d_SP_BReadsPerGene.out.tab'
PH1_3d_SP_C.f = 'PH1_3d_SP_CReadsPerGene.out.tab'
PH1_7d_R12A.f = 'PH1_7d_R1-2_AReadsPerGene.out.tab'
PH1_7d_R12C.f = 'PH1_7d_R1-2_CReadsPerGene.out.tab'
PH1_7d_R34A.f = 'PH1_7d_R3-4_AReadsPerGene.out.tab'
PH1_7d_R34B.f = 'PH1_7d_R3-4_BReadsPerGene.out.tab'
PH1_7d_R34C.f = 'PH1_7d_R3-4_CReadsPerGene.out.tab'
PH1_7d_R56A.f = 'PH1_7d_R5-6_AReadsPerGene.out.tab'
PH1_7d_R56B.f = 'PH1_7d_R5-6_BReadsPerGene.out.tab'
PH1_7d_R56C.f = 'PH1_7d_R5-6_CReadsPerGene.out.tab'
PH1_7d_R78A.f = 'PH1_7d_R7-8_AReadsPerGene.out.tab'
PH1_7d_R78B.f = 'PH1_7d_R7-8_BReadsPerGene.out.tab'
PH1_7d_R78C.f = 'PH1_7d_R7-8_CReadsPerGene.out.tab'
Mock_3d_SP_A.f = 'Mock_3d_SP_AReadsPerGene.out.tab'
Mock_3d_SP_B.f = 'Mock_3d_SP_BReadsPerGene.out.tab'
Mock_7d_SP_A.f ='Mock_7d_SP_AReadsPerGene.out.tab'
Mock_7d_SP_B.f = 'Mock_7d_SP_BReadsPerGene.out.tab'
Mock_7d_R12_A.f = 'Mock_7d_R1_AReadsPerGene.out.tab'
Mock_7d_R12_B.f = 'Mock_7d_R1_BReadsPerGene.out.tab'

# Read in our data

PH1_7d_SP_A = read.table(PH1_7d_SP_A.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_SP_B = read.table(PH1_7d_SP_B.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_SP_C = read.table(PH1_7d_SP_C.f, sep='\t', row.names = 1, skip = 4)
PH1_3d_SP_A = read.table(PH1_3d_SP_A.f, sep='\t', row.names = 1, skip = 4)
PH1_3d_SP_B = read.table(PH1_3d_SP_B.f, sep='\t', row.names = 1, skip = 4)
PH1_3d_SP_C = read.table(PH1_3d_SP_C.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R12A = read.table(PH1_7d_R12A.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R12C = read.table(PH1_7d_R12C.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R34A = read.table(PH1_7d_R34A.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R34B = read.table(PH1_7d_R34B.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R34C = read.table(PH1_7d_R34C.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R56A = read.table(PH1_7d_R56A.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R56B = read.table(PH1_7d_R56B.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R56C = read.table(PH1_7d_R56C.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R78A = read.table(PH1_7d_R78A.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R78B = read.table(PH1_7d_R78B.f, sep='\t', row.names = 1, skip = 4)
PH1_7d_R78C = read.table(PH1_7d_R78C.f, sep='\t', row.names = 1, skip = 4)
Mock_3d_SP_A = read.table(Mock_3d_SP_A.f, sep='\t', row.names = 1, skip = 4)
Mock_3d_SP_B = read.table(Mock_3d_SP_B.f, sep='\t', row.names = 1, skip = 4)
Mock_7d_R12_A = read.table(Mock_7d_R12_A.f, sep='\t', row.names = 1, skip = 4)
Mock_7d_R12_B = read.table(Mock_7d_R12_B.f, sep='\t', row.names = 1, skip = 4)
Mock_7d_SP_A = read.table(Mock_7d_SP_A.f, sep='\t', row.names = 1, skip = 4)
Mock_7d_SP_B = read.table(Mock_7d_SP_B.f, sep='\t', row.names = 1, skip = 4)

##########Overall MDS plot#########

###Sample info 
sample.info = read.csv('sampleinfo_mocks.csv')
sample.info$sample = factor(sample.info$sample)
sample.info$tissue = factor(sample.info$tissue)
sample.info$stage = factor(sample.info$tissue)
colnames(sample.info) = gsub("strain", "Treatment", colnames(sample.info))

counts= data.frame(PH1_7d_SP_A = PH1_7d_SP_A$V2, PH1_7d_SP_B = PH1_7d_SP_B$V2, 
                   PH1_7d_SP_C = PH1_7d_SP_C$V2, PH1_3d_SP_A = PH1_3d_SP_A$V2,
                   PH1_3d_SP_B = PH1_3d_SP_B$V2, PH1_3d_SP_C = PH1_3d_SP_C$V2,
                   PH1_7d_R12A = PH1_7d_R12A$V2, PH1_7d_R12C = PH1_7d_R12C$V2,
                   PH1_7d_R34A = PH1_7d_R34A$V2, PH1_7d_R34B = PH1_7d_R34B$V2,
                   PH1_7d_R34C = PH1_7d_R34C$V2, PH1_7d_R56A = PH1_7d_R56A$V2,
                   PH1_7d_R56B = PH1_7d_R56B$V2, PH1_7d_R56C = PH1_7d_R56C$V2,
                   PH1_7d_R78A = PH1_7d_R78A$V2, PH1_7d_R78B = PH1_7d_R78B$V2,
                   PH1_7d_R78C = PH1_7d_R78C$V2, Mock_3d_SP_A = Mock_3d_SP_A$V2,
                   Mock_3d_SP_B = Mock_3d_SP_B$V2,
                   Mock_7d_SP_A = Mock_7d_SP_A$V2, 
                   Mock_7d_SP_B = Mock_7d_SP_B$V2, 
                   Mock_7d_R12_A = Mock_7d_R12_A$V2,
                   Mock_7d_R12_B = Mock_7d_R12_B$V2, 
                   row.names = row.names(PH1_7d_R12A))

row.names(counts) <- gsub('transcript:', '', row.names(counts))

dds = DESeqDataSetFromMatrix(countData = counts, colData = sample.info, 
                             design = ~ tissue)
norm.counts = as.data.frame(assay(varianceStabilizingTransformation(dds)))

sampleDists = dist(t(norm.counts))
sampleDistMatrix = as.matrix(sampleDists)
mdsData = data.frame(cmdscale(sampleDistMatrix))
mds = cbind(mdsData, as.data.frame(colData(varianceStabilizingTransformation(dds))))


# Plot without arrows #
mds_p = ggplot(mds, aes(X1, X2, color=Treatment, shape = sample.1)) + geom_point(size=5)+ 
  scale_shape(name = "Sample") + scale_colour_hue(name = "Treatment") + 
  theme(legend.position="right")

setwd("Y://Fusarium-Wheat_Network//")
source("ggplot_theme_Publication.R")

mds_p + scale_colour_Publication()+ theme_Publication() 


################Filtering to remove PH-1 off-targets###########

###Mock only 

sample.info = read.csv('sampleinfo_onlymocks.csv')
sample.info$sample = factor(sample.info$sample)
sample.info$tissue = factor(sample.info$tissue)

counts_mocks= data.frame(Mock_3d_SP_A = Mock_3d_SP_A$V2,
                   Mock_3d_SP_B = Mock_3d_SP_B$V2, Mock_7d_SP_A = Mock_7d_SP_A$V2, 
                   Mock_7d_SP_B = Mock_7d_SP_B$V2, Mock_7d_R12_A = Mock_7d_R12_A$V2,
                   Mock_7d_R12_B = Mock_7d_R12_B$V2, row.names = row.names(PH1_7d_R12A))
row.names(counts_mocks) <- gsub('transcript:', '', row.names(counts_mocks))

counts_PH1_mock = counts_mocks[-c(1:132624), ]

#Group the samples together based on their conditions
# 1 - Mock 3d SP
# 2 - Mock 7d SP
# 3 - PH1 7d R1-2

group = factor(c(1, 1, 2, 2, 3, 3))

keep <- filterByExpr(counts_PH1_mock, group = group)
filt.counts <- counts_PH1_mock[keep, ]

##########The number of genes with F graminearum off targets############

print(nrow(filt.counts))

#########################################################################################################################
######################################### 02. Normalise/ quality control #################################################

######Generate counts that no longer have mock reads in them#########

sample.info_dilks = read.csv('sampleinfo_nomocks.csv')
sample.info_dilks$sample = factor(sample.info_dilks$sample)
sample.info_dilks$tissue = factor(sample.info_dilks$tissue)
sample.info_dilks$stage = factor(sample.info_dilks$stage)


counts_nomocks= data.frame(PH1_7d_SP_A = PH1_7d_SP_A$V2, PH1_7d_SP_B = PH1_7d_SP_B$V2, PH1_7d_SP_C = PH1_7d_SP_C$V2, PH1_3d_SP_A = PH1_3d_SP_A$V2,
                   PH1_3d_SP_B = PH1_3d_SP_B$V2, PH1_3d_SP_C = PH1_3d_SP_C$V2, PH1_7d_R12A = PH1_7d_R12A$V2, 
                   PH1_7d_R12C = PH1_7d_R12C$V2, PH1_7d_R34A = PH1_7d_R34A$V2, PH1_7d_R34B = PH1_7d_R34B$V2,  
                   PH1_7d_R34C = PH1_7d_R34C$V2, PH1_7d_R56A = PH1_7d_R56A$V2, PH1_7d_R56B = PH1_7d_R56B$V2, 
                   PH1_7d_R56C = PH1_7d_R56C$V2, PH1_7d_R78A = PH1_7d_R78A$V2, PH1_7d_R78B = PH1_7d_R78B$V2, 
                   PH1_7d_R78C = PH1_7d_R78C$V2, row.names = row.names(PH1_7d_R12A))
row.names(counts_nomocks) <- gsub('transcript:', '', row.names(counts_nomocks))



###Filter counts

#Group the samples together based on their conditions
# 1 - PH1 7d SP
# 2 - PH1 3d SP
# 3 - PH1 7d R1-2
# 4 - PH1 7d R3-4
# 5 - PH1 7d R5-6
# 6 - PH1 7d R7-8

group = factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6))

keep <- filterByExpr(counts_nomocks, group = group)
filt.counts <- counts_nomocks[keep, ]

#########################################################################################################################################
########################### 03. Separating Wheat Data for WGNCA ###############################################################################


######### separate wheat counts ##########

counts_wheat = filt.counts[grepl("Traes*", rownames(filt.counts)), ]

### Variance stabilising normalisation ###

dds_all = DESeqDataSetFromMatrix(countData = counts_wheat, colData = sample.info_dilks, design = ~ tissue + stage)
norm.counts_all = as.data.frame(assay(varianceStabilizingTransformation(dds_all)))

#####Write datExpr object#####

datExpr0 = as.data.frame(t(norm.counts_all))

#######Outlier Identification using WGCNA package#########

gsg = goodSamplesGenes(datExpr0, verbose = 5);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

########### sanity check that there are only wheat genes ####

PH1_counts <- datExpr0[ , grepl( "FGRAM" , names( datExpr0 ) ) ]
print(ncol(PH1_counts))

Wheat_counts <- datExpr0[ , grepl( "Traes" , names( datExpr0 ) ) ]
print(ncol(Wheat_counts))

####### Cluster samples to visualise outliers ##########

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.5);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     ylab = "",
     cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 1)

title(ylab= substitute(paste(italic("Height"))), line=2, cex.lab=1.2)

############Include Trait Data in Expr Vector###########

allTraits = read.csv("sampleinfo_nomocks.csv")
samples = rownames(datExpr0)
traitRows = match(samples, allTraits$sample)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1];

currentDate <- Sys.Date()
csvFileName <- 

save(datExpr0, datTraits, file = paste("EdgeR_Wheat_dataInput",currentDate,".RData",sep=""))

############ 04. Seperating PH1 data for WGCNA #################################################

counts_PH1 = filt.counts[grepl("FGRAM*", rownames(filt.counts)), ]

### Variance stabilising normalisation ###

dds_all = DESeqDataSetFromMatrix(countData = counts_PH1, colData = sample.info_dilks, design = ~ tissue + stage)
norm.counts_all = as.data.frame(assay(varianceStabilizingTransformation(dds_all)))

#Construct datExpr

datExpr0 = as.data.frame(t(norm.counts_all))

#######Outlier Identification using WGCNA filter ##########

gsg = goodSamplesGenes(datExpr0, verbose = 5);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

################Sanity check that no Wheat genes are present #####################

PH1_counts <- datExpr0[ , grepl( "FGRAM" , names( datExpr0 ) ) ]
print(ncol(PH1_counts))

Wheat_counts <- datExpr0[ , grepl( "Traes" , names( datExpr0 ) ) ]
print(ncol(Wheat_counts))

################ Cluster samples to visualise outliers ###########################

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.5);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", ylab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 1)

title(ylab= substitute(paste(italic("Height"))), line=2.5, cex.lab=1.2)

############Include Trait Data in Expr Vector###########

samples = rownames(datExpr0)
allTraits = read.csv("sampleinfo_nomocks.csv")
traitRows = match(samples, allTraits$sample)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1];

currentDate <- Sys.Date()
save(datExpr0, datTraits, file = paste("EdgeR_PH1_dataInput",currentDate,".RData",sep=""))


