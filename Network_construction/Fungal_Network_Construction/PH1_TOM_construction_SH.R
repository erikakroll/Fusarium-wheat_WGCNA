#########################################################################################################################
#
# SCRIPT CONTRIBUTORS     : E.Kroll <erika.kroll@rothamsted.ac.uk>; #
#
# DESCRIPTION : Script to construct PH1 signed hybrid TOM
#
# KEY STAGES  :
#        01. Calculate Soft Thresholding Power
#        02. Construct TOM 
#        
#       
#
# NOTES       :
#
###################################### 00 Setup ###########################################################################

### Load required libraries ###
library(WGCNA)

setwd("/home/data/fus_sysbio_ek2020/Fusarium-Wheat_Network/")
#setwd("Y:\\Fusarium-Wheat_Network\\")

# Extras #
options(stringsAsFactors = FALSE)
#enableWGCNAThreads() OMITTED HERE AS PER TUTORIAL RECOMMENDATION FOR WHEN USING R STUDIO
allowWGCNAThreads() # RStudio


###########################################################################################################################
######################################### 01 Calculate Soft Thresholding Power ############################################

load('EdgeR_PH1_dataInput2022-06-16.RData')
powers = c(1:20)
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType = "signed hybrid", verbose = 5)

# Plot the results:

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;

# Scale-free topology fit index as a function of the soft-thresholding power

x = sft$fitIndices[,1]

y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]

sftdf = data.frame(x,y)

source("ggplot_theme_Publication.R")



p <- ggplot(sftdf, aes(x, y, label = rownames(sftdf))) + 
  xlab("Soft Threshold (Power)") +
  ylab("Scale Free Topology Model Fit (R-squared)") +
  theme(text=element_text(size=16,  family="serif"))

p2 = p + geom_text() 

p2


plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h

abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Calculate adjacency matrix using appropriate soft thresholding power #
SoftPower = 9
adjacency = adjacency(datExpr0, type = "signed hybrid", power = SoftPower)

################################################################################################################################
################################# 03 Calculate TOM #############################################################################

# Create a Topological Overlap Matrix (TOM) from adjacency matrix #
TOM = TOMsimilarity(adjacency, TOMType = "signed")

currentDate <- Sys.Date()
save(TOM, file = paste("EdgeR_sh_9_PH1_TOM",currentDate,".RData",sep=""))

###################################################################################################################
