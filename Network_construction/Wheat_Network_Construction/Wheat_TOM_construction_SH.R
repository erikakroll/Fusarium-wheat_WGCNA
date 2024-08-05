#########################################################################################################################
#
# SCRIPT CONTRIBUTORS     : E.Kroll <erika.kroll@rothamsted.ac.uk>; #
#
# DESCRIPTION : Script to construct Wheat signed hybrid TOM
#
# KEY STAGES  :
#        
#        00. Setup  
#        01. Calculate Soft Thresholding Power
#        02. Construct TOM 
#        
#       
#
# NOTES       :
#
##########################################################################################################################
###################################### 00 Setup ###########################################################################


### Load required libraries ###
library(WGCNA)

setwd("/home/data/fus_sysbio_ek2020/Fusarium-Wheat_Network/")
#setwd("Y:\\Fusarium-Wheat_Network")

# Extras #
options(stringsAsFactors = FALSE)
#enableWGCNAThreads() OMITTED HERE AS PER TUTORIAL RECOMMENDATION FOR WHEN USING R STUDIO
allowWGCNAThreads() # RStudio

load('EdgeR_Wheat_dataInput2022-06-16.RData')

############################################################################################################################
########################################## 01. Calculate Soft Thresholding Power ##########################################

powers = c(1:30)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, networkType = "signed hybrid")

sft_wheat = sft

save(sft_wheat, file = "sft_wheat.RData")

############ Plot the results:

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

########### Scale-free topology fit index as a function of the soft-thresholding power

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
SoftPower = 17
adjacency = adjacency(datExpr0, type = "signed hybrid", power = SoftPower)


################################################################################################################################
##################################### 02. Construct TOM ########################################################################

# Create a Topological Overlap Matrix (TOM) from adjacency matrix #
TOM = TOMsimilarity(adjacency, TOMType = "signed")

currentDate <- Sys.Date()
save(TOM, file = paste("Wheat_sh_18_TOM",currentDate,".RData",sep=""))

