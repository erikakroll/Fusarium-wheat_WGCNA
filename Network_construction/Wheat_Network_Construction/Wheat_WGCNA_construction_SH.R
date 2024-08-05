#########################################################################################################################
#
# SCRIPT CONTRIBUTORS     : E.Kroll <erika.kroll@rothamsted.ac.uk>; #
#
# DESCRIPTION : Script to construct Wheat network  for signed dual F. graminearum PH1 
#                     and Wheat (Bobwhite) WGCNA Network. DeepSplit of 3 and merge cut height 0.25. Minmodule size 50.  
#
# KEY STAGES  :
#        00 Setup 
#        01 Construct Modules
#        02 Saves
#        
#       
#
# NOTES       :
#
###########################################################################################################################
######################################### 00 Setup ########################################################################

### Load required libraries ###
library(WGCNA)
library(ggplot2)



setwd("/home/data/fus_sysbio_ek2020/Fusarium-Wheat_Network/")
#setwd("Y://Fusarium-Wheat_Network")

# Extras #
options(stringsAsFactors = FALSE)
#enableWGCNAThreads() OMITTED HERE AS PER TUTORIAL RECOMMENDATION FOR WHEN USING R STUDIO
allowWGCNAThreads() # RStudio

load('EdgeR_Wheat_dataInput2022-06-16.RData')

load('Wheat_sh_18_TOM2022-06-27.RData')

##############################################################################################################################
#################################### 01 Construct Modules ####################################################################

dissTOM = 1-TOM


# Clustering using dissTOM #
geneTree = hclust(as.dist(dissTOM), method = "average")

# Identify modules using dynamic tree cut #
minModuleSize = 50
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 3, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric labels to colours #
dynamicColors = dynamicMods
table(dynamicColors)

# Calculate module eigengenes #
MEList = moduleEigengenes((datExpr0), colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes #
MEDiss = 1-cor(MEs)

# Cluster module eigengenes and plot result #
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red", lwd = 0.35)
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# Calculate the merged module colors Eigengenes of the new merged modules #
mergedColors = merge$colors
mergedMEs = merge$newMEs
table(mergedColors)

# Update module colours, labels and module eigengenes #

moduleColors = labels2colors(mergedColors)
colorOrder = c("grey", standardColors(100))
moduleLabels = as.character(mergedColors)
MEs = mergedMEs
row.names(MEs) = row.names(data.frame(datExpr0))
MEDiss = 1-cor(MEs)

sizeGrWindow(12, 9)
pdf(file = "geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Get gene-module assignments #
gene.module.assignments = data.frame("ID" = colnames(datExpr0),
                                     "Module.Assignment" = moduleColors,
                                     check.names = F)

gene.label.assignments = data.frame("ID" = colnames(datExpr0),
                                     "Module.Assignment" = moduleLabels,
                                     check.names = F)


##########################################################################################################
#################################### 02 Saves ############################################################

currentDate = Sys.Date()

# Module gene assignments #

write.table(gene.module.assignments, paste("Wheat_sh_GeneModule_Assignments_",currentDate,".txt"), sep = "\t", quote=F, row.names = F)
write.table(gene.label.assignments, paste("Wheat_sh_GeneLabel_Assignments_",currentDate,".txt"), sep = "\t", quote=F, row.names = F)

# Module eigengenes #
write.table(MEs, paste("Wheat_sh_ModuleEigengenes_",currentDate,".txt"), sep = "\t", quote=F, row.names = T)
write.table(MEDiss, paste("Wheat_sh_ModuleEigengenes_distances_",currentDate,".txt"), sep='\t', quote=F, row.names=T)

# Module size data #
ColourModules = data.frame(table(moduleColors))
write.table(ColourModules, paste("Wheat_sh_ModuleSizeDataColours_",currentDate,".txt"), sep="\t", quote=F, row.names = F)

# Save normalised counts, module eigengenes, colour/module labels and network data #
save(datExpr0, MEs, moduleLabels, moduleColors, geneTree, file = paste("Wheat_sh_NetworkData_",currentDate,".RData"))