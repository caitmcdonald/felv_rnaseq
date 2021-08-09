####

#### WGCNA, remote ####

### Goal: Examine patterns in gene expression using network analysis and variables of interest

# This code will do the following:
# 1. Load, transform data
# 2. ID and remove outliers
# 3. Perform stepwise gene network construction
# 3a. Pick soft thresholding power
# 3b. Hierarchical gene clustering
# 3c. Identify gene modules
# 4. Correlate gene modules with variables of interest
# 5. ID significant modules and genes

# Tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

# Preconditions: log-transformed gene counts, with or without filtered subset (e.g. subset by variancePartition)

###

library(WGCNA)
library(igraph)
library(tidyverse)
library(here)

#### 1. Load and transform data ####

# Read in counts data (change tissue as necessary)
# vP_dat <- read_delim("bpitanga_2018/bpit_varPar/varPar_0.66/top10_sep_nobatch_skin.txt", delim="\t")
# counts <- read_delim("results/edgeR/counts_fibrouninf_filt25_log.txt", delim="\t")
counts <- read_delim("results/edgeR/counts_fibros_all_filt25_log.txt", delim="\t")
names(counts)

dat0 <- counts
# dat0 <- semi_join(counts, vP_dat, by=c("gene_id"="gene")) %>% 
  # select(-transcript_id, -Blast_name, -Full_name, -GO_id) 


# Read in trait data
# traits <- read_delim("data/felv_metadata.tsv", delim="\t") %>% select(id_inf, num_LTR)
traits <- read_delim("data/felv_metadata.tsv", delim="\t") %>% select(id_inf, status, num_LTR)
traits0 <- filter(traits,id_inf %in% names(dat0)) %>% column_to_rownames("id_inf")
# traits1 <- traits0
traits1 <- binarizeCategoricalColumns(traits0, convertColumns = c("status"), includePairwise = T, includeLevelInformation = T, includeLevelVsAll = F)
# names(traits1) <- c("samples","PMBCvsfibroblast", "SPFvsoutbred","uninfectedvsinfected","LTR_number")


tdat <- as.data.frame(t(dat0[,-1]))
names(tdat) <- dat0$genes  
rownames(tdat) <- names(dat0)[-1]


#### 2. ID outliers and remove ####
## ID samples with too many missing values ##
gsg = goodSamplesGenes(tdat, verbose = 3);
gsg$allOK


## plot heirarchical clustering to ID outliers ##
sampleTree = hclust(dist(tdat), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

## Plot a line to show the (arbitrary) cut to remove outliers ##
abline(h = 130, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 130, minSize = 3)
table(clust)
# clust 1 contains the samples we want to keep.

## Remove outliers ##
keepSamples <- (clust==1)
dat <- tdat[keepSamples, ]
nGenes = ncol(dat)
nSamples = nrow(dat)


#### 3. Stepwise network construction

#### 3a. Pick soft thresholding power ####
## for large dataset, run on cluster

## Enable multi-threading
# enableWGCNAThreads(nThreads = 2)

## Load data saved from part I
#lnames = load(file = "Nvir_wgcna_input.RData");

## Choose soft thresholding powers 
powers <- c(c(1:10), seq(from=12, to=30, by=2))

## Call network topology analysis function
sft <- pickSoftThreshold(dat, powerVector=powers, verbose=5, networkType = "signed")

## Plot results:
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1=0.9

## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main=paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
## This line corresponds to using an R^2 cut-off of h
abline(h=0.80, col="red")
## Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")


#### 3b. Hierarchical gene clustering ####
# Calculate adjacencies
softPower <- 14 #recommended for signed networks
adjacency <- adjacency(dat, power=softPower, type = "signed")

# Transform adjacency to Topological Overlap Matrix and calculate corresponding dissimilarity
TOM=TOMsimilarity(adjacency, TOMType= "signed")
dissTOM <- 1-TOM

# Hierarchical clustering of genes
# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method="average")

# Plot the resulting dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity", labels=FALSE, hang=0.04)


#### 3c. Identify gene modules ####
# Set minimum module size
minModuleSize <- 30

# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)

# Plot module assignment under gene dendrogram:
# Convert numeric labels to colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Plot dendrogram and module colors
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, 
                    hang=0.03, addGuide = TRUE, guideHang = 0.05, 
                    main="Gene dendrogram and module colors")

# Merge modules with similar expression profiles
# Calculate eigengenes
MEList <- moduleEigengenes(dat, colors=dynamicColors)

MEs <- MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method="average")

# Plot result
sizeGrWindow(7,6)
plot(METree, main="Clustering of module eigengenes", xlab="", sub="")

# Choose dissimilarity threshold
MEDissThres <- 0.25

# Plot the cut line on the eigengene dendrogram
abline(h=MEDissThres, col="red")

# Call automatic merging function
merge <- mergeCloseModules(dat, dynamicColors, cutHeight = MEDissThres, verbose=3)

# Create merged module colors
mergedColors <- merge$colors
table(mergedColors)

# Create eigengenes for the new modules
mergedMEs <- merge$newMEs

# See effect of merging on module colors by plotting gene dendrogram with original and merged
sizeGrWindow(12,9)
#pdf(file="Plots/geneDendro.pdf", wi=9, he=6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang=0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Save relevant variables for subsequent analysis
# Rename moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs
# Save module colors and labels for subsequent analysis
# save(MEs, moduleLabels, moduleColors, geneTree, file="results/wgcna/networkConstruction_stepwise.RData")


# ### Export for Cytoscape
# # Recalculate topological overlap if needed
# TOM = TOMsimilarityFromExpr(dat, power = 14);
# # Select modules
# modules = c("green");
# # Select module probes
# probes = names(dat)
# inModule = is.finite(match(moduleColors, modules));
# modProbes = probes[inModule];
# modGenes = counts$Full_name[match(modProbes, counts$gene_id)]; #add alternate name (e.g. Blast name)
# # Select the corresponding Topological Overlap
# modTOM = TOM[inModule, inModule];
# dimnames(modTOM) = list(modProbes, modProbes)
# # Export the network into edge and node list files Cytoscape can read
# cyt = exportNetworkToCytoscape(modTOM,
#                                edgeFile = paste("./bpitanga_2018/bpit_wgcna/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
#                                nodeFile = paste("./bpitanga_2018/bpit_wgcna/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
#                                weighted = TRUE,
#                                threshold = 0.08, #default in vignette is 0.02
#                                nodeNames = modProbes,
#                                altNodeNames = modGenes,
#                                nodeAttr = moduleColors[inModule]);



#### 4. Relate modules to variables of interest ####

## Load saved data from module assignment
# Load expression and trait data
#lnames <- load(file="Nvir_wgcna_input.RData")

# Load network data
#lnames <- load(file="Nvir_networkConstruction_auto.RData")

## Quantify module-trait associations (correlate module eigengenes with external traits, look for most significant associations)

# Define numbers of genes and samples
nGenes <- ncol(dat)
nSamples <- nrow(dat)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(dat, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, traits1, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sigmodules0 <- as_tibble(moduleTraitPvalue, rownames="module") %>% 
  filter(num_LTR <= 0.05)

sigtraitcor <- as_tibble(moduleTraitCor, rownames="module") %>% 
  filter(module %in% sigmodules0$module) %>% 
  column_to_rownames("module") %>% as.matrix()

sigmodules1 <- sigmodules0 %>% column_to_rownames("module") %>% as.matrix()
sigMEs <- as_tibble(MEs) %>% select(rownames(sigmodules1))

#pdf("bpitanga_2018/bpit_wgcna/results/skin_sig_heatmap.pdf", height = 8, width=10)
sizeGrWindow(10,6)
textMatrix <- paste(signif(sigtraitcor, 2), "\n(",
                    signif(sigmodules1, 1), ")", sep="")
dim(textMatrix) <- dim(sigmodules1)
par(mar=c(6,8.5,3,3))
labeledHeatmap(Matrix = sigtraitcor,
               xLabels=names(traits1),
               yLabels=names(sigMEs),
               ySymbols=names(sigMEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=FALSE,
               cex.text=0.6,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
#dev.off()


# Examine associations visually
#sizeGrWindow(10,6)
# Display correlations and their p-values
# dev.off()
# textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
#                     signif(moduleTraitPvalue, 1), ")", sep="")
# dim(textMatrix) <- dim(moduleTraitCor)
# par(mar=c(6,8.5,3,3))
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels=names(traits2),
#                yLabels=names(MEs),
#                ySymbols=names(MEs),
#                colorLabels=FALSE,
#                colors=blueWhiteRed(50),
#                textMatrix=textMatrix,
#                setStdMargins=FALSE,
#                cex.text=0.5,
#                zlim=c(-1,1),
#                main=paste("Module-trait relationships"))



#### 5. ID significant modules and genes of interest ####

## Gene significance and module membership
# Quantify associations of individual genes with traits of interest by defining Gene Significance (GS) as the absolute value of the correlation between the gene and the trait. Also define quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes to every module.

# Define variable from traits2
variable <- as.data.frame(traits1$num_LTR)

#variable <- as.data.frame(traits2$time)
#names(variable)="time"

# names (colors) of the modules
modNames <- substring(names(MEs), 3)
modNames

geneModuleMembership <- as.data.frame(cor(dat, MEs, use="p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")

geneTraitsSignificance <- as.data.frame(cor(dat, variable, use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitsSignificance), nSamples))

names(geneTraitsSignificance) <- paste("GS", names(variable), sep="")
names(GSPvalue) <- paste("p.GS.", names(variable), sep="")


## Having ID'd modules with high association with the trait of interest and 'central players' by Module Membership, now merge stats information with gene annotation and write out a file that summarizes the most important results
# Match contig names in annotation file and dat
annot <- counts[,1]
contigs <- names(dat)
contigs
contigs2annot <- match(contigs, annot$genes)
sum(is.na(contigs2annot))

# Creat a data frame holding information for all contigs: ID, gene_id, module color, gene significance for TREATMENT, and module membership and p-values.
geneInfo0 <- data.frame(
  moduleColor=moduleColors,
  geneTraitsSignificance,
  GSPvalue)

# Order modules by their significance for treatment
modOrder <- order(-abs(cor(MEs, variable, use="p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames=names(geneInfo0)
  geneInfo0=data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                       MMPvalue[, modOrder[mod]])
  names(geneInfo0)=c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                     paste("p.MM", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitsSignificance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GStraits1.num_LTR))
geneInfo <- geneInfo0[geneOrder, ]

sigmodules0 <- as_tibble(moduleTraitPvalue, rownames="module") %>% 
  filter(num_LTR<=0.05) %>%
  #filter(time<0.05) %>% 
  #filter(load_at_death<0.05) %>%  
  select(module) %>% 
  mutate(
    mod=substring(module,3)
  )

sigmodules <- sigmodules0$mod

sigmods_topGO <- as_tibble(geneInfo, rownames="gene_id") %>% 
  filter(moduleColor %in% sigmodules) %>%
  filter(moduleColor=="black") %>%
  # select(genes, moduleColor)
  write_csv("results/wgcna/wgcna_fibro_all_black.csv")


## ID genes with high GS and MM: using GS and MM measures, ID genes with high significance for 'treatment' as well as high module membership in interesting modules. For example, look at the module with the highest association with treatment (e.g. skyblue), and plot a scatterplot of GS vs. MM in that module:
module <- "black"
column <- match(module, modNames)
moduleGenes <- moduleColors==module

#sizeGrWindow(7,7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitsSignificance[moduleGenes, 1]),
                   xlab=paste("Module Membership in", module, "module"),
                   ylab="Gene significanc for treatment",
                   main=paste("Module membership vs. gene significance\n"),
                   cex.main=1.2, cex.lab=1.2, cex.axis=1.2, col=module)

count(geneModuleMembership[moduleGenes,])

