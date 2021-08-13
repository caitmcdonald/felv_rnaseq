###

#### edgR: differential expression tests ####

# Dataset: endogenous/exogenous FeLV in vitro gene expression (fibroblasts infected and uninfected, PMBCs uninfected only)

# This code will do the following:
# 1. Generate tissue-specific filtered and normalized DGElists
# 2. Write tables of log-transformed, filtered counts
# 3. Generate a tissue-specific PCA
# 4. Design GLMs, estimate dispersion, and fit models
# 5. Generate DGE contrasts to be made
# 6. Identify sig DEGs in each comparison based on likelihood ratio tests (or quasi-likelihood F-tests which are more conservative), order genes by logFC and FDR<0.05
# 7. Write output tables of DEGs

# Preconditions: gene-level counts data, sample metadata

###


library(edgeR)
library(tidyverse)
library(matrixStats)
library(pheatmap)
library(here)


##counts and metadata; either full gene counts or only cat-puma orthologs
counts0 <- read_delim("results/star_quant_bam/readcountsmatrix.txt", delim = "\t") %>% select(-DC2Pool, -DC3Pool) %>% slice(-c(1:4, 31498))#all counts, not just orthologs
# counts0 <- read_delim("results/orthofinder/genecounts_orthologs.txt", delim="\t")

metadat <- read_tsv("data/felv_metadata.tsv") %>% filter(id_inf != c("DC2Pool", "DC3Pool"))

# Puma vs cat fibroblasts
# counts <- counts0 %>% select(!contains(c("44", "45")))
# meta <- metadat %>% filter(cell_type=="fibroblast")

#uninfected/baseline only: fibroblasts vs. PMBCs
# counts <- counts0 %>% select(!contains(c("Mischief", "PLUS")))
# meta <- metadat %>% filter(cat_id!="Mischief") %>% filter(status=="uninfected")

#uninfected vs infected fibroblasts
# counts <- counts0 %>% select(!contains(c("Mischief","44", "45")))
# meta <- metadat %>% filter(cat_id!="Mischief") %>% filter(cell_type=="fibroblast")

#PMBCs only
counts <- counts0 %>% select(!contains(c("Mischief", "DC", "X")))
meta <- metadat %>% filter(cell_type=="PBMC")

# group <- meta$status
# group <- meta$cell_type
# group <- meta$LTR_pmbc1
group <- meta$LTR_pmbc2
# group <- meta$LTR_pmbc3


#### 1. Create DGElist, filter, and calcNormFactors ####
dat <- DGEList(counts=counts[,-(1)], group=group, genes=counts[,1])
dat$samples

#head(dat$samples$group)
keep_counts<-rowSums(cpm(counts[,-(1)])>1) >= 0.25 * ncol(dat) #filter >1 cpm in >= 1/4 of samples
dat.filt<-dat[which(keep_counts==T), , keep.lib.sizes=FALSE]
dim(dat.filt)
dat.norm<-calcNormFactors(dat.filt)
dat.norm$samples

#### 2. Log-transform, and write tables for downstream analyses ####
#log-transform, all counts
dat.logCPM <- cpm(dat.norm, log=TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

# write table for downstream (e.g. WGCNA)
dat.logCPM.tbl <- cbind(dat.filt$genes, dat.logCPM)
# write_delim(dat.logCPM.tbl, "results/edgeR/counts_fibros_all_filt25_log.txt", delim="\t")


#### 3. Generate PCA ####

# PCA:
pca<-prcomp(t(dat.logCPM))
s<-summary(pca)
s$importance
scores<-as.data.frame(pca$x)
# PCA by Treatment
pca <- ggplot(data=scores, aes(x=PC1, y=PC2, group=group)) +
  geom_point(aes(color=group)) +
  scale_color_manual(values=c("#00b6bd", "lightgrey")) +
  # scale_color_manual(values=c("#008080", "#ef6079")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  xlab ("PC1 (45.6%)") +
  ylab ("PC2 (21.4%)") +
  # guides(colour = guide_legend(override.aes = list(size=2), title="FeLV status")) 
  guides(colour = guide_legend(override.aes = list(size=2), title="LTR integration"))
pca + ggrepel::geom_text_repel(aes(label=meta$cat_id))
# ggsave(file="plots/pca_lTR_pmbc2.jpg")
# ggsave(file="plots/pca_fibroblasts_all_orthologs.jpg")
# ggsave(file="plots/pca_celltype_nopuma.jpg")


#### 4. Design model, estimate dispersion, fit model ####
# design matrix
design.mat <- model.matrix(~ 0 + dat.norm$samples$group)

colnames(design.mat) <- levels(dat.norm$samples$group)
design.mat

dat.norm$samples$group

#estimate dispersion
d1 <- estimateGLMCommonDisp(dat.norm, design.mat, verbose = TRUE)
d1 <- estimateGLMTrendedDisp(d1,design.mat, method="power")
d1 <- estimateGLMTagwiseDisp(d1,design.mat)
plotBCV(d1)

# fit model
fit.glm <- glmFit(d1, design.mat)


#### 5. Set up contrasts ####
my.contrasts <- makeContrasts(
  # inf_uninf = uninfected - infected,
  # celltype = PBMC - fibroblast,
  # population = outbred - puma,
  # ltr_present1 = present - absent,
  # ltr_present2 = present - absent,
  ltr_present3 = present - absent,
  levels=design.mat
)


#### 6. Run LR tests on all contrasts ####

## Individually:
# test <- glmLRT(fit.glm, contrast=my.contrasts[,"celltype"])
# test <- glmLRT(fit.glm, contrast=my.contrasts[,"inf_uninf"])
# test <- glmLRT(fit.glm, contrast=my.contrasts[,"ltr_present1"])
# test <- glmLRT(fit.glm, contrast=my.contrasts[,"ltr_present2"])
test <- glmLRT(fit.glm, contrast=my.contrasts[,"ltr_present3"])

topDE <- topTags(test, adjust.method="BH", n=NULL)
topDE_sig <- topDE[topDE$table$FDR<0.05,]
# topDE_sig$table
summary(decideTests(test, p.value = 0.05))


## If interested in LTR-adjacent subset, specify here:
# pmbc_ltr <- read_delim("data/ltr_data/LTR_pmbc1_genids.txt", delim="\t") %>% distinct()
# pmbc_ltr <- read_delim("data/ltr_data/LTR_pmbc2_genids.txt", delim="\t") %>% distinct()
pmbc_ltr <- read_delim("data/ltr_data/LTR_pmbc3_genids.txt", delim="\t") %>% distinct()
# over10_ltr <- read_delim("data/ltr_data/LTR_over10_genids.txt", delim="\t") %>% distinct()

# Figure out rownumbers of LTR genes
rownumbers <- test$genes %>% mutate(
  rownumber=seq(1:15308)
)
rownumbers <- subset(rownumbers, gene_id %in% pmbc_ltr$gene_id)
rownum <- as_vector(rownumbers$rownumber)

# Restrict test to only LTR genes
test2 <- test[rownum,]
topDE <- topTags(test2, adjust.method="BH", n=NULL)
topDE_sig <- topDE[topDE$table$FDR<0.05,]
# topDE_sig$table
summary(decideTests(test2, p.value = 0.05)) #no significant genes, any PMBC LTR sets


## as loop:
contrast_list <- as.list(colnames(my.contrasts))

topDE_sig <- list()

for (f in 1:length(contrast_list))
{
  test <- glmLRT(fit.glm, contrast=my.contrasts[,f])
  topDE <- topTags(test, adjust.method="BH", n=NULL)
  topDE_order <- topDE[order(topDE$table$logFC),]
  topDE_sig[[f]] <- topDE_order[topDE_order$table$FDR<0.05,]
}
names(topDE_sig) <- unlist(contrast_list)
# view(topDE_sig$inf_uninf$table)


#### 7. Write tables ####
# setwd("results/edgeR")
# lapply(1:length(topDE_sig),
#        function(f) write.table(topDE_sig[[f]], file = paste0("sigDGE_",names(topDE_sig[f]),"_nopuma.txt"), row.names=F, quote=F, sep="\t"))
