###

# This code is used to explore gene expression data and conduct expression tests in edgeR

# Dataset: endogenous/exogenous FeLV in vitro gene expression (fibroblasts infected and uninfected, PMBCs uninfected only)

#### This code will do the following ####
# 1. Generate a filtered and normalized DGElist of all sample data
# 2. Generate a logCPM counts matrix
# 3. Use logCPM counts matrix to generate sample correlation pheatmap
# 4. Use logCPM counts matrix to generate PCAs to test for batch effects and outliers
# 5. Subset counts and metadata dfs by tissue
# 6. Write tissue-specific, raw counts and metadata files

# Preconditions: raw gene-level counts matrix (e.g. from STAR quantMode), metadata file with sample names that correspond to sample names in counts matrix

###

library(here)
library(tidyverse)
library(edgeR)
library(pheatmap)


#### 0. Get data ####
## counts
dat <- read_delim("readcountsmatrix.txt", delim = "\t") %>% select(-DC2Pool, -DC3Pool)
mapping_dat <- dat[c(1:4, 31498),]
counts <- dat[-c(1:4, 31498),]

## metadata
meta <- read_tsv("data/felv_metadata.tsv") %>% filter(id_inf != c("DC2Pool", "DC3Pool")) %>% filter(id_inf != "DC1PLUS")
Cell_type=meta$cell_type
Infection_status=meta$status
Population=meta$population
cell_inf=paste(meta$cell_type, meta$status, sep="_")


#### 1. Create DGElist, filter by CPM, and calcNormFactors ####
## Create DGEList
dat.full <- DGEList(counts=counts[,2:16], group=cell_inf, genes = counts[,1])
# dat.full$samples$group

## Filter
keep_counts<-rowSums(cpm(counts[,2:16])>1) >= 0.25*ncol(dat) #filter >1 cpm in >= 1/2 of samples

dat.filt<-dat.full[which(keep_counts==T), , keep.lib.sizes=FALSE]
# dim(dat.full)
# dim(dat.filt)

dat.norm<-calcNormFactors(dat.filt)
dat.norm$samples #normalization factors and library size
# plotMDS(dat.norm)


#### 2. Generate moderated, logCPM counts matrix for sample correlation and PCA ####
dat.logCPM <- cpm(dat.norm, log=TRUE, prior.count = 1, normalized.lib.sizes = TRUE)


#### 3. Generate sample correlation matrix and pheatmap ####
## annotate columns by tissue, treatment, and time
annotation_col1 <- data.frame(
  Cell_type,
  Infection_status,
  Population
)
rownames(annotation_col1) <- colnames(dat[,c(2:16)])

ann_colors1 = list(
  Cell_type =c(fibroblast="#008080",PBMC="#ef6079"),
  Infection_status = c(infected = "#00b6bd", uninfected = "lightgrey"),
  Population = c(outbred = "#065b9b", puma = "#c29a2b", SPF = "#83308c")
  )


## Generate pheatmap
# pdf(file="plots/sample_correlation.pdf", width=9,height=7)
pheat_samplecor <- pheatmap(cor(dat.logCPM),
                            annotation_col=annotation_col1,
                            annotation_names_col=F,
                            annotation_colors=ann_colors1,
                            annotation_row=annotation_col1,
                            annotation_names_row=F,
                            angle_col= "45",
                            border_color = NA,
                            show_rownames = F,
                            show_colnames = F
                            #cutree_cols = 3,
                            #cutree_rows=3
)
# dev.off()


#### 4. Generate PCAs, look for outliers and batch effects ####
## Generate PCA
counts.pca<-prcomp(t(dat.logCPM))
s<-summary(counts.pca)
s$importance
scores<-as.data.frame(counts.pca$x)

Mischief <- scores[12,]
meta_mischief <- meta[12,]

## By cell type:
pca <- ggplot(data=scores, aes(x=PC1, y=PC2, group=Cell_type)) + 
  geom_point(aes(color=Infection_status, shape=Cell_type, size=Cell_type)) +
  scale_color_manual(values=c("#00b6bd", "lightgrey")) +
  scale_size_manual(values=c(3,3)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  xlab ("PC1 (69.9%)") +
  ylab ("PC2 (21.2%)") +
  guides(colour = guide_legend(override.aes = list(size=3)))
pca
# pca + annotate("text", label=meta_mischief$population, x=Mischief$PC1, y=Mischief$PC2, size=3)
pca + annotate("text", label=meta_mischief$population, x=(-105), y=258, size=4)
# ggsave(file="plots/pca_all.eps")


### 5. Subset counts and metadata dataframes by cell type and infection status ####
## Subset metadata df

## No pooled samples
meta_nopool <- meta

## fibroblasts only
meta_fibro <- meta[meta$cell_type=="fibroblast",]

## PMBCs only
meta_pmbc <- meta[meta$cell_type=="PMBC",]

## uninfected only
meta_uninf <- meta[meta$status=="uninfected",]

## Subset counts df using dplyr
## fibroblasts
counts_genes <- counts[,1]
counts_f <- counts[,c(8:16)]
counts_fibro <- cbind(counts_genes, counts_f)

## PMBCs
counts_pmbc <- counts[,c(1:7)]

## uninfected
counts_uninf <- counts[,c(1:9,11,13,15)]


#### 6. Write tissue-specific, outlier-removed raw counts and metadata files ####
## Metadata
# write.table(meta,"results/edgeR/meta_all.txt", sep="\t", quote=F, row.names=F)
# write.table(meta_fibro,"results/edgeR/meta_fibro.txt", sep="\t", quote=F, row.names=F)
# write.table(meta_pmbc,"results/edgeR/meta_pmbc.txt", sep="\t", quote=F, row.names=F)
# write.table(meta_uninf,"results/edgeR/meta_uninf.txt", sep="\t", quote=F, row.names=F)

## Raw gene counts
# write.table(counts_fibro,"results/edgeR/counts_fibro.txt", sep="\t", quote=F, row.names=F)
# write.table(counts_pmbc,"results/edgeR/counts_pmbc.txt", sep="\t", quote=F, row.names=F)
# write.table(counts_uninf,"results/edgeR/counts_uninf.txt", sep="\t", quote=F, row.names=F)

# Full table, all tissues
# write.table(counts, "results/edgeR/counts_all.txt", sep="\t", quote=F, row.names = F)
