###

# This code is used to explore gene expression data and conduct expression tests in edgeR

# Dataset: endogenous/exogenous FeLV in vitro gene expression (fibroblasts infected and uninfected, PMBCs uninfected only)

#### This code will do the following ####
# 1. Generate a filtered and normalized DGElist of all sample data
# 2. Generate a logCPM counts matrix
# 3. Use logCPM counts matrix to generate sample correlation pheatmap
# 4. Use logCPM counts matrix to generate PCAs to test for batch effects and outliers

# Preconditions: raw gene-level counts matrix (e.g. from STAR quantMode), metadata file with sample names that correspond to sample names in counts matrix

###

library(here)
library(tidyverse)
library(edgeR)
library(pheatmap)


#### 0. Get data ####
## counts: either full counts dataset, or restricted to only cat-puma orthologs
dat <- read_delim("results/star_quant_bam/readcountsmatrix.txt", delim = "\t") %>% select(-DC2Pool, -DC3Pool)
counts <- dat %>% slice(-c(1:4, 31498)) #all counts, not just orthologs; remove sample totals rows
# counts <- read_delim("results/orthofinder/genecounts_orthologs.txt", delim = "\t") #all counts, not just orthologs; remove sample totals rows

## metadata
meta <- read_tsv("data/felv_metadata.tsv") %>% filter(id_inf != c("DC2Pool", "DC3Pool"))
cell_inf=paste(meta$cell_type, meta$status, sep="_") #grouping


#### 1. Create DGElist, filter by CPM, and calcNormFactors ####
## Create DGEList
dat.full <- DGEList(counts=counts[,2:17], group=cell_inf, genes = counts[,1])
dat.full$samples

## Filter
keep_counts<-rowSums(cpm(counts[,2:17])>1) >= 0.25*ncol(dat) #filter >1 cpm in >= 1/4 of samples
dat.filt<-dat.full[which(keep_counts==T), , keep.lib.sizes=FALSE]
dim(dat.full)
dim(dat.filt)

dat.norm<-calcNormFactors(dat.filt)
# dat.norm$samples #normalization factors and library size
plotMDS(dat.norm)


#### 2. Generate moderated, logCPM counts matrix for sample correlation and PCA ####
dat.logCPM <- cpm(dat.norm, log=TRUE, prior.count = 1, normalized.lib.sizes = TRUE)


#### 3. Generate sample correlation matrix and pheatmap ####
## annotate columns by tissue, treatment, and time
Cell_type=meta$cell_type
Infection_status=meta$status
Population=meta$population

annotation_col1 <- data.frame(
  Cell_type,
  Infection_status,
  Population
)
rownames(annotation_col1) <- colnames(dat[,c(2:17)])

ann_colors1 = list(
  Cell_type =c(fibroblast="#008080",PBMC="#ef6079"),
  Infection_status = c(infected = "#00b6bd", uninfected = "lightgrey"),
  Population = c(outbred = "#065b9b", puma = "#c29a2b", SPF = "#83308c")
  )


## Generate pheatmap
# pdf(file="plots/sample_correlation.pdf", width=9,height=7)
# pdf(file="plots/sample_correlation_orthologs.pdf", width=9,height=7)
pheatmap(cor(dat.logCPM),
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

meta_mischief <- meta[13,] #extract puma sample to label pca

## By cell type:
pca <- ggplot(data=scores, aes(x=PC1, y=PC2, group=Cell_type)) + 
  geom_point(aes(color=Infection_status, shape=Cell_type, size=Cell_type)) +
  scale_color_manual(values=c("#00b6bd", "lightgrey")) +
  scale_size_manual(values=c(3,3)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  xlab ("PC1 (71.9%)") +
  ylab ("PC2 (19.8%)") +
  guides(colour = guide_legend(override.aes = list(size=3)))
pca
pca + annotate("text", label=meta_mischief$population, x=(-85), y=246, size=4)
# ggsave(file="plots/pca_all.eps")
# ggsave(file="plots/pca_all_orthologs.eps")
