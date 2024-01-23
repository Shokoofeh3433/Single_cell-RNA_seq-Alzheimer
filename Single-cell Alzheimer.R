setwd("E:/Single Cell/ALZ_Blood")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(metap)
library(celldex)
library(scRNAseq)
library(ccImpute)
library(SingleR)
library(ExperimentHub)
library(scuttle)
library(SingleCellExperiment)
library(DESeq2)
library(scran)
library(scDblFinder)
library(BiocSingular)


# get data location
dirs <- list.dirs(path = "Data/", recursive = F, full.names = F)

for(x in dirs){
  name <- gsub('_f_bc_mx', '', x)
  
  cts <- ReadMtx(mtx = paste0('Data/',x,'/matrix.mtx.gz'),
                 features = paste0('Data/',x,'/features.tsv.gz'),
                 cells = paste0('Data/',x,'/barcodes.tsv.gz'))
  
  # create Seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}


# merge data sets
merged_seurat <- merge(AD1, y = c(AD2, AD3, NC1, NC2),
                       add.cell.ids = ls()[1:5],
                       project = 'Alzheimer-Covid-PBMC')
saveRDS(merged_seurat, file = 'merged_seurat_Breast_cancer')


# Reading merged RDS file
view(merged_seurat@meta.data)
merged_seurat$sample <- rownames(merged_seurat@meta.data)
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = "sample", into = c("sample", "barcode"), sep = "_")

# Mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
View(merged_seurat@meta.data)

# Quality Control
p1 <- VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
p2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 700 &
                                   nFeature_RNA < 2500 &
                                   nCount_RNA > 1500 & nCount_RNA < 9000 &
                                   mitoPercent > 2 & mitoPercent < 6)

# Removing doublets
dbl <- GetAssayData(merged_seurat_filtered, slot = "counts", assay = "RNA")
scores <- computeDoubletDensity(dbl)
merged_seurat_filtered$dbl.score <- scores
dbl.calls <- doubletThresholding(data.frame(score=merged_seurat_filtered$dbl.score),
                                 method="griffiths", returnType="call")
summary(dbl.calls)

merged_seurat_filtered2$doublets <- dbl.calls

merged_seurat_filtered <- subset(merged_seurat_filtered, subset = doublets == "singlet")


# Normalization and Scaling
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered, normalization.method = "LogNormalize", scale.factor = 1e4)
merged_seurat_filtered <- FindVariableFeatures(object =  merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object =  merged_seurat_filtered, vars.to.regress = "mitoPercent")
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
p3 <- ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered, resolution = 0.1)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)
saveRDS(merged_seurat_filtered, file = "merged_seurat_filtered")

p4 <- DimPlot(object = merged_seurat_filtered, reduction = "umap", group.by = "sample",
              cols = c("red", "green", "blue", "purple", "orange"), raster = FALSE)

# Finding doublets
dbl <- GetAssayData(merged_seurat_filtered, slot = "counts", assay = "RNA")
scores <- computeDoubletDensity(dbl)
merged_seurat_filtered$dbl.score <- scores
plot(merged_seurat_filtered$dbl.score, pch = 19, cex = 1, col = "blue")

dbl.calls <- doubletThresholding(data.frame(score=merged_seurat_filtered$dbl.score),                                 method="griffiths", returnType="call")
summary(dbl.calls)
merged_seurat_filtered$doublets <- dbl.calls
merged_seurat_filtered <- subset(merged_seurat_filtered, subset = doublets == "singlet")

# Correcting batch effect
obj.list <- SplitObject(merged_seurat_filtered, split.by = "sample")
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]]) 
}

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors)
saveRDS(seurat.integrated, file = "seurat_integrated")
seurat.integrated <- ScaleData(object =  merged_seurat_filtered, vars.to.regress = "mitoPercent")
seurat.integrated <- RunPCA(object = seurat.integrated)
ElbowPlot(seurat.integrated)
seurat.integrated <- FindNeighbors(object = seurat.integrated, dims = 1:15)
seurat.integrated <- FindClusters(object = seurat.integrated, resolution = 0.1)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:15)

p5 <- DimPlot(object = seurat.integrated, reduction = "umap", group.by = "sample",
              cols = c("red", "green", "blue", "purple", "orange"), raster = FALSE)

p6 <- DimPlot(object = seurat.integrated, reduction = "umap", group.by = "seurat_clusters",
              raster = FALSE, label = TRUE)

n <- DimPlot(object = seurat.integrated, reduction = "umap", group.by = "seurat_clusters",
             raster = FALSE, label = TRUE)

grid.arrange(p4, p5, ncol = 2)


# annotation cell clusters
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))
ALZ_count <- GetAssayData(seurat.integrated, slot = "counts", assay = "RNA")

pred <- SingleR(ALZ_count, ref = ref, labels = ref$label.main)
seurat.integrated$cell_type <- pred$labels[match(rownames(seurat.integrated@meta.data), rownames(pred))]

p7 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "seurat_clusters", raster = FALSE, label = TRUE)
p8 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "SingleR", raster = FALSE, label = TRUE)

grid.arrange(p7, p8, ncol = 2)
grid.arrange(p6, n, ncol = 2)

clusters <- DimPlot(seurat.integrated, reduction = "umap", group.by = "SingleR", raster = FALSE, label = TRUE)
condition <- DimPlot(seurat.integrated, reduction = "umap", group.by = "sample", raster = FALSE)
clusters|condition

seurat.integrated$condition <- "NA"
for (i in 1:20325) {
  if(seurat.integrated$sample[i] == "AD1" | seurat.integrated$sample[i] == "AD2" | seurat.integrated$sample[i] == "AD3")
    seurat.integrated$condition[i] <- "AD"
  else
    seurat.integrated$condition[i] <- "CTRL"
}

seurat.integrated$clus.cnd <- paste0(seurat.integrated$cell_type, '_', seurat.integrated$condition)
Idents(seurat.integrated) <- seurat.integrated$clus.cnd

DEG_T <- FindMarkers(seurat.integrated, ident.1 = "T_cells_AD", ident.2 = "T_cells_CTRL")
dim(DEG_T)

DEG_T_f <- subset(DEG_T, subset = p_val_adj < 0.05 &
                    (avg_log2FC < -0.01 | avg_log2FC > 0.01))

dim(DEG_T_f)
write.csv(DEG_T_f, file = "Result mohammad/T-cell.csv")
