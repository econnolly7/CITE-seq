library(Seurat)
library(cowplot
        )
setwd("/Users/erinconnolly/Downloads/")

test<-readRDS('seuratObj_10x_sc5p_v2_hs_PBMC.rds')
GSM4712907<-Read10X("Sample2")

CR<- CreateAssayObject(
  counts = GSM4712895$"Antibody Capture")

GSM4712907 <- CreateSeuratObject(
  counts = GSE157007_RAW$"Gene Expression",project="GSM4712907",
  assay = "RNA")

GSM4982334[["ADT"]] <- CR

# Validate that the object now contains multiple assays
Assays(dn_2811)

# Extract a list of features measured in the ADT assay
rownames(dn_2811[["ADT"]])

# Note that we can easily switch back and forth between the two assays to specify the default
# for visualization and analysis

# List the current default assay
DefaultAssay(dn_2811)

# Switch the default to ADT
DefaultAssay(dn_2811) <- "ADT"
DefaultAssay(dn_2811)

#normalize RNA-seq
DefaultAssay(dn_2811) <- "RNA"

# perform visualization and clustering steps
dn_2811 <- NormalizeData(dn_2811)
dn_2811 <- FindVariableFeatures(dn_2811)
dn_2811 <- ScaleData(dn_2811)
dn_2811 <- RunPCA(dn_2811, verbose = FALSE)
dn_2811 <- FindNeighbors(dn_2811, dims = 1:30)
dn_2811 <- FindClusters(dn_2811, resolution = 0.8, verbose = FALSE)
dn_2811 <- RunUMAP(dn_2811, dims = 1:30)
DimPlot(dn_2811, label = TRUE)

#normalize ADT
DefaultAssay(dn_2811) <- "ADT"
dn_2811 <- NormalizeData(dn_2811, normalization.method = "CLR", margin = 2)

# Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(dn_2811) <- "ADT"
p1 <- FeaturePlot(dn_2811, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(dn_2811) <- "RNA"
p2 <- FeaturePlot(dn_2811, "CD19") + ggtitle("CD19 RNA")

# place plots side-by-side
p1 | p2

