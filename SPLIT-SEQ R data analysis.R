###This assumes you ran split-pipe twice as previously explained. Outputs from split-pipe are stored in folder called DGE_filtered.

#Loading librairies
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

#Reading in data
mat_path <- "DGE_filtered/"

#The Seurat function ReadParseBio() provides a convenient way to read your expression matrix into R using the DGE folder path as input.
mat <- ReadParseBio(mat_path)

#Check to see if empty gene names are present, add name if so
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

#Read in cell meta data
meta_csv <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

#Create Seurat object
seur_obj <- CreateSeuratObject(mat, min_genes = 100, min_cells = 100,
                               names.feild = 0, meta.data = meta_csv)

#Check to see if empty gene names are present, add name if so
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

#Assign the meta slot to a separate variable. This so we can check the meta data object to make sure it is correct before assigning it back to our Seurat object.
meta <- seur_obj@meta.data
#Initialize an empty column where we will assign new values based on sample category.
meta$sample_category <- ""
#Here we are filling in our new category ONLY where HC and SOR are present in the original "sample" column that was in the original Seurat object. Note that "sample_category" and "sample" are two different columns. We are using the information from "sample" so we can assign new values to the column "sample_category".
meta$sample_category[grep("HC", meta$sample)] <- "HC"
meta$sample_category[grep("SOR", meta$sample)] <- "SOR"
#Verify that the new samples lineup before assigning back.
head(meta$sample_category)
seur_obj@meta.data <- meta

#Normalizing the data
seur_obj <- NormalizeData(seur_obj)

#Identification of highly variable features
seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)

#Scaling the data
seur_obj <- ScaleData(seur_obj, verbose = TRUE)

#Perform linear dimensional reduction
seur_obj <- RunPCA(seur_obj, npcs = 50, verbose = FALSE)

#Cluster the cells
seur_obj <- FindNeighbors(seur_obj, reduction = "pca", dims = 1:20)
seur_obj <- FindClusters(seur_obj, resolution = 0.3)

Idents(seur_obj) <- seur_obj@meta.data$integrated_snn_res.0.3
seur_obj@meta.data$
  seur_obj@meta.data
Idents(seur_obj) <- seur_obj@meta.data$RNA_snn_res.0.3
#Reorder clusters according to their similarity
seur_obj <- BuildClusterTree(seur_obj, reorder.numeric = TRUE)

#Run non-linear dimensional reduction (UMAP/tSNE)
seur_obj <- RunUMAP(seur_obj, reduction = "pca", dims = 1:20)
DimPlot(seur_obj, reduction = "umap", group.by = "orig.ident")

#Assessing the expression of bmarkers in specific clusters
markers <- c("Ahcyl2", "Gm10754", "Nectin3", "Rgs20", "Gad2", "Mog", "Tshz2",
"Pdgfra", "Il16", "Cp", "Inpp5d", "Filip1l", "Trp73", "Slc6a13")
FeaturePlot(pbmc, features = markers)

#Generating a list where each cell type is assigned multiple clusters so they can be merged
new_ids <- c("Dentate Gyrus", "CA1", "CA3", "Astrocytes", "Interneurons",
  "Oligodendrocytes", "Excitatory cortical neurons", "CA1-ProS", "OPC", "CA2", "Endothelial", "Microglia", "SMC-Peri", "Cajal–Retzius", "VLMC")

new_id_list <- list(Dentate Gyrus = c(0,2,16), CA1 = 1,
  CA3 = 3, Astrocytes = c(4,10), Interneurons = c(5,13), Oligodendrocytes = 6, Excitatory cortical neurons = c(7,8),
  CA1-ProS = 9, OPC = 11, CA2 = 12, Endothelial = 14, Microglia = 15, SMC-Peri = 17, Cajal–Retzius = 18, VLMC = 19)

for (i in 1:length(new_id_list)) {
  ind <- which(seur_obj@meta.data$RNA_snn_res.0.3 %in% new_id_list[[i]])
  seur_obj@meta.data$collapsed[ind] <- names(new_id_list)[i]
}

seur_obj@meta.data$collapsed <- factor(
  seur_obj@meta.data$collapsed, levels = names(new_id_list), ordered = TRUE)
Idents(seur_obj) <- pbmc@meta.data$collapsed

names(new_ids) <- levels(seur_obj)
seur_obj <- RenameIdents(seur_obj, new_ids)

#Differential gene expression analysis in Oligodendrocytes
seur_obj$celltype.sample_category <- paste(Idents(seur_obj), seur_obj$sample_category, sep = "_")
seur_obj$celltype <- Idents(seur_obj)
Idents(seur_obj) <- "celltype.sample_category"

oligo.DEG <- FindMarkers(seur_obj, ident.1 = "Oligodendrocytes_SOR", ident.2 = "Oligodendrocytes_HC", min.pct=0, logfc.threshold=0)
