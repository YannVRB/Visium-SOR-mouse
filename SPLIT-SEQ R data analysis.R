###This assumes you split-pipe twice as previously explained. Outputs from split-pipe are stored in folder called DGE_filtered

#Loading librairies
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

#Reading in data
mat_path <- "DGE_filtered/"

#The Seurat function ReadParseBio() provides a convenient way to read your expression matrix into R using the DGE folder path as input.
mat <- ReadParseBio(mat_path)
fig_path <- "Figures/"
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", "png"),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", "pdf"),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}
SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}
ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}
cleanUMAP <- function(plot_obj, dark = FALSE, axis_title_size = 16,
                      plot_title = "") {
  if (dark == TRUE) {
    plot_obj <- plot_obj + labs(title = element_text(plot_title), x = "UMAP 1", y = "UMAP 2") +
      theme(plot.background = element_rect(fill = "black",
                                           colour = "black", size = 0.5, linetype = "solid"),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.line.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), axis.line.y = element_blank(),
            axis.title = element_text(size = axis_title_size, color = "white"),
            legend.text = element_text(color = c("white"), size = 10)) +
      ggtitle(plot_title)
  } else {
    plot_obj <- plot_obj + labs(title = element_text(plot_title), x = "UMAP 1", y = "UMAP 2") +
      theme(plot.title = element_text(plot_title), axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), axis.line.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title = element_text(size = axis_title_size))
  }
}

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

