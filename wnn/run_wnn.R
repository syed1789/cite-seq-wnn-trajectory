#!/usr/bin/env Rscript
# CITE-seq WNN integration (RNA + ADT) on the public 10x "10k PBMC" dataset.
# Cell Ranger's *output* matrix is used as-provided -- see README for why
# Cell Ranger itself isn't re-run here.
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(1)
out_dir <- "/work/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data <- Read10X_h5("/work/data/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5")
cat("Modalities in h5:", paste(names(data), collapse = ", "), "\n")

obj <- CreateSeuratObject(counts = data[["Gene Expression"]], min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

adt_counts <- data[["Antibody Capture"]][, colnames(obj)]
obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)

cat(sprintf("Before QC: %d cells, %d genes, %d ADT features\n",
            ncol(obj), nrow(obj), nrow(obj[["ADT"]])))

obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
cat(sprintf("After QC: %d cells\n", ncol(obj)))

# RNA
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA(npcs = 30, verbose = FALSE)

# ADT (CLR-normalized, margin=2 -- across cells, the standard CITE-seq convention)
DefaultAssay(obj) <- "ADT"
VariableFeatures(obj) <- rownames(obj[["ADT"]])
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2) |>
  ScaleData() |>
  RunPCA(reduction.name = "apca", npcs = min(18, nrow(obj[["ADT"]]) - 1), verbose = FALSE)

# Weighted Nearest Neighbor integration
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:30, 1:min(18, nrow(obj[["ADT"]]) - 1)),
  modality.weight.name = "RNA.weight"
)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, resolution = 0.6, verbose = FALSE)

cat(sprintf("WNN clusters: %d\n", length(unique(Idents(obj)))))

p1 <- DimPlot(obj, reduction = "wnn.umap", label = TRUE, repel = TRUE) + ggtitle("WNN clusters (RNA + ADT)")
p2 <- FeaturePlot(obj, features = c("adt_CD3-TotalSeqB", "adt_CD19-TotalSeqB", "adt_CD14-TotalSeqB"), reduction = "wnn.umap", ncol = 3)
ggsave(file.path(out_dir, "wnn_umap_clusters.png"), p1, width = 7, height = 6, dpi = 150)
ggsave(file.path(out_dir, "wnn_umap_adt_markers.png"), p2, width = 12, height = 4, dpi = 150)

# --- Export for the trajectory (Monocle3) container ---
# Sparse RNA counts, restricted to variable features, to keep the handoff small.
vf <- VariableFeatures(obj, assay = "RNA")
counts_sub <- GetAssayData(obj, assay = "RNA", layer = "counts")[vf, ]
writeMM(counts_sub, file.path(out_dir, "rna_counts_hvg.mtx"))
fwrite(data.table(gene_short_name = rownames(counts_sub)), file.path(out_dir, "genes_hvg.csv"))
fwrite(data.table(barcode = colnames(counts_sub)), file.path(out_dir, "barcodes.csv"))

cell_meta <- data.table(
  barcode = colnames(obj),
  wnn_cluster = as.character(Idents(obj)),
  wnnUMAP_1 = Embeddings(obj, "wnn.umap")[, 1],
  wnnUMAP_2 = Embeddings(obj, "wnn.umap")[, 2],
  adt_CD3 = FetchData(obj, "adt_CD3-TotalSeqB")[, 1],
  adt_CD4 = FetchData(obj, "adt_CD4-TotalSeqB")[, 1],
  adt_CD8a = FetchData(obj, "adt_CD8a-TotalSeqB")[, 1],
  adt_CD19 = FetchData(obj, "adt_CD19-TotalSeqB")[, 1],
  adt_CD14 = FetchData(obj, "adt_CD14-TotalSeqB")[, 1],
  adt_CD45RA = FetchData(obj, "adt_CD45RA-TotalSeqB")[, 1]
)
fwrite(cell_meta, file.path(out_dir, "cell_metadata.csv"))

cat("Exported rna_counts_hvg.mtx, genes_hvg.csv, barcodes.csv, cell_metadata.csv for the trajectory step\n")
cat("Done.\n")
