#!/usr/bin/env Rscript
# Monocle3 pseudotime/trajectory on top of the WNN (RNA+ADT) integration from
# ../wnn/run_wnn.R. Runs in a separate R 3.6 container because bioconda's
# r-monocle3 build is pinned to R 3.6, while a modern Seurat (needed for
# FindMultiModalNeighbors/WNN) needs R >= 4.x -- a real dependency conflict,
# handled here by handing off flat files between two containers rather than
# forcing an unstable single environment.
suppressPackageStartupMessages({
  library(monocle3)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

in_dir <- "/work/results"
out_dir <- "/work/results"

# run_wnn.R writes counts_sub as genes (rows) x cells (cols) -- the same
# orientation Seurat and Monocle3 both use -- so no transpose here.
counts <- readMM(file.path(in_dir, "rna_counts_hvg.mtx"))
genes <- fread(file.path(in_dir, "genes_hvg.csv"))
barcodes <- fread(file.path(in_dir, "barcodes.csv"))
cell_meta <- fread(file.path(in_dir, "cell_metadata.csv"))
stopifnot(identical(barcodes$barcode, cell_meta$barcode))

rownames(counts) <- genes$gene_short_name
colnames(counts) <- barcodes$barcode

gene_metadata <- data.frame(gene_short_name = genes$gene_short_name, row.names = genes$gene_short_name)
cell_metadata <- as.data.frame(cell_meta)
rownames(cell_metadata) <- cell_metadata$barcode

cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# Use the WNN UMAP already computed from RNA+ADT (Seurat) as Monocle3's
# reduced-dimension space, instead of recomputing UMAP from RNA alone --
# the whole point of this repo is that the protein modality sharpens the
# embedding, so the trajectory should be learned on that embedding too.
reducedDims(cds)$UMAP <- as.matrix(cell_metadata[, c("wnnUMAP_1", "wnnUMAP_2")])

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# Root at the cell with the highest CD45RA (canonical naive-lymphocyte
# surface marker) -- a defensible, non-interactive, objective root-selection
# rule for order_cells() rather than a manually clicked point.
root_cell <- cell_metadata$barcode[which.max(cell_metadata$adt_CD45RA)]
cds <- order_cells(cds, root_cells = root_cell)

p <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 3
) + ggtitle("Monocle3 pseudotime on the WNN (RNA+ADT) embedding")
ggsave(file.path(out_dir, "monocle3_pseudotime.png"), p, width = 7, height = 6, dpi = 150)

pt <- data.table(
  barcode = colnames(cds),
  pseudotime = pseudotime(cds),
  wnn_cluster = cell_metadata$wnn_cluster
)
fwrite(pt, file.path(out_dir, "pseudotime.csv"))

cat(sprintf("Root cell: %s (adt_CD45RA = %.2f)\n", root_cell, max(cell_metadata$adt_CD45RA)))
cat("Pseudotime range:", paste(round(range(pt$pseudotime[is.finite(pt$pseudotime)]), 2), collapse = " to "), "\n")
cat("Mean pseudotime by WNN cluster:\n")
print(pt[, .(mean_pseudotime = mean(pseudotime[is.finite(pseudotime)])), by = wnn_cluster][order(mean_pseudotime)])

cat("Done.\n")
