# CITE-seq WNN Integration + Monocle3 Trajectory

**Author:** Syed Nurul Hasan, Ph.D.
**Contact:** nayanchoton1789@gmail.com | [ORCID: 0000-0002-4564-872X](https://orcid.org/0000-0002-4564-872X)

---

## Why this project

My `scATAC-WNN-multimodal` repo already demonstrates Weighted Nearest Neighbor
integration for RNA + ATAC. Several roles I've applied to specifically named
**CITE-seq** (RNA + surface protein) and **Monocle** for trajectory inference
— related techniques, no public repo. This one integrates RNA and antibody
(ADT) modalities from a real public CITE-seq dataset via WNN, then runs
Monocle3 pseudotime on top of that integrated embedding.

**Scope note, stated plainly:** Cell Ranger's *output* matrix is used
as-provided — Cell Ranger itself was not re-run here. There's no legally
redistributable Cell Ranger Docker image (10x's EULA doesn't permit
third-party redistribution), so re-running it isn't practical in a public
repo; the WNN integration, clustering, and trajectory analysis downstream of
it are all real, executed work.

## Pipeline

```
10x "10k PBMC" CITE-seq (RNA + ADT), Cell Ranger output used as-provided
        │
        ▼
1. QC + per-modality normalization  (Seurat, R 4.4 container)
   ├── RNA: NormalizeData, HVG selection, PCA
   └── ADT: CLR normalization, PCA ("apca")
        │
        ▼
2. Weighted Nearest Neighbor integration (FindMultiModalNeighbors)
   ├── Per-cell modality weights (RNA vs. ADT)
   ├── WNN graph → UMAP
   └── Louvain clustering on the WNN graph
        │
        ▼
3. Export: HVG counts (.mtx) + cell metadata (WNN UMAP coords, ADT levels)
        │
        ▼
4. Monocle3 (separate R 3.6-vs-4.x-conflict-free container, from-source build)
   ├── Rebuild a CellDataSet from the exported counts
   ├── Inject the precomputed WNN UMAP as Monocle3's reduced-dim space
   ├── learn_graph + order_cells (root = highest-CD45RA cell, naive marker)
   └── Pseudotime plot + per-cluster mean pseudotime
```

## Dataset

`pbmc_10k_protein_v3` — 10x Genomics' public "10k PBMCs from a Healthy
Donor" CITE-seq dataset (Cell Ranger 3.0.0, 17 TotalSeq-B antibodies),
downloaded directly from 10x's public sample repository (no login/EULA
click-through required for this demo dataset).

## Results (this run)

- 7,856 cells before QC → **7,722 cells** after (`nFeature_RNA` 200–6000, `percent.mt` < 15%).
- WNN clustering (Louvain, resolution 0.6) → **19 clusters**, cleanly separable
  by ADT surface markers with no RNA-only ambiguity:
  - Clusters 0, 11 — high CD14 (monocytes)
  - Clusters 1, 2, 5, 6, 7, 9, 15 — high CD3, low CD19/CD14 (T cells)
  - Clusters 4, 16, 17 — high CD19, low CD3 (B cells)
  - See [`results/cell_metadata.csv`](results/cell_metadata.csv) for full per-cluster ADT means.
- [`results/wnn_umap_clusters.png`](results/wnn_umap_clusters.png),
  [`results/wnn_umap_adt_markers.png`](results/wnn_umap_adt_markers.png)
- **Monocle3 pseudotime**, rooted at the highest-CD45RA cell (naive-lymphocyte
  marker): 666 of 7,722 cells (8.6%) received finite pseudotime, range 0–1.81.
  `learn_graph` partitions the WNN UMAP into multiple disconnected principal
  graphs, and `order_cells` only propagates pseudotime within the root's
  partition — everything else is legitimately undefined (`Inf`), not a bug.
  **What actually happened, reported precisely rather than assumed:** the
  root's partition turned out to be almost entirely the B-cell lineage —
  clusters 4/16/17 (CD19-high) contribute 645 of the 666 finite-pseudotime
  cells — plus a sliver of cluster 8 and just the single root cell itself
  from cluster 1 (T cells). CD45RA marks naive T *and* B lymphocytes, and
  this pure "highest single-cell value" root heuristic happened to land in
  the B-cell partition rather than the T-cell one. A more robust
  root-selection strategy would pick the target *partition* first (e.g. by
  mean CD45RA per partition), then the highest-marker cell within it — noted
  here rather than silently re-picking a root to make the story cleaner.
  [`results/monocle3_pseudotime.png`](results/monocle3_pseudotime.png),
  [`results/pseudotime.csv`](results/pseudotime.csv).

## Why two containers (and the real debugging trail to get Monocle3 working)

`wnn/` runs modern Seurat (needs R ≥ 4.x for `FindMultiModalNeighbors`).
`trajectory/` runs Monocle3, which took five real attempts to get installed
cleanly — recorded here instead of smoothed over, since it's a legitimate
environment-debugging trail:

1. bioconda's `r-monocle3` package (any version/build pin) turned out to be
   permanently unsolvable — its recipe pulls in `r-spelling` → `r-hunspell`
   version constraints that no longer resolve against current channel
   snapshots.
2. Installing from source (`remotes::install_github`, `dependencies=TRUE`)
   instead pulled in monocle3's full *Suggests* list, several of which
   (`sf`, `BPCells`) don't build in a minimal container.
3. `dependencies=FALSE` swung too far the other way — it also skipped
   genuinely required *Imports* (`RANN`, `uwot`, `Rtsne`, etc.), so the
   install failed on missing hard dependencies instead.
4. Pre-installing a manually-curated dependency list via conda got further,
   but master/HEAD monocle3 (v1.4+) added `BPCells`, `ggdist`, and `sf` as
   **hard** Imports — `BPCells` in particular is fetched from a separate
   GitHub fork with its own native-library build requirements, not
   available as a conda package.
5. **What worked:** pinning to monocle3 `v1.3.1` (the last tag before
   `BPCells` became a hard dependency) and installing every package in its
   *actual* `Imports:` list — read directly from that tag's `DESCRIPTION`
   file, not guessed — via conda, then `remotes::install_github(...,
   ref='v1.3.1', dependencies=FALSE)` for monocle3 and leidenbase themselves.

Both `wnn/Dockerfile` and `trajectory/Dockerfile` also needed
`ARG MAMBA_DOCKERFILE_ACTIVATE=1` added after `FROM` — without it, the
`mambaorg/micromamba` base image's installed R isn't on `PATH` for
subsequent `RUN` layers, which silently breaks the CellChat-install pattern
already sitting in my `heart-nextflow-r/Dockerfile` too (flagged there
separately, not fixed in this repo).

The two containers hand off through flat files (`results/*.mtx`,
`results/*.csv`) instead of a shared R session, since Seurat's WNN needs
modern R and Monocle3 v1.3.1 does not require it (and pinning both to the
same R version wasn't worth re-solving on top of everything above).

## Running it

```bash
# 1. Fetch the public CITE-seq dataset
curl -L --fail -o data/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5 \
  https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5

# 2. WNN integration
docker build -t cite-seq-wnn:latest wnn/
docker run --rm -v "$(pwd)":/work cite-seq-wnn:latest Rscript /work/wnn/run_wnn.R

# 3. Monocle3 trajectory on the WNN embedding
docker build -t monocle3-trajectory:latest trajectory/
docker run --rm -v "$(pwd)":/work monocle3-trajectory:latest Rscript /work/trajectory/run_trajectory.R
```

## Key Features

- Real 10x CITE-seq public dataset, real WNN integration (RNA + ADT)
- Clusters that separate cleanly by surface-protein identity, not just RNA
- Monocle3 pseudotime learned directly on the WNN embedding, not RNA-only UMAP
- Honest handling of a real bioconda packaging failure (documented, not hidden)
- Honest disclosure of the Cell Ranger scope boundary

## References

1. Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal
   single-cell data. *Cell*. 2021;184(13):3573-3587.e29. (Seurat WNN)
2. Stoeckius M, Hafemeister C, Stephenson W, et al. Simultaneous epitope and
   transcriptome measurement in single cells. *Nature Methods*.
   2017;14(9):865-868. (CITE-seq)
3. Cao J, Spielmann M, Qiu X, et al. The single-cell transcriptional
   landscape of mammalian organogenesis. *Nature*. 2019;566(7745):496-502.
   (Monocle3)

## License

MIT — see [LICENSE](LICENSE).

## Running on HPC (SLURM) and AWS

Includes `run_analysis.sh` (Docker or Apptainer), a SLURM job script, and a one-shot EC2 launcher. See [HPC_AWS.md](HPC_AWS.md) for steps and what was verified.
