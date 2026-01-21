# ==============================================================================
# Peak Calling Pipeline
# ==============================================================================
# 
# This script performs pseudobulk peak calling using MACS2, filters cell types 
# by minimum representation, and creates integrated peak-based embeddings.
#
# Input:  Annotated ArchR project with unified cell type labels
# Output: Project with reproducible peak set and PeakMatrix


#### 1) Setup ####
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1)
threads <- parallel::detectCores() - 1
addArchRThreads(threads = threads)
addArchRLocking(FALSE)
addArchRGenome("hg38")
macs2path <- path.expand("~/miniconda3/bin/macs2")

# 2) Load annotated project
proj <- loadArchRProject("ArchR_Collab_vs_Pub_Annotated")

#--- CONFIG ---------------------------------------------------------------
# Cell types with fewer than minCells are excluded to ensure robust pseudobulk 
# peak calling. Excluded populations: B lineage, T lineage, monocytes, MDP.

minCells <- 250
oldLabel <- "CellType_Unified"
newProj  <- "ArchR_Collab_vs_Pub_Peaks"
#-------------------------------------------------------------------------

# 3) Count cells per cell type
ct_vec <- getCellColData(proj)[[oldLabel]]
celltype_n <- table(ct_vec, useNA = "no")
cat("\n=== Cell counts per unified cell type ===\n")
print(celltype_n)

# 4) Identify which cell types meet the threshold
good_celltypes <- names(celltype_n)[celltype_n >= minCells]
cat("\n=== Cell types with >=", minCells, "cells ===\n")
print(good_celltypes)

# 5) Keep only cells whose cell type is in good_celltypes
keep_cells <- getCellNames(proj)[ct_vec %in% good_celltypes]
cat("\nCells before filtering:", length(getCellNames(proj)), "\n")
cat("Cells after filtering:", length(keep_cells), "\n")

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = keep_cells,
  dropCells = TRUE,
  outputDirectory = "ArchR_Collab_vs_Pub_MinCells_Temp",
  force = TRUE
)

# 6) Create tissue_celltype label for peak calling
meta <- getCellColData(proj)
tissue_celltype <- paste0(meta$Tissue, "_", meta$CellType_Unified)
proj$Tissue_CellType <- tissue_celltype

# Check the groupings
cat("\n=== Tissue_CellType groups ===\n")
print(sort(table(proj$Tissue_CellType), decreasing = TRUE))

# 7) Create pseudobulks
cat("\n=== Adding group coverages ===\n")
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Tissue_CellType",
  force = TRUE
)

# 8) Call peaks using tissue-specific cell type labels
cat("\n=== Calling peaks with MACS2 ===\n")
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Tissue_CellType",
  pathToMacs2 = macs2path,
  force = TRUE
)

# 9) Add PeakMatrix
cat("\n=== Creating peak matrix ===\n")
proj <- addPeakMatrix(proj, force = TRUE)



# 10) Save final project
cat("\n=== Saving project ===\n")
saveArchRProject(
  proj,
  outputDirectory = newProj,
  load = FALSE
)

cat("\nProject saved to:", newProj, "\n")
cat("Total cells in final project:", length(getCellNames(proj)), "\n")