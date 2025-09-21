# For efficiency reasons this is handled via an R script - a jupyter notebook with visual representation of filtering is available

#### 1) Setup ####
# Imports
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1)
threads <- parallel::detectCores() - 1
addArchRThreads(threads = threads)
addArchRLocking(FALSE)
addArchRGenome("hg38")

#Paths
collab_dir <- "../../../data/data_raw/collaborators"
pub_dir   <- "../../../data/data_raw/granja/hg38"

#### 2) Set up samples ####
pub_frags    <- list.files(pub_dir,
                           pattern = "\\.fragments\\.tsv\\.gz$",
                           full.names = TRUE)
pub_raw <- sub("_hg38\\.fragments\\.tsv\\.gz$", "", basename(pub_frags))
pub_suffixes <- sub(".*(_D\\d+T\\d+)$", "\\1", pub_raw)
pub_names <- paste0("granja", pub_suffixes)

# Combine fragments files from data sources for arrow creation
all_frags  <- c(
  file.path(collab_dir, "collab_donor0_fragments.tsv.gz"),
  file.path(collab_dir, "collab_donor1_fragments.tsv.gz"),
  pub_frags
)
all_names  <- c("collab_0", "collab_1", pub_names)


#### 3) Create Arrow Files needed for ArchR
arrow_files <- createArrowFiles(
  inputFiles      = all_frags,
  sampleNames     = all_names,
  minFrags        = 1000,
  minTSS          = 8,
  addTileMat      = TRUE,
  addGeneScoreMat = TRUE,
  force           = FALSE
)

#### 4) Initialize Project & QC ####
proj <- ArchRProject(
  ArrowFiles      = arrow_files,
  outputDirectory = "ArchR_Collab_vs_Pub",
  copyArrows      = FALSE
)


# Quick sanity checks
table(proj$Sample)
length(proj$cellNames)

# Pull out per‐cell Sample & TSSEnrichment
df <- getCellColData(proj, select = c("Sample","TSSEnrichment"))

# Compute median and mean per sample
tss_median <- aggregate(TSSEnrichment ~ Sample, data = df, FUN = median)
tss_mean   <- aggregate(TSSEnrichment ~ Sample, data = df, FUN = mean)

# Merge into one data.frame
tss_stats <- merge(tss_median, tss_mean, by="Sample",
                   suffixes = c("_median","_mean"))
print(tss_stats)


# Filter all cells not passing QC. This is required since ArchR does not do this during arrow creation
minTSS   <- 8
minFrags <- 1000
goodCells <- proj$cellNames[
  proj$TSSEnrichment >= minTSS &
  proj$nFrags        >= minFrags
]

# This is required to 
proj <- subsetArchRProject(
  proj,
  cells          = goodCells,
  dropCells      = TRUE,
  outputDirectory = "ArchR_Collab_vs_Pub_QCFiltered",
  force          = TRUE
)



#### 4a) Doublet detection & filtering ####
proj <- addDoubletScores(proj)
proj <- filterDoublets(proj)

proj <- subsetArchRProject(
  proj,
  cells          = proj$cellNames,
  dropCells      = TRUE,
  outputDirectory = "ArchR_Collab_vs_Pub_QC1",
  force          = TRUE
)

#### Add Sample-level metadata ####
# Named vector mapping Sample → Tissue
all_names <- proj$Sample
tissue_lookup <- setNames(
  ifelse(grepl("^collab_", all_names), "MobilizedBlood", "BoneMarrow"),
  all_names
)

proj <- addCellColData(
  proj,
  data      = tissue_lookup[proj$Sample],
  cells     = proj$cellNames,
  name      = "Tissue",
  force     = TRUE
)



### Clustering ###

# LSI
proj <- addIterativeLSI(
  ArchRProj  = proj,
  useMatrix  = "TileMatrix",
  name       = "IterLSI1",
  force      = TRUE
)

proj <- addIterativeLSI(
  ArchRProj  = proj,
  useMatrix  = "TileMatrix",
  name       = "IterLSI",
  iterations = 3,
  dimsToUse  = keep_dims,
  force      = TRUE
)
                     

table(proj$Tissue)

proj <- addHarmony(
  proj,
  reducedDims = "IterLSI",
  name        = "Harmony",
  groupBy     = "Sample",
  force       = TRUE
)

proj <- addImputeWeights(proj, reducedDims = "Harmony")

proj <- addClusters(
  proj,
  reducedDims = "Harmony",
  name        = "Clusters",
  force       = TRUE
)

# with batch correction
proj <- addUMAP(
  proj,
  reducedDims = "Harmony",
  name        = "UMAP",
  force       = TRUE
)

# w/O batch correction
proj <- addUMAP(
  proj,
  reducedDims = "IterLSI",
  name        = "UMAPX",
  force       = TRUE
)


#### Save ####
saveArchRProject(proj, load = FALSE)