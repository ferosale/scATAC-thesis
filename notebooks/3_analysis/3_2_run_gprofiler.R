#===========================================================================
# GO Enrichment Analysis
#===========================================================================

# Set paths
go_dir <- "/work/project/ladcol_013/bt-atac/notebooks/1_processing/test/ArchR_HSC_Only/Plots"

cat("=== Running gprofiler2 GO Analysis (Headless) ===\n")
cat("Directory:", go_dir, "\n\n")

# Load library
library(gprofiler2)

#===========================================================================
# 1. READ GENE LISTS & BACKGROUND
#===========================================================================

# Read the significant gene lists
genes_mb <- readLines(file.path(go_dir, "GeneList_MobilizedBlood.txt"))
genes_bm <- readLines(file.path(go_dir, "GeneList_BoneMarrow.txt"))

# Read the Background Universe
bg_file <- file.path(go_dir, "GeneList_Background_Universe_CLEAN.txt")

if(!file.exists(bg_file)) {
  stop("ERROR: Could not find 'GeneList_Background_Universe_CLEAN.txt'. Please export this from your ArchR script first.")
}
background_genes <- readLines(bg_file)

cat("Loaded gene lists:\n")
cat("  MobilizedBlood Hits: ", length(genes_mb), "\n")
cat("  BoneMarrow Hits:     ", length(genes_bm), "\n")
cat("  HSC Background:      ", length(background_genes), "(Used for correction)\n\n")

#===========================================================================
# 2. GO ENRICHMENT (With Custom Background)
#===========================================================================

run_hsc_gost <- function(query_genes, bg_genes, label) {
  cat(paste0("Running GO enrichment for ", label, "...\n"))
  
  tryCatch({
    res <- gost(query = query_genes, 
                organism = "hsapiens",
                ordered_query = FALSE,
                multi_query = FALSE,
                significant = TRUE,
                exclude_iea = TRUE,
                user_threshold = 0.05,
                correction_method = "fdr",
                domain_scope = "custom",
                custom_bg = bg_genes,
                sources = c("GO:BP", "GO:MF", "KEGG", "REAC"),
                evcodes = TRUE)
    
    if(!is.null(res) && !is.null(res$result)) {
      cat("  Total significant terms:", nrow(res$result), "\n")
      return(res)
    } else {
      cat("  No significant terms found.\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("  Error running gost:", e$message, "\n")
    return(NULL)
  })
}

# Run the analysis
go_mb <- run_hsc_gost(genes_mb, background_genes, "MobilizedBlood")
go_bm <- run_hsc_gost(genes_bm, background_genes, "BoneMarrow")

#===========================================================================
# 3. PRINT TOP RESULTS TO LOG
#===========================================================================

print_top_hits <- function(gost_obj, label) {
  cat(paste0("\n=== TOP 10 TERMS - ", label, " ===\n"))
  if(!is.null(gost_obj) && !is.null(gost_obj$result)) {
    # Sort by p-value
    top <- head(gost_obj$result[order(gost_obj$result$p_value), 
                                c("source", "term_name", "p_value", "intersection_size")], 10)
    # Format p-value for readability
    top$p_value <- format(top$p_value, scientific = TRUE, digits = 2)
    # Truncate long names
    top$term_name <- substr(top$term_name, 1, 50)
    print(top, row.names = FALSE)
  } else {
    cat("No results to display.\n")
  }
}

print_top_hits(go_mb, "MOBILIZED BLOOD")
print_top_hits(go_bm, "BONE MARROW")

#===========================================================================
# 4. EXPORT RESULTS TO CSV
#===========================================================================

flatten_for_csv <- function(df) {
  # Convert list columns to character
  df_out <- df
  for(col in names(df_out)) {
    if(is.list(df_out[[col]])) {
      df_out[[col]] <- sapply(df_out[[col]], function(x) paste(x, collapse = ", "))
    }
  }
  return(df_out)
}

export_gost <- function(gost_obj, prefix) {
  if(!is.null(gost_obj) && !is.null(gost_obj$result)) {
    # Export ALL results
    write.csv(flatten_for_csv(gost_obj$result), 
              file.path(go_dir, paste0("GO_", prefix, "_ALL.csv")), row.names = FALSE)
  }
}

cat("\nExporting results...\n")
export_gost(go_mb, "MobilizedBlood")
export_gost(go_bm, "BoneMarrow")

cat("\n=== DONE ===\n")