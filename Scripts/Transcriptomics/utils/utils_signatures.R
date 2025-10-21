# =============================================================================
# utils_signatures.R — Gene signature and scoring utilities
# =============================================================================
#
# Utility functions for working with gene signatures, immune cell deconvolution,
# and pathway scoring in the DGAT immunology analysis pipeline.
#
# =============================================================================

#' Load gene signatures from GMT file
#' @param gmt_file Path to GMT file
#' @return Named list of gene sets
load_gmt_signatures <- function(gmt_file) {
  
  if (!file.exists(gmt_file)) {
    stop("GMT file does not exist: ", gmt_file)
  }
  
  cat("📖 Loading gene signatures from:", basename(gmt_file), "\n")
  
  signatures <- list()
  lines <- readLines(gmt_file)
  
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      sig_name <- parts[1]
      sig_desc <- parts[2]
      genes <- parts[3:length(parts)]
      genes <- genes[genes != ""]  # Remove empty entries
      signatures[[sig_name]] <- genes
    }
  }
  
  cat("✅ Loaded", length(signatures), "gene signatures\n")
  return(signatures)
}

#' Calculate GSVA scores for gene signatures
#' @param expr_matrix Expression matrix (genes x samples)
#' @param gene_sets Named list of gene sets
#' @param method GSVA method ("gsva", "ssgsea", "zscore", "plage")
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @return GSVA score matrix (signatures x samples)
calculate_gsva_scores <- function(expr_matrix, gene_sets, method = "ssgsea", 
                                  min_size = 10, max_size = 500) {
  
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package required. Install with: BiocManager::install('GSVA')")
  }
  
  cat("🧮 Calculating GSVA scores using", method, "method...\n")
  
  # Filter gene sets by size
  filtered_sets <- lapply(gene_sets, function(genes) {
    genes[genes %in% rownames(expr_matrix)]
  })
  
  # Remove gene sets that are too small or too large
  set_sizes <- sapply(filtered_sets, length)
  valid_sets <- set_sizes >= min_size & set_sizes <= max_size
  filtered_sets <- filtered_sets[valid_sets]
  
  cat("📊 Using", length(filtered_sets), "gene sets (filtered from", length(gene_sets), ")\n")
  cat("📊 Gene set sizes: min =", min(set_sizes[valid_sets]), ", max =", max(set_sizes[valid_sets]), "\n")
  
  # Calculate GSVA scores
  gsva_scores <- GSVA::gsva(expr_matrix, filtered_sets, 
                           method = method,
                           min.sz = min_size,
                           max.sz = max_size,
                           verbose = TRUE)
  
  cat("✅ Calculated GSVA scores for", nrow(gsva_scores), "signatures x", ncol(gsva_scores), "samples\n")
  return(gsva_scores)
}

#' Calculate AddModuleScore for Seurat objects
#' @param seurat_obj Seurat object
#' @param gene_sets Named list of gene sets
#' @param name_prefix Prefix for score column names
#' @param ctrl Number of control genes
#' @param seed Random seed for reproducibility
#' @return Seurat object with added module scores
calculate_module_scores <- function(seurat_obj, gene_sets, name_prefix = "Module", 
                                   ctrl = 100, seed = 42) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package required")
  }
  
  cat("🧮 Calculating module scores for", length(gene_sets), "gene sets...\n")
  
  set.seed(seed)
  
  for (i in seq_along(gene_sets)) {
    set_name <- names(gene_sets)[i]
    genes <- gene_sets[[i]]
    
    # Filter genes that exist in the Seurat object
    available_genes <- genes[genes %in% rownames(seurat_obj)]
    
    if (length(available_genes) < 5) {
      cat("⚠️  Skipping", set_name, "- too few genes available (", length(available_genes), ")\n")
      next
    }
    
    cat("  📊", set_name, ":", length(available_genes), "genes\n")
    
    # Calculate module score
    seurat_obj <- Seurat::AddModuleScore(
      object = seurat_obj,
      features = list(available_genes),
      name = paste0(name_prefix, "_", set_name),
      ctrl = ctrl,
      seed = seed
    )
  }
  
  cat("✅ Module scores calculated and added to Seurat object\n")
  return(seurat_obj)
}

#' Load immune cell signatures
#' @param source Source of signatures ("IOBR", "custom")
#' @param custom_file Path to custom signature file (if source = "custom")
#' @return Named list of immune cell signatures
load_immune_signatures <- function(source = "IOBR", custom_file = NULL) {
  
  if (source == "IOBR") {
    if (!requireNamespace("IOBR", quietly = TRUE)) {
      stop("IOBR package required. Install with: remotes::install_github('IOBR/IOBR')")
    }
    
    cat("📖 Loading immune signatures from IOBR...\n")
    
    # Load IOBR signatures
    signatures <- list(
      "T_cells" = IOBR::sig_tcell,
      "B_cells" = IOBR::sig_bcell,
      "NK_cells" = IOBR::sig_nk,
      "Macrophages" = IOBR::sig_macrophage,
      "Dendritic_cells" = IOBR::sig_dc,
      "Neutrophils" = IOBR::sig_neutrophil,
      "Mast_cells" = IOBR::sig_mast,
      "Eosinophils" = IOBR::sig_eosinophil,
      "Monocytes" = IOBR::sig_monocyte
    )
    
    # Remove NULL entries
    signatures <- signatures[!sapply(signatures, is.null)]
    
    cat("✅ Loaded", length(signatures), "immune signatures from IOBR\n")
    
  } else if (source == "custom") {
    if (is.null(custom_file)) {
      stop("custom_file must be specified when source = 'custom'")
    }
    signatures <- load_gmt_signatures(custom_file)
  } else {
    stop("Unknown signature source: ", source)
  }
  
  return(signatures)
}

#' Calculate DGAT pathway scores
#' @param expr_matrix Expression matrix (genes x samples)
#' @param dgat_genes Vector of DGAT-related genes
#' @return Named list with DGAT scores
calculate_dgat_scores <- function(expr_matrix, dgat_genes = NULL) {
  
  if (is.null(dgat_genes)) {
    # Default DGAT genes
    dgat_genes <- c("DGAT1", "DGAT2", "ACAT1", "ACAT2", "SCD", "FASN", "ACC1")
  }
  
  cat("🧮 Calculating DGAT pathway scores...\n")
  
  # Filter genes that exist in expression matrix
  available_genes <- dgat_genes[dgat_genes %in% rownames(expr_matrix)]
  cat("📊 Using", length(available_genes), "DGAT genes:", paste(available_genes, collapse = ", "), "\n")
  
  if (length(available_genes) == 0) {
    stop("No DGAT genes found in expression matrix")
  }
  
  # Calculate mean expression of DGAT genes
  dgat_scores <- colMeans(expr_matrix[available_genes, , drop = FALSE], na.rm = TRUE)
  
  # Calculate individual gene scores
  individual_scores <- expr_matrix[available_genes, , drop = FALSE]
  
  # Calculate pathway activity score (mean z-score)
  z_scores <- scale(t(individual_scores))
  pathway_scores <- rowMeans(z_scores, na.rm = TRUE)
  names(pathway_scores) <- colnames(expr_matrix)
  
  result <- list(
    "DGAT_mean" = dgat_scores,
    "DGAT_pathway" = pathway_scores,
    "DGAT_individual" = individual_scores,
    "genes_used" = available_genes
  )
  
  cat("✅ Calculated DGAT scores for", length(dgat_scores), "samples\n")
  return(result)
}

#' Create gene signature summary table
#' @param signatures Named list of gene signatures
#' @return Data.frame with signature information
summarize_signatures <- function(signatures) {
  
  summary_df <- data.frame(
    signature_name = names(signatures),
    gene_count = sapply(signatures, length),
    stringsAsFactors = FALSE
  )
  
  summary_df$gene_list <- sapply(signatures, function(genes) {
    paste(genes, collapse = ";")
  })
  
  cat("📊 Signature summary:\n")
  print(summary_df[, c("signature_name", "gene_count")])
  
  return(summary_df)
}

cat("✅ Signature utilities loaded\n")
