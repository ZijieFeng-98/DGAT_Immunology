# =============================================================================
# utils_io.R — Input/Output helper functions
# =============================================================================
#
# Utility functions for reading and writing data files in the DGAT immunology
# analysis pipeline. Handles common file formats and provides consistent
# error handling and logging.
#
# =============================================================================

#' Read expression matrix with error handling
#' @param file_path Path to the expression matrix file
#' @param format File format ("csv", "tsv", "gct", "h5")
#' @param gene_col Column name for gene identifiers (default: "gene_id")
#' @param sample_col Column name for sample identifiers (default: "sample_id")
#' @return Expression matrix with genes as rows, samples as columns
read_expression_matrix <- function(file_path, format = "auto", gene_col = "gene_id", sample_col = "sample_id") {
  
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  # Auto-detect format from extension
  if (format == "auto") {
    ext <- tools::file_ext(file_path)
    format <- switch(ext,
      "csv" = "csv",
      "tsv" = "tsv", 
      "txt" = "tsv",
      "gct" = "gct",
      "h5" = "h5",
      "csv"  # default
    )
  }
  
  cat("📖 Reading", format, "file:", basename(file_path), "\n")
  
  switch(format,
    "csv" = {
      data <- data.table::fread(file_path, data.table = FALSE)
      rownames(data) <- data[[gene_col]]
      data[[gene_col]] <- NULL
      return(as.matrix(data))
    },
    "tsv" = {
      data <- data.table::fread(file_path, data.table = FALSE, sep = "\t")
      rownames(data) <- data[[gene_col]]
      data[[gene_col]] <- NULL
      return(as.matrix(data))
    },
    "gct" = {
      # GCT format: first line is version, second is dimensions, data starts at line 3
      lines <- readLines(file_path)
      data_start <- 3
      data <- data.table::fread(file_path, skip = 2, data.table = FALSE)
      rownames(data) <- data[[1]]  # First column is gene ID
      data[[1]] <- NULL
      data[[1]] <- NULL  # Second column is gene description
      return(as.matrix(data))
    },
    "h5" = {
      if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("hdf5r package required for reading H5 files. Install with: install.packages('hdf5r')")
      }
      # H5 file reading would go here
      stop("H5 format reading not yet implemented")
    },
    stop("Unsupported format: ", format)
  )
}

#' Write expression matrix with metadata
#' @param expr_matrix Expression matrix (genes x samples)
#' @param file_path Output file path
#' @param format Output format ("csv", "tsv")
#' @param gene_metadata Optional data.frame with gene annotations
write_expression_matrix <- function(expr_matrix, file_path, format = "auto", gene_metadata = NULL) {
  
  # Auto-detect format from extension
  if (format == "auto") {
    ext <- tools::file_ext(file_path)
    format <- ifelse(ext %in% c("csv"), "csv", "tsv")
  }
  
  cat("💾 Writing", format, "file:", basename(file_path), "\n")
  
  # Create output directory if it doesn't exist
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  
  # Prepare data for writing
  if (!is.null(gene_metadata)) {
    # Combine expression data with gene metadata
    output_data <- cbind(gene_metadata, expr_matrix)
  } else {
    # Add gene IDs as first column
    output_data <- data.frame(
      gene_id = rownames(expr_matrix),
      expr_matrix,
      stringsAsFactors = FALSE
    )
  }
  
  # Write file
  if (format == "csv") {
    data.table::fwrite(output_data, file_path, sep = ",")
  } else {
    data.table::fwrite(output_data, file_path, sep = "\t")
  }
  
  cat("✅ Written", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples to", basename(file_path), "\n")
}

#' Read clinical data
#' @param file_path Path to clinical data file
#' @param format File format ("csv", "tsv")
#' @param sample_col Column name containing sample IDs
#' @return Clinical data data.frame
read_clinical_data <- function(file_path, format = "auto", sample_col = "sample_id") {
  
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  # Auto-detect format
  if (format == "auto") {
    ext <- tools::file_ext(file_path)
    format <- ifelse(ext == "csv", "csv", "tsv")
  }
  
  cat("📖 Reading clinical data:", basename(file_path), "\n")
  
  if (format == "csv") {
    clinical <- data.table::fread(file_path, data.table = FALSE)
  } else {
    clinical <- data.table::fread(file_path, data.table = FALSE, sep = "\t")
  }
  
  # Set sample IDs as rownames if specified column exists
  if (sample_col %in% colnames(clinical)) {
    rownames(clinical) <- clinical[[sample_col]]
  }
  
  cat("✅ Read clinical data for", nrow(clinical), "samples\n")
  return(clinical)
}

#' Write results table with timestamp
#' @param results Data.frame or matrix to write
#' @param file_path Output file path (without extension)
#' @param format Output format ("csv", "tsv")
#' @param add_timestamp Whether to add timestamp to filename
write_results <- function(results, file_path, format = "csv", add_timestamp = TRUE) {
  
  # Add timestamp if requested
  if (add_timestamp) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    base_name <- tools::file_path_sans_ext(file_path)
    ext <- tools::file_ext(file_path)
    if (ext == "") ext <- format
    file_path <- paste0(base_name, "_", timestamp, ".", ext)
  }
  
  # Create output directory
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  
  # Write file
  if (format == "csv") {
    data.table::fwrite(results, file_path, sep = ",")
  } else {
    data.table::fwrite(results, file_path, sep = "\t")
  }
  
  cat("💾 Results written to:", basename(file_path), "\n")
}

#' Check if required columns exist in data
#' @param data Data.frame to check
#' @param required_cols Vector of required column names
#' @param data_name Name of dataset for error messages
check_required_columns <- function(data, required_cols, data_name = "data") {
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", data_name, ": ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

#' Create output directory if it doesn't exist
#' @param dir_path Directory path to create
ensure_output_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    cat("📁 Created directory:", dir_path, "\n")
  }
  invisible(dir_path)
}

cat("✅ I/O utilities loaded\n")
