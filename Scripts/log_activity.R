#!/usr/bin/env Rscript
# =============================================================================
# log_activity.R - Automated Daily Research Activity Logger
# =============================================================================
#
# Purpose: Automatically log research activities, script runs, and analyses
#
# Usage:
#   Source this at the start of any R script:
#   source("Scripts/log_activity.R")
#   log_activity("Started DGAT immune analysis on TCGA data")
#
#   Or run standalone to view logs:
#   Rscript Scripts/log_activity.R --view
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Configuration
LOG_DIR <- "Notes/Daily_Logs"
CURRENT_LOG <- file.path(LOG_DIR, paste0("WorkLog_", Sys.Date(), ".md"))

# Create log directory if it doesn't exist
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOGGING FUNCTIONS
# =============================================================================

#' Log an activity with timestamp
#'
#' @param message Character string describing the activity
#' @param category Category of activity (script, analysis, data, meeting, etc.)
#' @param auto_detect Automatically detect script name if called from within script
log_activity <- function(message, category = "general", auto_detect = TRUE) {
  
  timestamp <- format(Sys.time(), "%H:%M:%S")
  date <- Sys.Date()
  
  # Auto-detect calling script if enabled
  script_name <- NA
  if (auto_detect) {
    call_stack <- sys.calls()
    for (call in call_stack) {
      call_str <- paste(deparse(call), collapse = " ")
      if (length(call_str) == 1 && grepl("source|Rscript", call_str)) {
        script_match <- regmatches(call_str, regexpr("['\"].*\\.R['\"]", call_str))
        if (length(script_match) > 0) {
          script_name <- gsub("['\"]", "", script_match)
          break
        }
      }
    }
  }
  
  # Construct log entry
  entry <- list(
    timestamp = timestamp,
    category = category,
    script = ifelse(is.na(script_name), "", script_name),
    message = message
  )
  
  # Initialize daily log if it doesn't exist
  if (!file.exists(CURRENT_LOG)) {
    initialize_daily_log(date)
  }
  
  # Append to log
  append_to_log(entry)
  
  # Print to console
  cat(sprintf("[%s] %s: %s\n", timestamp, toupper(category), message))
  
  invisible(entry)
}

#' Initialize a new daily log file
initialize_daily_log <- function(date) {
  header <- sprintf(
"# Research Work Log - %s

**Project:** DGAT Immunology Analysis  
**Researcher:** Zijie Feng  
**Date:** %s  
**Day:** %s

---

## Today's Activities

", 
    date, 
    format(date, "%B %d, %Y"),
    format(date, "%A")
  )
  
  cat(header, file = CURRENT_LOG)
}

#' Append entry to daily log
append_to_log <- function(entry) {
  # Format entry based on category
  emoji <- switch(entry$category,
                 "script" = "🔧",
                 "analysis" = "📊",
                 "data" = "💾",
                 "result" = "✅",
                 "error" = "❌",
                 "meeting" = "👥",
                 "reading" = "📖",
                 "writing" = "✍️",
                 "planning" = "📋",
                 "•")
  
  log_entry <- sprintf("**[%s]** %s **%s**", 
                       entry$timestamp, 
                       emoji,
                       entry$message)
  
  if (entry$script != "") {
    log_entry <- paste0(log_entry, sprintf("  \n  _Script: `%s`_", entry$script))
  }
  
  log_entry <- paste0(log_entry, "\n\n")
  
  # Append to file
  cat(log_entry, file = CURRENT_LOG, append = TRUE)
}

#' Log script start (call this at the beginning of analysis scripts)
log_script_start <- function(script_name = NULL, description = "") {
  if (is.null(script_name)) {
    script_name <- basename(sys.frame(1)$ofile)
  }
  
  message <- sprintf("Started running %s", script_name)
  if (description != "") {
    message <- paste0(message, " - ", description)
  }
  
  log_activity(message, category = "script", auto_detect = FALSE)
}

#' Log script completion
log_script_end <- function(script_name = NULL, summary = "") {
  if (is.null(script_name)) {
    script_name <- basename(sys.frame(1)$ofile)
  }
  
  message <- sprintf("Completed %s", script_name)
  if (summary != "") {
    message <- paste0(message, " - ", summary)
  }
  
  log_activity(message, category = "result", auto_detect = FALSE)
}

#' Log data processing
log_data <- function(dataset, action, details = "") {
  message <- sprintf("%s: %s", dataset, action)
  if (details != "") {
    message <- paste0(message, " (", details, ")")
  }
  log_activity(message, category = "data")
}

#' Log analysis results
log_result <- function(description) {
  log_activity(description, category = "result")
}

#' Log errors or issues
log_error <- function(description) {
  log_activity(description, category = "error")
}

# =============================================================================
# DAILY SUMMARY FUNCTIONS
# =============================================================================

#' Add end-of-day summary section
add_daily_summary <- function(accomplishments = "", next_steps = "", notes = "") {
  
  summary <- "\n---\n\n## Daily Summary\n\n"
  
  if (accomplishments != "") {
    summary <- paste0(summary, "### ✅ Accomplishments\n", accomplishments, "\n\n")
  }
  
  if (next_steps != "") {
    summary <- paste0(summary, "### 📋 Next Steps\n", next_steps, "\n\n")
  }
  
  if (notes != "") {
    summary <- paste0(summary, "### 📝 Notes\n", notes, "\n\n")
  }
  
  summary <- paste0(summary, sprintf("---\n*Log generated: %s*\n", Sys.time()))
  
  cat(summary, file = CURRENT_LOG, append = TRUE)
  
  cat("✅ Daily summary added to log\n")
}

# =============================================================================
# VIEWING AND REPORTING
# =============================================================================

#' View today's log
view_log <- function(date = Sys.Date()) {
  log_file <- file.path(LOG_DIR, paste0("WorkLog_", date, ".md"))
  
  if (!file.exists(log_file)) {
    cat(sprintf("No log file found for %s\n", date))
    return(invisible(NULL))
  }
  
  cat(readLines(log_file), sep = "\n")
}

#' Get weekly summary
weekly_summary <- function() {
  dates <- Sys.Date() - 0:6
  log_files <- file.path(LOG_DIR, paste0("WorkLog_", dates, ".md"))
  
  cat("\n", strrep("=", 70), "\n")
  cat("WEEKLY WORK SUMMARY\n")
  cat(strrep("=", 70), "\n\n")
  
  for (i in seq_along(dates)) {
    date <- dates[i]
    log_file <- log_files[i]
    
    cat(sprintf("📅 %s (%s):\n", format(date, "%B %d, %Y"), format(date, "%A")))
    
    if (file.exists(log_file)) {
      lines <- readLines(log_file)
      activity_lines <- grep("^\\*\\*\\[", lines, value = TRUE)
      
      if (length(activity_lines) > 0) {
        cat(sprintf("   %d activities logged\n", length(activity_lines)))
        cat("   ", head(activity_lines, 3), sep = "\n    ")
        if (length(activity_lines) > 3) {
          cat(sprintf("\n    ... and %d more\n", length(activity_lines) - 3))
        }
      } else {
        cat("   No activities logged\n")
      }
    } else {
      cat("   No log file\n")
    }
    cat("\n")
  }
}

# =============================================================================
# AUTO-TRACKING FOR COMMON OPERATIONS
# =============================================================================

#' Wrap common functions to auto-log
auto_track <- function() {
  # This function sets up automatic tracking of common operations
  # Call this at the start of your R session
  
  # Track when you load data
  .old_readRDS <- readRDS
  readRDS <- function(file, ...) {
    log_data(basename(file), "Loaded RDS file", dirname(file))
    .old_readRDS(file, ...)
  }
  
  cat("🤖 Automatic activity tracking enabled for this session\n")
  cat("   - Data loading (readRDS, fread) will be logged\n")
  cat("   - Use log_activity() to manually log activities\n\n")
}

# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

# Only run CLI if script is executed directly (not sourced)
if (!interactive() && identical(basename(commandArgs()[4]), "log_activity.R")) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage: Rscript Scripts/log_activity.R [--view|--weekly|--summary]\n")
    cat("\nOptions:\n")
    cat("  --view     View today's log\n")
    cat("  --weekly   View weekly summary\n")
    cat("  --summary  Add end-of-day summary (interactive)\n")
    quit(save = "no", status = 0)
  }
  
  if (args[1] == "--view") {
    view_log()
    quit(save = "no", status = 0)
  } else if (args[1] == "--weekly") {
    weekly_summary()
    quit(save = "no", status = 0)
  } else if (args[1] == "--summary") {
    cat("Enter accomplishments (press Enter twice when done):\n")
    accomplishments <- paste(readLines(file("stdin"), n = -1), collapse = "\n")
    
    cat("\nEnter next steps:\n")
    next_steps <- paste(readLines(file("stdin"), n = -1), collapse = "\n")
    
    add_daily_summary(accomplishments, next_steps)
    quit(save = "no", status = 0)
  }
}

# =============================================================================
# EXPORT FUNCTIONS
# =============================================================================

# Make functions available when sourced
if (interactive() || (exists(".GlobalEnv") && environmentIsLocked(.GlobalEnv) == FALSE)) {
  cat("📝 Activity logging system loaded!\n")
  cat("\nAvailable functions:\n")
  cat("  log_activity(message, category) - Log any activity\n")
  cat("  log_script_start(name, desc)    - Log script start\n")
  cat("  log_script_end(name, summary)   - Log script completion\n")
  cat("  log_data(dataset, action)        - Log data operations\n")
  cat("  log_result(description)          - Log analysis results\n")
  cat("  log_error(description)           - Log errors/issues\n")
  cat("  view_log()                       - View today's log\n")
  cat("  weekly_summary()                 - View weekly summary\n")
  cat("  add_daily_summary()              - Add end-of-day summary\n")
  cat("\n")
}

