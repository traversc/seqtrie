# Ensure R6 method subsections have unique titles to avoid duplicate HTML anchors.
# Run after roxygen has generated the Rd files.
# This is a bug with R-devel that will likely be fixed in the future

files <- c("man/RadixForest.Rd", "man/RadixTree.Rd")

sanitize_method <- function(method_call) {
  name <- sub("\\(.*$", "", method_call)
  name <- sub(".*\\$", "", name)
  trimws(name)
}

process_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  current_method <- NA_character_
  changed <- FALSE

  method_pattern <- "^\\\\subsection\\{Method \\\\code\\{([^}]*)\\}\\}\\{$"

  for (i in seq_along(lines)) {
    line <- lines[[i]]

    if (grepl(method_pattern, line, perl = TRUE)) {
      method_call <- sub(method_pattern, "\\1", line, perl = TRUE)
      current_method <- sanitize_method(method_call)
      next
    }

    if (is.na(current_method) || !nzchar(current_method)) {
      next
    }

    if (grepl("^\\\\subsection\\{Usage", line)) {
      if (!grepl(" - ", line, fixed = TRUE)) {
        lines[[i]] <- sprintf("\\subsection{Usage - %s}{", current_method)
        changed <- TRUE
      }
      next
    }

    if (grepl("^\\\\subsection\\{Arguments", line)) {
      if (!grepl(" - ", line, fixed = TRUE)) {
        lines[[i]] <- sprintf("\\subsection{Arguments - %s}{", current_method)
        changed <- TRUE
      }
      next
    }

    if (grepl("^\\\\subsection\\{Returns", line)) {
      if (!grepl(" - ", line, fixed = TRUE)) {
        lines[[i]] <- sprintf("\\subsection{Returns - %s}{", current_method)
        changed <- TRUE
      }
    }
  }

  if (changed) {
    writeLines(lines, path)
    message("Updated ", path)
  } else {
    message("No changes made to ", path)
  }
}

invisible(lapply(files, process_file))
