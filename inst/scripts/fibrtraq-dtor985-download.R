#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Usage: fibrtraq-dtor985-download.R <output-dir> [--article-id <id>] [--file-name <name>] [--download]\n",
    "\n",
    "Resolves the public Figshare metadata for the dTOR-985 full tractogram and\n",
    "optionally downloads the selected .trk file into <output-dir>.\n",
    "\n",
    "Examples:\n",
    "  Rscript fibrtraq-dtor985-download.R /path/to/dtor985\n",
    "  Rscript fibrtraq-dtor985-download.R /path/to/dtor985 --download\n",
    sep = ""
  )
  quit(status = 1L)
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[1] %in% c("-h", "--help")) usage()

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[1])) y else x

out_dir <- path.expand(args[1])
article_id <- "25209947"
file_name <- NULL
download <- FALSE

i <- 2L
while (i <= length(args)) {
  arg <- args[[i]]
  if (arg == "--article-id") {
    i <- i + 1L
    if (i > length(args)) usage()
    article_id <- args[[i]]
  } else if (arg == "--file-name") {
    i <- i + 1L
    if (i > length(args)) usage()
    file_name <- args[[i]]
  } else if (arg == "--download") {
    download <- TRUE
  } else {
    stop("Unknown argument: ", arg)
  }
  i <- i + 1L
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
api_url <- paste0("https://api.figshare.com/v2/articles/", article_id)
message("Fetching Figshare metadata: ", api_url)
meta <- jsonlite::fromJSON(api_url, simplifyVector = FALSE)
if (is.null(meta$files) || !length(meta$files)) {
  stop("No files were listed for Figshare article ", article_id)
}

files <- meta$files
names_available <- vapply(files, function(x) as.character(x$name), character(1))
if (is.null(file_name)) {
  trk_files <- files[grepl("\\.trk$", names_available, ignore.case = TRUE)]
  if (!length(trk_files)) stop("No .trk files were listed for article ", article_id)
  sizes <- vapply(trk_files, function(x) as.numeric(x$size), numeric(1))
  selected <- trk_files[[which.max(sizes)]]
} else {
  idx <- match(file_name, names_available)
  if (is.na(idx)) {
    stop("Requested file was not found. Available files:\n", paste(names_available, collapse = "\n"))
  }
  selected <- files[[idx]]
}

metadata_path <- file.path(out_dir, "dtor985-figshare-metadata.json")
selected_path <- file.path(out_dir, "dtor985-selected-file.json")
jsonlite::write_json(meta, metadata_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
jsonlite::write_json(selected, selected_path, auto_unbox = TRUE, pretty = TRUE, null = "null")

cat("Selected dTOR-985 file\n")
cat("  article:  ", meta$title, "\n", sep = "")
cat("  file:     ", selected$name, "\n", sep = "")
cat("  size:     ", format(as.numeric(selected$size), big.mark = ","), " bytes\n", sep = "")
cat("  md5:      ", selected$computed_md5 %||% selected$supplied_md5 %||% "", "\n", sep = "")
cat("  url:      ", selected$download_url, "\n", sep = "")
cat("  metadata: ", metadata_path, "\n", sep = "")
cat("  selected: ", selected_path, "\n", sep = "")

if (isTRUE(download)) {
  raw_dir <- file.path(out_dir, "raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(raw_dir, selected$name)
  if (file.exists(dest) && file.info(dest)$size == as.numeric(selected$size)) {
    cat("Download already complete: ", dest, "\n", sep = "")
  } else {
    tmp <- paste0(dest, ".part")
    if (file.exists(dest) && file.info(dest)$size != as.numeric(selected$size) && !file.exists(tmp)) {
      file.rename(dest, tmp)
    }
    cat("Downloading to: ", dest, "\n", sep = "")
    curl <- Sys.which("curl")
    if (nzchar(curl)) {
      status <- system2(
        curl,
        c(
          "-L",
          "--fail",
          "--retry", "5",
          "--continue-at", "-",
          "--output", tmp,
          selected$download_url
        )
      )
      if (!identical(status, 0L)) stop("curl download failed with status ", status)
    } else {
      utils::download.file(selected$download_url, tmp, mode = "wb", quiet = FALSE)
    }
    file.rename(tmp, dest)
    if (file.info(dest)$size != as.numeric(selected$size)) {
      warning("Downloaded file size differs from Figshare metadata.")
    }
    cat("Download complete: ", dest, "\n", sep = "")
  }
}
