#' Shared tractogram-reader dispatch
#'
#' @keywords internal
#' @noRd

ftq_tractogram_format <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("tractogram must be one nonempty file path.", call. = FALSE)
  }
  extension <- tolower(tools::file_ext(path))
  if (!extension %in% c("trk", "tck")) {
    stop("Supported tractogram files are uncompressed .trk and .tck files.", call. = FALSE)
  }
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  if (isTRUE(file.info(path)$isdir)) {
    stop("tractogram must be a file, not a directory.", call. = FALSE)
  }
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", n = 14L)
  expected <- if (extension == "trk") c(charToRaw("TRACK"), as.raw(0)) else {
    charToRaw("mrtrix tracks")
  }
  if (length(magic) < length(expected) ||
      !identical(magic[seq_along(expected)], expected)) {
    stop("Tractogram signature does not match its .", extension, " extension.", call. = FALSE)
  }
  extension
}

ftq_read_tractogram <- function(path, reader = "auto", chunk_size = 50000L) {
  reader <- match.arg(reader, c("auto", "stream", "memory"))
  format <- ftq_tractogram_format(path)
  if (reader == "auto") reader <- "stream"
  if (format == "trk") {
    if (reader == "stream") return(ftq_trk_stream(path, chunk_size = chunk_size))
    read_trk <- freesurferformats::read.dti.trk
    read_args <- list(filepath = path)
    # FibrTraQ applies the selected mapping itself. Newer dependency versions
    # expose a coordinate choice; request the unchanged stored coordinates.
    if ("coords" %in% names(formals(read_trk))) read_args$coords <- "native"
    trk <- do.call(read_trk, read_args)
    # fs.tracts is a compact S3 container: for() and C++ see its storage fields,
    # while as.list() dispatches to the dependency's per-streamline accessor.
    trk$tracks <- as.list(trk$tracks)
    return(trk)
  }
  stream <- ftq_tck_stream(path, chunk_size = chunk_size)
  if (reader == "stream") return(stream)
  tracks <- list()
  ftq_tck_stream_iterate(stream, function(chunk, start, end) {
    tracks <<- c(tracks, chunk)
  })
  structure(list(header = stream$header, tracks = tracks), class = "ftq_tck")
}

ftq_is_track_stream <- function(x) {
  inherits(x, "ftq_trk_stream") || inherits(x, "ftq_tck_stream")
}

ftq_track_format <- function(x) {
  if (inherits(x, "ftq_tck_stream") || inherits(x, "ftq_tck")) return("tck")
  if (inherits(x, "ftq_trk_stream")) return("trk")
  # A bare compact track collection has no tractogram header; its [[ method
  # addresses streamlines by index rather than named metadata fields.
  if (inherits(x, "fs.tracts")) return(NULL)
  if (!is.list(x) || !is.list(x[["header", exact = TRUE]])) return(NULL)
  header <- x[["header", exact = TRUE]]
  if (identical(header[["id", exact = TRUE]], "mrtrix tracks")) return("tck")
  if (identical(header[["id_string", exact = TRUE]], "TRACK")) return("trk")
  NULL
}

ftq_track_count <- function(x) {
  if (ftq_is_track_stream(x)) {
    if (identical(ftq_track_format(x), "tck")) {
      n <- x$header$derived$declared_count
      if (is.null(n) || length(n) != 1L || is.na(n)) return(NA_integer_)
    } else {
      n <- x$header$n_count
      if (is.null(n) || length(n) != 1L || is.na(n) || n <= 0) return(NA_integer_)
    }
  } else if (inherits(x, "fs.tracts")) {
    n <- length(x)
  } else if (is.list(x) && is.list(x[["tracks", exact = TRUE]])) {
    n <- length(x[["tracks", exact = TRUE]])
  } else if (is.list(x)) {
    n <- length(x)
  } else return(NA_integer_)
  if (n <= .Machine$integer.max) as.integer(n) else n
}

ftq_track_stream_iterate <- function(stream, callback, max_tracks = Inf,
                                     chunk_size = NULL) {
  if (inherits(stream, "ftq_tck_stream")) {
    return(ftq_tck_stream_iterate(stream, callback, max_tracks, chunk_size))
  }
  if (inherits(stream, "ftq_trk_stream")) {
    return(ftq_trk_stream_iterate(stream, callback, max_tracks, chunk_size))
  }
  stop("stream must be a TCK or TRK stream reader.", call. = FALSE)
}

ftq_track_stream_sample <- function(stream, n) {
  if (inherits(stream, "ftq_tck_stream")) return(ftq_tck_stream_sample(stream, n))
  if (inherits(stream, "ftq_trk_stream")) return(ftq_trk_stream_sample(stream, n))
  stop("stream must be a TCK or TRK stream reader.", call. = FALSE)
}
