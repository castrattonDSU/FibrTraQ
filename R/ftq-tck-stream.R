#' MRtrix streamline reading
#'
#' TCK coordinates are read unchanged in RAS+ world millimetres. Image geometry
#' and any coordinate conversion belong to the analysis layer, not the reader.
#'
#' @keywords internal
#' @noRd

ftq_tck_stream <- function(path, chunk_size = 50000L, block_points = 8192L,
                           chunk_points = 1000000L,
                           max_points_per_track = 10000000L) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  header <- ftq_read_tck_header(path)
  structure(list(
    path = path,
    header = header,
    chunk_size = ftq_tck_positive_integer(chunk_size, "track_chunk_size"),
    block_points = ftq_tck_positive_integer(block_points, "block_points"),
    chunk_points = ftq_tck_positive_integer(chunk_points, "chunk_points"),
    max_points_per_track = ftq_tck_positive_integer(
      max_points_per_track, "max_points_per_track"
    )
  ), class = "ftq_tck_stream")
}

ftq_tck_positive_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value != floor(value) || value > .Machine$integer.max) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
}

ftq_tck_header_integer <- function(value, name) {
  if (length(value) != 1L || is.na(value) ||
      !grepl("^[0-9]+$", value)) {
    stop("Invalid TCK ", name, ": expected a nonnegative integer.", call. = FALSE)
  }
  number <- suppressWarnings(as.numeric(value))
  if (!is.finite(number) || number > 2^53 - 1) {
    stop("Invalid TCK ", name, ": integer is too large.", call. = FALSE)
  }
  number
}

ftq_read_tck_header <- function(path, max_header_bytes = 1048576L) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  size <- file.info(path)$size
  if (isTRUE(file.info(path)$isdir) || length(size) != 1L || !is.finite(size)) {
    stop("TCK input must be a readable regular file.", call. = FALSE)
  }
  max_header_bytes <- ftq_tck_positive_integer(max_header_bytes, "max_header_bytes")
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  buffered <- raw(0)
  consumed <- 0
  lines <- character(0)
  ended <- FALSE
  repeat {
    # Read a bounded header block; do not decode binary bytes after END.
    remaining <- max_header_bytes - consumed - length(buffered)
    if (remaining <= 0) {
      stop("TCK header exceeds the supported byte limit or has no END line.", call. = FALSE)
    }
    bytes <- readBin(con, "raw", n = min(4096L, remaining))
    if (!length(bytes)) {
      stop("Unexpected end of TCK header: missing newline-terminated END.", call. = FALSE)
    }
    buffered <- c(buffered, bytes)
    repeat {
      newline <- match(as.raw(10), buffered)
      if (is.na(newline)) break
      line_bytes <- if (newline > 1L) buffered[seq_len(newline - 1L)] else raw(0)
      if (length(line_bytes) && utils::tail(line_bytes, 1L) == as.raw(13)) {
        line_bytes <- utils::head(line_bytes, -1L)
      }
      codes <- as.integer(line_bytes)
      if (any(!(codes == 9L | (codes >= 32L & codes <= 126L)))) {
        stop("TCK header contains invalid non-ASCII or control bytes.", call. = FALSE)
      }
      line <- rawToChar(line_bytes)
      lines <- c(lines, line)
      consumed <- consumed + newline
      buffered <- if (newline < length(buffered)) {
        buffered[seq.int(newline + 1L, length(buffered))]
      } else raw(0)
      if (identical(line, "END")) {
        ended <- TRUE
        break
      }
    }
    if (ended) break
  }
  # MRtrix may pad the magic line with spaces. Accept only trailing horizontal
  # ASCII whitespace; keep the original line in raw_header for provenance.
  if (!length(lines) || !grepl("^mrtrix tracks[ \t]*$", lines[[1L]])) {
    stop("File is not a valid MRtrix TCK file: expected 'mrtrix tracks'.", call. = FALSE)
  }
  fields <- list()
  if (length(lines) > 2L) {
    for (line in lines[seq.int(2L, length(lines) - 1L)]) {
      colon <- regexpr(":", line, fixed = TRUE)[[1L]]
      if (colon < 2L) {
        stop("Invalid TCK header entry: expected key: value.", call. = FALSE)
      }
      key <- trimws(substr(line, 1L, colon - 1L))
      value <- trimws(substr(line, colon + 1L, nchar(line)))
      if (!nzchar(key)) stop("Invalid empty TCK header key.", call. = FALSE)
      if (key %in% c("file", "datatype", "count") && !is.null(fields[[key]])) {
        stop("Duplicate TCK header field: ", key, ".", call. = FALSE)
      }
      fields[[key]] <- c(fields[[key]], value)
    }
  }
  file_field <- fields[["file", exact = TRUE]]
  datatype <- fields[["datatype", exact = TRUE]]
  count_field <- fields[["count", exact = TRUE]]
  if (is.null(file_field) || is.null(datatype)) {
    stop("TCK header requires 'file' and 'datatype' fields.", call. = FALSE)
  }
  if (!grepl("^\\.[ \t]+[0-9]+$", file_field)) {
    stop("TCK supports only a single-file payload: 'file: . offset'.", call. = FALSE)
  }
  offset <- ftq_tck_header_integer(
    sub("^\\.[ \t]+", "", file_field), "data offset"
  )
  if (offset < consumed || offset > size) {
    stop("TCK data offset is outside the payload or overlaps the header.", call. = FALSE)
  }
  types <- c("Float32LE", "Float32BE", "Float64LE", "Float64BE")
  if (!datatype %in% types) {
    stop("Unsupported TCK datatype; expected Float32LE/BE or Float64LE/BE.", call. = FALSE)
  }
  dsize <- if (startsWith(datatype, "Float32")) 4L else 8L
  payload_size <- size - offset
  if (payload_size < 3L * dsize || payload_size %% (3L * dsize) != 0) {
    stop("Truncated TCK payload: expected complete coordinate triplets and an EOF marker.", call. = FALSE)
  }
  count <- if (is.null(count_field)) NA_real_ else {
    ftq_tck_header_integer(count_field, "count")
  }
  header <- fields
  # The original text retains repeated/custom fields even if their names
  # coincide with the derived fields used internally by this reader.
  header$id <- "mrtrix tracks"
  header$raw_header <- lines
  header$derived <- list(
    filename_part = ".", data_offset = offset,
    endian = if (endsWith(datatype, "LE")) "little" else "big",
    dsize = dsize, header_bytes = consumed, file_bytes = size,
    declared_count = count, coordinate_space = "RAS+ world millimetres"
  )
  header
}

ftq_tck_stream_iterate <- function(stream, callback, max_tracks = Inf,
                                   chunk_size = NULL) {
  if (!inherits(stream, "ftq_tck_stream")) {
    stop("stream must be an ftq_tck_stream object.", call. = FALSE)
  }
  if (!is.function(callback)) stop("callback must be a function.", call. = FALSE)
  if (!is.numeric(max_tracks) || length(max_tracks) != 1L || is.na(max_tracks) ||
      max_tracks < 0 || (is.finite(max_tracks) && max_tracks != floor(max_tracks))) {
    stop("max_tracks must be a nonnegative integer or Inf.", call. = FALSE)
  }
  chunk_size <- ftq_tck_positive_integer(
    if (is.null(chunk_size)) stream$chunk_size else chunk_size, "chunk_size"
  )
  # Revalidate the header each pass so stale offsets/counts are never trusted.
  header <- ftq_read_tck_header(stream$path)
  if (!identical(header, stream$header)) {
    stop("TCK header or file size changed after the reader was created.", call. = FALSE)
  }
  if (max_tracks == 0) return(invisible(0L))
  con <- file(stream$path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = header$derived$data_offset, origin = "start")
  chunk <- list()
  chunk_npoints <- 0
  chunk_start <- 1
  track_index <- 0
  fragments <- list()
  npoints <- 0
  finished <- FALSE
  stopped <- FALSE
  flush_chunk <- function() {
    if (!length(chunk)) return(FALSE)
    start <- if (chunk_start <= .Machine$integer.max) as.integer(chunk_start) else chunk_start
    end <- if (track_index <= .Machine$integer.max) as.integer(track_index) else track_index
    requested_stop <- identical(callback(chunk, start, end), FALSE)
    chunk <<- list()
    chunk_npoints <<- 0
    chunk_start <<- track_index + 1
    requested_stop
  }
  append_points <- function(points) {
    if (!nrow(points)) return(invisible(NULL))
    npoints <<- npoints + nrow(points)
    if (npoints > stream$max_points_per_track) {
      stop("TCK streamline ", track_index + 1,
           " exceeds max_points_per_track; input was not truncated.", call. = FALSE)
    }
    fragments[[length(fragments) + 1L]] <<- points
    invisible(NULL)
  }
  repeat {
    values <- readBin(con, numeric(), n = 3 * stream$block_points,
                      size = header$derived$dsize, endian = header$derived$endian)
    if (!length(values)) {
      stop("Unexpected end of TCK payload: missing Inf EOF marker.", call. = FALSE)
    }
    if (length(values) %% 3L != 0L) {
      stop("Truncated TCK coordinate triplet.", call. = FALSE)
    }
    points <- matrix(values, ncol = 3L, byrow = TRUE)
    markers <- which(rowSums(is.finite(points)) != 3L)
    next_point <- 1L
    for (marker in markers) {
      if (marker > next_point) {
        append_points(points[seq.int(next_point, marker - 1L), , drop = FALSE])
      }
      value <- points[marker, ]
      if (all(is.na(value))) {
        if (npoints == 0) {
          stop("Empty TCK streamline at record ", track_index + 1, ".", call. = FALSE)
        }
        track_index <- track_index + 1
        coords <- if (length(fragments) == 1L) fragments[[1L]] else do.call(rbind, fragments)
        rownames(coords) <- NULL
        chunk[[length(chunk) + 1L]] <- list(
          num_points = nrow(coords), coords = coords, scalars = NULL, properties = NULL
        )
        chunk_npoints <- chunk_npoints + npoints
        fragments <- list()
        npoints <- 0
        if (length(chunk) >= chunk_size || chunk_npoints >= stream$chunk_points ||
            track_index >= max_tracks) {
          stopped <- flush_chunk()
        }
        if (stopped || track_index >= max_tracks) {
          stopped <- TRUE
          break
        }
      } else if (all(is.infinite(value) & value == Inf)) {
        if (npoints != 0) {
          stop("TCK EOF reached before the final streamline's NaN delimiter.", call. = FALSE)
        }
        if (marker != nrow(points) || seek(con) != header$derived$file_bytes) {
          stop("Unexpected trailing data after the TCK EOF marker.", call. = FALSE)
        }
        finished <- TRUE
        break
      } else {
        stop("Invalid non-finite TCK coordinate triplet at streamline ",
             track_index + 1, ".", call. = FALSE)
      }
      next_point <- marker + 1L
    }
    if (stopped || finished) break
    if (next_point <= nrow(points)) {
      append_points(points[seq.int(next_point, nrow(points)), , drop = FALSE])
    }
  }
  if (finished) {
    declared <- header$derived$declared_count
    if (!is.na(declared) && track_index != declared) {
      stop("TCK count mismatch: header declares ", declared,
           " streamlines but payload contains ", track_index, ".", call. = FALSE)
    }
    flush_chunk()
  }
  invisible(if (track_index <= .Machine$integer.max) as.integer(track_index) else track_index)
}

ftq_tck_stream_sample <- function(stream, n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n)) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }
  if (n == 0) return(list())
  out <- list()
  ftq_tck_stream_iterate(stream, function(chunk, start, end) {
    out <<- c(out, chunk)
  }, max_tracks = n, chunk_size = min(n, stream$chunk_size))
  out
}
