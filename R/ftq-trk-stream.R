#' TrackVis streaming utilities
#'
#' Minimal exact reader for TrackVis `.trk` files. It reads the header once and
#' then yields streamlines in bounded chunks so large tractograms do not need to
#' be materialized as one large R list.
#'
#' @keywords internal
#' @noRd

ftq_trk_stream <- function(path, chunk_size = 50000L) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  header <- ftq_read_trk_header(path)
  chunk_size <- as.integer(chunk_size)[1]
  if (is.na(chunk_size) || chunk_size <= 0L) {
    stop("track_chunk_size must be a positive integer.")
  }
  structure(
    list(path = path, header = header, endian = header$endian, chunk_size = chunk_size),
    class = "ftq_trk_stream"
  )
}

ftq_read_trk_header <- function(path) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  endian <- ftq_trk_endian(path)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)

  header <- list()
  header$id_string <- ftq_read_fixed_char(con, 6L)
  header$dim <- readBin(con, integer(), n = 3L, size = 2L, endian = endian)
  header$voxel_size <- readBin(con, numeric(), n = 3L, size = 4L, endian = endian)
  header$origin <- readBin(con, numeric(), n = 3L, size = 4L, endian = endian)
  header$n_scalars <- readBin(con, integer(), n = 1L, size = 2L, endian = endian)
  header$scalar_names <- ftq_read_fixed_char(con, 200L)
  header$n_properties <- readBin(con, integer(), n = 1L, size = 2L, endian = endian)
  header$property_names <- ftq_read_fixed_char(con, 200L)
  header$vox2ras <- matrix(
    readBin(con, numeric(), n = 16L, size = 4L, endian = endian),
    ncol = 4L,
    byrow = TRUE
  )
  header$reserved <- ftq_read_fixed_char(con, 444L)
  header$voxel_order <- ftq_read_fixed_char(con, 4L)
  header$pad2 <- ftq_read_fixed_char(con, 4L)
  header$image_orientation_patient <- readBin(con, numeric(), n = 6L, size = 4L, endian = endian)
  header$pad1 <- ftq_read_fixed_char(con, 2L)
  header$invert_x <- readBin(con, integer(), n = 1L, size = 1L, signed = FALSE, endian = endian)
  header$invert_y <- readBin(con, integer(), n = 1L, size = 1L, signed = FALSE, endian = endian)
  header$invert_z <- readBin(con, integer(), n = 1L, size = 1L, signed = FALSE, endian = endian)
  header$swap_xy <- readBin(con, integer(), n = 1L, size = 1L, signed = FALSE, endian = endian)
  header$swap_yz <- readBin(con, integer(), n = 1L, size = 1L, signed = FALSE, endian = endian)
  header$swap_zx <- readBin(con, integer(), n = 1L, size = 1L, signed = FALSE, endian = endian)
  header$n_count <- readBin(con, integer(), n = 1L, size = 4L, endian = endian)
  header$version <- readBin(con, integer(), n = 1L, size = 4L, endian = endian)
  header$hdr_size <- readBin(con, integer(), n = 1L, size = 4L, endian = endian)
  header$endian <- endian

  if (!identical(header$hdr_size, 1000L)) {
    stop("Invalid TrackVis header size in ", path, ": ", header$hdr_size)
  }
  if (!identical(header$version, 2L)) {
    warning("TRK file uses TrackVis version ", header$version, "; version 2 is expected.")
  }
  header
}

ftq_trk_stream_iterate <- function(stream, callback, max_tracks = Inf, chunk_size = NULL) {
  if (!inherits(stream, "ftq_trk_stream")) {
    stop("stream must be an ftq_trk_stream object.")
  }
  if (!is.function(callback)) stop("callback must be a function.")
  chunk_size <- as.integer(if (is.null(chunk_size)) stream$chunk_size else chunk_size)[1]
  if (is.na(chunk_size) || chunk_size <= 0L) stop("chunk_size must be a positive integer.")

  header <- stream$header
  endian <- stream$endian
  n_scalars <- as.integer(header$n_scalars)
  n_properties <- as.integer(header$n_properties)
  declared_n <- as.integer(header$n_count)
  max_tracks <- if (is.infinite(max_tracks)) .Machine$integer.max else as.integer(max_tracks)

  con <- file(stream$path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = 1000L, origin = "start")

  chunk <- vector("list", chunk_size)
  chunk_n <- 0L
  chunk_start <- 1L
  track_idx <- 0L
  stop_requested <- FALSE

  flush_chunk <- function() {
    if (!chunk_n) return(invisible(NULL))
    res <- callback(chunk[seq_len(chunk_n)], chunk_start, chunk_start + chunk_n - 1L)
    if (identical(res, FALSE)) stop_requested <<- TRUE
    invisible(NULL)
  }

  repeat {
    if (stop_requested) break
    if (track_idx >= max_tracks) break
    if (declared_n > 0L && track_idx >= declared_n) break

    n_points <- readBin(con, integer(), n = 1L, size = 4L, endian = endian)
    if (!length(n_points)) break
    if (is.na(n_points) || n_points < 0L) stop("Invalid track point count at streamline ", track_idx + 1L)

    values_per_point <- 3L + n_scalars
    point_values <- if (n_points > 0L) {
      readBin(con, numeric(), n = n_points * values_per_point, size = 4L, endian = endian)
    } else {
      numeric(0)
    }
    expected <- n_points * values_per_point
    if (length(point_values) != expected) {
      stop("Unexpected end of TRK file while reading streamline ", track_idx + 1L)
    }

    coords <- if (n_points > 0L) {
      point_matrix <- matrix(point_values, ncol = values_per_point, byrow = TRUE)
      point_matrix[, 1:3, drop = FALSE]
    } else {
      matrix(numeric(0), ncol = 3L)
    }
    scalars <- NULL
    if (n_scalars > 0L && n_points > 0L) {
      scalars <- point_matrix[, 4:values_per_point, drop = FALSE]
    }
    properties <- if (n_properties > 0L) {
      readBin(con, numeric(), n = n_properties, size = 4L, endian = endian)
    } else {
      NULL
    }
    if (n_properties > 0L && length(properties) != n_properties) {
      stop("Unexpected end of TRK file while reading properties for streamline ", track_idx + 1L)
    }

    track_idx <- track_idx + 1L
    chunk_n <- chunk_n + 1L
    chunk[[chunk_n]] <- list(
      num_points = n_points,
      coords = coords,
      scalars = scalars,
      properties = properties
    )

    if (chunk_n == chunk_size) {
      flush_chunk()
      if (stop_requested) break
      chunk <- vector("list", chunk_size)
      chunk_n <- 0L
      chunk_start <- track_idx + 1L
    }
  }

  if (!stop_requested) flush_chunk()
  invisible(track_idx)
}

ftq_trk_stream_sample <- function(stream, n) {
  n <- as.integer(n)[1]
  if (is.na(n) || n <= 0L) return(list())
  out <- list()
  ftq_trk_stream_iterate(
    stream,
    callback = function(chunk, start, end) {
      out <<- c(out, chunk)
    },
    max_tracks = n,
    chunk_size = min(n, stream$chunk_size)
  )
  out
}

ftq_trk_endian <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = 996L, origin = "start")
  little <- readBin(con, integer(), n = 1L, size = 4L, endian = "little")
  if (identical(little, 1000L)) return("little")
  seek(con, where = 996L, origin = "start")
  big <- readBin(con, integer(), n = 1L, size = 4L, endian = "big")
  if (identical(big, 1000L)) return("big")
  stop("File is not a valid TrackVis .trk file; header size was ",
       little, " little-endian and ", big, " big-endian.")
}

ftq_read_fixed_char <- function(con, n) {
  bytes <- readBin(con, "raw", n = n)
  if (length(bytes) != n) stop("Unexpected end of file while reading TrackVis header.")
  nul <- which(bytes == as.raw(0))
  if (length(nul)) bytes <- bytes[seq_len(nul[1L] - 1L)]
  if (!length(bytes)) return("")
  rawToChar(bytes, multiple = FALSE)
}
