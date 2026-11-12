ftq_reader_test_write_tck <- function(path, tracks = list(), datatype = "Float32LE",
                                      count = length(tracks), extra = character(),
                                      payload = NULL, eof = TRUE, offset = 512L,
                                      newline = "\n", magic = "mrtrix tracks") {
  lines <- c(magic, paste0("datatype: ", datatype),
             if (!is.null(count)) paste0("count: ", format(count, scientific = FALSE, trim = TRUE)), extra,
             paste0("file: . ", offset), "END")
  header <- charToRaw(paste0(paste(lines, collapse = newline), newline))
  stopifnot(length(header) <= offset)
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(c(header, raw(offset - length(header))), con)
  if (is.null(payload)) {
    payload <- unlist(lapply(tracks, function(x) c(as.numeric(t(x)), rep(NaN, 3L))),
                      use.names = FALSE)
    if (eof) payload <- c(payload, rep(Inf, 3L))
  }
  size <- if (startsWith(datatype, "Float64")) 8L else 4L
  endian <- if (endsWith(datatype, "BE")) "big" else "little"
  writeBin(as.numeric(payload), con, size = size, endian = endian)
  invisible(path)
}

ftq_reader_test_tracks <- function() {
  list(
    rbind(c(-91.5, 44.25, 0.125), c(-90.5, 43.25, 1.125)),
    matrix(c(2, -3, 5), nrow = 1L),
    cbind(seq_len(31L) / 8, -seq_len(31L), seq_len(31L) * 1.5),
    rbind(c(1, 2, 3), c(1, 2, 3), c(4, 5, 6))
  )
}

test_that("all standard TCK numeric encodings preserve coordinate/order data", {
  tracks <- ftq_reader_test_tracks()
  for (datatype in c("Float32LE", "Float32BE", "Float64LE", "Float64BE")) {
    path <- tempfile(fileext = ".tck")
    ftq_reader_test_write_tck(path, tracks, datatype,
                              extra = c("total_count: 47", "command_history: one: two",
                                        "command_history: three"))
    before <- tools::md5sum(path)
    stream <- ftq_read_tractogram(path)
    expect_s3_class(stream, "ftq_tck_stream")
    expect_identical(ftq_track_format(stream), "tck")
    expect_identical(ftq_track_count(stream), 4L)
    expect_identical(stream$header$total_count, "47")
    expect_identical(stream$header$command_history, c("one: two", "three"))
    expect_identical(stream$header$derived$coordinate_space, "RAS+ world millimetres")
    memory <- ftq_read_tractogram(path, reader = "memory")
    expect_s3_class(memory, "ftq_tck")
    expect_equal(lapply(memory$tracks, `[[`, "coords"), tracks)
    expect_identical(ftq_track_count(memory), 4L)
    expect_identical(ftq_track_format(memory), "tck")
    expect_identical(tools::md5sum(path), before)
  }
})

test_that("TCK readers accept MRtrix padded magic lines without altering provenance", {
  tracks <- ftq_reader_test_tracks()
  for (padding in c("    ", "\t", " \t  ")) {
    for (newline in c("\n", "\r\n")) {
      path <- tempfile(fileext = ".tck")
      magic <- paste0("mrtrix tracks", padding)
      ftq_reader_test_write_tck(path, tracks, magic = magic, newline = newline)
      before <- tools::md5sum(path)
      stream <- ftq_read_tractogram(path, reader = "stream", chunk_size = 2L)
      expect_identical(stream$header$id, "mrtrix tracks")
      expect_identical(stream$header$raw_header[[1L]], magic)
      expect_identical(ftq_track_count(stream), 4L)
      captured <- list()
      count <- ftq_track_stream_iterate(stream, function(chunk, start, end) {
        captured <<- c(captured, chunk)
      })
      expect_identical(count, 4L)
      expect_equal(lapply(captured, `[[`, "coords"), tracks)
      memory <- ftq_read_tractogram(path, reader = "memory")
      expect_identical(memory$header$raw_header[[1L]], magic)
      expect_identical(ftq_track_count(memory), 4L)
      expect_equal(lapply(memory$tracks, `[[`, "coords"), tracks)
      expect_identical(tools::md5sum(path), before)
    }
  }
})

test_that("TCK magic lines reject extra text and malformed prefixes", {
  for (magic in c("mrtrix tracksX", "mrtrix tracks x", "mrtrix tracks    x",
                  " mrtrix tracks", "mrtrix  tracks", "mrtrix track",
                  "MRtrix tracks")) {
    path <- tempfile(fileext = ".tck")
    ftq_reader_test_write_tck(path, ftq_reader_test_tracks(), magic = magic)
    expect_error(ftq_read_tck_header(path), "expected 'mrtrix tracks'", fixed = TRUE)
    expect_error(ftq_read_tractogram(path, reader = "memory"),
                 "expected 'mrtrix tracks'|signature does not match")
  }
})

test_that("TCK streaming crosses numeric blocks and callback chunks exactly", {
  path <- tempfile(fileext = ".tck")
  tracks <- ftq_reader_test_tracks()
  ftq_reader_test_write_tck(path, tracks)
  for (block in c(1L, 2L, 3L, 7L, 32L, 80L)) {
    for (chunk_size in c(1L, 3L, 10L)) {
      stream <- ftq_tck_stream(path, block_points = block, chunk_size = chunk_size)
      captured <- list()
      bounds <- list()
      n <- ftq_track_stream_iterate(stream, function(chunk, start, end) {
        captured <<- c(captured, chunk)
        bounds[[length(bounds) + 1L]] <<- c(start, end)
      })
      expect_identical(n, 4L)
      expect_equal(lapply(captured, `[[`, "coords"), tracks)
      expect_equal(unlist(lapply(bounds, function(x) seq.int(x[[1]], x[[2]]))), seq_len(4L))
      expect_true(all(vapply(bounds, function(x) x[[2]] - x[[1]] + 1 <= chunk_size, logical(1))))
    }
  }
})

test_that("TCK chunks also honor a coordinate budget", {
  path <- tempfile(fileext = ".tck")
  tracks <- rep(list(matrix(seq_len(24) / 4, ncol = 3L)), 7L)
  ftq_reader_test_write_tck(path, tracks)
  stream <- ftq_tck_stream(path, chunk_size = 100L, chunk_points = 16L, block_points = 5L)
  sizes <- integer()
  ftq_tck_stream_iterate(stream, function(chunk, start, end) {
    sizes <<- c(sizes, length(chunk))
  })
  expect_identical(sizes, c(2L, 2L, 2L, 1L))
})

test_that("TCK samples, limited iterations and callback stops preserve prefixes", {
  path <- tempfile(fileext = ".tck")
  tracks <- ftq_reader_test_tracks()
  ftq_reader_test_write_tck(path, tracks)
  stream <- ftq_tck_stream(path, block_points = 2L, chunk_size = 2L)
  expect_equal(lapply(ftq_track_stream_sample(stream, 2L), `[[`, "coords"), tracks[1:2])
  expect_equal(lapply(ftq_track_stream_sample(stream, 20L), `[[`, "coords"), tracks)
  expect_identical(ftq_track_stream_sample(stream, 0L), list())
  calls <- 0L
  n <- ftq_track_stream_iterate(stream, function(chunk, start, end) {
    calls <<- calls + 1L
    expect_equal(start, 1)
    expect_equal(end, 2)
    FALSE
  })
  expect_identical(n, 2L)
  expect_identical(calls, 1L)
  expect_identical(ftq_tck_stream_iterate(stream, function(...) stop("not called"),
                                         max_tracks = 0L), 0L)
  for (bad in list(-1, NA_real_, NaN, c(1, 2), "1", 1.5, -Inf)) {
    expect_error(ftq_tck_stream_iterate(stream, function(...) NULL, max_tracks = bad),
                 "max_tracks")
  }
})

test_that("TCK count is optional and total_count is not the stored record count", {
  path <- tempfile(fileext = ".tck")
  tracks <- ftq_reader_test_tracks()
  ftq_reader_test_write_tck(path, tracks, count = NULL, extra = "total_count: 10000",
                            newline = "\r\n")
  stream <- ftq_read_tractogram(path)
  expect_true(is.na(ftq_track_count(stream)))
  expect_equal(lapply(ftq_track_stream_sample(stream, 100L), `[[`, "coords"), tracks)
  expect_identical(ftq_track_count(ftq_read_tractogram(path, "memory")), 4L)
  for (count in c(0, 1, 3, 5, 5000000000)) {
    ftq_reader_test_write_tck(path, tracks, count = count)
    expect_error(ftq_read_tractogram(path, "memory"), "TCK count mismatch")
  }
  # A deliberately limited read is a prefix, not a validation of unread bytes.
  ftq_reader_test_write_tck(path, tracks, count = 999)
  expect_length(ftq_track_stream_sample(ftq_read_tractogram(path), 1L), 1L)
})

test_that("an empty TCK file parses without inventing a streamline", {
  path <- tempfile(fileext = ".tck")
  ftq_reader_test_write_tck(path)
  stream <- ftq_read_tractogram(path)
  expect_identical(ftq_track_count(stream), 0L)
  expect_identical(ftq_track_stream_sample(stream, 1L), list())
  expect_identical(ftq_read_tractogram(path, "memory")$tracks, list())
})

test_that("custom TCK fields cannot partially match the optional count field", {
  path <- tempfile(fileext = ".tck")
  tracks <- ftq_reader_test_tracks()
  for (custom in c("countdown: 99", "count_extra: not a record count")) {
    ftq_reader_test_write_tck(path, tracks, count = NULL, extra = custom)
    stream <- ftq_read_tractogram(path)
    expect_null(stream$header[["count", exact = TRUE]])
    expect_true(is.na(ftq_track_count(stream)))
    memory <- ftq_read_tractogram(path, reader = "memory")
    expect_identical(ftq_track_count(memory), 4L)
    expect_equal(lapply(memory$tracks, `[[`, "coords"), tracks)
  }
})

test_that("malformed TCK payloads fail rather than silently dropping coordinates", {
  path <- tempfile(fileext = ".tck")
  cases <- list(
    list(c(1, 2, 3, NaN, NaN, NaN), "missing Inf EOF"),
    list(c(1, 2, 3, Inf, Inf, Inf), "final streamline's NaN delimiter"),
    list(c(NaN, NaN, NaN, Inf, Inf, Inf), "Empty TCK streamline"),
    list(c(1, 2, 3, NaN, NaN, NaN, NaN, NaN, NaN, Inf, Inf, Inf), "Empty TCK streamline"),
    list(c(1, 2, NaN, Inf, Inf, Inf), "Invalid non-finite"),
    list(c(Inf, NaN, Inf, Inf, Inf, Inf), "Invalid non-finite"),
    list(c(-Inf, -Inf, -Inf, Inf, Inf, Inf), "Invalid non-finite"),
    list(c(Inf, Inf, Inf, 1, 2, 3), "trailing data"),
    list(c(Inf, Inf, Inf, Inf, Inf, Inf), "trailing data")
  )
  for (case in cases) {
    ftq_reader_test_write_tck(path, count = NULL, payload = case[[1L]])
    expect_error(ftq_read_tractogram(path, "memory"), case[[2L]])
  }
  ftq_reader_test_write_tck(path, count = NULL, payload = c(1, 2, 3, NaN, NaN, NaN, Inf))
  expect_error(ftq_read_tractogram(path), "Truncated TCK payload")
  ftq_reader_test_write_tck(path, count = NULL, payload = c(Inf, Inf, Inf))
  con <- file(path, "ab")
  writeBin(as.raw(1), con)
  close(con)
  expect_error(ftq_read_tractogram(path), "Truncated TCK payload")
})

test_that("TCK EOF must be final even when it ends a read block", {
  path <- tempfile(fileext = ".tck")
  ftq_reader_test_write_tck(path, count = NULL,
                            payload = c(1, 2, 3, rep(NaN, 3L), rep(Inf, 3L), 4, 5, 6))
  stream <- ftq_tck_stream(path, block_points = 3L)
  expect_error(ftq_tck_stream_iterate(stream, function(...) NULL), "trailing data")
})

test_that("TCK header parsing rejects invalid geometry-free storage metadata", {
  path <- tempfile(fileext = ".tck")
  make_raw <- function(text, pad = 512L, payload = rep(Inf, 3L)) {
    bytes <- charToRaw(text)
    con <- file(path, "wb")
    on.exit(close(con), add = TRUE)
    writeBin(c(bytes, raw(max(0, pad - length(bytes)))), con)
    writeBin(payload, con, size = 4L, endian = "little")
  }
  cases <- list(
    c("mrtrix tracks\ndatatype: Float32LE\nfile: ../payload.bin 0\nEND\n", "single-file"),
    c("mrtrix tracks\ndatatype: Float32LE\nfile: /private/payload.bin 0\nEND\n", "single-file"),
    c("mrtrix tracks\ndatatype: Float32LE\nfile: . -1\nEND\n", "single-file"),
    c("mrtrix tracks\ndatatype: Float32LE\nfile: . 1\nEND\n", "overlaps the header"),
    c("mrtrix tracks\ndatatype: Float32LE\nfile: . 999999\nEND\n", "outside the payload"),
    c("mrtrix tracks\ndatatype: UInt32LE\nfile: . 512\nEND\n", "Unsupported TCK datatype"),
    c("mrtrix tracks\nfile: . 512\nEND\n", "requires 'file' and 'datatype'"),
    c("mrtrix tracks\ndatatype: Float32LE\nEND\n", "requires 'file' and 'datatype'"),
    c("mrtrix tracks\ndatatype_extra: Float32LE\nfile: . 512\nEND\n", "requires 'file' and 'datatype'"),
    c("mrtrix tracks\ndatatype: Float32LE\nfile_extra: . 512\nEND\n", "requires 'file' and 'datatype'"),
    c("mrtrix tracks\ndatatype: Float32LE\ndatatype: Float64LE\nfile: . 512\nEND\n", "Duplicate"),
    c("mrtrix tracks\ndatatype: Float32LE\ncount: 1\ncount: 2\nfile: . 512\nEND\n", "Duplicate"),
    c("mrtrix tracks\ndatatype: Float32LE\ncount: -1\nfile: . 512\nEND\n", "Invalid TCK count"),
    c("mrtrix tracks\ndatatype: Float32LE\ncount: 1.5\nfile: . 512\nEND\n", "Invalid TCK count"),
    c("mrtrix tracks\ndatatype: Float32LE\ncount: 9007199254740992\nfile: . 512\nEND\n", "too large"),
    c("mrtrix tracks\nnot-a-header-entry\nfile: . 512\ndatatype: Float32LE\nEND\n", "key: value"),
    c("mrtrix tracks\n: value\nfile: . 512\ndatatype: Float32LE\nEND\n", "key: value")
  )
  for (case in cases) {
    make_raw(case[[1L]])
    expect_error(ftq_read_tck_header(path), case[[2L]])
  }
  writeBin(charToRaw("mrtrix tracks\nEND"), path)
  expect_error(ftq_read_tck_header(path), "missing newline-terminated END")
  writeBin(charToRaw(paste0("mrtrix tracks\ncomment: ", strrep("x", 100L))), path)
  expect_error(ftq_read_tck_header(path, max_header_bytes = 64L), "header exceeds")
  make_raw("mrtrix tracks\ncomment: \001\nfile: . 512\ndatatype: Float32LE\nEND\n")
  expect_error(ftq_read_tck_header(path), "non-ASCII or control")
})

test_that("reader signatures and options cannot silently select a wrong file parser", {
  path <- tempfile(fileext = ".tck")
  ftq_reader_test_write_tck(path, ftq_reader_test_tracks())
  renamed <- tempfile(fileext = ".trk")
  expect_true(file.copy(path, renamed))
  expect_error(ftq_read_tractogram(renamed), "signature does not match")
  expect_error(ftq_read_tractogram(paste0(path, ".gz")), "uncompressed")
  expect_error(ftq_read_tractogram(path, reader = "unsupported"), "arg")
  for (bad in list(0, -1, 1.2, NA_real_, Inf, c(1, 2), "2")) {
    expect_error(ftq_tck_stream(path, chunk_size = bad), "positive integer")
  }
  expect_false(ftq_is_track_stream(list()))
  expect_null(ftq_track_format(list()))
  expect_null(ftq_track_format(list(header = list(vox2ras = diag(4L)))))
  expect_null(ftq_track_format(list(header_other = list(id = "mrtrix tracks"))))
  expect_identical(ftq_track_format(list(header = list(id = "mrtrix tracks"))), "tck")
  expect_identical(ftq_track_count(list(header = list(count = 999), tracks = list(matrix(1, 2, 3)))), 1L)
  expect_identical(ftq_track_count(list(matrix(1, 2, 3), matrix(2, 3, 3))), 2L)
  expect_true(is.na(ftq_track_count(NULL)))
  expect_error(ftq_track_stream_sample(list(), 1L), "TCK or TRK")
  expect_error(ftq_track_stream_iterate(list(), function(...) NULL), "TCK or TRK")
})

test_that("TCK readers detect changed headers and exceptionally long records", {
  path <- tempfile(fileext = ".tck")
  ftq_reader_test_write_tck(path, ftq_reader_test_tracks())
  stream <- ftq_tck_stream(path, max_points_per_track = 20L, block_points = 7L)
  expect_error(ftq_tck_stream_iterate(stream, function(...) NULL), "max_points_per_track")
  stream <- ftq_tck_stream(path)
  ftq_reader_test_write_tck(path, ftq_reader_test_tracks(), count = 99L)
  expect_error(ftq_tck_stream_iterate(stream, function(...) NULL), "changed after")
})

test_that("shared dispatch preserves the existing TRK reader results", {
  path <- system.file("extdata", "fibrtraq-known-answer", "input", "synthetic.trk",
                      package = "FibrTraQ")
  if (!nzchar(path)) {
    fixtures <- list.files(system.file("extdata", "fibrtraq-known-answer", package = "FibrTraQ"),
                           pattern = "\\.trk$", recursive = TRUE, full.names = TRUE)
    expect_length(fixtures, 1L)
    path <- fixtures[[1L]]
  }
  direct <- ftq_trk_stream(path, chunk_size = 2L)
  generic <- ftq_read_tractogram(path, chunk_size = 2L)
  expect_identical(generic, direct)
  expect_identical(ftq_track_format(generic), "trk")
  read_trk <- freesurferformats::read.dti.trk
  read_args <- list(filepath = path)
  if ("coords" %in% names(formals(read_trk))) read_args$coords <- "native"
  memory_expected <- do.call(read_trk, read_args)
  memory_expected$tracks <- as.list(memory_expected$tracks)
  expect_identical(ftq_track_format(memory_expected), "trk")
  expect_equal(ftq_track_stream_sample(generic, 100L), ftq_trk_stream_sample(direct, 100L))
  expect_equal(ftq_read_tractogram(path, "memory"), memory_expected)
  expect_identical(ftq_track_count(generic), as.integer(generic$header$n_count))
})
