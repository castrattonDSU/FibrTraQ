#' FibrTraQ utilities (internal)
#'
#' Helpers for package provenance, world-to-voxel conversion, and key-based
#' membership.
#' @keywords internal
#' @noRd

ftq_package_identity <- function(description = NULL,
                                 package_path = NULL,
                                 package = "FibrTraQ") {
  package <- as.character(package)[1L]
  if (is.na(package) || !nzchar(package)) {
    stop("package must be a non-empty package name.")
  }

  if (is.null(package_path)) {
    package_path <- tryCatch(
      getNamespaceInfo(asNamespace(package), "path"),
      error = function(e) ""
    )
    if (!length(package_path) || is.na(package_path[1L]) ||
        !nzchar(package_path[1L])) {
      package_path <- system.file(package = package)
    }
  }
  package_path <- as.character(package_path)[1L]
  if (is.na(package_path) || !nzchar(package_path)) {
    package_path <- NA_character_
  } else {
    package_path <- normalizePath(package_path, mustWork = FALSE)
  }

  if (is.null(description)) {
    description_path <- if (is.na(package_path)) {
      NA_character_
    } else {
      file.path(package_path, "DESCRIPTION")
    }
    description <- if (!is.na(description_path) &&
                       file.exists(description_path)) {
      as.list(read.dcf(description_path)[1L, , drop = TRUE])
    } else {
      tryCatch(
        as.list(utils::packageDescription(package)),
        error = function(e) list()
      )
    }
  }

  field <- function(name, default = NA_character_) {
    value <- if (is.matrix(description)) {
      if (name %in% colnames(description)) description[1L, name] else NULL
    } else {
      description[[name]]
    }
    if (is.null(value) || !length(value)) return(default)
    value <- trimws(as.character(value)[1L])
    if (is.na(value) || !nzchar(value)) default else value
  }
  sanitize_remote_url <- function(url) {
    if (is.na(url)) return(url)
    url <- sub(
      "^([[:alpha:]][[:alnum:]+.-]*://)[^/@[:space:]]+@",
      "\\1",
      url,
      perl = TRUE
    )
    url <- sub("^[^/@[:space:]]+@", "", url, perl = TRUE)
    sub("[?#].*$", "", url)
  }

  field_any <- function(names, default = NA_character_) {
    for (name in names) {
      value <- field(name)
      if (!is.na(value)) return(value)
    }
    default
  }

  remote_type <- field("RemoteType")
  if (is.na(remote_type) && !is.na(field_any(c("GithubRepo", "GithubUsername")))) {
    remote_type <- "github"
  }
  remote_url <- field("RemoteUrl")
  if (!is.na(remote_url)) {
    remote_url <- sanitize_remote_url(remote_url)
  } else if (identical(tolower(remote_type), "github")) {
    remote_host <- field_any(c("RemoteHost", "GithubHost"), "github.com")
    remote_username <- field_any(c("RemoteUsername", "GithubUsername"))
    remote_repo <- field_any(c("RemoteRepo", "GithubRepo"))

    remote_host <- sub(
      "^[[:alpha:]][[:alnum:]+.-]*://",
      "",
      remote_host,
      perl = TRUE
    )
    remote_host <- sub("/.*$", "", remote_host)
    if (identical(tolower(remote_host), "api.github.com")) {
      remote_host <- "github.com"
    }
    safe_host <- grepl(
      "^[A-Za-z0-9.-]+(?::[0-9]+)?$",
      remote_host,
      perl = TRUE
    )
    safe_component <- function(value) {
      !is.na(value) && grepl("^[A-Za-z0-9._-]+$", value, perl = TRUE)
    }
    if (safe_host && safe_component(remote_username) &&
        safe_component(remote_repo)) {
      remote_url <- paste0(
        "https://",
        remote_host,
        "/",
        remote_username,
        "/",
        remote_repo
      )
    }
  }

  list(
    name = field("Package", package),
    version = field("Version"),
    remote_type = remote_type,
    remote_url = remote_url,
    remote_ref = field_any(c("RemoteRef", "GithubRef")),
    remote_sha = field_any(c("RemoteSha", "GithubSHA1")),
    built = field("Built"),
    loaded_path = package_path
  )
}

# Convert all streamline points from world (mm) to voxel indices using a reference image.
# Requires RNifti; respects qform/sform.
# Returns a deep copy of 'trk' with coords replaced.
#' @importFrom RNifti worldToVoxel
ftq_trk_world_to_voxel <- function(trk, ref_image, add_one = FALSE) {
  if (is.null(trk) || is.null(trk$tracks)) stop("Invalid TRK object.")
  trk2 <- trk
  for (i in seq_along(trk2$tracks)) {
    pts_mm <- trk2$tracks[[i]]$coords
    if (!is.matrix(pts_mm) || ncol(pts_mm) < 3L) next
    vox <- ftq_world_mm_to_voxel(pts_mm, ref_image)
    if (add_one) vox <- vox + 1
    trk2$tracks[[i]]$coords <- vox
  }
  trk2
}

ftq_trackvis_affine <- function(header) {
  if (is.null(header)) return(NULL)
  vox2ras <- ftq_header_matrix(header, c("vox2ras", "voxel_to_rasmm"))
  voxel_size <- ftq_header_vector(header, c("voxel_size", "voxel_sizes"), 3L)
  dims <- ftq_header_vector(header, c("dim", "dimensions"), 3L)
  voxel_order <- ftq_header_value(header, c("voxel_order", "voxel_order_original"))
  if (is.null(vox2ras)) return(NULL)
  if (is.null(voxel_size) || is.null(dims) || is.null(voxel_order)) {
    return(vox2ras)
  }

  voxel_order <- toupper(gsub("[^A-Za-z]", "", as.character(voxel_order)[1]))
  voxel_order <- substr(voxel_order, 1L, 3L)
  if (nchar(voxel_order) != 3L) return(vox2ras)

  scale <- diag(4)
  scale[cbind(1:3, 1:3)] <- 1 / as.numeric(voxel_size)
  offset <- diag(4)
  offset[1:3, 4] <- -0.5

  start_ornt <- ftq_axcodes_to_ornt(strsplit(voxel_order, "", fixed = TRUE)[[1]])
  end_ornt <- ftq_affine_ornt(vox2ras)
  orient <- ftq_ornt_transform(start_ornt, end_ornt)
  orient_aff <- ftq_inv_ornt_aff(orient, as.numeric(dims))

  vox2ras %*% orient_aff %*% offset %*% scale
}

ftq_trackvis_voxel_size <- function(header) {
  ftq_header_vector(header, c("voxel_size", "voxel_sizes"), 3L)
}

ftq_validate_track_voxel_size <- function(voxel_size, caller = "FibrTraQ") {
  voxel_size <- as.numeric(voxel_size)
  if (length(voxel_size) < 3L || anyNA(voxel_size[1:3]) ||
      any(!is.finite(voxel_size[1:3])) || any(voxel_size[1:3] <= 0)) {
    stop(caller, ": track_voxel_size must contain three positive finite values in mm.")
  }
  voxel_size[1:3]
}

ftq_trackvis_voxelmm_to_voxel <- function(coords, voxel_size, caller = "FibrTraQ") {
  voxel_size <- ftq_validate_track_voxel_size(voxel_size, caller = caller)
  coords <- as.matrix(coords)[, 1:3, drop = FALSE]
  sweep(coords, 2L, voxel_size, "/") + 0.5
}

ftq_densify_voxel_path <- function(coords, max_axis_step = 0.5) {
  coords <- as.matrix(coords)[, 1:3, drop = FALSE]
  if (nrow(coords) < 2L) return(coords)
  max_axis_step <- as.numeric(max_axis_step)[1]
  if (!is.finite(max_axis_step) || max_axis_step <= 0) {
    stop("max_axis_step must be a positive finite number.")
  }

  pieces <- vector("list", nrow(coords))
  for (i in seq_len(nrow(coords) - 1L)) {
    start <- coords[i, ]
    delta <- coords[i + 1L, ] - start
    if (any(!is.finite(c(start, delta)))) {
      pieces[[i]] <- matrix(start, nrow = 1L)
      next
    }
    n_steps <- max(1L, as.integer(ceiling(max(abs(delta)) / max_axis_step)))
    fraction <- seq.int(0L, n_steps - 1L) / n_steps
    pieces[[i]] <- sweep(fraction %o% delta, 2L, start, "+")
  }
  pieces[[nrow(coords)]] <- matrix(coords[nrow(coords), ], nrow = 1L)
  out <- do.call(rbind, pieces)
  colnames(out) <- colnames(coords)
  out
}

ftq_header_value <- function(header, names) {
  for (name in names) {
    if (!is.null(header[[name]])) return(header[[name]])
  }
  NULL
}

ftq_header_vector <- function(header, names, n) {
  value <- ftq_header_value(header, names)
  if (is.null(value)) return(NULL)
  value <- as.numeric(value)
  if (length(value) < n || anyNA(value[seq_len(n)])) return(NULL)
  value[seq_len(n)]
}

ftq_header_matrix <- function(header, names) {
  value <- ftq_header_value(header, names)
  if (is.null(value)) return(NULL)
  value <- as.matrix(value)
  if (!identical(dim(value), c(4L, 4L)) || anyNA(value)) return(NULL)
  storage.mode(value) <- "numeric"
  value
}

ftq_axcodes_to_ornt <- function(codes) {
  axes <- c(R = 1, L = 1, A = 2, P = 2, S = 3, I = 3)
  signs <- c(R = 1, L = -1, A = 1, P = -1, S = 1, I = -1)
  out <- cbind(axis = axes[codes], sign = signs[codes])
  if (anyNA(out)) stop("Unsupported voxel_order in TrackVis header: ", paste(codes, collapse = ""))
  out
}

ftq_affine_ornt <- function(affine) {
  codes <- character(3)
  for (j in 1:3) {
    axis <- which.max(abs(affine[1:3, j]))
    sign <- sign(affine[axis, j])
    codes[j] <- switch(axis,
      if (sign >= 0) "R" else "L",
      if (sign >= 0) "A" else "P",
      if (sign >= 0) "S" else "I"
    )
  }
  ftq_axcodes_to_ornt(codes)
}

ftq_ornt_transform <- function(start_ornt, end_ornt) {
  out <- matrix(NA_real_, nrow = 3L, ncol = 2L)
  for (i in 1:3) {
    src <- which(start_ornt[, 1] == end_ornt[i, 1])
    if (length(src) != 1L) stop("Invalid TrackVis orientation transform.")
    out[i, 1] <- src
    out[i, 2] <- start_ornt[src, 2] * end_ornt[i, 2]
  }
  out
}

ftq_inv_ornt_aff <- function(ornt, dims) {
  aff <- matrix(0, nrow = 4L, ncol = 4L)
  aff[4, 4] <- 1
  for (i in 1:3) {
    src <- as.integer(ornt[i, 1])
    direction <- as.numeric(ornt[i, 2])
    aff[i, src] <- direction
    if (direction < 0) aff[i, 4] <- dims[src] - 1
  }
  aff
}

# Build a key vector "x,y,z" for integer coords
ftq_keys <- function(m) {
  if (!is.matrix(m) || ncol(m) < 3L) stop("coords matrix with 3 columns required")
  paste(m[,1], m[,2], m[,3], sep = ",")
}
