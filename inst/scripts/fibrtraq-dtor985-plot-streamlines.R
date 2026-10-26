# Plot representative real streamlines through the dTOR-985 left-thalamus ROI.
#
# Run from RStudio with:
# source(system.file("scripts", "fibrtraq-dtor985-plot-streamlines.R", package = "FibrTraQ"))
#
# The quantitative validation remains the FibrTraQ per-voxel count table. This
# script is a visual companion: it samples actual trajectories from the full TRK
# and renders their axial/coronal/sagittal projections.

repo_root <- Sys.getenv("FIBRTRAQ_REPO_ROOT", unset = normalizePath(getwd(), mustWork = FALSE))
data_root <- Sys.getenv(
  "FIBRTRAQ_DATA_ROOT",
  unset = file.path(dirname(repo_root), "FibrTraQ-data")
)
dtor_root <- Sys.getenv(
  "FIBRTRAQ_DTOR_ROOT",
  unset = file.path(data_root, "public-validation", "dtor985")
)
left_root <- file.path(
  dtor_root,
  "outputs",
  "harvard-oxford-left-thalamus"
)
report_dir <- file.path(left_root, "full")
legacy_report <- file.path(left_root, "full-exploratory-report",
  "sub-dtor985_desc-left-thalamus-harvardoxford-full_fibrtraq_exploratory-report.html")
fig_dir <- file.path(report_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  tractogram = file.path(dtor_root, "raw", "dTOR_full_tractogram.trk"),
  focal = file.path(dtor_root, "rois", "harvard-oxford-thalamus", "focal", "Left_Thalamus.nii.gz"),
  background = file.path(dtor_root, "anat", "MNI152_T1_2mm_on_HarvardOxford.nii.gz"),
  target_table = file.path(dtor_root, "rois", "harvard-oxford-thalamus", "harvard-oxford-targets.csv"),
  report_html = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fibrtraq_target-report.html"),
  figure = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines.png"),
  figure_3d = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines-3d.png"),
  figure_3d_interactive = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines-3d-interactive.html"),
  anatomical_figure = file.path(fig_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_fig_sampled-streamlines-anatomical-overlay.png"),
  streamlines_csv = file.path(report_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_sampled-streamline-points.csv"),
  summary_csv = file.path(report_dir, "sub-dtor985_desc-left-thalamus-harvardoxford-full_sampled-streamline-summary.csv")
)

selected_targets <- c(
  "Left_Frontal_Pole",
  "Left_Superior_Frontal_Gyrus",
  "Left_Precentral_Gyrus",
  "Left_Postcentral_Gyrus",
  "Left_Putamen",
  "Left_Brain_Stem"
)

sample_per_target <- 12L
max_streamlines <- 3000000L
max_points_per_streamline <- 110L
track_chunk_size <- 50000L
overwrite <- TRUE
show_target_labels <- FALSE

required <- unlist(paths[c("tractogram", "focal", "background", "target_table")], use.names = TRUE)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing required input file(s):\n", paste(missing, collapse = "\n"))
}
if (file.exists(paths$figure) && !isTRUE(overwrite)) {
  stop("Output exists. Set overwrite <- TRUE to replace: ", paths$figure)
}

load_fibrtraq <- function() {
  if (requireNamespace("pkgload", quietly = TRUE) &&
      file.exists(file.path(getwd(), "DESCRIPTION")) &&
      identical(basename(normalizePath(getwd())), "FibrTraQ")) {
    pkgload::load_all(getwd())
  } else {
    library(FibrTraQ)
  }
  asNamespace("FibrTraQ")
}

open_output <- function(path) {
  FibrTraQ:::ftq_qc_open_output(path)
}
ns <- load_fibrtraq()
ftq_trk_stream <- get("ftq_trk_stream", envir = ns)
ftq_trk_stream_iterate <- get("ftq_trk_stream_iterate", envir = ns)
ftq_trackvis_affine <- get("ftq_trackvis_affine", envir = ns)
as_track <- get("as_track", envir = ns)
iso_coords <- get("iso_coords", envir = ns)

if (!requireNamespace("RNifti", quietly = TRUE)) stop("RNifti is required.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")

target_table <- utils::read.csv(paths$target_table, stringsAsFactors = FALSE)
missing_targets <- setdiff(selected_targets, target_table$label)
if (length(missing_targets)) {
  stop("Selected target(s) are missing from the target table: ", paste(missing_targets, collapse = ", "))
}
target_table <- target_table[match(selected_targets, target_table$label), , drop = FALSE]
target_paths <- file.path(dtor_root, target_table$path)

load_roi <- function(path) {
  img <- RNifti::readNifti(path, internal = FALSE)
  list(
    image = img,
    mask = as.array(img) > 0,
    coords = iso_coords(img)
  )
}

focal <- load_roi(paths$focal)
background <- RNifti::readNifti(paths$background, internal = FALSE)
targets <- stats::setNames(lapply(target_paths, load_roi), selected_targets)
vol_dim <- dim(focal$mask)[1:3]
ref_image <- focal$image

encode_voxels <- function(coords) {
  mm <- as.matrix(coords)[, 1:3, drop = FALSE]
  storage.mode(mm) <- "integer"
  x <- mm[, 1]
  y <- mm[, 2]
  z <- mm[, 3]
  in_bounds <- x >= 1L & x <= vol_dim[1] &
    y >= 1L & y <= vol_dim[2] &
    z >= 1L & z <= vol_dim[3]
  key <- x + (y - 1L) * vol_dim[1] + (z - 1L) * (vol_dim[1] * vol_dim[2])
  key[!in_bounds] <- NA_integer_
  key
}

make_lut <- function(coords) {
  lut <- rep(FALSE, prod(vol_dim))
  key <- encode_voxels(coords)
  lut[key[!is.na(key)]] <- TRUE
  lut
}

focal_lut <- make_lut(focal$coords)
target_luts <- lapply(targets, function(x) make_lut(x$coords))
methods <- c("round", "floor", "ceil")

discretize <- function(mapped, method) {
  out <- switch(
    method,
    round = round(mapped),
    floor = floor(mapped),
    ceil = ceiling(mapped)
  )
  storage.mode(out) <- "integer"
  out
}

hits_lut <- function(mapped, lut) {
  for (method in methods) {
    key <- encode_voxels(discretize(mapped, method))
    if (any(!is.na(key) & lut[key])) return(TRUE)
  }
  FALSE
}

decimate_points <- function(points, n) {
  if (nrow(points) <= n) return(points)
  keep <- unique(round(seq(1, nrow(points), length.out = n)))
  points[keep, , drop = FALSE]
}

world_points_for_roi <- function(roi, max_points = Inf) {
  coords <- roi$coords
  if (is.finite(max_points) && nrow(coords) > max_points) {
    keep <- unique(round(seq(1L, nrow(coords), length.out = max_points)))
    coords <- coords[keep, , drop = FALSE]
  }
  RNifti::voxelToWorld(coords, roi$image, simple = FALSE)
}

normalize_intensity <- function(x) {
  q <- stats::quantile(x[is.finite(x)], c(0.01, 0.995), na.rm = TRUE)
  out <- (x - q[1]) / max(q[2] - q[1], .Machine$double.eps)
  pmin(pmax(out, 0), 1)
}

anatomical_projection <- function(img, projection) {
  arr <- as.array(img)
  dims <- dim(arr)[1:3]
  if (identical(projection, "Axial (x-y)")) {
    mat <- apply(arr, c(1, 2), max, na.rm = TRUE)
    vox <- expand.grid(i = seq_len(dims[1]), j = seq_len(dims[2]))
    world <- RNifti::voxelToWorld(cbind(vox$i, vox$j, round(dims[3] / 2)), img, simple = FALSE)
    out <- data.frame(projection = projection, px = world[, 1], py = world[, 2], intensity = as.vector(mat))
  } else if (identical(projection, "Coronal (x-z)")) {
    mat <- apply(arr, c(1, 3), max, na.rm = TRUE)
    vox <- expand.grid(i = seq_len(dims[1]), k = seq_len(dims[3]))
    world <- RNifti::voxelToWorld(cbind(vox$i, round(dims[2] / 2), vox$k), img, simple = FALSE)
    out <- data.frame(projection = projection, px = world[, 1], py = world[, 3], intensity = as.vector(mat))
  } else if (identical(projection, "Sagittal (y-z)")) {
    mat <- apply(arr, c(2, 3), max, na.rm = TRUE)
    vox <- expand.grid(j = seq_len(dims[2]), k = seq_len(dims[3]))
    world <- RNifti::voxelToWorld(cbind(round(dims[1] / 2), vox$j, vox$k), img, simple = FALSE)
    out <- data.frame(projection = projection, px = world[, 2], py = world[, 3], intensity = as.vector(mat))
  } else {
    stop("Unsupported projection: ", projection)
  }
  out$projection <- factor(out$projection, levels = c("Axial (x-y)", "Coronal (x-z)", "Sagittal (y-z)"))
  out
}

make_anatomical_projection_df <- function(img) {
  out <- do.call(
    rbind,
    lapply(c("Axial (x-y)", "Coronal (x-z)", "Sagittal (y-z)"), anatomical_projection, img = img)
  )
  out$intensity <- normalize_intensity(out$intensity)
  out
}

project_3d <- function(points, azimuth_degrees = -42, elevation_degrees = 24) {
  p <- as.matrix(points)[, 1:3, drop = FALSE]
  az <- azimuth_degrees * pi / 180
  el <- elevation_degrees * pi / 180
  rz <- matrix(
    c(cos(az), -sin(az), 0,
      sin(az),  cos(az), 0,
      0,        0,       1),
    nrow = 3,
    byrow = TRUE
  )
  rx <- matrix(
    c(1, 0,       0,
      0, cos(el), -sin(el),
      0, sin(el),  cos(el)),
    nrow = 3,
    byrow = TRUE
  )
  rotated <- p %*% t(rx %*% rz)
  data.frame(u = rotated[, 1], v = rotated[, 2], depth = rotated[, 3])
}

brain_cloud_world <- function(img, max_points = 12000L) {
  arr <- as.array(img)
  threshold <- stats::quantile(arr[arr > 0], 0.20, na.rm = TRUE)
  ijk <- which(arr > threshold, arr.ind = TRUE)
  if (nrow(ijk) > max_points) {
    set.seed(20260605)
    ijk <- ijk[sort(sample.int(nrow(ijk), max_points)), , drop = FALSE]
  }
  world <- RNifti::voxelToWorld(ijk, img, simple = FALSE)
  data.frame(x = world[, 1], y = world[, 2], z = world[, 3])
}

brain_cloud <- function(img, max_points = 12000L) {
  world <- brain_cloud_world(img, max_points = max_points)
  projected <- project_3d(world)
  depth <- normalize_intensity(projected$depth)
  data.frame(
    u = projected$u,
    v = projected$v,
    depth = projected$depth,
    alpha = 0.035 + 0.16 * depth
  )
}

rows_3d <- function(x) {
  m <- as.matrix(x[, c("x", "y", "z"), drop = FALSE])
  lapply(seq_len(nrow(m)), function(i) unname(as.numeric(m[i, ])))
}

write_interactive_3d_html <- function(path,
                                      streamline_points,
                                      focal_world,
                                      target_centers,
                                      brain_world,
                                      colors,
                                      scanned,
                                      focal_hits) {
  target_labels <- gsub("_", " ", names(colors))
  names(colors) <- target_labels

  streamlines <- unname(lapply(split(streamline_points, streamline_points$streamline_id), function(x) {
    label <- gsub("_", " ", x$target[1])
    list(
      id = x$streamline_id[1],
      target = label,
      color = unname(colors[[label]]),
      points = rows_3d(x)
    )
  }))

  focal_df <- data.frame(x = focal_world[, 1], y = focal_world[, 2], z = focal_world[, 3])
  target_df <- target_centers
  target_df$target_label <- gsub("_", " ", target_df$target)
  target_payload <- lapply(seq_len(nrow(target_df)), function(i) {
    label <- target_df$target_label[i]
    list(
      target = label,
      color = unname(colors[[label]]),
      point = unname(as.numeric(target_df[i, c("x", "y", "z")]))
    )
  })

  all_world <- rbind(
    streamline_points[, c("x", "y", "z")],
    focal_df,
    target_centers[, c("x", "y", "z")],
    brain_world[, c("x", "y", "z")]
  )
  bounds <- apply(all_world, 2, range, finite = TRUE)
  center <- colMeans(bounds)
  radius <- max(bounds[2, ] - bounds[1, ]) / 2

  payload <- list(
    metadata = list(
      title = "Interactive 3D QC: Sampled Streamlines Through Left Thalamus",
      subtitle = paste0(
        "Sampled from full dTOR-985 tractogram; scanned ",
        format(scanned, big.mark = ","),
        " streamlines and observed ",
        format(focal_hits, big.mark = ","),
        " focal hits while sampling."
      ),
      note = paste0(
        "Drag to rotate. Scroll or pinch to zoom. All ",
        format(nrow(focal_df), big.mark = ","),
        " focal ROI voxel centers are shown as a translucent blue volume. ",
        "This viewer is for qualitative QC; use FibrTraQ tables for exact counts."
      )
    ),
    bounds = list(center = unname(as.numeric(center)), radius = unname(as.numeric(radius))),
    brain = rows_3d(brain_world),
    focal = rows_3d(focal_df),
    focalCount = nrow(focal_df),
    targets = target_payload,
    streamlines = streamlines,
    legend = lapply(target_labels, function(label) list(target = label, color = unname(colors[[label]])))
  )
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = 7, pretty = FALSE)

  html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Interactive 3D FibrTraQ QC</title>
<style>
  :root { color-scheme: light; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }
  body { margin: 0; background: #f8fafc; color: #111827; }
  main { width: min(1500px, calc(100vw - 28px)); margin: 14px auto; }
  h1 { margin: 0 0 4px; font-size: 24px; line-height: 1.15; }
  p { margin: 0; color: #475569; font-size: 14px; }
  .toolbar { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin: 12px 0; }
  button { border: 1px solid #cbd5e1; background: white; border-radius: 6px; padding: 7px 10px; color: #111827; cursor: pointer; }
  button[aria-pressed="true"] { background: #e0f2fe; border-color: #0284c7; }
  .viewer { position: relative; height: min(78vh, 900px); min-height: 560px; border: 1px solid #cbd5e1; border-radius: 8px; overflow: hidden; background: white; }
  canvas { display: block; width: 100%; height: 100%; cursor: grab; touch-action: none; }
  canvas.dragging { cursor: grabbing; }
  .hud { position: absolute; left: 12px; top: 12px; max-width: 330px; padding: 9px 10px; border-radius: 7px; background: rgba(255,255,255,0.88); border: 1px solid rgba(203,213,225,0.9); font-size: 12px; color: #334155; }
  .legend { display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); gap: 6px 14px; margin-top: 10px; font-size: 13px; color: #334155; }
  .legend span { display: inline-flex; align-items: center; gap: 7px; min-width: 0; }
  .swatch { width: 23px; height: 4px; border-radius: 99px; display: inline-block; flex: 0 0 auto; }
</style>
</head>
<body>
<main>
  <h1 id="title"></h1>
  <p id="subtitle"></p>
  <div class="toolbar">
    <button id="reset" type="button">Reset View</button>
    <button id="toggleBrain" type="button" aria-pressed="true">Brain Context</button>
    <button id="toggleFocal" type="button" aria-pressed="true">Focal ROI Volume</button>
    <button id="toggleLabels" type="button" aria-pressed="true">Target Labels</button>
  </div>
  <section class="viewer" aria-label="Interactive 3D streamline quality-control viewer">
    <canvas id="scene"></canvas>
    <div class="hud" id="hud"></div>
  </section>
  <div class="legend" id="legend"></div>
</main>
<script id="scene-data" type="application/json">', json, '</script>
<script>
const data = JSON.parse(document.getElementById("scene-data").textContent);
const canvas = document.getElementById("scene");
const ctx = canvas.getContext("2d");
const state = { yaw: -0.75, pitch: 0.42, zoom: 1.0, showBrain: true, showFocal: true, showLabels: true };
let pointer = null;

document.getElementById("title").textContent = data.metadata.title;
document.getElementById("subtitle").textContent = data.metadata.subtitle;
document.getElementById("hud").textContent = data.metadata.note;
document.getElementById("legend").innerHTML = data.legend.map(d => `<span><i class="swatch" style="background:${d.color}"></i>${d.target}</span>`).join("");

function setPressed(id, value) {
  document.getElementById(id).setAttribute("aria-pressed", value ? "true" : "false");
}

document.getElementById("reset").addEventListener("click", () => {
  state.yaw = -0.75; state.pitch = 0.42; state.zoom = 1.0; render();
});
document.getElementById("toggleBrain").addEventListener("click", () => {
  state.showBrain = !state.showBrain; setPressed("toggleBrain", state.showBrain); render();
});
document.getElementById("toggleFocal").addEventListener("click", () => {
  state.showFocal = !state.showFocal; setPressed("toggleFocal", state.showFocal); render();
});
document.getElementById("toggleLabels").addEventListener("click", () => {
  state.showLabels = !state.showLabels; setPressed("toggleLabels", state.showLabels); render();
});

function resize() {
  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.max(1, Math.round(rect.width * dpr));
  canvas.height = Math.max(1, Math.round(rect.height * dpr));
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  render();
}

window.addEventListener("resize", resize);
canvas.addEventListener("pointerdown", event => {
  pointer = { x: event.clientX, y: event.clientY };
  canvas.setPointerCapture(event.pointerId);
  canvas.classList.add("dragging");
});
canvas.addEventListener("pointermove", event => {
  if (!pointer) return;
  const dx = event.clientX - pointer.x;
  const dy = event.clientY - pointer.y;
  pointer = { x: event.clientX, y: event.clientY };
  state.yaw += dx * 0.008;
  state.pitch = Math.max(-1.45, Math.min(1.45, state.pitch + dy * 0.008));
  render();
});
canvas.addEventListener("pointerup", event => {
  pointer = null;
  canvas.classList.remove("dragging");
});
canvas.addEventListener("pointercancel", () => {
  pointer = null;
  canvas.classList.remove("dragging");
});
canvas.addEventListener("wheel", event => {
  event.preventDefault();
  const factor = Math.exp(-event.deltaY * 0.001);
  state.zoom = Math.max(0.35, Math.min(5.0, state.zoom * factor));
  render();
}, { passive: false });

function project(point) {
  const c = data.bounds.center;
  let x = point[0] - c[0];
  let y = point[1] - c[1];
  let z = point[2] - c[2];
  const cy = Math.cos(state.yaw), sy = Math.sin(state.yaw);
  const cp = Math.cos(state.pitch), sp = Math.sin(state.pitch);
  const x1 = cy * x + sy * z;
  const z1 = -sy * x + cy * z;
  const y1 = cp * y - sp * z1;
  const z2 = sp * y + cp * z1;
  const rect = canvas.getBoundingClientRect();
  const scale = Math.min(rect.width, rect.height) / (data.bounds.radius * 2.25) * state.zoom;
  return { x: rect.width / 2 + x1 * scale, y: rect.height / 2 - y1 * scale, z: z2 };
}

function rgba(hex, alpha) {
  const h = hex.replace("#", "");
  const r = parseInt(h.slice(0, 2), 16);
  const g = parseInt(h.slice(2, 4), 16);
  const b = parseInt(h.slice(4, 6), 16);
  return `rgba(${r},${g},${b},${alpha})`;
}

function renderPointCloud(points, color, alpha, radius) {
  const projected = points.map(project).sort((a, b) => a.z - b.z);
  ctx.fillStyle = rgba(color, alpha);
  for (const p of projected) {
    ctx.beginPath();
    ctx.arc(p.x, p.y, radius, 0, Math.PI * 2);
    ctx.fill();
  }
}

function renderFocalVolume(points) {
  const projected = points.map(project).sort((a, b) => a.z - b.z);
  ctx.fillStyle = "rgba(0, 166, 214, 0.14)";
  ctx.strokeStyle = "rgba(0, 109, 149, 0.28)";
  ctx.lineWidth = 0.55;
  for (const p of projected) {
    ctx.beginPath();
    ctx.arc(p.x, p.y, 3.1, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
  }
  ctx.fillStyle = "rgba(0, 109, 149, 0.45)";
  for (const p of projected) {
    ctx.beginPath();
    ctx.arc(p.x, p.y, 1.15, 0, Math.PI * 2);
    ctx.fill();
  }
}

function renderStreamlines() {
  for (const s of data.streamlines) {
    ctx.strokeStyle = rgba(s.color, 0.72);
    ctx.lineWidth = 1.8;
    ctx.lineJoin = "round";
    ctx.lineCap = "round";
    ctx.beginPath();
    s.points.forEach((pt, i) => {
      const p = project(pt);
      if (i === 0) ctx.moveTo(p.x, p.y); else ctx.lineTo(p.x, p.y);
    });
    ctx.stroke();
  }
}

function renderTargets() {
  ctx.font = "12px -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif";
  ctx.textBaseline = "middle";
  for (const t of data.targets) {
    const p = project(t.point);
    ctx.fillStyle = "white";
    ctx.strokeStyle = t.color;
    ctx.lineWidth = 3;
    ctx.beginPath();
    ctx.arc(p.x, p.y, 7, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
    if (state.showLabels) {
      ctx.fillStyle = "rgba(255,255,255,0.86)";
      const width = ctx.measureText(t.target).width + 10;
      ctx.fillRect(p.x + 10, p.y - 10, width, 20);
      ctx.fillStyle = "#111827";
      ctx.fillText(t.target, p.x + 15, p.y);
    }
  }
}

function renderAxes() {
  const c = data.bounds.center;
  const r = data.bounds.radius * 0.42;
  const axes = [
    { label: "x", color: "#64748b", pts: [[c[0]-r, c[1]-r, c[2]-r], [c[0], c[1]-r, c[2]-r]] },
    { label: "y", color: "#64748b", pts: [[c[0]-r, c[1]-r, c[2]-r], [c[0]-r, c[1], c[2]-r]] },
    { label: "z", color: "#64748b", pts: [[c[0]-r, c[1]-r, c[2]-r], [c[0]-r, c[1]-r, c[2]]] }
  ];
  ctx.font = "12px -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif";
  for (const axis of axes) {
    const a = project(axis.pts[0]);
    const b = project(axis.pts[1]);
    ctx.strokeStyle = axis.color;
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    ctx.moveTo(a.x, a.y);
    ctx.lineTo(b.x, b.y);
    ctx.stroke();
    ctx.fillStyle = "#334155";
    ctx.fillText(axis.label, b.x + 5, b.y);
  }
}

function render() {
  const rect = canvas.getBoundingClientRect();
  ctx.clearRect(0, 0, rect.width, rect.height);
  if (state.showBrain) renderPointCloud(data.brain, "#64748b", 0.15, 1.15);
  renderStreamlines();
  if (state.showFocal) renderFocalVolume(data.focal);
  renderTargets();
  renderAxes();
}

resize();
</script>
</body>
</html>'
  )

  writeLines(html, path)
  invisible(path)
}

sampling_done <- function() {
  stop(structure(
    list(message = "Requested streamline samples collected.", call = NULL),
    class = c("fibrtraq_sampling_done", "error", "condition")
  ))
}

stream <- ftq_trk_stream(paths$tractogram, chunk_size = track_chunk_size)
track_affine <- ftq_trackvis_affine(stream$header)
if (is.null(track_affine)) stop("The TRK header did not provide a usable TrackVis affine.")

map_track <- function(track) {
  trk <- as_track(track)
  raw <- as.matrix(trk$coords)[, 1:3, drop = FALSE]
  world <- cbind(raw, 1) %*% t(track_affine)
  world <- world[, 1:3, drop = FALSE]
  mapped <- RNifti::worldToVoxel(world, ref_image, simple = FALSE)
  list(world = world, voxel = mapped)
}

samples <- stats::setNames(vector("list", length(selected_targets)), selected_targets)
for (target in selected_targets) samples[[target]] <- list()
hits_seen <- stats::setNames(integer(length(selected_targets)), selected_targets)
scanned <- 0L
focal_hits <- 0L
done <- FALSE

message("Sampling real streamlines from ", paths$tractogram)
message("Targets: ", paste(gsub("_", " ", selected_targets), collapse = "; "))

tryCatch(
  ftq_trk_stream_iterate(
    stream,
    callback = function(chunk, start_index, end_index) {
      if (done) sampling_done()
    for (track in chunk) {
      scanned <<- scanned + 1L
      if (scanned > max_streamlines) {
        done <<- TRUE
        break
      }

      mapped <- tryCatch(map_track(track), error = function(e) NULL)
      if (is.null(mapped)) next
      if (!hits_lut(mapped$voxel, focal_lut)) next

      focal_hits <<- focal_hits + 1L
      target_hit <- vapply(target_luts, function(lut) hits_lut(mapped$voxel, lut), logical(1))
      available <- names(target_hit)[target_hit & lengths(samples) < sample_per_target]
      if (!length(available)) next

      n_by_target <- lengths(samples)[available]
      chosen <- available[which.min(n_by_target)]
      hits_seen[chosen] <<- hits_seen[chosen] + 1L
      samples[[chosen]][[length(samples[[chosen]]) + 1L]] <<- decimate_points(
        mapped$world,
        max_points_per_streamline
      )

      if (scanned %% 10000L == 0L) {
        message(
          "Scanned ", format(scanned, big.mark = ","),
          "; focal hits ", format(focal_hits, big.mark = ","),
          "; sampled ", paste(lengths(samples), collapse = "/")
        )
      }

      if (all(lengths(samples) >= sample_per_target)) {
        done <<- TRUE
        sampling_done()
      }
    }
    },
    max_tracks = max_streamlines,
    chunk_size = track_chunk_size
  ),
  fibrtraq_sampling_done = function(e) invisible(NULL)
)

sample_counts <- lengths(samples)
if (!any(sample_counts)) {
  stop("No streamlines matching the selected target set were found.")
}

streamline_points <- do.call(
  rbind,
  unlist(
    lapply(names(samples), function(target) {
      target_samples <- samples[[target]]
      if (!length(target_samples)) return(list())
      lapply(seq_along(target_samples), function(i) {
        pts <- target_samples[[i]]
        data.frame(
          target = target,
          streamline_id = paste(target, i, sep = "__"),
          point_index = seq_len(nrow(pts)),
          x = pts[, 1],
          y = pts[, 2],
          z = pts[, 3],
          stringsAsFactors = FALSE
        )
      })
    }),
    recursive = FALSE
  )
)

summary <- data.frame(
  target = selected_targets,
  sampled_streamlines = as.integer(sample_counts[selected_targets]),
  target_hits_seen_while_sampling = as.integer(hits_seen[selected_targets]),
  stringsAsFactors = FALSE
)
summary <- rbind(
  data.frame(
    target = c("streamlines_scanned", "focal_hits_seen"),
    sampled_streamlines = c(scanned, focal_hits),
    target_hits_seen_while_sampling = NA_integer_,
    stringsAsFactors = FALSE
  ),
  summary
)
utils::write.csv(streamline_points, paths$streamlines_csv, row.names = FALSE)
utils::write.csv(summary, paths$summary_csv, row.names = FALSE)

projection_df <- rbind(
  transform(streamline_points, projection = "Axial (x-y)", px = x, py = y),
  transform(streamline_points, projection = "Coronal (x-z)", px = x, py = z),
  transform(streamline_points, projection = "Sagittal (y-z)", px = y, py = z)
)
projection_df$target_label <- gsub("_", " ", projection_df$target)
projection_df$projection <- factor(
  projection_df$projection,
  levels = c("Axial (x-y)", "Coronal (x-z)", "Sagittal (y-z)")
)

focal_world <- world_points_for_roi(focal)
focal_projection <- rbind(
  data.frame(projection = "Axial (x-y)", px = focal_world[, 1], py = focal_world[, 2]),
  data.frame(projection = "Coronal (x-z)", px = focal_world[, 1], py = focal_world[, 3]),
  data.frame(projection = "Sagittal (y-z)", px = focal_world[, 2], py = focal_world[, 3])
)
focal_projection$projection <- factor(focal_projection$projection, levels = levels(projection_df$projection))

target_centers <- do.call(
  rbind,
  lapply(names(targets), function(target) {
    pts <- RNifti::voxelToWorld(targets[[target]]$coords, targets[[target]]$image, simple = FALSE)
    center <- colMeans(pts)
    data.frame(target = target, x = center[1], y = center[2], z = center[3], stringsAsFactors = FALSE)
  })
)
center_projection <- rbind(
  transform(target_centers, projection = "Axial (x-y)", px = x, py = y),
  transform(target_centers, projection = "Coronal (x-z)", px = x, py = z),
  transform(target_centers, projection = "Sagittal (y-z)", px = y, py = z)
)
center_projection$target_label <- gsub("_", " ", center_projection$target)
center_projection$projection <- factor(center_projection$projection, levels = levels(projection_df$projection))

colors <- stats::setNames(
  grDevices::hcl.colors(length(selected_targets), "Dark 3"),
  gsub("_", " ", selected_targets)
)

target_label_layer <- if (!isTRUE(show_target_labels)) {
  NULL
} else if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(
    data = center_projection,
    ggplot2::aes(px, py, label = target_label),
    size = 2.0,
    min.segment.length = 0,
    segment.size = 0.2,
    box.padding = 0.18,
    point.padding = 0.08,
    max.overlaps = Inf,
    show.legend = FALSE
  )
} else {
  ggplot2::geom_text(
    data = center_projection,
    ggplot2::aes(px, py, label = target_label),
    size = 1.9,
    hjust = -0.04,
    vjust = -0.35,
    check_overlap = TRUE,
    show.legend = FALSE
  )
}

p <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = focal_projection,
    ggplot2::aes(px, py),
    color = "#1f9bd1",
    alpha = 0.11,
    size = 0.7
  ) +
  ggplot2::geom_path(
    data = projection_df,
    ggplot2::aes(px, py, group = streamline_id, color = target_label),
    linewidth = 0.34,
    alpha = 0.42,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    data = center_projection,
    ggplot2::aes(px, py, fill = target_label),
    shape = 21,
    color = "black",
    stroke = 0.35,
    size = 2.3
  ) +
  target_label_layer +
  ggplot2::facet_wrap(~ projection, nrow = 1) +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::scale_color_manual(values = colors, name = "Selected target") +
  ggplot2::scale_fill_manual(values = colors, guide = "none") +
  ggplot2::labs(
    title = "Sampled Real Streamlines Through Left Thalamus",
    subtitle = paste0(
      "Sampled from the full dTOR-985 tractogram; ",
      sample_per_target,
      " trajectories requested per selected target"
    ),
    x = "MNI coordinate (mm)",
    y = "MNI coordinate (mm)",
    caption = paste0(
      "Blue points show the focal left-thalamus ROI. Lines follow target quotas in scan order, not random sampling or population proportions. ",
      "Scanned ", format(scanned, big.mark = ","), " streamlines; focal hits seen while sampling: ",
      format(focal_hits, big.mark = ","), "."
    )
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 9),
    plot.caption = ggplot2::element_text(size = 7.5, color = "#374151"),
    strip.background = ggplot2::element_rect(fill = "#f3f4f6", color = "#d1d5db"),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold"),
    legend.text = ggplot2::element_text(size = 7.5),
    legend.key.width = grid::unit(18, "pt"),
    plot.margin = ggplot2::margin(10, 20, 10, 10)
  )

ggplot2::ggsave(paths$figure, p, width = 14.5, height = 6.1, dpi = 320)

streamline_3d <- cbind(
  streamline_points,
  project_3d(streamline_points[, c("x", "y", "z")])
)
streamline_3d$target_label <- gsub("_", " ", streamline_3d$target)

focal_3d <- data.frame(
  project_3d(focal_world),
  type = "Left thalamus focal ROI"
)

center_3d <- cbind(
  target_centers,
  project_3d(target_centers[, c("x", "y", "z")])
)
center_3d$target_label <- gsub("_", " ", center_3d$target)

brain_world <- brain_cloud_world(background)
brain_3d <- brain_cloud(background)

axis_origin <- c(x = -82, y = -112, z = -65)
axes_world <- rbind(
  data.frame(axis = "x", x = c(axis_origin[["x"]], axis_origin[["x"]] + 45), y = axis_origin[["y"]], z = axis_origin[["z"]]),
  data.frame(axis = "y", x = axis_origin[["x"]], y = c(axis_origin[["y"]], axis_origin[["y"]] + 45), z = axis_origin[["z"]]),
  data.frame(axis = "z", x = axis_origin[["x"]], y = axis_origin[["y"]], z = c(axis_origin[["z"]], axis_origin[["z"]] + 45))
)
axes_3d <- cbind(axes_world, project_3d(axes_world[, c("x", "y", "z")]))

p3d <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = brain_3d,
    ggplot2::aes(u, v, alpha = alpha),
    color = "#6b7280",
    size = 0.24,
    show.legend = FALSE
  ) +
  ggplot2::geom_path(
    data = axes_3d,
    ggplot2::aes(u, v, group = axis),
    color = "#6b7280",
    linewidth = 0.35,
    alpha = 0.8
  ) +
  ggplot2::geom_path(
    data = streamline_3d,
    ggplot2::aes(u, v, group = streamline_id, color = target_label),
    linewidth = 0.52,
    alpha = 0.70,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    data = focal_3d,
    ggplot2::aes(u, v),
    color = "#00a6d6",
    alpha = 0.32,
    size = 0.75
  ) +
  ggplot2::geom_point(
    data = center_3d,
    ggplot2::aes(u, v, color = target_label),
    shape = 21,
    fill = "white",
    stroke = 0.8,
    size = 3.1
  ) +
  ggplot2::geom_text(
    data = axes_3d[c(2, 4, 6), , drop = FALSE],
    ggplot2::aes(u, v, label = paste0(axis, " axis")),
    color = "#374151",
    size = 2.7,
    hjust = -0.15,
    vjust = -0.15
  ) +
  ggplot2::coord_equal() +
  ggplot2::scale_color_manual(values = colors, name = "Selected target") +
  ggplot2::scale_alpha_identity() +
  ggplot2::labs(
    title = "3D Perspective Rendering of Sampled Streamlines",
    subtitle = "MNI-space streamline trajectories through the left thalamus, shown with target centroids and a light MNI brain context cloud.",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Static 3D perspective view: azimuth -42 degrees, elevation 24 degrees. ",
      "Blue points mark the focal ROI; colored lines are sampled real streamlines."
    )
  ) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.title = ggplot2::element_text(face = "bold", size = 14, color = "#111827"),
    plot.subtitle = ggplot2::element_text(size = 9, color = "#374151"),
    plot.caption = ggplot2::element_text(size = 7.6, color = "#374151"),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", color = "#111827"),
    legend.text = ggplot2::element_text(size = 7.5, color = "#111827"),
    legend.key.width = grid::unit(18, "pt"),
    plot.margin = ggplot2::margin(10, 12, 10, 12)
  )

ggplot2::ggsave(paths$figure_3d, p3d, width = 11.5, height = 8.0, dpi = 320, bg = "white")

qc_manifest <- list(
  dataset_root = dtor_root,
  subject = "sub-dtor985",
  description = "left-thalamus-harvardoxford-streamline-qc",
  fnode = paths$focal,
  tractogram = paths$tractogram,
  targets = stats::setNames(as.list(target_paths), selected_targets),
  background_image = paths$background,
  coordinate = list(
    track_space = "trackvis",
    discretization = "any",
    axis_flip = c(FALSE, FALSE, FALSE)
  ),
  execution = list(
    track_reader = "stream",
    track_chunk_size = track_chunk_size
  ),
  outputs = list(
    output_dir = report_dir,
    write_maps = FALSE,
    write_report = FALSE,
    overwrite = TRUE
  )
)
qc_viewer <- ftq_interactive_qc_viewer(
  qc_manifest,
  output_html = paths$figure_3d_interactive,
  target_labels = selected_targets,
  sample_per_target = sample_per_target,
  max_streamlines = max_streamlines,
  max_points_per_streamline = max_points_per_streamline,
  brain_points = 12000L,
  overwrite = TRUE
)

background_projection <- make_anatomical_projection_df(background)
anatomical <- ggplot2::ggplot() +
  ggplot2::geom_raster(
    data = background_projection,
    ggplot2::aes(px, py, fill = intensity),
    interpolate = FALSE
  ) +
  ggplot2::geom_point(
    data = focal_projection,
    ggplot2::aes(px, py),
    color = "#00a6d6",
    alpha = 0.18,
    size = 0.62
  ) +
  ggplot2::geom_path(
    data = projection_df,
    ggplot2::aes(px, py, group = streamline_id, color = target_label),
    linewidth = 0.32,
    alpha = 0.58,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    data = center_projection,
    ggplot2::aes(px, py, color = target_label),
    shape = 21,
    fill = "white",
    stroke = 0.55,
    size = 2.4
  ) +
  ggplot2::facet_wrap(~ projection, nrow = 1) +
  ggplot2::coord_equal(clip = "on") +
  ggplot2::scale_fill_gradient(low = "black", high = "white", guide = "none") +
  ggplot2::scale_color_manual(values = colors, name = "Selected target") +
  ggplot2::labs(
    title = "Sampled Streamlines Overlaid on MNI T1 Anatomy",
    subtitle = "Background is the actual MNI152 T1 image projected through the same NIfTI affine; streamline points are in the same MNI world coordinates.",
    x = "MNI coordinate (mm)",
    y = "MNI coordinate (mm)",
    caption = paste0(
      "Blue points mark the focal left-thalamus ROI. Colored lines follow target quotas in scan order, not random sampling or population proportions. ",
      "This is an affine-registered anatomical projection, not a single histological slice."
    )
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 8.8),
    plot.caption = ggplot2::element_text(size = 7.4, color = "#374151"),
    strip.background = ggplot2::element_rect(fill = "#f3f4f6", color = "#d1d5db"),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "black", color = NA),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold"),
    legend.text = ggplot2::element_text(size = 7.5),
    legend.key.width = grid::unit(18, "pt")
  )

ggplot2::ggsave(paths$anatomical_figure, anatomical, width = 14.5, height = 6.1, dpi = 320)

append_report_figure <- function(html_path, figure_path, section_title, paragraph, caption) {
  if (!file.exists(html_path)) return(invisible(FALSE))
  html <- readLines(html_path, warn = FALSE)
  rel <- FibrTraQ:::ftq_report_href(figure_path, dirname(html_path))
  if (any(grepl(rel, html, fixed = TRUE) | grepl(basename(figure_path), html, fixed = TRUE))) return(invisible(TRUE))
  section <- paste0(
    "<h2>", section_title, "</h2>",
    "<p>", paragraph, "</p>",
    "<figure><img src='", rel, "' alt='", caption, "'>",
    "<figcaption>", caption, "</figcaption></figure>"
  )
  marker <- "<h2>Output Files</h2>"
  pos <- grep(marker, html, fixed = TRUE)
  if (length(pos)) {
    html <- append(html, section, after = pos[1] - 1L)
  } else {
    html <- append(html, section, after = length(html) - 1L)
  }
  writeLines(html, html_path)
  invisible(TRUE)
}

append_report_interactive_3d <- function(html_path, viewer_path, section_title, paragraph) {
  if (!file.exists(html_path)) return(invisible(FALSE))
  html <- readLines(html_path, warn = FALSE)
  rel <- FibrTraQ:::ftq_report_href(viewer_path, dirname(html_path))
  if (any(grepl(rel, html, fixed = TRUE) | grepl(basename(viewer_path), html, fixed = TRUE))) return(invisible(TRUE))
  section <- paste0(
    "<h2>", section_title, "</h2>",
    "<p>", paragraph, "</p>",
    "<p><a href='", rel, "'>Open the interactive 3D QC viewer in a separate browser tab.</a></p>",
    "<iframe src='", rel, "' title='", section_title, "' ",
    "style='width:100%;height:760px;border:1px solid #d1d5db;border-radius:8px;background:#fff;'></iframe>"
  )
  marker <- "<h2>Output Files</h2>"
  pos <- grep(marker, html, fixed = TRUE)
  if (length(pos)) {
    html <- append(html, section, after = pos[1] - 1L)
  } else {
    html <- append(html, section, after = length(html) - 1L)
  }
  writeLines(html, html_path)
  invisible(TRUE)
}

refresh_static_3d_report_section <- function(html_path, figure_path) {
  if (!file.exists(html_path)) return(invisible(FALSE))
  html <- readLines(html_path, warn = FALSE)
  rel <- FibrTraQ:::ftq_report_href(figure_path, dirname(html_path))
  hit <- which(grepl(rel, html, fixed = TRUE) | grepl(basename(figure_path), html, fixed = TRUE))
  if (!length(hit)) return(invisible(FALSE))
  replacement <- paste0(
    "<h2>Static 3D Snapshot</h2>",
    "<p>This PNG is only a static snapshot of the interactive 3D QC view. ",
    "Use the interactive viewer above when checking spatial geometry.</p>",
    "<figure><img src='", rel,
    "' alt='Static snapshot of sampled streamlines through the left thalamus.'>",
    "<figcaption>Static snapshot of sampled streamlines through the left thalamus.</figcaption></figure>"
  )
  html[hit[1]] <- replacement
  writeLines(html, html_path)
  invisible(TRUE)
}

order_report_3d_sections <- function(html_path, interactive_path, static_path) {
  if (!file.exists(html_path)) return(invisible(FALSE))
  html <- readLines(html_path, warn = FALSE)
  interactive_href <- FibrTraQ:::ftq_report_href(interactive_path, dirname(html_path))
  static_href <- FibrTraQ:::ftq_report_href(static_path, dirname(html_path))
  interactive <- which(grepl(interactive_href, html, fixed = TRUE) | grepl(basename(interactive_path), html, fixed = TRUE))
  static <- which(grepl(static_href, html, fixed = TRUE) | grepl(basename(static_path), html, fixed = TRUE))
  if (!length(interactive) || !length(static) || interactive[1] < static[1]) {
    return(invisible(TRUE))
  }
  interactive_line <- html[interactive[1]]
  html <- html[-interactive[1]]
  static <- which(grepl(static_href, html, fixed = TRUE) | grepl(basename(static_path), html, fixed = TRUE))
  html <- append(html, interactive_line, after = static[1] - 1L)
  writeLines(html, html_path)
  invisible(TRUE)
}

append_report_interactive_3d(
  paths$report_html,
  paths$figure_3d_interactive,
  "Interactive 3D QC Viewer",
  "This local HTML viewer renders sampled real streamlines in rotatable MNI-space 3D. Drag to rotate and scroll to zoom. It is intended for quality control and spatial interpretation, while the exact quantitative evidence remains in the FibrTraQ tables."
)
append_report_figure(
  paths$report_html,
  paths$figure_3d,
  "Static 3D Snapshot",
  "This figure is a static snapshot of the interactive 3D QC view. Use the interactive viewer above when checking spatial geometry.",
  "Static snapshot of sampled streamlines through the left thalamus."
)
refresh_static_3d_report_section(paths$report_html, paths$figure_3d)
order_report_3d_sections(paths$report_html, paths$figure_3d_interactive, paths$figure_3d)
append_report_figure(
  paths$report_html,
  paths$anatomical_figure,
  "Anatomical Streamline Overlay",
  "This figure overlays sampled real trajectories on the actual MNI T1 background image. The underlay and streamlines share the same MNI world-coordinate space, so the projection is anatomically registered rather than manually aligned.",
  "Representative sampled streamlines projected onto the MNI T1 anatomy."
)
append_report_figure(
  paths$report_html,
  paths$figure,
  "Representative Streamlines",
  "This figure shows sampled actual trajectories from the full tractogram. It is qualitative visual evidence; use the FibrTraQ tables for exact quantitative counts.",
  "Representative sampled streamlines through left thalamus to selected targets."
)

cat("\nRepresentative streamline sampling complete\n")
print(summary, row.names = FALSE)
cat("\nFigure:\n", paths$figure, "\n", sep = "")
cat("Interactive 3D QC viewer:\n", paths$figure_3d_interactive, "\n", sep = "")
cat("Static 3D snapshot:\n", paths$figure_3d, "\n", sep = "")
cat("Anatomical overlay figure:\n", paths$anatomical_figure, "\n", sep = "")
cat("Streamline point CSV:\n", paths$streamlines_csv, "\n", sep = "")
cat("Sampling summary CSV:\n", paths$summary_csv, "\n", sep = "")

if (!exists("OPEN_OUTPUTS", inherits = TRUE) || isTRUE(get("OPEN_OUTPUTS", inherits = TRUE))) {
  entrypoints <- c(
    file.path(left_root, "fibrtraq-report-index-portable.html"),
    file.path(left_root, "fibrtraq-report-index.html"),
    paths$report_html,
    legacy_report,
    paths$figure_3d_interactive
  )
  open_output(entrypoints[file.exists(entrypoints)][1L])
}
