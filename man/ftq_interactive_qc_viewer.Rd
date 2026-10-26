\name{ftq_interactive_qc_viewer}
\alias{ftq_interactive_qc_viewer}
\title{Build an interactive 3D streamline QC viewer}
\usage{
ftq_interactive_qc_viewer(
  manifest,
  output_html = NULL,
  target_labels = NULL,
  sample_per_target = 25L,
  max_streamlines = 3000000L,
  max_points_per_streamline = 140L,
  brain_points = 12000L,
  target_points = 2000L,
  focal_window = 3L,
  overwrite = FALSE,
  open = FALSE,
  tractogram_context = NULL
)
}
\arguments{
\item{manifest}{Manifest path/list, or a completed object returned by \code{\link[=ftq_run]{ftq_run()}}.}
\item{output_html}{Optional output HTML path.}
\item{target_labels}{Optional subset of manifest target labels to sample.}
\item{sample_per_target}{Number of streamlines to sample per target. Supply
a single value for all targets, or a named integer vector with per-target
quotas.}
\item{max_streamlines}{Maximum streamlines to scan while sampling.}
\item{max_points_per_streamline}{Maximum points retained per rendered streamline segment.}
\item{brain_points}{Maximum anatomical context points sampled from \code{background_image} or tractogram-derived density context.}
\item{target_points}{Maximum target ROI boundary points rendered per target.}
\item{focal_window}{Number of points on each side of a focal hit to include in focal-segment mode.}
\item{overwrite}{Whether to replace an existing viewer file.}
\item{open}{Whether to open the viewer after writing.}
\item{tractogram_context}{Optional tractogram-derived context output or density image path used when no anatomical background is available.}
}
\value{
A list with viewer path and sampling summary.
}
\description{
Sample real streamlines from a FibrTraQ manifest and write a standalone HTML
viewer for spatial quality control. The viewer renders a boundary
representation of the focal ROI, a target-balanced quota subset of positive
streamlines, target ROI extents, end voxels, and context points from an
anatomical image or tractogram-derived density context. Positive streamlines
are accepted in tractogram scan order until each target quota or the scan
limit is reached; the subset is not random or prevalence weighted. Sampling
follows the manifest's primary \code{target_intersection_mode}. Endpoint-enabled
manifests also add exact versus one-voxel sensitivity controls,
endpoint-class filters, and optional focal/target hit markers without
changing the selected-method target quotas. The initial display uses an
anatomical RAS basis with superior up and an orientation compass. Endpoint
markers use continuous streamline coordinates, with colours determined by
voxel-based classifications. Circle/square shapes identify the first/last
stored point, not biological direction. Sampling and method details are
available below the canvas; the nonrandom-subset warning remains visible.
Drag to rotate freely, drag near an edge to roll, and scroll or pinch with
two fingers to zoom. Reset Orientation restores the initial view and zoom.
}
