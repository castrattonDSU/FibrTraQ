# bench_fibrtraq.R
# Usage:
#   Rscript bench_fibrtraq.R <ROOT> <SUB> <SES> <FNODE_DESC> [DISCRETIZATION]
# Example:
#   Rscript bench_fibrtraq.R "$HOME/Projects/FibrTraQ-data" sub-001 ses-01 roi95

saferead <- function(path, gz_ok = TRUE) {
  if (file.exists(path)) return(path)
  if (gz_ok && file.exists(paste0(path, ".gz"))) return(paste0(path, ".gz"))
  alt <- sub("\\.nii\\.gz$", ".nii", path)
  if (file.exists(alt)) return(alt)
  return(NA_character_)
}

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else file.path(Sys.getenv("HOME"), "Projects", "FibrTraQ-data")
subj <- if (length(args) >= 2) args[2] else "sub-001"
sess <- if (length(args) >= 3) args[3] else "ses-01"
desc <- if (length(args) >= 4) args[4] else "roi95"
discret <- if (length(args) >= 5) args[5] else "any"
if (!discret %in% c("any", "round", "floor", "ceil")) {
  stop("DISCRETIZATION must be one of: any, round, floor, ceil")
}

roip <- file.path(root, "derivatives", "rois")
trkp <- file.path(root, "derivatives", "tractography")
outp <- file.path(root, "derivatives", "fibrtraq-v0.2.0")
dir.create(outp, showWarnings = FALSE, recursive = TRUE)

# candidate files
fnode_file <- saferead(file.path(roip, sprintf("%s_%s_label-fnode_desc-%s_space-vox_mask.nii.gz", subj, sess, desc)))
local_file <- saferead(file.path(roip, sprintf("%s_%s_label-local_space-vox_mask.nii.gz", subj, sess)), gz_ok = TRUE)

cnode_labels <- c("IPT","MPT","SPT")
cnode_files <- vapply(cnode_labels, function(L)
  saferead(file.path(roip, sprintf("%s_%s_label-%s_space-vox_mask.nii.gz", subj, sess, L))), character(1))

trk_file <- file.path(trkp, sprintf("%s_%s_space-vox_desc-wholebrain_tractogram.trk", subj, sess))

message("Inputs:")
message("  fNode:  ", fnode_file)
message("  local:  ", ifelse(is.na(local_file), "(none; will use fNode)", local_file))
for (i in seq_along(cnode_labels)) message("  cNode[", cnode_labels[i], "]: ", cnode_files[i])
message("  tracks: ", trk_file)
message("  discretization: ", discret)

# fail fast on missing files
stopifnot(!is.na(fnode_file), file.exists(fnode_file))
stopifnot(file.exists(trk_file))

suppressPackageStartupMessages({
  library(FibrTraQ)
  library(freesurferformats)
})

# read ROIs and tracks
fnode <- read.nifti1.data(fnode_file)
local <- if (!is.na(local_file) && file.exists(local_file)) read.nifti1.data(local_file) else fnode

cnodes <- list()
labels <- character(0)
for (i in seq_along(cnode_labels)) {
  if (!is.na(cnode_files[i]) && file.exists(cnode_files[i])) {
    cnodes[[ cnode_labels[i] ]] <- read.nifti1.data(cnode_files[i])
    labels <- c(labels, cnode_labels[i])
  }
}
if (!length(cnodes)) stop("No cNode ROI files found.")

# basic dimension checks
dims <- lapply(c(list(fnode=fnode, local=local), cnodes), dim)
ok_dims <- all(vapply(dims, function(x) length(x) == 3, logical(1))) &&
           all(vapply(dims, function(x) all(x == dim(fnode)), logical(1)))
if (!ok_dims) stop("ROI array dimensions are not all identical.")

trk <- read.dti.trk(trk_file)

# run
t <- system.time({
  res <- fibrtraq(fNode = fnode, localROI = local,
                  cNode_list = cnodes, cNode_labels = labels,
                  all_tracks = trk,
                  discretization = discret)
  rsc <- rsc_summary(res, percent = TRUE)
})

prefix <- sprintf("%s_%s_desc-%s", subj, sess, desc)
out_csv1 <- file.path(outp, sprintf("%s_fibrtraq_per-voxel-hits.csv", prefix))
out_csv2 <- file.path(outp, sprintf("%s_fibrtraq_rsc-summary.csv",     prefix))
write.csv(res, out_csv1, row.names = FALSE)
write.csv(rsc, out_csv2, row.names = FALSE)

cat("\nTiming (sec):\n"); print(t)
cat("\nTop RSC rows by 'tracks':\n")
print(head(rsc[order(-rsc$tracks), , drop = FALSE], 5))
cat("\nWrote:\n  ", out_csv1, "\n  ", out_csv2, "\n")
