#!/usr/bin/env python3
"""Prepare a dTOR-985 + Harvard-Oxford exploratory FibrTraQ validation set.

This script intentionally keeps Python neuroimaging dependencies outside the R
package. It fetches MNI Harvard-Oxford atlases through Nilearn, writes one
left-thalamus focal mask, writes cortical/subcortical target masks, and creates
JSON FibrTraQ manifests for pilot and full exploratory runs.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Iterable

import nibabel as nib
import numpy as np
from nilearn.datasets import fetch_atlas_harvard_oxford, load_mni152_template
from nilearn.image import resample_to_img


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare dTOR-985 Harvard-Oxford ROIs and FibrTraQ manifests."
    )
    parser.add_argument("--data-root", required=True, help="Root of FibrTraQ-data or another validation data root.")
    parser.add_argument(
        "--out-root",
        help="Output root for the prepared dTOR validation set. Defaults to <data-root>/public-validation/dtor985.",
    )
    parser.add_argument(
        "--tractogram",
        help="Path to dTOR_full_tractogram.trk. Defaults to <out-root>/raw/dTOR_full_tractogram.trk.",
    )
    parser.add_argument(
        "--atlas-data-dir",
        help="Nilearn atlas cache directory. Defaults to <out-root>/nilearn-data.",
    )
    parser.add_argument(
        "--focal-regex",
        default=r"left.*thalamus|thalamus.*left",
        help="Case-insensitive regex used to select the focal subcortical atlas label.",
    )
    parser.add_argument(
        "--pilot-targets",
        type=int,
        default=12,
        help="Number of target ROIs included in the pilot manifest.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing generated masks/manifests.")
    return parser.parse_args()


def safe_label(label: str) -> str:
    label = label.strip().replace("&", "and")
    label = re.sub(r"[^A-Za-z0-9]+", "_", label)
    label = re.sub(r"_+", "_", label).strip("_")
    return label or "label"


def load_img(path_or_img) -> nib.Nifti1Image:
    if isinstance(path_or_img, nib.spatialimages.SpatialImage):
        return path_or_img
    return nib.load(str(path_or_img))


def atlas_labels(atlas) -> list[tuple[int, str]]:
    if hasattr(atlas, "lut") and atlas.lut is not None:
        lut = atlas.lut
        if "index" in lut and "name" in lut:
            return [(int(idx), str(name)) for idx, name in zip(lut["index"], lut["name"])]
    labels = getattr(atlas, "labels", None)
    if labels is None:
        raise RuntimeError("Atlas object did not include labels or lut.")
    return [(idx, str(label)) for idx, label in enumerate(labels)]


def write_mask(atlas_img: nib.Nifti1Image, label_index: int, output_path: Path, force: bool) -> int:
    if output_path.exists() and not force:
        data = np.asanyarray(nib.load(str(output_path)).dataobj)
        return int(np.count_nonzero(data))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = np.asanyarray(atlas_img.dataobj)
    mask = (data == label_index).astype(np.uint8)
    out = nib.Nifti1Image(mask, atlas_img.affine, atlas_img.header)
    out.set_data_dtype(np.uint8)
    nib.save(out, str(output_path))
    return int(np.count_nonzero(mask))


def rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def write_json(path: Path, obj: dict, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"{path} exists. Use --force to overwrite.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def find_focal_label(labels: Iterable[tuple[int, str]], pattern: str) -> tuple[int, str]:
    rx = re.compile(pattern, flags=re.IGNORECASE)
    matches = [(idx, name) for idx, name in labels if idx != 0 and rx.search(name)]
    if len(matches) != 1:
        detail = "\n".join(f"{idx}: {name}" for idx, name in labels)
        raise RuntimeError(
            f"Expected exactly one focal label match for regex {pattern!r}, found {len(matches)}.\n"
            f"Available labels:\n{detail}"
        )
    return matches[0]


def manifest(
    *,
    out_root: Path,
    tractogram: Path,
    fnode: Path,
    targets: list[dict],
    background: Path,
    output_dir: Path,
    description: str,
    write_maps: bool,
    write_report: bool,
    overwrite: bool,
) -> dict:
    target_map = {item["label"]: rel(item["path"], out_root) for item in targets}
    return {
        "dataset_root": str(out_root),
        "subject": "sub-dtor985",
        "description": description,
        "tractogram": rel(tractogram, out_root),
        "fnode": rel(fnode, out_root),
        "local": None,
        "ref_image": rel(fnode, out_root),
        "targets": target_map,
        "coordinate": {
            "track_space": "trackvis",
            "discretization": "any",
            "axis_flip": [False, False, False],
        },
        "execution": {
            "track_reader": "stream",
            "track_chunk_size": 50000,
        },
        "outputs": {
            "derivative_name": "fibrtraq-dtor985",
            "output_dir": str(output_dir),
            "write_maps": write_maps,
            "write_report": write_report,
            "overwrite": overwrite,
        },
        "report": {
            "background_image": rel(background, out_root),
        },
    }


def main() -> None:
    args = parse_args()
    data_root = Path(args.data_root).expanduser().resolve()
    out_root = Path(args.out_root).expanduser().resolve() if args.out_root else data_root / "public-validation" / "dtor985"
    tractogram = Path(args.tractogram).expanduser().resolve() if args.tractogram else out_root / "raw" / "dTOR_full_tractogram.trk"
    atlas_data_dir = Path(args.atlas_data_dir).expanduser().resolve() if args.atlas_data_dir else out_root / "nilearn-data"

    roi_root = out_root / "rois" / "harvard-oxford-thalamus"
    focal_dir = roi_root / "focal"
    target_dir = roi_root / "targets"
    anat_dir = out_root / "anat"
    manifest_dir = out_root / "manifests"
    output_root = out_root / "outputs" / "harvard-oxford-left-thalamus"

    print("Fetching Harvard-Oxford subcortical atlas...")
    sub_atlas = fetch_atlas_harvard_oxford(
        "sub-maxprob-thr25-2mm",
        data_dir=str(atlas_data_dir),
        symmetric_split=True,
    )
    print("Fetching Harvard-Oxford lateralized cortical atlas...")
    cort_atlas = fetch_atlas_harvard_oxford(
        "cortl-maxprob-thr25-2mm",
        data_dir=str(atlas_data_dir),
    )

    sub_img = load_img(sub_atlas.maps)
    cort_img = load_img(cort_atlas.maps)
    if sub_img.shape[:3] != cort_img.shape[:3] or not np.allclose(sub_img.affine, cort_img.affine):
        raise RuntimeError("Harvard-Oxford cortical and subcortical atlases do not share the same grid.")

    sub_labels = atlas_labels(sub_atlas)
    cort_labels = atlas_labels(cort_atlas)
    focal_index, focal_name = find_focal_label(sub_labels, args.focal_regex)
    focal_safe = safe_label(focal_name)
    focal_path = focal_dir / f"{focal_safe}.nii.gz"
    focal_voxels = write_mask(sub_img, focal_index, focal_path, args.force)

    print("Preparing anatomical background on the atlas grid...")
    anat_dir.mkdir(parents=True, exist_ok=True)
    bg_path = anat_dir / "MNI152_T1_2mm_on_HarvardOxford.nii.gz"
    if not bg_path.exists() or args.force:
        template = load_mni152_template(resolution=2)
        bg = resample_to_img(template, sub_img, interpolation="continuous")
        nib.save(bg, str(bg_path))

    targets: list[dict] = []
    for source, atlas_img, labels in (
        ("cortical", cort_img, cort_labels),
        ("subcortical", sub_img, sub_labels),
    ):
        for idx, name in labels:
            if idx == 0:
                continue
            if source == "subcortical" and idx == focal_index:
                continue
            label = safe_label(name)
            path = target_dir / source / f"{label}.nii.gz"
            voxels = write_mask(atlas_img, idx, path, args.force)
            if voxels <= 0:
                continue
            targets.append(
                {
                    "source": source,
                    "atlas_index": idx,
                    "name": name,
                    "label": label,
                    "path": path,
                    "voxels": voxels,
                }
            )

    target_table = roi_root / "harvard-oxford-targets.csv"
    if target_table.exists() and not args.force:
        raise FileExistsError(f"{target_table} exists. Use --force to overwrite.")
    target_table.parent.mkdir(parents=True, exist_ok=True)
    with target_table.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["label", "name", "source", "atlas_index", "voxels", "path"],
        )
        writer.writeheader()
        for item in targets:
            writer.writerow(
                {
                    "label": item["label"],
                    "name": item["name"],
                    "source": item["source"],
                    "atlas_index": item["atlas_index"],
                    "voxels": item["voxels"],
                    "path": rel(item["path"], out_root),
                }
            )

    pilot_targets = targets[: max(1, min(args.pilot_targets, len(targets)))]
    full_manifest = manifest(
        out_root=out_root,
        tractogram=tractogram,
        fnode=focal_path,
        targets=targets,
        background=bg_path,
        output_dir=output_root / "full",
        description="dtor985-left-thalamus-harvardoxford-exploratory-full",
        write_maps=False,
        write_report=False,
        overwrite=False,
    )
    pilot_manifest = manifest(
        out_root=out_root,
        tractogram=tractogram,
        fnode=focal_path,
        targets=pilot_targets,
        background=bg_path,
        output_dir=output_root / "pilot",
        description="dtor985-left-thalamus-harvardoxford-exploratory-pilot",
        write_maps=True,
        write_report=True,
        overwrite=False,
    )
    full_manifest_path = manifest_dir / "dtor985-left-thalamus-harvardoxford-full.json"
    pilot_manifest_path = manifest_dir / "dtor985-left-thalamus-harvardoxford-pilot.json"
    write_json(full_manifest_path, full_manifest, args.force)
    write_json(pilot_manifest_path, pilot_manifest, args.force)

    readme = out_root / "README.txt"
    readme.write_text(
        "\n".join(
            [
                "dTOR-985 Harvard-Oxford exploratory FibrTraQ validation set",
                "",
                "Inputs:",
                f"- Tractogram: {tractogram}",
                f"- Focal ROI: {focal_path} ({focal_name}, {focal_voxels} voxels)",
                f"- Target ROIs: {len(targets)} Harvard-Oxford cortical/subcortical masks",
                f"- Target table: {target_table}",
                f"- Anatomical background: {bg_path}",
                "",
                "Manifests:",
                f"- Pilot: {pilot_manifest_path}",
                f"- Full: {full_manifest_path}",
                "",
                "The full manifest writes exact FibrTraQ CSV summaries without the standard target-map report.",
                "After the full run, build the many-target exploratory report with ftq_exploratory_report().",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("Prepared dTOR-985 exploratory validation inputs")
    print(f"  out_root:         {out_root}")
    print(f"  tractogram:       {tractogram}")
    print(f"  focal:            {focal_path} ({focal_voxels} voxels)")
    print(f"  targets:          {len(targets)}")
    print(f"  target table:     {target_table}")
    print(f"  pilot manifest:   {pilot_manifest_path}")
    print(f"  full manifest:    {full_manifest_path}")
    if not tractogram.exists():
        print("  warning: tractogram file does not exist yet; download it before running preflight/full analysis.")


if __name__ == "__main__":
    main()
