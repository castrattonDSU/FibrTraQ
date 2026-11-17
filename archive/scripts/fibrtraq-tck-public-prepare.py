#!/usr/bin/env python3
"""Prepare the pinned Neurodesk native-TCK interoperability demonstration.

Requires NiBabel, NumPy and an existing MRtrix mrconvert executable only for
this optional validation workflow; these are not R-package dependencies.
No registration is inferred, no streamlines are filtered, and raw inputs are
never changed. The x-partitioned tissue masks are not anatomical atlas parcels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import nibabel as nib
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def write_json(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")
    path.chmod(0o600)


def acquire(item: dict, raw_dir: Path) -> Path:
    path = raw_dir / item["filename"]
    require(path.name == item["filename"], "Source filename must be a basename")
    if path.exists():
        require(path.is_file() and not path.is_symlink(), f"Unsafe cached source: {path}")
        require(path.stat().st_size == item["bytes"] and sha256(path) == item["sha256"],
                f"Cached source differs from its pinned checksum: {path}")
        return path
    require(item["url"].startswith("https://"), "Source URLs must use HTTPS")
    fd, temporary = tempfile.mkstemp(prefix=".download-", dir=raw_dir)
    temp = Path(temporary)
    try:
        with os.fdopen(fd, "wb") as out, urllib.request.urlopen(item["url"], timeout=120) as response:
            shutil.copyfileobj(response, out)
        require(temp.stat().st_size == item["bytes"] and sha256(temp) == item["sha256"],
                f"Downloaded source failed size/SHA-256 verification: {item['filename']}")
        os.link(temp, path)
        path.chmod(0o444)
    finally:
        temp.unlink(missing_ok=True)
    return path


def validate_image(image: nib.Nifti1Image, label: str) -> dict:
    require(len(image.shape) == 3, f"{label}: expected a three-dimensional image")
    require(image.header.get_xyzt_units()[0] == "mm", f"{label}: expected declared mm units")
    require(np.isfinite(image.affine).all() and abs(np.linalg.det(image.affine[:3, :3])) > 1e-12,
            f"{label}: invalid or singular affine")
    qform, qcode = image.get_qform(coded=True)
    sform, scode = image.get_sform(coded=True)
    require(qcode > 0 or scode > 0, f"{label}: no declared qform or sform")
    disagreement = None
    if qcode > 0 and scode > 0:
        corners = np.array(np.meshgrid(*[(0, n - 1) for n in image.shape], indexing="ij")).reshape(3, -1).T
        disagreement = float(np.max(np.abs(nib.affines.apply_affine(qform, corners) -
                                           nib.affines.apply_affine(sform, corners))))
        require(disagreement <= 0.001, f"{label}: qform/sform disagree by {disagreement} mm")
    return {"shape": list(image.shape), "affine_mm": image.affine.tolist(),
            "voxel_sizes_mm": nib.affines.voxel_sizes(image.affine).tolist(),
            "qform_code": int(qcode), "sform_code": int(scode),
            "qform_sform_corner_difference_mm": disagreement}


def remove_trailing_singletons(path: Path) -> dict:
    image = nib.load(str(path))
    original_shape = image.shape
    require(len(original_shape) >= 3, f"{path.name}: fewer than three spatial dimensions")
    require(all(n == 1 for n in original_shape[3:]),
            f"{path.name}: non-singleton extra volumes require an explicit scientific selection")
    if len(original_shape) > 3:
        values = np.asarray(image.dataobj).reshape(original_shape[:3])
        # affine=None preserves the existing coded qform/sform header fields.
        # This removes storage-only singleton axes, without resampling voxels.
        normalized = nib.Nifti1Image(values, None, image.header.copy())
        nib.save(normalized, str(path))
        reloaded = nib.load(str(path))
        require(np.array_equal(np.asarray(reloaded.dataobj), values),
                f"{path.name}: voxel values changed when dropping singleton axes")
        require(np.array_equal(reloaded.affine, image.affine),
                f"{path.name}: affine changed when dropping singleton axes")
        for name in ("qform", "sform"):
            before, before_code = getattr(image, "get_" + name)(coded=True)
            after, after_code = getattr(reloaded, "get_" + name)(coded=True)
            require(before_code == after_code and
                    (before_code == 0 or np.array_equal(before, after)),
                    f"{path.name}: {name} changed when dropping singleton axes")
    return {"file": path.name, "original_shape": list(original_shape),
            "prepared_shape": list(original_shape[:3]),
            "operation": "drop trailing singleton axes only" if len(original_shape) > 3 else "none",
            "resampled": False, "spatial_transforms_changed": False}


def save_mask(mask: np.ndarray, reference: nib.Nifti1Image, path: Path) -> int:
    require(bool(mask.any()), f"Demonstration mask is empty: {path.name}")
    image = nib.Nifti1Image(mask.astype(np.uint8), reference.affine, reference.header.copy())
    image.set_data_dtype(np.uint8)
    for name in ("qform", "sform"):
        transform, code = getattr(reference, "get_" + name)(coded=True)
        getattr(image, "set_" + name)(transform, int(code))
    nib.save(image, str(path))
    reloaded = nib.load(str(path))
    require(np.array_equal(np.asarray(reloaded.dataobj), mask.astype(np.uint8)),
            f"Mask values changed on NIfTI roundtrip: {path.name}")
    require(np.array_equal(reloaded.affine, reference.affine),
            f"Mask affine changed on NIfTI roundtrip: {path.name}")
    validate_image(reloaded, path.name)
    return int(np.count_nonzero(mask))


def reader_oracle(tracks) -> dict:
    """Hash every decoded point and record length in an explicit byte order."""
    coordinate_hashes = {"md5": hashlib.md5(), "sha256": hashlib.sha256()}
    count_hashes = {"md5": hashlib.md5(), "sha256": hashlib.sha256()}
    streamlines = points = 0
    for track in tracks:
        coordinates = np.asarray(track, dtype="<f8", order="C")
        require(coordinates.ndim == 2 and coordinates.shape[1] == 3 and
                0 < len(coordinates) <= np.iinfo(np.int32).max and
                np.isfinite(coordinates).all(), "Invalid track for coordinate oracle")
        coordinate_bytes = coordinates.tobytes(order="C")
        count_bytes = np.asarray([len(coordinates)], dtype="<i4").tobytes()
        for digest in coordinate_hashes.values():
            digest.update(coordinate_bytes)
        for digest in count_hashes.values():
            digest.update(count_bytes)
        streamlines += 1
        points += len(coordinates)
    return {
        "schema_version": 1, "streamlines": streamlines, "points": points,
        "coordinate_encoding": "Float64LE; C-order xyz triples, concatenated in stored streamline/point order; no delimiters",
        "point_count_encoding": "Int32LE; one point count per stored streamline",
        "coordinates": {"bytes": points * 3 * 8,
                        **{name: digest.hexdigest() for name, digest in coordinate_hashes.items()}},
        "points_per_streamline": {"bytes": streamlines * 4,
                                  **{name: digest.hexdigest() for name, digest in count_hashes.items()}}
    }


def tractogram_check(source: Path, reference: nib.Nifti1Image, destination: Path,
                     expected_count: int) -> dict:
    native = nib.streamlines.load(str(source), lazy_load=False)
    tracks = native.streamlines
    require(len(tracks) == expected_count, "Published TCK streamline count changed")
    require(np.array_equal(native.tractogram.affine_to_rasmm, np.eye(4)),
            "NiBabel did not expose the TCK as RAS+ world millimetres")
    require(all(len(track) > 0 and np.isfinite(track).all() for track in tracks),
            "TCK contains an empty or nonfinite streamline")
    header = nib.streamlines.TrkFile.create_empty_header()
    header["dimensions"] = reference.shape
    header["voxel_sizes"] = nib.affines.voxel_sizes(reference.affine)
    header["voxel_to_rasmm"] = reference.affine
    header["voxel_order"] = "".join(nib.aff2axcodes(reference.affine)).encode("ascii")
    tractogram = nib.streamlines.Tractogram(tracks, affine_to_rasmm=np.eye(4))
    nib.streamlines.save(nib.streamlines.TrkFile(tractogram, header), str(destination))
    equivalent = nib.streamlines.load(str(destination), lazy_load=False).streamlines
    require(len(equivalent) == len(tracks), "TRK roundtrip changed streamline count")
    maximum = 0.0
    max_length_difference = 0.0
    lengths = []
    minimum = np.full(3, np.inf)
    maximum_bound = np.full(3, -np.inf)
    point_count = 0
    for index, (original, converted) in enumerate(zip(tracks, equivalent)):
        require(original.shape == converted.shape, f"TRK point-count change at streamline {index}")
        require(np.isfinite(converted).all(), f"TRK contains nonfinite points at streamline {index}")
        original64 = np.asarray(original, dtype=np.float64)
        converted64 = np.asarray(converted, dtype=np.float64)
        maximum = max(maximum, float(np.max(np.abs(original64 - converted64))))
        length = float(np.linalg.norm(np.diff(original64, axis=0), axis=1).sum())
        converted_length = float(np.linalg.norm(np.diff(converted64, axis=0), axis=1).sum())
        max_length_difference = max(max_length_difference, abs(length - converted_length))
        lengths.append(length)
        minimum = np.minimum(minimum, original64.min(axis=0))
        maximum_bound = np.maximum(maximum_bound, original64.max(axis=0))
        point_count += len(original)
    require(maximum <= 0.0001, f"TRK coordinate roundtrip exceeds 0.0001 mm: {maximum}")
    return {"streamlines": len(tracks), "points": point_count,
            "reader_oracle": reader_oracle(tracks),
            "native_world_bounds_mm": [minimum.tolist(), maximum_bound.tolist()],
            "native_length_mm": dict(zip(["min", "q25", "median", "q75", "max"],
                                          np.quantile(lengths, [0, .25, .5, .75, 1]).tolist())),
            "native_length_mean_mm": float(np.mean(lengths)),
            "trk_coordinate_max_abs_difference_mm": maximum,
            "trk_coordinate_tolerance_mm": 0.0001,
            "trk_length_max_abs_difference_mm": max_length_difference,
            "streamline_and_point_order_preserved": True}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--mrconvert", required=True, help="Existing MRtrix mrconvert executable")
    parser.add_argument("--sources", type=Path, default=Path(__file__).resolve().parent.parent /
                        "validation" / "tck-public-sources.json")
    args = parser.parse_args()
    sources_path = args.sources.resolve(strict=True)
    specification = json.loads(sources_path.read_text(encoding="utf-8"))
    root = args.output_root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    raw = root / "raw"
    raw.mkdir(exist_ok=True, mode=0o700)
    require(not raw.is_symlink(), "Raw cache must not be a symbolic link")
    prepared = root / "prepared"
    require(not prepared.exists(), "Prepared output already exists; use a new output root. Raw cache can be reused.")
    mrconvert = shutil.which(args.mrconvert)
    require(mrconvert is not None, "mrconvert was not found")
    version_run = subprocess.run([mrconvert, "-version"], check=True, capture_output=True, text=True)
    paths = {item["role"]: acquire(item, raw) for item in specification["sources"]}
    stage = Path(tempfile.mkdtemp(prefix=".prepared-", dir=root))
    commands = []
    try:
        for role, name in (("cortical_gm", "reference.nii.gz"),
                           ("gm_wm_interface", "published-interface.nii.gz")):
            command = [mrconvert, str(paths[role]), str(stage / name), "-quiet"]
            converted = subprocess.run(command, check=True, capture_output=True, text=True)
            commands.append({"argv": command, "final_output": str(prepared / name),
                             "stdout": converted.stdout, "stderr": converted.stderr})
            for message in (converted.stdout, converted.stderr):
                if message.strip():
                    print(message.strip(), file=sys.stderr)
        normalization = [remove_trailing_singletons(stage / name)
                         for name in ("reference.nii.gz", "published-interface.nii.gz")]
        reference = nib.load(str(stage / "reference.nii.gz"))
        interface = nib.load(str(stage / "published-interface.nii.gz"))
        geometry = validate_image(reference, "cortical GM")
        validate_image(interface, "GM/WM interface")
        require(reference.shape == interface.shape and np.allclose(reference.affine, interface.affine,
                                                                  atol=1e-6, rtol=0),
                "Published paired images do not share one grid")
        gm = np.asarray(reference.dataobj)
        boundary = np.asarray(interface.dataobj)
        require(np.isfinite(gm).all() and np.isfinite(boundary).all(), "Tissue volumes contain nonfinite values")
        require(float(gm.min()) >= 0 and float(gm.max()) <= 1 + 1e-6 and
                float(boundary.min()) >= 0 and float(boundary.max()) <= 1 + 1e-6,
                "Tissue volumes are outside the expected [0,1] range")
        grid = np.ogrid[tuple(slice(0, n) for n in reference.shape)]
        world_x = reference.affine[0, 3] + sum(reference.affine[0, axis] * grid[axis] for axis in range(3))
        focal = (boundary > .5) & (world_x < 0)
        masks = {"focal-left-interface.nii.gz": focal,
                 "target-right-gm.nii.gz": (gm > .5) & (world_x >= 0),
                 "target-left-gm.nii.gz": (gm > .5) & (world_x < 0) & ~focal}
        counts = {name: save_mask(mask, reference, stage / name) for name, mask in masks.items()}
        require(not np.any(masks["target-right-gm.nii.gz"] & focal), "Focal/right target overlap")
        require(not np.any(masks["target-left-gm.nii.gz"] & focal), "Focal/left target overlap")
        copied = stage / paths["native_tractogram"].name
        shutil.copyfile(paths["native_tractogram"], copied)
        require(sha256(copied) == sha256(paths["native_tractogram"]), "Native TCK copy changed bytes")
        comparison = tractogram_check(copied, reference, stage / "native-equivalent.trk",
                                      specification["expected_streamlines"])
        write_json(stage / "reader-oracle.json", {
            **comparison["reader_oracle"], "source_filename": copied.name,
            "source_sha256": sha256(copied), "nibabel_version": nib.__version__})
        provenance = {"schema_version": 1, "status": "prepared_and_independently_checked",
                      "created_utc": datetime.now(timezone.utc).isoformat(),
                      "source_specification": specification,
                      "source_manifest_sha256": sha256(sources_path),
                      "preparation_script_sha256": sha256(Path(__file__).resolve()),
                      "software": {"python": platform.python_version(), "platform": platform.platform(),
                                   "nibabel": nib.__version__, "numpy": np.__version__,
                                   "mrconvert": (version_run.stdout + version_run.stderr).strip()},
                      "invocation": sys.argv, "conversion_commands": commands,
                      "dimensional_normalization": normalization,
                      "reference_kind": "Published coregistered cortical GM probability volume, not T1 intensity",
                      "geometry": geometry, "mask_voxels": counts,
                      "independent_tractogram_check": comparison,
                      "derived_files": {p.name: {"bytes": p.stat().st_size, "sha256": sha256(p)}
                                        for p in sorted(stage.iterdir()) if p.is_file()}}
        for item in specification["sources"]:
            require(sha256(raw / item["filename"]) == item["sha256"], "Raw source changed during preparation")
        write_json(stage / "provenance.json", provenance)
        os.rename(stage, prepared)
        print(json.dumps({"prepared": str(prepared), "mask_voxels": counts,
                          "tractogram_check": comparison}, indent=2))
    except Exception:
        shutil.rmtree(stage)
        raise


if __name__ == "__main__":
    main()
