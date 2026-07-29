#!/usr/bin/env python3
"""Run identity and provenance capture.

Nothing here is hardcoded except the mapping from a tool to the command that reports its
version. If a version cannot be determined we raise at the *start* of the run rather than
discovering it eight hours later while writing the manifest.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from schemas import PIPELINE_VERSION


class ProvenanceError(RuntimeError):
    """Raised when provenance required by the manifest contract cannot be determined."""


# tool name in the manifest -> (argv, regex capturing the version)
TOOL_PROBES: dict[str, tuple[list[str], str]] = {
    "STAR": (["STAR", "--version"], r"(\d[\w.]*)"),
    "kallisto": (["kallisto", "version"], r"version\s+(\S+)"),
    "rMATS": (["rmats.py", "--version"], r"v?(\d[\w.]*)"),
    "fastp": (["fastp", "--version"], r"fastp\s+(\S+)"),
    "fastqc": (["fastqc", "--version"], r"v?(\d[\w.]*)"),
    "RSeQC": (["infer_experiment.py", "--version"], r"(\d[\w.]*)"),
    "samtools": (["samtools", "--version"], r"samtools\s+(\S+)"),
    "bedtools": (["bedtools", "--version"], r"v?(\d[\w.]*)"),
    "wigToBigWig": (["wigToBigWig"], r"v(\d[\w.]*)"),
}

# Fallback when a tool has no usable --version flag (GSEA) or is not on PATH in a dev checkout.
CONDA_PACKAGES: dict[str, str] = {
    "STAR": "star",
    "kallisto": "kallisto",
    "rMATS": "rmats",
    "fastp": "fastp",
    "fastqc": "fastqc",
    "RSeQC": "rseqc",
    "samtools": "samtools",
    "bedtools": "bedtools",
    "GSEA": "gsea",
    "wigToBigWig": "ucsc-wigtobigwig",
    "bigWigAverageOverBed": "ucsc-bigwigaverageoverbed",
}

R_PACKAGES = {"DESeq2": "DESeq2", "sleuth": "sleuth"}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _run(argv: list[str], timeout: int = 60) -> str:
    try:
        completed = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.SubprocessError):
        return ""
    return f"{completed.stdout}\n{completed.stderr}"


def _conda_version(package: str) -> str | None:
    """Read a version out of the active conda environment's package metadata."""
    prefix = os.environ.get("CONDA_PREFIX")
    candidates = [Path(prefix)] if prefix else []
    candidates.append(Path("/opt/conda/envs/rnaseqpipeline"))
    for env_dir in candidates:
        meta_dir = env_dir / "conda-meta"
        if not meta_dir.is_dir():
            continue
        for meta in sorted(meta_dir.glob(f"{package}-*.json")):
            # Filenames are <name>-<version>-<build>.json.
            stem = meta.stem[len(package) + 1:]
            version = stem.rsplit("-", 1)[0]
            if re.match(r"^\d", version):
                return version
    return None


def _r_package_version(package: str) -> str | None:
    output = _run(["Rscript", "-e", f'cat(as.character(packageVersion("{package}")))'])
    match = re.search(r"(\d+\.\d+[\w.\-]*)", output)
    return match.group(1) if match else None


def collect_tool_versions() -> dict[str, str]:
    """Version of every tool that can influence a result, obtained by asking the tool."""
    versions: dict[str, str] = {}
    unresolved: list[str] = []

    for tool, (argv, pattern) in TOOL_PROBES.items():
        match = re.search(pattern, _run(argv))
        version = match.group(1) if match else _conda_version(CONDA_PACKAGES.get(tool, tool))
        if version:
            versions[tool] = version
        else:
            unresolved.append(tool)

    for tool, package in R_PACKAGES.items():
        version = _r_package_version(package)
        if version:
            versions[tool] = version
        else:
            unresolved.append(tool)

    # GSEA's CLI has no version flag; its conda package is the only reliable source.
    gsea_version = _conda_version("gsea")
    if gsea_version:
        versions["GSEA"] = gsea_version
    else:
        unresolved.append("GSEA")

    if unresolved:
        raise ProvenanceError(
            "could not determine version(s) for: " + ", ".join(sorted(unresolved)) + "\n"
            "Every result must be attributable to a specific tool version, so the run is "
            "stopping now rather than producing an unattributable manifest. Check that the "
            "'rnaseqpipeline' environment is active and that these tools are on PATH."
        )
    return dict(sorted(versions.items()))


def repo_root(scripts_dir: Path | None = None) -> Path:
    """Repository root: the parent of rnaseq_helper_scripts."""
    if scripts_dir is not None:
        return Path(scripts_dir).resolve().parent
    return Path(__file__).resolve().parent.parent


def git_commit(root: Path) -> tuple[str, bool]:
    """``(commit, dirty)``. Falls back to the VERSION file baked in by the Dockerfile."""
    try:
        completed = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        completed = None

    if completed is not None and completed.returncode == 0:
        commit = completed.stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            capture_output=True, text=True, timeout=60,
        )
        return commit, bool(status.stdout.strip())

    for version_file in (root / "VERSION", Path("/work/VERSION")):
        if version_file.is_file():
            commit = version_file.read_text().strip()
            if commit:
                # An image is immutable, so a commit read from VERSION is by definition clean.
                return commit, False

    raise ProvenanceError(
        f"no git repository at {root} and no VERSION file; cannot record git_commit. "
        "Build the image with the provided Dockerfile (it writes /work/VERSION) or run from "
        "a git checkout."
    )


def image_digest(root: Path, commit: str) -> tuple[str, str]:
    """``(digest, source)`` for the container image.

    A true registry digest is only knowable outside the image, so it is passed in when
    available. Otherwise we fingerprint the inputs that determine the image contents, which is
    what the digest is actually used for: detecting that two runs used different software.
    """
    env_digest = os.environ.get("RNASEQ_DOCKER_IMAGE_DIGEST", "").strip()
    if env_digest:
        return env_digest, "env:RNASEQ_DOCKER_IMAGE_DIGEST"

    for build_info in (root / "BUILD_INFO.json", Path("/work/BUILD_INFO.json")):
        if build_info.is_file():
            try:
                data = json.loads(build_info.read_text())
            except json.JSONDecodeError:
                continue
            digest = str(data.get("docker_image_digest", "")).strip()
            if digest:
                return digest, f"file:{build_info}"

    digest = hashlib.sha256()
    digest.update(commit.encode())
    for name in ("Dockerfile", "rnaseqpipeline.yml"):
        path = root / name
        if path.is_file():
            digest.update(path.read_bytes())
    return f"sha256:{digest.hexdigest()}", "computed:dockerfile+env+commit"


def reference_dir_sha256(reference_dir: Path) -> str:
    """Cheap fingerprint of the reference tree that detects a swapped reference.

    Hashes the sorted ``(relative_path, size, mtime)`` triples rather than file contents, so it
    costs a directory walk instead of reading a 30 GB STAR index.
    """
    reference_dir = Path(reference_dir)
    if not reference_dir.is_dir():
        raise ProvenanceError(f"reference_dir does not exist: {reference_dir}")
    entries = []
    for path in sorted(reference_dir.rglob("*")):
        if path.is_file():
            stat = path.stat()
            rel = path.relative_to(reference_dir).as_posix()
            entries.append(f"{rel}\t{stat.st_size}\t{int(stat.st_mtime)}")
    digest = hashlib.sha256("\n".join(entries).encode())
    return digest.hexdigest()


def annotation_version(reference_dir: Path, genome_build: str) -> str:
    """Annotation release backing ``genome_build``.

    Read from ``annotations/ANNOTATION_VERSION`` if present, otherwise recovered from the GTF
    header. Custom builds with viral contigs appended usually lose the header, hence the file.
    """
    annotations = Path(reference_dir) / genome_build / "annotations"
    marker = annotations / "ANNOTATION_VERSION"
    if marker.is_file():
        value = marker.read_text().strip()
        if value:
            return value

    gtf = annotations / f"{genome_build}.gtf"
    if gtf.is_file():
        with open(gtf, "r", errors="replace") as handle:
            for _ in range(50):
                line = handle.readline()
                if not line or not line.startswith("#"):
                    break
                gencode = re.search(r"version\s+(\d+)", line)
                if "gencode" in line.lower() and gencode:
                    return f"gencode_v{gencode.group(1)}"
                ensembl = re.search(r"genebuild-last-updated\s+(\S+)", line)
                if ensembl:
                    return f"ensembl_{ensembl.group(1)}"

    raise ProvenanceError(
        f"cannot determine the annotation version for {genome_build}. Record it once with:\n"
        f"    echo 'gencode_v44' > {marker}\n"
        f"or pass --annotation-version to rnaseq.py."
    )


def make_run_id(experiment_id: str, start_utc: str, commit: str) -> str:
    """``R`` + 12 hex chars. Deterministic, collision-safe, no central counter."""
    digest = hashlib.sha256(f"{experiment_id}{start_utc}{commit}".encode()).hexdigest()
    return f"R{digest[:12]}"


class Manifest:
    """Accumulates provenance during the run and writes ``run_manifest.json`` at the end."""

    def __init__(self, results_dir: Path, data: dict):
        self.results_dir = Path(results_dir)
        self.data = data

    @property
    def path(self) -> Path:
        return self.results_dir / "run_manifest.json"

    @property
    def run_id(self) -> str:
        return self.data["run_id"]

    def set(self, key: str, value) -> None:
        self.data[key] = value

    def record_stage(self, stage: str, status: str) -> None:
        self.data.setdefault("stage_status", {})[stage] = status

    def write(self, exit_status: str, failed_stage: str | None = None) -> Path:
        self.data["run_end_utc"] = utc_now()
        self.data["exit_status"] = exit_status
        if failed_stage:
            self.data["failed_stage"] = failed_stage
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.data, indent=2, sort_keys=False) + "\n")
        return self.path


# Documented in docs/SCHEMA.md, which is generated from this list. start_run() asserts that the
# manifest it builds has exactly these keys, so a new field cannot be added without documenting it.
MANIFEST_FIELDS = (
    ("run_id", "`R` + 12 hex characters, derived from experiment, start time and commit. "
               "Joins to the `run_id` column of every table."),
    ("experiment_id", "Input directory name, `Model_Experiment`."),
    ("pipeline_version", "Version of this pipeline's output contract."),
    ("git_commit", "Commit the code was run from."),
    ("git_dirty", "`true` if the working tree had uncommitted changes. Results from a dirty "
                  "tree are not reproducible from the commit alone."),
    ("docker_image_digest", "Image the run executed in."),
    ("docker_image_digest_source", "How the digest was obtained: the runner's environment, "
                                   "`/work/VERSION`, or the local git commit for a bare run."),
    ("genome_build", "Reference build name, matching a directory under `reference_dir`."),
    ("annotation_version", "GTF version, e.g. `gencode_v44`. Read from an "
                           "`ANNOTATION_VERSION` marker, the GTF header, or `--annotation-version`."),
    ("reference_dir", "Path to the reference directory as mounted."),
    ("reference_dir_sha256", "Digest over the reference file inventory, so two runs against "
                             "different references are distinguishable."),
    ("investigator", "Who ran it."),
    ("library_layout", "`PE` or `SE`."),
    ("run_start_utc", "ISO 8601, UTC."),
    ("run_end_utc", "ISO 8601, UTC. `null` while the run is in progress."),
    ("exit_status", "`running`, `success` or `failed`."),
    ("tool_versions", "Object mapping tool name to version string, captured at runtime rather "
                      "than hardcoded."),
    ("parameters", "Object recording the non-default arguments passed to each tool."),
)

# Added as the run progresses; absent from a manifest written at startup.
OPTIONAL_MANIFEST_FIELDS = (
    ("stage_status", "Object mapping stage name to `ok`, `skipped` or `partial`."),
    ("failed_stage", "Stage that raised, present only when `exit_status` is `failed`."),
    ("error", "Exception type and message, present only when `exit_status` is `failed`."),
    ("validation", "Result of `validate_outputs.py`: `status`, `failures` and `warnings`."),
)


def start_run(*, results_dir: Path, reference_dir: Path, scripts_dir: Path, experiment_id: str,
              genome_build: str, investigator: str, library_layout: str,
              parameters: dict, annotation_version_override: str | None = None) -> Manifest:
    """Collect everything knowable at the start of a run and return the open manifest."""
    root = repo_root(scripts_dir)
    start_utc = utc_now()
    commit, dirty = git_commit(root)
    digest, digest_source = image_digest(root, commit)

    data = {
        "run_id": make_run_id(experiment_id, start_utc, commit),
        "experiment_id": experiment_id,
        "pipeline_version": PIPELINE_VERSION,
        "git_commit": commit,
        "git_dirty": dirty,
        "docker_image_digest": digest,
        "docker_image_digest_source": digest_source,
        "genome_build": genome_build,
        "annotation_version": annotation_version_override
                              or annotation_version(reference_dir, genome_build),
        "reference_dir": str(reference_dir),
        "reference_dir_sha256": reference_dir_sha256(reference_dir),
        "investigator": investigator,
        "library_layout": library_layout,
        "run_start_utc": start_utc,
        "run_end_utc": None,
        "exit_status": "running",
        "tool_versions": collect_tool_versions(),
        "parameters": parameters,
    }

    documented = [name for name, _ in MANIFEST_FIELDS]
    if list(data) != documented:
        raise ProvenanceError(
            "the manifest and MANIFEST_FIELDS disagree; add the field to MANIFEST_FIELDS and "
            f"regenerate docs/SCHEMA.md.\n  built: {list(data)}\n  documented: {documented}"
        )
    return Manifest(results_dir, data)
