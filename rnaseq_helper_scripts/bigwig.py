#!/usr/bin/env python3
"""Tasks 12 and 13: deterministic bigWig naming, uniform normalization, gene-level signal.

Two things the legacy archive got wrong are fixed here.

*Strand*: STAR's ``str1``/``str2`` means opposite genomic strands depending on library
chemistry, so it is resolved through ``strandedness.star_strand_label`` and written into the
filename and the manifest. A library RSeQC could not classify arrives here as unstranded
(strandedness.analysed_library_type) and gets a single unstranded track.

*Normalization*: STAR is run with ``--outWigNorm None``, so raw tracks are not comparable
between samples. Every track is scaled to CPM and the factor actually applied is recorded per
file. Minus-strand tracks hold negative values (the lab's genome-browser convention); plus and
unstranded tracks are positive. Strand is also carried by the filename and the manifest.
signal_over_gene flips the minus-strand sign back, so its coverage values are positive on both
strands -- compare coverage between tracks by absolute value.

*Layout*: tracks are split by content into ``artifacts/bigwig/unique/`` and
``artifacts/bigwig/uniquemulti/`` so a unique track cannot be picked up by mistake for a
unique+multi one (the legacy pipeline's 1_unique / 2_unique_multiple split).

*Raw counts*: alongside each CPM track the pipeline writes the same signal unscaled -- STAR's
own read counts, minus strand negative -- to ``raw_count_bigwigs/<content>/*.raw.bw`` with
``raw_count_bigwigs/raw_count_manifest.tsv``. That folder sits outside ``tables/`` and
``artifacts/`` on purpose: the update kit uploads it to the lab share's
``raw_count_bigwigs/<experiment>/`` rather than into ``pipeline_output/``. Unique tracks hold
whole reads; unique+multi tracks hold fractions, because STAR splits a read that maps to N loci
as 1/N per locus.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

import strandedness
import vocab
from annotation import Annotation
from outputs import RunContext, sha256_file, write_table
from schemas import BIGWIG_MANIFEST, NA, SIGNAL_OVER_GENE

CONTENT_BY_STAR_LABEL = {"Unique": "unique", "UniqueMultiple": "uniquemulti"}
WIG_NAME_RE = re.compile(r"Signal\.(Unique|UniqueMultiple)\.(str\d)\.out\.wig$")
DEFAULT_NORMALIZATION = "cpm"
# signal_over_gene aggregates the uniquely-mapping tracks; multi-mapping signal is ambiguous per
# locus and would double-count repeats. The uniquemulti tracks are still in the manifest.
SIGNAL_CONTENT = "unique"
DEFAULT_FEATURE_SET = "gene_body"
RAW_COUNT_DIR = "raw_count_bigwigs"
RAW_COUNT_MANIFEST = "raw_count_manifest.tsv"
# Same columns as the manifests of the raw copies made from the 70 experiments converted in
# Oct 2026. Those were derived by dividing CPM back out and rounding, hence the rounding
# columns; tracks written here come straight from STAR, so both are always 0.
RAW_COUNT_COLUMNS = ["run_id", "sample_id", "file_path", "source_file", "content", "strand",
                     "genome_build", "scale_factor", "values", "max_rounding_residual",
                     "ambiguous_positions", "sha256", "bytes"]


class BigWigError(RuntimeError):
    pass


def star_library_sizes(results_dir: Path, sample_name: str) -> dict[str, int]:
    """Uniquely-mapped and multi-mapped read counts from STAR's final log."""
    path = Path(results_dir) / "star" / sample_name / f"{sample_name}_Log.final.out"
    if not path.is_file():
        raise BigWigError(f"STAR log missing for {sample_name}: {path}")
    text = path.read_text()

    def value(label: str) -> int:
        match = re.search(re.escape(label) + r"\s*\|\s*(\d+)", text)
        if not match:
            raise BigWigError(f"{path} has no '{label}' line; cannot compute a scale factor.")
        return int(match.group(1))

    unique = value("Uniquely mapped reads number")
    multi = value("Number of reads mapped to multiple loci")
    return {"unique": unique, "uniquemulti": unique + multi}


def scale_factor(library_size: int, normalization: str = DEFAULT_NORMALIZATION) -> float:
    if normalization != "cpm":
        raise BigWigError(f"unsupported normalization {normalization!r}")
    if library_size <= 0:
        raise BigWigError("library size is zero; cannot normalize to CPM")
    return 1e6 / library_size


def strand_sign(strand: str) -> int:
    """-1 for minus-strand tracks (stored negative), +1 for plus and unstranded."""
    return -1 if strand == "minus" else 1


def _rescale_wig(source: Path, target: Path, factor: float,
                 raw_target: Path | None = None, raw_sign: int = 1) -> None:
    """Write ``target`` with every value multiplied by ``factor``. With ``raw_target``, also
    write the unscaled values (only ``raw_sign`` applied) in the same pass."""
    with open(source, "r") as handle_in, open(target, "w") as handle_out, \
            (open(raw_target, "w") if raw_target else open(os.devnull, "w")) as handle_raw:
        write, write_raw = handle_out.write, handle_raw.write
        for line in handle_in:
            if line[:1].isdigit():
                position, _, value = line.partition("\t")
                number = float(value)
                write(f"{position}\t{number * factor:.6g}\n")
                write_raw(f"{position}\t{number * raw_sign:.9g}\n")
            else:
                # variableStep / fixedStep / track declarations pass through untouched.
                write(line)
                write_raw(line)


def _has_stranded_wigs(star_dir: Path) -> bool:
    """True if STAR left a per-strand signal pair, i.e. anything beyond ``str1``."""
    for wig in star_dir.glob("*Signal.*.out.wig"):
        match = WIG_NAME_RE.search(wig.name)
        if match and match.group(2) != "str1":
            return True
    return False


def regenerate_unstranded_wigs(results_dir: Path, sample_name: str) -> None:
    """Rebuild one sample's signal as a single unstranded track, straight from its BAM.

    STAR's wiggle strandedness is fixed at alignment time, but the library type is only known
    afterwards -- RSeQC needs the BAM to measure it. So STAR always aligns with the default
    ``--outWigStrand Stranded`` and emits ``str1``/``str2``. For a genuinely unstranded library
    that pair is an arbitrary split of one signal: neither track is a genomic strand, and
    dropping either would silently discard half the reads.

    ``--runMode inputAlignmentsFromBAM`` re-derives the signal from the finished BAM with
    ``--outWigStrand Unstranded``, producing one ``str1`` track per content type that holds all
    the reads. No realignment happens, so this costs one pass over the BAM.
    """
    star_dir = Path(results_dir) / "star" / sample_name
    bam = star_dir / f"{sample_name}_Aligned.sortedByCoord.out.bam"
    if not bam.is_file():
        raise BigWigError(
            f"cannot regenerate unstranded signal for {sample_name}: BAM missing at {bam}"
        )
    if not shutil.which("STAR"):
        raise BigWigError("STAR is not on PATH; cannot regenerate unstranded signal")

    # Build into a scratch directory so a failure leaves the existing wiggles untouched.
    with tempfile.TemporaryDirectory(dir=star_dir) as scratch:
        prefix = Path(scratch) / "regen_"
        command = [
            "STAR",
            "--runMode", "inputAlignmentsFromBAM",
            "--inputBAMfile", str(bam),
            "--outWigType", "wiggle",
            "--outWigStrand", "Unstranded",
            "--outWigNorm", "None",
            "--outFileNamePrefix", str(prefix),
        ]
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as error:
            raise BigWigError(
                f"STAR failed to regenerate unstranded signal for {sample_name}: {error}"
            ) from error

        rebuilt = sorted(Path(scratch).glob("*Signal.*.out.wig"))
        if not rebuilt:
            raise BigWigError(
                f"STAR produced no signal wiggles for {sample_name} in unstranded mode"
            )

        # Only now discard the stranded pair, then move the unstranded tracks into place.
        for stale in star_dir.glob("*Signal.*.out.wig"):
            stale.unlink()
        for wig in rebuilt:
            match = WIG_NAME_RE.search(wig.name)
            if not match:
                raise BigWigError(f"unexpected STAR signal filename: {wig.name}")
            shutil.move(str(wig), str(star_dir / f"{sample_name}_{match.group(0)}"))


def convert_sample_bigwigs(ctx: RunContext, results_dir: Path, sample_name: str,
                           library_type: str, chr_sizes: Path,
                           normalization: str = DEFAULT_NORMALIZATION) -> list[dict]:
    """Convert one sample's STAR wiggles to normalized, canonically named bigWigs, plus an
    unscaled raw-count bigWig for each. Raw-file details ride along in each row under
    ``raw_*`` keys for :func:`write_raw_count_manifest`; the bigWig manifest ignores them."""
    results_dir = Path(results_dir)
    star_dir = results_dir / "star" / sample_name
    bigwig_dir = results_dir / "artifacts" / "bigwig"
    raw_dir = results_dir / RAW_COUNT_DIR

    if not shutil.which("wigToBigWig"):
        raise BigWigError("wigToBigWig is not on PATH")

    # An unstranded library cannot be described by STAR's str1/str2 pair, so rebuild the signal
    # as a single track before anything tries to give those two files a genomic strand.
    if library_type == strandedness.UNSTRANDED and _has_stranded_wigs(star_dir):
        regenerate_unstranded_wigs(results_dir, sample_name)

    sample_id = ctx.sample_id(sample_name)
    sizes = star_library_sizes(results_dir, sample_name)
    rows: list[dict] = []

    for wig in sorted(star_dir.glob("*Signal.*.out.wig")):
        match = WIG_NAME_RE.search(wig.name)
        if not match:
            raise BigWigError(f"unexpected STAR signal filename: {wig.name}")
        star_content, star_strand = match.group(1), match.group(2)
        content = CONTENT_BY_STAR_LABEL[star_content]
        strand = strandedness.star_strand_label(library_type, star_strand)

        factor = scale_factor(sizes[content], normalization)
        scaled_wig = wig.with_suffix(".normalized.wig")
        raw_wig = wig.with_suffix(".raw.wig")
        # scale_factor in the manifest stays the (positive) CPM factor; the sign is the strand.
        _rescale_wig(wig, scaled_wig, factor * strand_sign(strand),
                     raw_target=raw_wig, raw_sign=strand_sign(strand))

        name = f"{sample_id}.{content}.{ctx.genome_build}.{strand}"
        target = bigwig_dir / content / f"{name}.bw"
        raw_target = raw_dir / content / f"{name}.raw.bw"
        target.parent.mkdir(parents=True, exist_ok=True)
        raw_target.parent.mkdir(parents=True, exist_ok=True)
        try:
            for source, output in ((scaled_wig, target), (raw_wig, raw_target)):
                command = ["wigToBigWig", str(source), str(chr_sizes), str(output)]
                try:
                    subprocess.run(command, check=True)
                except subprocess.CalledProcessError as error:
                    raise BigWigError(f"wigToBigWig failed for {wig.name}: {error}") from error
        finally:
            scaled_wig.unlink(missing_ok=True)
            raw_wig.unlink(missing_ok=True)

        rows.append({
            "run_id": ctx.run_id,
            "sample_id": sample_id,
            "file_path": target.relative_to(results_dir).as_posix(),
            "content": content,
            "strand": strand,
            "genome_build": ctx.genome_build,
            "normalization": normalization,
            "scale_factor": factor,
            "sha256": sha256_file(target),
            "bytes": target.stat().st_size,
            "raw_file_path": raw_target.relative_to(raw_dir).as_posix(),
            "raw_sha256": sha256_file(raw_target),
            "raw_bytes": raw_target.stat().st_size,
        })
    return rows


def write_bigwig_manifest(results_dir: Path, rows: list[dict]) -> Path:
    frame = pd.DataFrame(rows, columns=BIGWIG_MANIFEST.column_names)
    return write_table(results_dir, "bigwig_manifest", frame)


def write_raw_count_manifest(ctx: RunContext, results_dir: Path, rows: list[dict]) -> Path:
    """``raw_count_bigwigs/raw_count_manifest.tsv``: one row per raw-count bigWig.

    ``file_path`` is relative to the raw-count folder (it becomes
    ``raw_count_bigwigs/<experiment>/`` on the lab share); ``source_file`` names the CPM track
    it pairs with, as it is laid out under the share's ``pipeline_output/``.
    """
    records = [{
        "run_id": row["run_id"],
        "sample_id": row["sample_id"],
        "file_path": row["raw_file_path"],
        "source_file": f"pipeline_output/{ctx.experiment_id}/{row['file_path']}",
        "content": row["content"],
        "strand": row["strand"],
        "genome_build": row["genome_build"],
        "scale_factor": row["scale_factor"],
        "values": "integer_reads" if row["content"] == "unique" else "fractional_reads",
        "max_rounding_residual": 0,
        "ambiguous_positions": 0,
        "sha256": row["raw_sha256"],
        "bytes": row["raw_bytes"],
    } for row in rows]
    target = Path(results_dir) / RAW_COUNT_DIR / RAW_COUNT_MANIFEST
    target.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records, columns=RAW_COUNT_COLUMNS).to_csv(
        target, sep="\t", index=False, lineterminator="\n")
    return target


def feature_bed_files(reference_dir: Path, genome_build: str) -> dict[str, Path]:
    """Feature sets to aggregate over, keyed by ``feature_set``.

    Any BED dropped into ``annotations/features/`` becomes a new feature set without a schema
    change, which is how ``promoter_tss_500``, ``exon`` and ``dog_10kb`` get added later.
    """
    annotations = Path(reference_dir) / genome_build / "annotations"
    features_dir = annotations / "features"
    if features_dir.is_dir():
        beds = {bed.stem: bed for bed in sorted(features_dir.glob("*.bed"))}
        if beds:
            return beds
    default = annotations / f"{genome_build}.bed"
    return {DEFAULT_FEATURE_SET: default} if default.is_file() else {}


AVERAGE_OVER_BED_COLUMNS = ["name", "size", "covered", "sum", "mean0", "mean", "min", "max"]


def _bed_with_unique_names(bed: Path, scratch: Path) -> Path:
    """Write a BED that ``bigWigAverageOverBed`` will accept.

    UCSC's tool rejects duplicate names in column 4. The combined
    ``hg38plusAkataInverted.bed`` has exact-duplicate ERCC rows (same interval, same
    name twice), which aborted every ``signal_over_gene`` call. Collapse identical
    lines and, if a name still appears more than once on distinct intervals, keep the
    first and drop the rest — gene-level rollup does not need every isoform copy when
    the name collides.
    """
    seen_lines: set[str] = set()
    seen_names: set[str] = set()
    kept = 0
    dropped = 0
    with open(bed, "r") as handle_in, open(scratch, "w") as handle_out:
        for line in handle_in:
            if not line.strip() or line.startswith(("track", "browser", "#")):
                handle_out.write(line)
                continue
            if line in seen_lines:
                dropped += 1
                continue
            seen_lines.add(line)
            parts = line.rstrip("\n").split("\t")
            name = parts[3] if len(parts) > 3 else ""
            if name in seen_names:
                dropped += 1
                continue
            seen_names.add(name)
            handle_out.write(line)
            kept += 1
    if dropped:
        print(f"WARNING: dropped {dropped} duplicate BED row(s) from {bed.name} "
              f"({kept} unique names kept) so bigWigAverageOverBed can run.")
    return scratch


def _average_over_bed(bigwig: Path, bed: Path, output: Path) -> pd.DataFrame:
    command = ["bigWigAverageOverBed", "-minMax", str(bigwig), str(bed), str(output)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return pd.read_csv(output, sep="\t", header=None, names=AVERAGE_OVER_BED_COLUMNS)


def write_signal_over_gene(ctx: RunContext, results_dir: Path, reference_dir: Path,
                           annotation: Annotation,
                           manifest_rows: list[dict]) -> tuple[int, str]:
    """Task 13. Returns ``(row_count, status)``."""
    results_dir = Path(results_dir)
    beds = feature_bed_files(reference_dir, ctx.genome_build)
    tool = shutil.which("bigWigAverageOverBed")

    if not tool or not beds:
        reason = ("bigWigAverageOverBed is not on PATH" if not tool
                  else f"no feature BED found for {ctx.genome_build}")
        print(f"WARNING: {reason}; writing signal_over_gene.tsv with headers only.")
        write_table(results_dir, "signal_over_gene",
                    pd.DataFrame(columns=SIGNAL_OVER_GENE.column_names))
        return 0, "skipped"

    scratch = results_dir / "logs" / "signal_over_gene.tmp"
    # Deduplicate each feature BED once; bigWigAverageOverBed is then called per track.
    unique_beds: dict[str, Path] = {}
    for feature_set, bed in beds.items():
        deduped = results_dir / "logs" / f"signal_over_gene.{feature_set}.dedup.bed"
        unique_beds[feature_set] = _bed_with_unique_names(bed, deduped)

    frames = []
    try:
        for row in manifest_rows:
            if row["content"] != SIGNAL_CONTENT:
                continue
            bigwig_path = results_dir / row["file_path"]
            for feature_set, bed in unique_beds.items():
                try:
                    raw = _average_over_bed(bigwig_path, bed, scratch)
                except (subprocess.CalledProcessError, OSError) as error:
                    detail = ""
                    if isinstance(error, subprocess.CalledProcessError) and error.stderr:
                        detail = f" ({error.stderr.strip()[:200]})"
                    print(f"WARNING: bigWigAverageOverBed failed for {bigwig_path.name} / "
                          f"{feature_set}: {error}{detail}")
                    continue
                # Minus-strand tracks are stored negative; report coverage as a positive
                # magnitude on every strand (the most negative value is the peak).
                if strand_sign(row["strand"]) < 0:
                    raw = raw.assign(sum=-raw["sum"], max=-raw["min"])
                # BED12 annotations are keyed by transcript; roll them up so the table is keyed
                # by gene like every other table.
                gene_ids = annotation.annotate_targets(raw["name"].astype(str))["gene_id"]
                raw = raw.assign(
                    gene_id=gene_ids.where(gene_ids != NA, raw["name"].astype(str)))
                grouped = raw.groupby("gene_id", as_index=False).agg(
                    sum_coverage=("sum", "sum"),
                    covered_bases=("covered", "sum"),
                    feature_length=("size", "sum"),
                    max_coverage=("max", "max"),
                )
                grouped["mean_coverage"] = (
                    grouped["sum_coverage"] / grouped["feature_length"].where(
                        grouped["feature_length"] > 0)
                )
                grouped["run_id"] = ctx.run_id
                grouped["sample_id"] = row["sample_id"]
                grouped["strand"] = row["strand"]
                grouped["feature_set"] = feature_set
                frames.append(grouped[SIGNAL_OVER_GENE.column_names])
    finally:
        scratch.unlink(missing_ok=True)
        for bed in unique_beds.values():
            bed.unlink(missing_ok=True)

    if not frames:
        write_table(results_dir, "signal_over_gene",
                    pd.DataFrame(columns=SIGNAL_OVER_GENE.column_names))
        return 0, "skipped"

    combined = pd.concat(frames, ignore_index=True)
    write_table(results_dir, "signal_over_gene", combined)
    return len(combined), "ok"


def assert_strand_assignment(results_dir: Path, organism: str, library_type: str,
                             annotation: Annotation) -> str:
    """Check the str1/str2 -> strand mapping against known-strand control genes.

    A silently inverted mapping is the single failure mode that made the legacy bigWig archive
    unusable, so it is checked rather than trusted.
    """
    if library_type in (strandedness.UNSTRANDED, strandedness.UNDETERMINED):
        return "not_applicable"

    path = Path(results_dir) / SIGNAL_OVER_GENE.relpath
    controls = vocab.STRAND_CONTROL_GENES.get(organism)
    if not path.is_file() or controls is None:
        return "skipped"

    signal = pd.read_csv(path, sep="\t")
    if signal.empty:
        return "skipped"

    checked = 0
    for expected_strand, symbols in controls.items():
        other = "minus" if expected_strand == "plus" else "plus"
        for symbol in symbols:
            gene_id = annotation.gene_id_for_symbol(symbol)
            rows = signal[signal["gene_id"].isin({gene_id, symbol})]
            on = rows.loc[rows["strand"] == expected_strand, "sum_coverage"].sum()
            off = rows.loc[rows["strand"] == other, "sum_coverage"].sum()
            if on == 0 and off == 0:
                continue
            checked += 1
            if off > on:
                raise BigWigError(
                    f"bigWig strand assignment looks inverted: {symbol} is on the "
                    f"{expected_strand} strand but carries more signal on the {other} track "
                    f"({off:.1f} vs {on:.1f}). Library type was {library_type}; check "
                    f"strandedness.STAR_STRAND_LABELS before trusting these tracks."
                )
    if not checked:
        print("WARNING: no strand control genes found in signal_over_gene.tsv; strand "
              "assignment was not verified.")
        return "skipped"
    return "pass"
