"""bigWig naming, the strand mapping, and the manifest-to-disk correspondence."""

import subprocess
from pathlib import Path

import pandas as pd
import pytest

import bigwig
import strandedness
from outputs import sha256_file
from schemas import table as get_table

from conftest import GENOME_BUILD

FIRSTSTRAND_RSEQC = (
    'This is PairEnd Data\n'
    'Fraction of reads failed to determine: 0.0300\n'
    'Fraction of reads explained by "1++,1--,2+-,2-+": 0.0400\n'
    'Fraction of reads explained by "1+-,1-+,2++,2--": 0.9300\n'
)
SECONDSTRAND_RSEQC = (
    'This is PairEnd Data\n'
    'Fraction of reads failed to determine: 0.0300\n'
    'Fraction of reads explained by "1++,1--,2+-,2-+": 0.9300\n'
    'Fraction of reads explained by "1+-,1-+,2++,2--": 0.0400\n'
)
UNSTRANDED_RSEQC = (
    'This is SingleEnd Data\n'
    'Fraction of reads failed to determine: 0.0500\n'
    'Fraction of reads explained by "++,--": 0.4800\n'
    'Fraction of reads explained by "+-,-+": 0.4700\n'
)


def test_rseqc_parsing_maps_patterns_to_library_types():
    assert strandedness.parse_rseqc(FIRSTSTRAND_RSEQC)["library_type"] == \
        strandedness.FR_FIRSTSTRAND
    assert strandedness.parse_rseqc(SECONDSTRAND_RSEQC)["library_type"] == \
        strandedness.FR_SECONDSTRAND
    assert strandedness.parse_rseqc(UNSTRANDED_RSEQC)["library_type"] == \
        strandedness.UNSTRANDED
    assert strandedness.parse_rseqc("nothing useful")["library_type"] == \
        strandedness.UNDETERMINED


def test_strand_assignment_mapping_is_explicit_and_inverts_with_chemistry():
    """str1 is the minus strand for dUTP libraries and the plus strand otherwise."""
    assert strandedness.star_strand_label(strandedness.FR_FIRSTSTRAND, "str1") == "minus"
    assert strandedness.star_strand_label(strandedness.FR_FIRSTSTRAND, "str2") == "plus"
    assert strandedness.star_strand_label(strandedness.FR_SECONDSTRAND, "str1") == "plus"
    assert strandedness.star_strand_label(strandedness.FR_SECONDSTRAND, "str2") == "minus"
    assert strandedness.star_strand_label(strandedness.UNSTRANDED, "str1") == "unstranded"


def test_strand_assignment_refuses_to_guess():
    with pytest.raises(ValueError):
        strandedness.star_strand_label(strandedness.UNDETERMINED, "str1")
    with pytest.raises(ValueError):
        # An unstranded library must not produce a second track.
        strandedness.star_strand_label(strandedness.UNSTRANDED, "str2")


def test_strand_assignment_verified_against_control_genes(results_dir, built_results):
    status = bigwig.assert_strand_assignment(
        results_dir, "human", strandedness.FR_FIRSTSTRAND, built_results["annotation"])
    assert status == "pass"


def test_strand_assignment_detects_inversion(results_dir, built_results):
    """Swapping the strand labels must be caught, not averaged away."""
    path = results_dir / get_table("signal_over_gene").relpath
    frame = pd.read_csv(path, sep="\t")
    frame["strand"] = frame["strand"].map({"plus": "minus", "minus": "plus"})
    frame.to_csv(path, sep="\t", index=False)

    with pytest.raises(bigwig.BigWigError, match="inverted"):
        bigwig.assert_strand_assignment(
            results_dir, "human", strandedness.FR_FIRSTSTRAND, built_results["annotation"])


def test_bigwig_manifest_matches_files(results_dir):
    """Every manifest row points at a real file whose size and checksum agree."""
    manifest = pd.read_csv(results_dir / get_table("bigwig_manifest").relpath, sep="\t")
    assert not manifest.empty

    for row in manifest.itertuples(index=False):
        path = results_dir / row.file_path
        assert path.is_file(), f"{row.file_path} is in the manifest but not on disk"
        assert path.stat().st_size == row.bytes
        assert sha256_file(path) == row.sha256

    on_disk = {p.relative_to(results_dir).as_posix()
               for p in (results_dir / "artifacts" / "bigwig").glob("*.bw")}
    assert on_disk == set(manifest["file_path"])


def test_bigwig_filenames_carry_strand_and_content(results_dir):
    manifest = pd.read_csv(results_dir / get_table("bigwig_manifest").relpath, sep="\t")
    for row in manifest.itertuples(index=False):
        name = row.file_path.rsplit("/", 1)[-1]
        assert name == f"{row.sample_id}.{row.content}.{GENOME_BUILD}.{row.strand}.bw"
        assert "str1" not in name and "str2" not in name


def test_bigwig_manifest_records_normalization(results_dir):
    manifest = pd.read_csv(results_dir / get_table("bigwig_manifest").relpath, sep="\t")
    assert set(manifest["normalization"]) == {"cpm"}
    assert (manifest["scale_factor"] > 0).all()


def test_cpm_scale_factor():
    assert bigwig.scale_factor(1_600_000) == pytest.approx(0.625)
    with pytest.raises(bigwig.BigWigError):
        bigwig.scale_factor(0)


def test_bed_dedupe_drops_exact_and_name_collisions(tmp_path):
    """bigWigAverageOverBed rejects duplicate names; ERCC rows in the lab BED are exact dups."""
    bed = tmp_path / "features.bed"
    bed.write_text(
        "chr1\t0\t10\tGENE_A\t0\t+\n"
        "chr1\t0\t10\tGENE_A\t0\t+\n"          # exact duplicate
        "chrERCC\t0\t100\tDQ459430\t100\t+\n"
        "chrERCC\t0\t100\tDQ459430\t100\t+\n"  # exact duplicate (the smoke-run failure)
        "chr2\t0\t5\tGENE_B\t0\t-\n"
        "chr2\t10\t20\tGENE_B\t0\t-\n"         # same name, different interval -> keep first
    )
    out = tmp_path / "dedup.bed"
    result = bigwig._bed_with_unique_names(bed, out)
    lines = [line for line in result.read_text().splitlines() if line.strip()]
    names = [line.split("\t")[3] for line in lines]
    assert names == ["GENE_A", "DQ459430", "GENE_B"]
    assert len(lines) == 3


def test_wig_rescaling_preserves_declarations(tmp_path):
    source = tmp_path / "in.wig"
    source.write_text("variableStep chrom=chr1\n1\t10\n2\t20\n")
    target = tmp_path / "out.wig"
    bigwig._rescale_wig(source, target, 0.5)
    assert target.read_text() == "variableStep chrom=chr1\n1\t5\n2\t10\n"


def test_unstranded_signal_is_rebuilt_from_bam(tmp_path, monkeypatch):
    """An unstranded library must have its str1/str2 pair rebuilt as one track.

    STAR aligns before RSeQC can measure the library type, so it always writes a stranded
    pair. Regenerating from the BAM is what keeps an unstranded sample from either crashing on
    str2 or silently losing the reads in it.
    """
    sample = "S1"
    star_dir = tmp_path / "star" / sample
    star_dir.mkdir(parents=True)
    (star_dir / f"{sample}_Aligned.sortedByCoord.out.bam").write_text("bam")
    for strand in ("str1", "str2"):
        for content in ("Unique", "UniqueMultiple"):
            (star_dir / f"{sample}_Signal.{content}.{strand}.out.wig").write_text("old\n")

    def fake_run(command, *args, **kwargs):
        assert "--outWigStrand" in command
        assert command[command.index("--outWigStrand") + 1] == "Unstranded"
        assert "inputAlignmentsFromBAM" in command
        prefix = Path(command[command.index("--outFileNamePrefix") + 1])
        for content in ("Unique", "UniqueMultiple"):
            path = prefix.parent / f"{prefix.name}Signal.{content}.str1.out.wig"
            path.write_text("rebuilt\n")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(bigwig.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(bigwig.subprocess, "run", fake_run)

    bigwig.regenerate_unstranded_wigs(tmp_path, sample)

    remaining = sorted(p.name for p in star_dir.glob("*Signal.*.out.wig"))
    assert remaining == [
        f"{sample}_Signal.Unique.str1.out.wig",
        f"{sample}_Signal.UniqueMultiple.str1.out.wig",
    ]
    for name in remaining:
        assert (star_dir / name).read_text() == "rebuilt\n"
    # Every surviving track must now be labellable for an unstranded library.
    for name in remaining:
        star_strand = bigwig.WIG_NAME_RE.search(name).group(2)
        assert strandedness.star_strand_label(
            strandedness.UNSTRANDED, star_strand) == "unstranded"


def test_unstranded_regeneration_keeps_wigs_when_star_fails(tmp_path, monkeypatch):
    """A failed rebuild must not destroy the existing signal."""
    sample = "S1"
    star_dir = tmp_path / "star" / sample
    star_dir.mkdir(parents=True)
    (star_dir / f"{sample}_Aligned.sortedByCoord.out.bam").write_text("bam")
    original = star_dir / f"{sample}_Signal.Unique.str1.out.wig"
    original.write_text("old\n")

    def failing_run(command, *args, **kwargs):
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(bigwig.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(bigwig.subprocess, "run", failing_run)

    with pytest.raises(bigwig.BigWigError, match="regenerate"):
        bigwig.regenerate_unstranded_wigs(tmp_path, sample)
    assert original.read_text() == "old\n"


def test_minus_strand_is_stored_negative():
    assert bigwig.strand_sign("minus") == -1
    assert bigwig.strand_sign("plus") == 1
    assert bigwig.strand_sign("unstranded") == 1


def test_wig_rescaling_with_negative_factor(tmp_path):
    source = tmp_path / "in.wig"
    source.write_text("variableStep chrom=chr1\n1\t10\n2\t20\n")
    target = tmp_path / "out.wig"
    bigwig._rescale_wig(source, target, -0.5)
    assert target.read_text() == "variableStep chrom=chr1\n1\t-5\n2\t-10\n"


def test_signal_over_gene_reports_minus_strand_as_positive_magnitude(tmp_path, monkeypatch):
    """A minus track stored negative must give the same signal_over_gene values as its magnitude."""
    from types import SimpleNamespace

    def fake_average(bigwig_path, bed, output):
        sign = -1 if ".minus." in bigwig_path.name else 1
        lo, hi = sorted((sign * 1.0, sign * 9.0))
        return pd.DataFrame([{"name": "GENE1", "size": 100, "covered": 50, "sum": sign * 400.0,
                              "mean0": sign * 4.0, "mean": sign * 8.0, "min": lo, "max": hi}])

    captured = {}
    monkeypatch.setattr(bigwig.shutil, "which", lambda tool: "/usr/bin/" + tool)
    monkeypatch.setattr(bigwig, "feature_bed_files", lambda ref, build: {"gene_body": tmp_path / "g.bed"})
    monkeypatch.setattr(bigwig, "_bed_with_unique_names", lambda bed, out: bed)
    monkeypatch.setattr(bigwig, "_average_over_bed", fake_average)
    monkeypatch.setattr(bigwig, "write_table", lambda d, name, frame: captured.setdefault(name, frame))
    ctx = SimpleNamespace(run_id="R1", genome_build="hg38")
    annotation = SimpleNamespace(annotate_targets=lambda names: pd.DataFrame({"gene_id": names.values}))
    rows = [{"content": "unique", "file_path": f"artifacts/bigwig/S1.unique.hg38.{s}.bw",
             "sample_id": "S1", "strand": s} for s in ("plus", "minus")]
    (tmp_path / "logs").mkdir()

    n, status = bigwig.write_signal_over_gene(ctx, tmp_path, tmp_path, annotation, rows)
    out = captured["signal_over_gene"].set_index("strand")
    assert (n, status) == (2, "ok")
    for col in ("sum_coverage", "max_coverage", "mean_coverage"):
        assert out.loc["minus", col] == out.loc["plus", col] > 0
