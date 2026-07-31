"""bigWig naming, the strand mapping, and the manifest-to-disk correspondence."""

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
