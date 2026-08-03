"""Fixtures that build a complete, contract-satisfying results directory.

The tables are produced by the *real* writers reading realistic tool output, so these tests
cover the parsing and schema code rather than a hand-written copy of it. Tool output is
synthesized because STAR, kallisto and GSEA are not available on a CI runner; the end-to-end
run against real FASTQs is described in tests/fixtures/README.md.
"""

import json
import sys
import zipfile
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "rnaseq_helper_scripts"))

import checksums  # noqa: E402
import tables_de  # noqa: E402
import tables_expression  # noqa: E402
import tables_gsea  # noqa: E402
import tables_junction  # noqa: E402
import tables_qc  # noqa: E402
import tables_rmats  # noqa: E402
from annotation import Annotation  # noqa: E402
from outputs import RunContext, create_layout, sha256_file, write_table  # noqa: E402
from schemas import BIGWIG_MANIFEST, METADATA, METADATA_LEGACY_COLUMNS, SIGNAL_OVER_GENE  # noqa: E402

GENES = [
    # (transcript_id, symbol, ensembl gene id, strand)
    ("ENST00000000001.1", "GAPDH", "ENSG00000111640.15", "plus"),
    ("ENST00000000002.3", "ACTB", "ENSG00000075624.17", "minus"),
    ("ENST00000000003.2", "MYC", "ENSG00000136997.18", "plus"),
    ("ENST00000000004.1", "RPMS1_1", "RPMS1", "plus"),
]
SAMPLES = [("cntl_a", "cntl", 1), ("test_a", "test", 1)]
EXPERIMENT_ID = "SNU719_Zta-plus-Rta"
GENOME_BUILD = "hg38plusAkataInverted"
RUN_ID = "Rabcdef123456"


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _tsv(path: Path, header: list[str], rows: list[list]) -> Path:
    lines = ["\t".join(header)]
    lines += ["\t".join(str(value) for value in row) for row in rows]
    return _write(path, "\n".join(lines) + "\n")


def write_reference(reference_dir: Path) -> Path:
    annotations = reference_dir / GENOME_BUILD / "annotations"
    _tsv(reference_dir / GENOME_BUILD / "biomart" / f"{GENOME_BUILD}.mart_export.txt",
         ["target_id", "gene", "gene_id"],
         [[transcript, symbol, gene_id] for transcript, symbol, gene_id, _ in GENES])
    _write(annotations / "ANNOTATION_VERSION", "gencode_v44\n")
    return reference_dir


def write_metadata(results_dir: Path) -> pd.DataFrame:
    rows = []
    for sample_name, condition, replicate in SAMPLES:
        sample_id = f"{EXPERIMENT_ID}_{condition}{replicate}"
        rows.append({
            "sample_id": sample_id, "experiment_id": EXPERIMENT_ID, "condition": condition,
            "replicate_index": replicate, "cell_line": "SNU719", "organism": "human",
            "genome_build": GENOME_BUILD,
            "perturbation_type": "transfection" if condition == "test" else "none",
            "perturbation_agent": "Zta" if condition == "test" else "NA",
            "perturbation_target": "NA",
            "perturbation_dose": "NA",
            "co_treatment": "none", "co_treatment_target": "NA",
            "facs_purified": "yes", "facs_gfp_promoter": "pCMV",
            "timepoint_hours": "24.0", "library_layout": "PE",
            "library_selection": "polyA", "library_strandedness": "stranded",
            "strandedness": "NA", "fastq_r1": f"/data/{sample_name}_1.fq.gz",
            "fastq_r2": f"/data/{sample_name}_2.fq.gz", "fastq_r1_md5": "0" * 32,
            "fastq_r2_md5": "1" * 32, "investigator": "ethan",
            "sequencing_run_date": "2026-01-15", "notes": "NA",
            "Path Read 1": f"/data/{sample_name}_1.fq.gz",
            "Path Read 2": f"/data/{sample_name}_2.fq.gz",
            "Species": GENOME_BUILD, "Sample name": sample_name,
            "Condition": "test" if condition == "test" else "control",
            "Control?": "no" if condition == "test" else "yes",
            "ConditionReplicate": f"{condition}{replicate}",
            "Experiment name": EXPERIMENT_ID,
        })
    frame = pd.DataFrame(rows)[list(METADATA.column_names) + list(METADATA_LEGACY_COLUMNS)]
    frame.to_csv(results_dir / "metadata.tsv", sep="\t", index=False)
    return frame


def write_tool_output(results_dir: Path) -> None:
    for index, (sample_name, condition, _) in enumerate(SAMPLES):
        # kallisto
        kallisto_dir = results_dir / "kallisto" / sample_name / f"{sample_name}_{GENOME_BUILD}"
        tpms = [400000.0, 300000.0, 200000.0, 100000.0]
        _tsv(kallisto_dir / "abundance.tsv",
             ["target_id", "length", "eff_length", "est_counts", "tpm"],
             [[transcript, 1000, 800.5, 100.0 * (position + 1), tpms[position]]
              for position, (transcript, *_rest) in enumerate(GENES)])
        _write(kallisto_dir / "run_info.json", json.dumps(
            {"n_processed": 1000000, "n_pseudoaligned": 900000, "p_pseudoaligned": 90.0}))

        # fastp
        _write(results_dir / "fastp" / sample_name / f"{sample_name}.json", json.dumps({
            "summary": {
                "before_filtering": {"total_reads": 2000000, "total_bases": 200000000,
                                     "q30_rate": 0.95, "gc_content": 0.48,
                                     "read1_mean_length": 100},
                "after_filtering": {"total_reads": 1990000},
            },
            "duplication": {"rate": 0.12},
        }))

        # FastQC
        fastqc_dir = results_dir / "fastqc" / sample_name
        fastqc_dir.mkdir(parents=True, exist_ok=True)
        for read in ("1", "2"):
            archive = fastqc_dir / f"{sample_name}_{read}_fastqc.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                handle.writestr(
                    f"{sample_name}_{read}_fastqc/summary.txt",
                    f"PASS\tBasic Statistics\t{sample_name}_{read}.fq.gz\n"
                    f"WARN\tPer base sequence content\t{sample_name}_{read}.fq.gz\n",
                )

        # STAR
        star_dir = results_dir / "star" / sample_name
        _write(star_dir / f"{sample_name}_Log.final.out", "\n".join([
            "                          Number of input reads |\t2000000",
            "                      Average input read length |\t200",
            "                        Uniquely mapped reads number |\t1600000",
            "                             Uniquely mapped reads % |\t80.00%",
            "                           Average mapped length |\t198.50",
            "       Number of reads mapped to multiple loci |\t200000",
            "            % of reads mapped to multiple loci |\t10.00%",
            "  Number of reads mapped to too many loci |\t10000",
            "                 % of reads unmapped: too short |\t9.00%",
            "                    Mismatch rate per base, % |\t0.30%",
            "                       Number of splices: Total |\t500000",
            "            Number of splices: Annotated (sjdb) |\t490000",
        ]) + "\n")
        _tsv(star_dir / f"{sample_name}_SJ.out.tab", [], [
            ["chr1", 1000 + index, 2000 + index, 1, 1, 1, 30, 3, 45],
            ["chr1", 5000 + index, 6000 + index, 2, 2, 0, 12, 1, 20],
        ])
        # RSeQC: a first-strand (dUTP) library.
        _write(results_dir / "rseqc" /
               f"{sample_name}_Aligned.sortedByCoord.out_strandedness.txt",
               'This is PairEnd Data\n'
               'Fraction of reads failed to determine: 0.0300\n'
               'Fraction of reads explained by "1++,1--,2+-,2-+": 0.0400\n'
               'Fraction of reads explained by "1+-,1-+,2++,2--": 0.9300\n'
               '\nLibrary Type: fr-firststrand\n')

    # DESeq2, keyed by biomart symbol like the real counts matrix.
    _write(results_dir / "deseq2" / "deseq2_results.csv", "\n".join([
        ",baseMean,log2FoldChange,lfcSE,stat,pvalue,padj",
        "GAPDH,1000.5,0.12,0.05,2.4,0.016,0.05",
        "ACTB,2000.25,-1.5,0.2,-7.5,1e-13,1e-11",
        "MYC,500.0,3.2,0.4,8.0,1e-15,NA",
        "RPMS1_1,25.0,0.0,1.0,0.0,1.0,1.0",
    ]) + "\n")

    sleuth_dir = results_dir / "sleuth"
    _write(sleuth_dir / "sleuth_DETranscripts_wt.csv", "\n".join([
        '"",target_id,gene,pval,qval,b,se_b,mean_obs',
        '"1",ENST00000000002.3,ACTB,1e-10,1e-8,-1.4,0.2,6.5',
        '"2",ENST00000000001.1,GAPDH,0.02,0.09,0.11,0.05,8.1',
    ]) + "\n")
    _write(sleuth_dir / "sleuth_DETranscripts_lrt.csv", "\n".join([
        '"",target_id,gene,test_stat,pval,qval,mean_obs',
        '"1",ENST00000000002.3,ACTB,42.5,1e-11,1e-9,6.5',
        '"2",ENST00000000001.1,GAPDH,3.1,0.08,0.2,8.1',
    ]) + "\n")

    gsea_dir = results_dir / "gsea" / "h.all.v2023.2.Hs.symbols.Gsea.1700000000"
    _tsv(gsea_dir / "gsea_report_for_test_1700000000.tsv",
         ["NAME", "GS DETAILS", "SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val",
          "RANK AT MAX", "LEADING EDGE"],
         [["HALLMARK_APOPTOSIS", "details", 161, 0.45, 1.9, 0.001, 0.01, 0.02, 1200,
           "tags=25%"]])
    _tsv(gsea_dir / "HALLMARK_APOPTOSIS.tsv",
         ["PROBE", "GENE SYMBOL", "RANK IN GENE LIST", "RUNNING ES", "CORE ENRICHMENT"],
         [["MYC", "MYC", 10, 0.4, "Yes"], ["GAPDH", "GAPDH", 900, 0.1, "No"]])

    write_rmats(results_dir / "rmats" / "rmats_out")


SE_HEADER = ["ID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
             "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "ID.1",
             "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncFormLen",
             "SkipFormLen", "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference"]
MXE_HEADER = ["ID", "GeneID", "geneSymbol", "chr", "strand", "1stExonStart_0base",
              "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE",
              "downstreamES", "downstreamEE", "ID.1", "IJC_SAMPLE_1", "SJC_SAMPLE_1",
              "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncFormLen", "SkipFormLen", "PValue", "FDR",
              "IncLevel1", "IncLevel2", "IncLevelDifference"]


def write_rmats(rmats_out: Path, id_offset: int = 0) -> Path:
    """Write a minimal but structurally faithful rMATS output directory."""
    for mode in ("JC", "JCEC"):
        _tsv(rmats_out / f"SE.MATS.{mode}.txt", SE_HEADER, [
            [id_offset + 0, '"ENSG00000136997.18"', '"MYC"', "chr8", "+",
             1000, 1200, 500, 700, 1500, 1700, id_offset + 0,
             "50,60", "10,12", "20,25", "40,45", 300, 150, 0.001, 0.004,
             "0.71,0.72", "0.20,0.22", 0.505],
            [id_offset + 1, '"ENSG00000075624.17"', '"ACTB"', "chr7", "-",
             2000, 2200, 1500, 1700, 2500, 2700, id_offset + 1,
             "5,6", "80,82", "7,8", "75,77", 300, 150, 0.4, 0.6,
             "0.06,0.07", "0.09,0.10", -0.03],
        ])
        _tsv(rmats_out / f"MXE.MATS.{mode}.txt", MXE_HEADER, [
            [id_offset + 0, '"ENSG00000111640.15"', '"GAPDH"', "chr12", "+",
             100, 200, 300, 400, 50, 80, 500, 600, id_offset + 0,
             "30,31", "3,4", "10,11", "20,21", 300, 300, 0.002, 0.01,
             "0.90,0.89", "0.33,0.34", 0.56],
        ])
    _tsv(rmats_out / "summary.txt",
         ["EventType", "TotalEventsJC", "TotalEventsJCEC", "SignificantEventsJC"],
         [["SE", 2, 2, 1], ["MXE", 1, 1, 1]])
    return rmats_out


def _bigwig_rows(ctx: RunContext, results_dir: Path) -> list[dict]:
    """Stand-in bigWigs: wigToBigWig is not available off-cluster, but the manifest contract
    (path, strand, content, checksum, size) is exactly what the test needs to check."""
    rows = []
    for sample_name, condition, replicate in SAMPLES:
        sample_id = ctx.sample_id(sample_name)
        for content in ("unique", "uniquemulti"):
            for strand in ("plus", "minus"):
                path = (results_dir / "artifacts" / "bigwig" /
                        f"{sample_id}.{content}.{GENOME_BUILD}.{strand}.bw")
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(f"bigwig:{sample_id}:{content}:{strand}".encode())
                rows.append({
                    "run_id": ctx.run_id, "sample_id": sample_id,
                    "file_path": path.relative_to(results_dir).as_posix(),
                    "content": content, "strand": strand, "genome_build": GENOME_BUILD,
                    "normalization": "cpm", "scale_factor": 0.625,
                    "sha256": sha256_file(path), "bytes": path.stat().st_size,
                })
    return rows


def _signal_rows(ctx: RunContext) -> pd.DataFrame:
    """Signal consistent with the annotated strand of each control gene."""
    rows = []
    for sample_name, *_ in SAMPLES:
        sample_id = ctx.sample_id(sample_name)
        for _transcript, _symbol, gene_id, gene_strand in GENES:
            for strand in ("plus", "minus"):
                on_strand = strand == gene_strand
                total = 10000.0 if on_strand else 100.0
                rows.append({
                    "run_id": ctx.run_id, "sample_id": sample_id, "strand": strand,
                    "gene_id": gene_id.split(".")[0], "feature_set": "gene_body",
                    "mean_coverage": total / 1000, "max_coverage": total / 100,
                    "sum_coverage": total, "covered_bases": 900, "feature_length": 1000,
                })
    return pd.DataFrame(rows, columns=SIGNAL_OVER_GENE.column_names)


def build_manifest(results_dir: Path, run_id: str = RUN_ID) -> dict:
    manifest = {
        "run_id": run_id, "experiment_id": EXPERIMENT_ID, "pipeline_version": "1.4.0",
        "git_commit": "0" * 40, "git_dirty": False,
        "docker_image_digest": "sha256:" + "a" * 64,
        "docker_image_digest_source": "test", "genome_build": GENOME_BUILD,
        "annotation_version": "gencode_v44", "reference_dir": "/data/referenceFiles",
        "reference_dir_sha256": "b" * 64, "investigator": "ethan", "library_layout": "PE",
        "run_start_utc": "2026-01-15T10:00:00Z", "run_end_utc": "2026-01-15T18:00:00Z",
        "exit_status": "success", "gene_id_source": "ensembl",
        "tool_versions": {"STAR": "2.7.11b", "kallisto": "0.52.0"},
        "parameters": {"deseq2_alpha": 0.05},
    }
    (results_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


@pytest.fixture(scope="session")
def built_results(tmp_path_factory) -> dict:
    """A complete results directory produced by the real table writers."""
    root = tmp_path_factory.mktemp("run")
    results_dir = root / "output"
    reference_dir = write_reference(root / "referenceFiles")
    create_layout(results_dir)
    metadata = write_metadata(results_dir)
    write_tool_output(results_dir)

    annotation = Annotation(
        reference_dir / GENOME_BUILD / "biomart" / f"{GENOME_BUILD}.mart_export.txt")
    ctx = RunContext(
        run_id=RUN_ID, experiment_id=EXPERIMENT_ID, results_dir=results_dir,
        genome_build=GENOME_BUILD,
        sample_ids=dict(zip(metadata["Sample name"], metadata["sample_id"])),
    )
    sample_names = [name for name, *_ in SAMPLES]

    tables_qc.write_qc_metrics(ctx, results_dir, sample_names, GENOME_BUILD)
    tables_qc.update_metadata_strandedness(
        results_dir, {name: "fr-firststrand" for name in sample_names})
    tables_de.write_de_gene(ctx, results_dir, annotation)
    tables_de.write_transcript_de(ctx, results_dir, annotation)
    tables_expression.write_expression_tables(
        ctx, results_dir, sample_names, GENOME_BUILD, annotation)
    tables_gsea.write_gene_set_enrichment(ctx, results_dir)
    tables_rmats.write_splicing_tables(ctx, results_dir, results_dir / "rmats" / "rmats_out")
    tables_junction.write_junctions(ctx, results_dir, sample_names)

    bigwig_rows = _bigwig_rows(ctx, results_dir)
    write_table(results_dir, "bigwig_manifest",
                pd.DataFrame(bigwig_rows, columns=BIGWIG_MANIFEST.column_names))
    write_table(results_dir, "signal_over_gene", _signal_rows(ctx))

    # Raw tool output belongs under artifacts/, as organize_artifacts() would leave it.
    for stage in ("kallisto", "deseq2", "sleuth", "gsea", "fastp", "fastqc", "rseqc"):
        source = results_dir / stage
        if source.is_dir():
            source.rename(results_dir / "artifacts" / source.name)

    checksums.write_checksums(results_dir,
                              known={row["file_path"]: row["sha256"] for row in bigwig_rows})
    build_manifest(results_dir)
    return {"results_dir": results_dir, "reference_dir": reference_dir, "ctx": ctx,
            "annotation": annotation, "bigwig_rows": bigwig_rows,
            "sample_names": sample_names}


@pytest.fixture
def results_dir(built_results, tmp_path) -> Path:
    """A writable copy, so a test that corrupts a table cannot affect its neighbours."""
    import shutil
    target = tmp_path / "output"
    shutil.copytree(built_results["results_dir"], target)
    return target
