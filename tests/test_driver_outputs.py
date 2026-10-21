"""Researcher-facing behaviour of rnaseq.py itself: which analyses run, and where their output
lands. Each test pins a behaviour of the original (pre-refactor) pipeline that the refactor had
dropped; the finding ID from the October 2026 drift audit is in each docstring.
"""

import rnaseq


def test_gsea_runs_every_gene_set_collection(tmp_path):
    """F1: the full gene_ontology GMT is run alongside the BP/MF splits. It is the only file
    with the cellular-component sets, and the original pipeline ran every *.gmt."""
    names = ["biological_process.v2024.1.Hs.symbols.gmt", "gene_ontology.v2024.1.Hs.symbols.gmt",
             "h.all.v2024.1.Hs.symbols.gmt", "molecular_function.v2024.1.Hs.symbols.gmt"]
    for name in names:
        (tmp_path / name).write_text("SET\tdesc\tGENE\n")
    (tmp_path / "README.txt").write_text("not a gene set\n")
    assert [path.name for path in rnaseq.gsea_gmt_files(tmp_path)] == names


SAMPLE = "CNTL-1_L02"
EXPERIMENT = "Mutu_DDX5-KO_2026-09-30"
RMATS_EVENTS = ("A3SS", "A5SS", "MXE", "RI", "SE")


def _touch(path, text="x\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _pre_organize_tree(results_dir):
    """The working tree rnaseq.py has just before organize_artifacts(): every stage's output
    where the stage wrote it."""
    star = results_dir / "star" / SAMPLE
    for name in ("Aligned.sortedByCoord.out.bam", "Aligned.sortedByCoord.out.bam.bai",
                 "SJ.out.tab", "SJ.out.tab.positive.bed", "SJ.out.tab.negative.bed",
                 "Chimeric.out.junction", "Log.final.out", "Log.out"):
        _touch(star / f"{SAMPLE}_{name}")
    _touch(results_dir / "kallisto" / SAMPLE / f"{SAMPLE}_hg38" / "abundance.tsv")
    for stage in ("fastp", "fastqc", "rseqc"):
        _touch(results_dir / stage / SAMPLE / "report.txt")
    for name in ("deseq2_results_genes.xlsx", "deseq2_results_genes.tsv", "deseq2_results.csv",
                 "deseq2_counts_matrix.tsv", "deseq2_metadata.tsv", "volcano_plot.png",
                 "PCA_plot.png"):
        _touch(results_dir / "deseq2" / name)
    for name in ("sleuth_results_wald_transcripts.xlsx", "sleuth_DETranscripts_wt.csv",
                 "sleuth_metadata.tsv"):
        _touch(results_dir / "sleuth" / name)
    _touch(results_dir / "gsea" / "test.cls")
    _touch(results_dir / "gsea" / "h.all.Gsea.1700000000" / "index.html")
    rmats_out = results_dir / "rmats" / "rmats_out"
    for event in RMATS_EVENTS:
        _touch(rmats_out / f"{event}.MATS.JCEC.txt")
        _touch(rmats_out / f"{event}.MATS.JC.txt")
    for name in ("summary.txt", "barChart_rmatsASE_fdr0.05.png"):
        _touch(rmats_out / name)
    _touch(results_dir / "rmats" / "testBAMs.txt")
    for name in ("topGenes.tsv", "alignmentSummary.tsv", "mycoplasma_report.tsv"):
        _touch(results_dir / name)
    _touch(results_dir / "reports" / "report.html")
    _touch(results_dir / "reports" / "de_gene.xlsx")


def _organized(tmp_path):
    results_dir = tmp_path / "output"
    _pre_organize_tree(results_dir)
    rnaseq.organize_artifacts(results_dir, [SAMPLE], EXPERIMENT)
    return results_dir


def test_rmats_summary_sits_with_the_jcec_results(tmp_path):
    """F2: summary.txt counts the events in the JCEC files, so it goes in JCEC/ with them, as
    it did in the original 1_JCEC folder. JC files and intermediates stay one level up."""
    jcec = _organized(tmp_path) / "artifacts" / "rmats_raw" / "JCEC"
    names = sorted(path.name for path in jcec.iterdir())
    expected = [f"{e}.MATS.JCEC.txt" for e in RMATS_EVENTS] + ["summary.txt"]
    assert names == sorted(f"{EXPERIMENT}_test_vs_cntl_{name}" for name in expected)
    assert (jcec.parent / "SE.MATS.JC.txt").is_file()


def test_kallisto_output_is_not_nested_twice(tmp_path):
    """F3: artifacts/kallisto/<sample>/, as in the original 04_abundance_estimation/kallisto/,
    not artifacts/kallisto/kallisto/<sample>/."""
    kallisto = _organized(tmp_path) / "artifacts" / "kallisto"
    assert (kallisto / SAMPLE / f"{SAMPLE}_hg38" / "abundance.tsv").is_file()
    assert not (kallisto / "kallisto").exists()


def test_de_spreadsheets_are_kept_apart_from_prep_files(tmp_path):
    """F5: the DESeq2/Sleuth XLSX files sit at the top of their folder and the intermediates in
    prep_files/, as in the original 05_differential_expression/<tool>/prep_files/."""
    artifacts = _organized(tmp_path) / "artifacts"
    for tool in ("deseq2", "sleuth"):
        top = sorted(path.suffix for path in (artifacts / tool).iterdir() if path.is_file())
        assert top == [".xlsx"]
        assert any((artifacts / tool / "prep_files").iterdir())
    assert (artifacts / "deseq2" / "prep_files" / f"{EXPERIMENT}_deseq2_results.csv").is_file()


def test_researcher_files_carry_the_experiment_name(tmp_path):
    """F4: the original rename passes' names, so files copied out of the run still say which
    experiment they are from."""
    results_dir = _organized(tmp_path)
    artifacts = results_dir / "artifacts"
    contrast = f"{EXPERIMENT}_test_vs_cntl_"
    assert sorted(path.name for path in (results_dir / "reports").glob("*.png")) == [
        f"{contrast}PCA_plot.png", f"{contrast}barChart_rmatsASE_fdr0.05.png",
        f"{contrast}volcano_plot.png",
    ]
    assert (artifacts / "deseq2" / f"{EXPERIMENT}_deseq2_results_genes.xlsx").is_file()
    assert (artifacts / "sleuth" / f"{EXPERIMENT}_sleuth_results_wald_transcripts.xlsx").is_file()
    prep = artifacts / "deseq2" / "prep_files"
    assert (prep / "deseq2_metadata.tsv").is_file()  # the original left this one unprefixed
    assert (prep / f"{EXPERIMENT}_deseq2_counts_matrix.tsv").is_file()
    assert (artifacts / "sleuth" / "prep_files" / f"{EXPERIMENT}_sleuth_metadata.tsv").is_file()
    assert (results_dir / "reports" / f"{EXPERIMENT}_de_gene.xlsx").is_file()
    # The original never renamed the summary tables or the HTML report.
    assert (results_dir / "reports" / "topGenes.tsv").is_file()
    assert (results_dir / "reports" / "report.html").is_file()


def test_prefixing_is_idempotent(tmp_path):
    """A resumed run must not stack a second prefix on files that already have one."""
    results_dir = _organized(tmp_path)
    before = sorted(p.relative_to(results_dir) for p in results_dir.rglob("*"))
    rnaseq.prefix_experiment_name(results_dir, EXPERIMENT)
    assert sorted(p.relative_to(results_dir) for p in results_dir.rglob("*")) == before


HUMAN_READABLE_STEPS = (
    "collate_tpms", "merge_deseq2_results", "reform_deseq2_results", "collate_transcript_tpms",
    "merge_sleuth_wald_results", "merge_sleuth_lrt_results", "reformat_sleuth_results",
    "check_mycoplasma_sequences", "order_mycoplasma_report", "prepare_volcano",
    "convert_SJtab_to_bed", "report_top_genes_and_transcripts", "report_top_gene_tpms",
    "summarize_alignments", "summarize_transcript_coverage", "generate_results_html",
)


def test_a_failed_researcher_output_is_returned_as_an_error(tmp_path, monkeypatch):
    """F6: the run carries on past a failed volcano plot, but the failure comes back to main()
    to be recorded, rather than ending as a WARNING line in pipeline.log."""
    for name in HUMAN_READABLE_STEPS:
        monkeypatch.setattr(rnaseq, name, lambda *args: None)
    monkeypatch.setattr(rnaseq.tables_de, "write_de_gene_xlsx", lambda *args: None)

    def broken_volcano(results_dir):
        raise RuntimeError("volcano.py exited 1")
    monkeypatch.setattr(rnaseq, "prepare_volcano", broken_volcano)

    later_steps_ran = []
    monkeypatch.setattr(rnaseq, "generate_results_html",
                        lambda *args: later_steps_ran.append("report"))
    errors = rnaseq.write_human_readable_outputs(
        tmp_path / "metadata.tsv", tmp_path, [SAMPLE], "hg38")
    assert errors == ["volcano plot: RuntimeError: volcano.py exited 1"]
    assert later_steps_ran == ["report"]


def test_output_errors_are_written_to_the_error_log(tmp_path):
    rnaseq.record_output_errors(tmp_path, ["volcano plot: RuntimeError: boom",
                                           "GSEA gene_ontology.gmt: gsea-cli failed"])
    lines = (tmp_path / "logs" / "errors.log").read_text().splitlines()
    assert [line.split("\t", 1)[1] for line in lines] == [
        "output\tvolcano plot: RuntimeError: boom",
        "output\tGSEA gene_ontology.gmt: gsea-cli failed",
    ]


def test_no_output_errors_writes_no_error_log(tmp_path):
    rnaseq.record_output_errors(tmp_path, [])
    assert not (tmp_path / "logs").exists()


UNCLASSIFIABLE_RSEQC = (
    "This is PairEnd Data\n"
    "Fraction of reads failed to determine: 0.9100\n"
    'Fraction of reads explained by "1++,1--,2+-,2-+": 0.0500\n'
    'Fraction of reads explained by "1+-,1-+,2++,2--": 0.0400\n'
)


def test_undetermined_strandedness_is_analysed_as_unstranded(tmp_path, monkeypatch, capsys):
    """F9: the original pipeline carried on as unstranded when RSeQC could not tell; the
    refactor stopped the run. Carry on again, but say so."""
    monkeypatch.setattr(rnaseq, "REFERENCE_DIR", tmp_path / "ref", raising=False)
    monkeypatch.setattr(rnaseq, "species_name", "hg38", raising=False)
    monkeypatch.setattr(rnaseq.subprocess, "run", lambda *args, **kwargs:
                        rnaseq.subprocess.CompletedProcess(args, 0, UNCLASSIFIABLE_RSEQC))
    bam = tmp_path / f"{SAMPLE}_Aligned.sortedByCoord.out.bam"

    assert rnaseq.infer_strandedness(bam, "hg38", tmp_path) == "unstranded"
    assert "could not determine strandedness" in capsys.readouterr().out

    strandedness = rnaseq.strandedness_module
    record = strandedness.rseqc_path(tmp_path, SAMPLE).read_text()
    assert record.startswith(UNCLASSIFIABLE_RSEQC)
    assert "Library Type: unstranded" in record and strandedness.FALLBACK_LINE in record
    # The rest of the run reads the analysed type back; metadata.tsv gets the measured one.
    assert strandedness.read_library_type(tmp_path, SAMPLE) == "unstranded"
    assert strandedness.read_measured_library_type(tmp_path, SAMPLE) == "undetermined"
    assert rnaseq.strandedness_label(record).startswith("unstranded (RSeQC could not")


def test_report_uses_the_original_strandedness_words():
    """F7: the original report said rf-stranded / fr-stranded / unstranded. The canonical names
    ride along in brackets; fr-firststrand alone is easily misread as the old fr-stranded."""
    def rseqc(forward, reverse):
        return (UNCLASSIFIABLE_RSEQC.replace("0.9100", "0.0100")
                .replace("0.0500", forward).replace("0.0400", reverse))
    assert rnaseq.strandedness_label(rseqc("0.0300", "0.9600")) == "rf-stranded (fr-firststrand)"
    assert rnaseq.strandedness_label(rseqc("0.9600", "0.0300")) == "fr-stranded (fr-secondstrand)"
    assert rnaseq.strandedness_label(rseqc("0.4900", "0.5000")) == "unstranded"


def test_manifest_deseq2_alpha_matches_what_r_uses():
    """F15: the manifest said 0.05 while the R script used DESeq2's default of 0.1. If the R
    call ever sets alpha, this has to change with it."""
    script = (rnaseq.Path(rnaseq.__file__).parent / "rnaseq_helper_scripts" /
              "deseq2_analysis_ercc.R").read_text()
    assert "res <- results(dds)\n" in script
    assert "alpha" not in script
    assert rnaseq.DESEQ2_ALPHA == 0.1


def test_manifest_records_which_kallisto_index_was_used(tmp_path, monkeypatch):
    """F10: local (macOS) runs use the pre-container ORIGINAL index and Docker runs the one
    rebuilt in April 2026. Abundances can differ slightly between the two, so the manifest
    names the index and its checksum."""
    monkeypatch.setattr(rnaseq, "REFERENCE_DIR", tmp_path, raising=False)
    index_dir = tmp_path / "hg38" / "kallisto"
    index_dir.mkdir(parents=True)
    (index_dir / "hg38.kallisto.idx").write_bytes(b"rebuilt index")

    monkeypatch.setattr(rnaseq.platform, "system", lambda: "Linux")
    record = rnaseq.kallisto_index_record("hg38")
    assert record["species"] == {"path": str(index_dir / "hg38.kallisto.idx"),
                                 "sha256": rnaseq.sha256_file(index_dir / "hg38.kallisto.idx")}
    assert record["contaminants"] == {
        "path": str(index_dir / "all_contaminants.kallisto.idx"), "sha256": "missing"}

    monkeypatch.setattr(rnaseq.platform, "system", lambda: "Darwin")
    assert rnaseq.kallisto_index_path("hg38") == index_dir / "hg38_ORIGINAL.kallisto.idx"
    assert rnaseq.kallisto_index_path("hg38", is_contaminant_run=True) == \
        index_dir / "all_contaminants_ORIGINAL.kallisto.idx"


def test_splicing_bar_chart_draws(tmp_path):
    """I1 moved the matplotlib imports into this function; draw it once where matplotlib is
    installed (the pipeline image) so a missing name cannot slip through again."""
    pytest = __import__("pytest")
    pytest.importorskip("matplotlib")
    rmats_out = tmp_path / "rmats" / "rmats_out"
    rmats_out.mkdir(parents=True)
    for event in RMATS_EVENTS:
        (rmats_out / f"{event}.MATS.JCEC.txt").write_text(
            "FDR\tIncLevelDifference\n0.001\t0.3\n0.001\t-0.2\n0.5\t0.1\n")
    rnaseq.generate_splicing_bar_chart(tmp_path)
    assert sorted(path.name for path in rmats_out.glob("*.png"))
