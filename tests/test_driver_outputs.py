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
