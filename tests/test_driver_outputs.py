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


def _organized(tmp_path):
    results_dir = tmp_path / "output"
    _pre_organize_tree(results_dir)
    rnaseq.organize_artifacts(results_dir, [SAMPLE])
    return results_dir


def test_rmats_summary_sits_with_the_jcec_results(tmp_path):
    """F2: summary.txt counts the events in the JCEC files, so it goes in JCEC/ with them, as
    it did in the original 1_JCEC folder. JC files and intermediates stay one level up."""
    jcec = _organized(tmp_path) / "artifacts" / "rmats_raw" / "JCEC"
    names = sorted(path.name for path in jcec.iterdir())
    assert names == sorted([f"{e}.MATS.JCEC.txt" for e in RMATS_EVENTS] + ["summary.txt"])
    assert (jcec.parent / "SE.MATS.JC.txt").is_file()


def test_kallisto_output_is_not_nested_twice(tmp_path):
    """F3: artifacts/kallisto/<sample>/, as in the original 04_abundance_estimation/kallisto/,
    not artifacts/kallisto/kallisto/<sample>/."""
    kallisto = _organized(tmp_path) / "artifacts" / "kallisto"
    assert (kallisto / SAMPLE / f"{SAMPLE}_hg38" / "abundance.tsv").is_file()
    assert not (kallisto / "kallisto").exists()
