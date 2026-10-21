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
